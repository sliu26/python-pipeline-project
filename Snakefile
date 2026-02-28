samples = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

rule all:
    input:
        "Liu_PipelineReport.txt"        #final output

#build index. ref folder stores all, assuming that the same reference index will be used and can be hard coded for index? 
rule build_index:
    input:
        "ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna"
    output:
        expand("ref/HCMV.{ext}.bt2", ext=["1","2","3","4","rev.1","rev.2"])
    shell:
        "bowtie2-build {input} ref/HCMV"

#count length before 
rule count_before:
    input:
        "data/{sample}_1.fastq"     #only need one of the paired files for counting
    output:
        "results/counts/{sample}_before.txt"
    shell:      #each read pair consists of 4 lines in FASTQ, so we divide total lines by 4
        """
        echo $(( $(wc -l < {input}) / 4 )) > {output}
        """
#maps sample reads to HCMV genome
rule map_reads:
    input:
        r1="data/{sample}_1.fastq",
        r2="data/{sample}_2.fastq",
        index=expand("ref/HCMV.{ext}.bt2", ext=["1","2","3","4","rev.1","rev.2"])   #require all Bowtie2 index files before mapping
    output:
        "results/mapped/{sample}.bam"
    shell:  # -F 4 removes unmapped reads (line 39), output is a filtered BAM containing only mapped reads
        """
        bowtie2 -x ref/HCMV \
            -1 {input.r1} \
            -2 {input.r2} | \
        samtools view -b -F 4 - > {output}
        """

#count reads after mapping
rule count_after:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/counts/{sample}_after.txt"
    shell:      #divide samtools count by 2 to get read pairs
        """
        echo $(( $(samtools view -c {input}) / 2 )) > {output}
        """
#convert mapped BAM to paired FASTQ
rule bam_to_fastq:
    input:
        "results/mapped/{sample}.bam"
    output:
        r1="results/filtered_fastq/{sample}_1.fastq",
        r2="results/filtered_fastq/{sample}_2.fastq"
    shell:
        """
        mkdir -p results/filtered_fastq
        samtools sort -n {input} -o results/mapped/{wildcards.sample}_sorted.bam
        samtools fastq \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null -n \
            results/mapped/{wildcards.sample}_sorted.bam
        """
#assemble mapped reads using SPAdes
rule assemble_spades:
    input:
        r1="results/filtered_fastq/{sample}_1.fastq",
        r2="results/filtered_fastq/{sample}_2.fastq"
    output:
        "results/assembly/{sample}/contigs.fasta"
    shell: #uses kmer size 99
        """
        mkdir -p results/assembly/{wildcards.sample}
        spades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -k 99 \
            -o results/assembly/{wildcards.sample}
        """

rule assembly_stats: #calculate assembly stat
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/assembly_stats/{sample}_assembly.txt"
    run:        #use biopython SeqIO to parse FASTA
        from Bio import SeqIO
        import os

        os.makedirs("results/assembly_stats", exist_ok=True)

        count = 0
        total_length = 0
        #iterate over all contigs and only pick out > 1000 bp
        for record in SeqIO.parse(input[0], "fasta"):
            length = len(record.seq)
            if length > 1000:
                count += 1
                total_length += length

        with open(output[0], "w") as out:
            out.write(f"{count}\t{total_length}")

rule longest_contig:
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/blast/{sample}_longest.fasta"
    run:
        from Bio import SeqIO
        import os
        os.makedirs("results/blast", exist_ok=True)
        longest = None
        max_len = 0
        for record in SeqIO.parse(input[0], "fasta"):
            if len(record.seq) > max_len:
                max_len = len(record.seq)
                longest = record
        if longest:     #used in blast query later
            SeqIO.write(longest, output[0], "fasta")

rule blast_longest:
    input:
        "results/blast/{sample}_longest.fasta"
    output:
        "results/blast/{sample}_blast.txt"
    params:
        db="data/blast_db/betaherpes_db"  
    shell:
        """
        blastn \
            -query {input} \
            -db {params.db} \
            -max_hsps 1 \
            -max_target_seqs 5 \
            > {output}
        """

rule build_report:
    input:
        before=expand("results/counts/{sample}_before.txt", sample=samples),
        after=expand("results/counts/{sample}_after.txt", sample=samples),
        assembly=expand("results/assembly_stats/{sample}_assembly.txt", sample=samples),
        blast=expand("results/blast/{sample}_blast.txt", sample=samples),

    output:
        "Liu_PipelineReport.txt"
    run:
        with open(output[0], "w") as out:
            for s in samples:
                with open(f"results/counts/{s}_before.txt") as b:
                    before = b.read().strip()
                
                with open(f"results/counts/{s}_after.txt") as a:
                    after = a.read().strip()

                with open(f"results/assembly_stats/{s}_assembly.txt") as f:
                    count, total = f.read().strip().split("\t")

                out.write(f"Sample {s} had {before} read pairs before and {after} read pairs after Bowtie2 filtering.\n")

                out.write(f"In the assembly of sample {s}, there are {count} contigs > 1000 bp and {total} total bp.\n\n")

                blast_file = f"results/blast/{s}_blast.txt"

                with open(blast_file) as bf:
                    lines = bf.readlines()

                if lines:
                    out.write("top hits for the longest contig:\n")
                    for line in lines:
                        parts = line.strip().split("\t")
                        if len(parts) != 10:
                            continue  # skip malformed lines
                        sacc, pident, length, qstart, qend, sstart, send, bitscore, evalue, stitle = parts
                        out.write(f"{sacc}\t{pident}\t{length}\t{qstart}\t{qend}\t{sstart}\t{send}\t{bitscore}\t{evalue}\t{stitle}\n")
                else:
                    out.write("no hits.\n")

                out.write("\n")

rule clean:
    shell:
        """
        rm Liu_PipelineReport.txt
        rm -rf results/*
        rm -rf ref/*
        """
