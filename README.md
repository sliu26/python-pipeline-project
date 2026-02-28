dependencies: Python, Snakemake, Bowtie2, Samtools, SPAdes, BLAST, Biopython

# python-pipeline-project
step 1: retrieved sra and used wget to retrieve the srr 
    wget [sra normalized data link from ncbi]

step 2: used fasterq-dump to get paired end reads
    (fasterq-dump ./[corresponding srr])

step 3: used head to reduce size of fastq files and created sample fastq and some renaming after to keep consistent 
    (head -n 40000 SRR56600*.fastq > SRR56600*.fastq)

bowtie: 
step 4: went to ncbi and retrieve HCMV genome (GCF_000845245.1), and run datasets; https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000845245.1/ 
    [download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report]

step 5: unzip ncbi dataset.zip and use bowtie to create index to map to: (later used snakemake to download for me)
- snakemake -c1 build_index
    [bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMV]

step 6: From here I started to do everything in snakemake. Instead of making a clean rule, I used forcerun instead: counted reads before mapping to HCMV genome to track initial read numbers
- snakemake -c4 --forcerun

step 7: mapped reads to HCMV genome using Bowtie2 and created BAM files

step 8: counted reads after mapping to see how many mapped

step 9: converted mapped BAMs to filtered paired-end FASTQ files for assembly

step 10: assembled reads using SPAdes with k=99

step 11: calculated assembly statistics (number of contigs >1000 bp and total bp)

step 12: made blast db manually on ncbi: searched 'Betaherpesvirinae[Organism] AND complete genome' and set filter source database to refseq to ensure complete db. used top 5 options and exported as FASTA file
https://www.ncbi.nlm.nih.gov/nuccore 

step 13: extracted longest contig for each sample for BLAST analysis

step 14: ran BLAST of longest contigs against Betaherpesvirinae database

step 15: built final report summarizing read counts, assembly stats, and BLAST results
