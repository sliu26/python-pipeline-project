dependencies: Python, Snakemake, Bowtie2, Samtools, SPAdes, BLAST, Biopython

# python-pipeline-project
step 1: retrieved sra and used wget to retrieve the srr 
    wget [sra normalized data link from ncbi]

step 2: used fasterq-dump to get paired end reads
    (fasterq-dump ./[corresponding srr])

step 3: used head to reduce size of fastq files and created sample fastq and some renaming after to keep consistent 
    (head -n 40000 SRR56600*.fastq > SRR56600*.fastq)

bowtie: 

4: went to ncbi and retrieve HCMV genome (GCF_000845245.1), and run datasets; https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000845245.1/ 
    [download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report]

5: unzip ncbi dataset.zip and use bowtie to create index to map to: (later used snakemake to download for me)
- snakemake -c1 build_index
    [bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMV]

6: 