# SNV Detection from FASTQ Files Using Docker: A Guide for Beginners

This guide will walk you through detecting Single Nucleotide Variants (SNVs) from FASTQ sequencing files using Docker containers. This approach keeps your system clean by installing all required bioinformatics tools inside containers rather than directly on your computer.

## Table of Contents
- [Docker Basics](#docker-basics)
- [Project Setup](#project-setup)
- [Getting the Reference Genome](#getting-the-reference-genome)
- [Aligning Reads](#aligning-reads)
- [Processing Alignments](#processing-alignments)
- [Calling Variants](#calling-variants)
- [Annotating Variants](#annotating-variants-optional)
- [Helpful Docker Commands](#helpful-docker-commands)
- [Complete Workflow Script](#complete-workflow-script)

## Docker Basics

Docker is a platform that uses containers to isolate applications and their dependencies from the host system. Think of containers as lightweight virtual machines.

**Key Concepts:**
- **Image**: A template with pre-installed software (like a snapshot of a computer)
- **Container**: A running instance of an image
- **Volume**: A way to share files between your computer and containers

### Why Use Docker for Bioinformatics?
- No need to install complex bioinformatics tools directly
- Ensures everyone runs the exact same software versions
- Keeps your system clean and organized
- Works consistently across different computers

## Project Setup

### 1. Install Docker

Download and install Docker Desktop from [docker.com](https://www.docker.com/products/docker-desktop/)

Verify installation:
```bash
docker --version
```

### 2. Create Project Directory

```bash
# Create project directory and data subdirectory
mkdir -p ~/snv_project/data
cd ~/snv_project
```

### 3. Copy FASTQ Files

```bash
# Copy your FASTQ files to the data directory
cp path/to/your/fastq/files/*.fastq.gz ~/snv_project/data/
```

### 4. Pull Required Docker Images

```bash
# Download the necessary bioinformatics tools as Docker images
docker pull biocontainers/bwa:v0.7.17_cv1
docker pull quay.io/biocontainers/samtools:1.15--h1170115_1
docker pull quay.io/biocontainers/bcftools:1.15--h0ea216a_2
```

## Getting the Reference Genome

### 1. Download Human Reference Genome

```bash
# Download and decompress the reference genome
docker run -v $(pwd)/data:/data -w /data --rm ubuntu:20.04 bash -c \
  "apt-get update && apt-get install -y wget && \
   wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && \
   gzip -d hg38.fa.gz"
```

**What this command does:**
- `-v $(pwd)/data:/data`: Links your data folder to /data in the container
- `-w /data`: Sets working directory inside container
- `--rm`: Removes container after it finishes
- The commands download and decompress the human genome

### 2. Index the Reference Genome

```bash
# Create BWA index
docker run -v $(pwd)/data:/data -w /data --rm biocontainers/bwa:v0.7.17_cv1 \
  bwa index /data/hg38.fa

# Create samtools index
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/samtools:1.15--h1170115_1 \
  samtools faidx /data/hg38.fa
```

These steps prepare the reference genome for fast alignment and variant calling.

## Aligning Reads

### 1. Align FASTQ Files to Reference

**For paired-end reads:**
```bash
docker run -v $(pwd)/data:/data -w /data --rm biocontainers/bwa:v0.7.17_cv1 \
  bwa mem -t 4 /data/hg38.fa /data/sample_R1.fastq.gz /data/sample_R2.fastq.gz > /data/aligned.sam
```

**For single-end reads:**
```bash
docker run -v $(pwd)/data:/data -w /data --rm biocontainers/bwa:v0.7.17_cv1 \
  bwa mem -t 4 /data/hg38.fa /data/sample.fastq.gz > /data/aligned.sam
```

**Note**: Replace `sample_R1.fastq.gz` and `sample_R2.fastq.gz` with your actual FASTQ filenames.

## Processing Alignments

### 1. Convert, Sort, and Index BAM Files

```bash
# Process the alignment file
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/samtools:1.15--h1170115_1 bash -c \
  "samtools view -b /data/aligned.sam > /data/aligned.bam && \
   samtools sort /data/aligned.bam -o /data/aligned.sorted.bam && \
   samtools index /data/aligned.sorted.bam"
```

This command:
1. Converts SAM to compressed BAM format
2. Sorts the reads by position in the genome
3. Creates an index for fast access

### 2. Generate Alignment Statistics

```bash
# Check alignment quality
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/samtools:1.15--h1170115_1 \
  samtools flagstat /data/aligned.sorted.bam > /data/alignment_stats.txt
```

This produces a summary of your alignment quality.

## Calling Variants

### 1. Detect Variants

```bash
# Identify variants (SNVs, indels)
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/bcftools:1.15--h0ea216a_2 bash -c \
  "bcftools mpileup -f /data/hg38.fa /data/aligned.sorted.bam | \
   bcftools call -mv -o /data/variants.vcf"
```

This command identifies positions where your reads differ from the reference genome.

### 2. Filter Variants

```bash
# Filter for high-quality variants
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/bcftools:1.15--h0ea216a_2 \
  bcftools filter -i 'QUAL>20 && DP>10' /data/variants.vcf > /data/filtered_variants.vcf
```

This removes low-quality variant calls, keeping only those with:
- Quality score greater than 20
- Read depth greater than 10

### 3. Examine Variant Statistics

```bash
# Get statistics on called variants
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/bcftools:1.15--h0ea216a_2 \
  bcftools stats /data/filtered_variants.vcf > /data/variant_stats.txt
```

## Annotating Variants (Optional)

### 1. Install and Run SnpEff

```bash
# Pull the SnpEff container
docker pull quay.io/biocontainers/snpeff:5.0--hdfd78af_1

# Download database and annotate variants
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/snpeff:5.0--hdfd78af_1 bash -c \
  "snpEff download -v GRCh38.99 && \
   snpEff ann GRCh38.99 /data/filtered_variants.vcf > /data/annotated_variants.vcf"
```

This adds information about the potential impact of each variant.

## Helpful Docker Commands

```bash
# List all downloaded Docker images
docker images

# See running containers
docker ps

# See all containers (including stopped ones)
docker ps -a

# Remove a container
docker rm [container-id]

# Remove an image
docker rmi [image-id]

# Clean up all unused containers
docker container prune

# Clean up all unused images
docker image prune
```

## Complete Workflow Script

You can save the following as `run_snv_analysis.sh` to automate the whole process:

```bash
#!/bin/bash
# SNV Detection Workflow

# Set file names - CHANGE THESE to match your files
FASTQ_R1="sample_R1.fastq.gz"
FASTQ_R2="sample_R2.fastq.gz"  
# Remove the FASTQ_R2 line if using single-end data

# Create data directory if needed
mkdir -p data

echo "=== Step 1: Downloading reference genome ==="
docker run -v $(pwd)/data:/data -w /data --rm ubuntu:20.04 bash -c \
  "if [ ! -f /data/hg38.fa ]; then 
     apt-get update && apt-get install -y wget && \
     wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && \
     gzip -d hg38.fa.gz; 
   fi"

echo "=== Step 2: Indexing reference genome ==="
docker run -v $(pwd)/data:/data -w /data --rm biocontainers/bwa:v0.7.17_cv1 \
  bwa index /data/hg38.fa

docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/samtools:1.15--h1170115_1 \
  samtools faidx /data/hg38.fa

echo "=== Step 3: Aligning reads to reference ==="
# For paired-end reads
docker run -v $(pwd)/data:/data -w /data --rm biocontainers/bwa:v0.7.17_cv1 \
  bwa mem -t 4 /data/hg38.fa /data/$FASTQ_R1 /data/$FASTQ_R2 > /data/aligned.sam

# For single-end reads (uncomment if needed)
# docker run -v $(pwd)/data:/data -w /data --rm biocontainers/bwa:v0.7.17_cv1 \
#   bwa mem -t 4 /data/hg38.fa /data/$FASTQ_R1 > /data/aligned.sam

echo "=== Step 4: Processing alignments ==="
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/samtools:1.15--h1170115_1 bash -c \
  "samtools view -b /data/aligned.sam > /data/aligned.bam && \
   samtools sort /data/aligned.bam -o /data/aligned.sorted.bam && \
   samtools index /data/aligned.sorted.bam && \
   samtools flagstat /data/aligned.sorted.bam > /data/alignment_stats.txt"

echo "=== Step 5: Calling variants ==="
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/bcftools:1.15--h0ea216a_2 bash -c \
  "bcftools mpileup -f /data/hg38.fa /data/aligned.sorted.bam | \
   bcftools call -mv -o /data/variants.vcf"

echo "=== Step 6: Filtering variants ==="
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/bcftools:1.15--h0ea216a_2 \
  bcftools filter -i 'QUAL>20 && DP>10' /data/variants.vcf > /data/filtered_variants.vcf

echo "=== Step 7: Getting variant statistics ==="
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/bcftools:1.15--h0ea216a_2 \
  bcftools stats /data/filtered_variants.vcf > /data/variant_stats.txt

echo "=== Process complete! ==="
echo "Results are in the data directory:"
echo "- Filtered variants: data/filtered_variants.vcf"
echo "- Alignment stats: data/alignment_stats.txt"
echo "- Variant stats: data/variant_stats.txt"
```

Make the script executable and run it:
```bash
chmod +x run_snv_analysis.sh
./run_snv_analysis.sh
```

This script will run the entire pipeline, from downloading the reference genome to calling and filtering variants.

## Final Notes

- The first run will take longer as Docker downloads and indexes the reference genome
- You may need to adjust memory settings in Docker Desktop for large genomes
- For production use, consider adding error checking to the script
- Results will be saved in your data directory and persist even after containers are removed

---

If you encounter any issues or have questions about the process, check Docker documentation or bioinformatics forums for additional guidance.