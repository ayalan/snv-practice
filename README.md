# SNV Detection from FASTQ Files Using Docker

This beginner's guide will walk you through detecting Single Nucleotide Variants (SNVs) from FASTQ sequencing files using Docker containers. This approach keeps your system clean by installing all required bioinformatics tools inside containers rather than directly on your computer.

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
mkdir -p ~/snv-practice/data
cd ~/snv-practice
```

### 3. Copy FASTQ Files

```bash
# Copy your FASTQ files to the data directory
cp path/to/your/fastq/files/*.fastq.gz ~/snv-practice/data/
```

### 4. Pull Required Docker Images

```bash
# Download the necessary bioinformatics tools as Docker images
docker pull biocontainers/bwa:v0.7.17_cv1
docker pull quay.io/biocontainers/samtools:1.15--h1170115_1
docker pull quay.io/biocontainers/bcftools:1.15--h0ea216a_2
```

**About these images:**

- **biocontainers/bwa**: BWA (Burrows-Wheeler Aligner) is the standard tool for aligning DNA sequences to a reference genome. It's efficient for short reads typical in next-generation sequencing.

- **quay.io/biocontainers/samtools**: SAMtools is essential for manipulating SAM/BAM files (sequence alignment formats). It handles conversion, sorting, indexing, and basic statistics.

- **quay.io/biocontainers/bcftools**: BCFtools works with VCF files (Variant Call Format) and is the standard for calling and filtering genetic variants.

**Finding alternative images:**

1. **BioContainers Registry**: [biocontainers.pro](https://biocontainers.pro) provides curated bioinformatics tools as Docker images

2. **Docker Hub**: Search at [hub.docker.com](https://hub.docker.com) for alternatives like:
   - `broadinstitute/gatk` - A popular alternative for variant discovery
   - `staphb/freebayes` - Another variant caller with different algorithms

3. **Version selection**: You can choose different versions by changing the tag (e.g., `biocontainers/bwa:v0.7.15_cv1`)

4. **Integrated solutions**: Some images combine multiple tools, like `nfcore/sarek` which includes a complete variant calling pipeline

When selecting alternative images, check documentation for compatibility and look for those with regular updates and community support.

## Getting the Reference Genome

### 1. Download Human Reference Genome

```bash
# Download and decompress the reference genome directly
cd ~/snv-practice/data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

### 2. Index the Reference Genome

```bash
# Create BWA index
docker run -v $(pwd)/data:/data -w /data --rm biocontainers/bwa:v0.7.17_cv1 \
  bwa index /data/hg38.fa

# Create samtools index
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/samtools:1.15--h1170115_1 \
  samtools faidx /data/hg38.fa
```

**What indexing does:**

- **BWA indexing**: Creates specialized data structures (suffix arrays and FM-index) that allow BWA to rapidly find where a short read sequence might match in the reference genome. Without this index, BWA would need to search through the entire reference genome for each read, which would be extremely slow. Think of it like creating an index in a book - it helps you quickly find specific content without reading the entire book.

- **Samtools indexing**: Creates a companion file (`.fai`) that contains the positions of each chromosome/contig in the FASTA file along with their lengths. This allows tools to quickly jump to a specific region in the genome without scanning the whole file. It's like adding bookmarks and page numbers to the reference genome.

## Aligning Reads

### 1. Align FASTQ Files to Reference

**Understanding Paired-End vs. Single-End Sequencing:**

- **Single-end sequencing**: The sequencer reads the DNA fragment from only one end. This produces one FASTQ file per sample.

- **Paired-end sequencing**: The sequencer reads the DNA fragment from both ends. This produces two FASTQ files per sample (typically named with R1 and R2 in the filename).

**How to tell what type of data you have:**

1. **Check the number of files per sample**: 
   - One file per sample → likely single-end
   - Two files per sample (often with R1/R2 or 1/2 in names) → paired-end

2. **Look at the filenames**:
   - Paired-end often follows patterns like:
     - `sample_R1.fastq.gz` and `sample_R2.fastq.gz`
     - `sample_1.fastq.gz` and `sample_2.fastq.gz`
     - `sample_forward.fastq.gz` and `sample_reverse.fastq.gz`

3. **Check file contents**:
   ```bash
   # Look at read IDs in the FASTQ file
   zcat your_file.fastq.gz | head -n 4
   ```
   - In paired-end data, corresponding reads in R1 and R2 files have identical IDs (except possibly for a /1 or /2 suffix)

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

**What these commands do:**

- **SAM to BAM conversion** (`samtools view -b`): 
  Converts the text-based SAM format to the binary BAM format. This is like converting a large text document to a compressed format. BAM files are typically 3-5x smaller than SAM files, much faster to process, and allow for random access when indexed. The conversion parses each alignment record and stores it in a specialized binary representation optimized for genomic data.

- **Sorting** (`samtools sort`):
  Rearranges the alignment records so they're in order by genomic position (chromosome and coordinate). This is required for most downstream processes and improves access efficiency. Think of it like sorting a filing cabinet by page number instead of having pages in random order.

- **Indexing** (`samtools index`):
  Creates a companion file (`.bai`) that works like a table of contents, allowing tools to quickly jump to specific genomic regions without scanning the entire file. This is essential for efficient visualization and analysis of specific regions.

### 2. Generate Alignment Statistics

```bash
# Check alignment quality
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/samtools:1.15--h1170115_1 \
  samtools flagstat /data/aligned.sorted.bam > /data/alignment_stats.txt
```

This produces a summary of your alignment quality, which is crucial for quality control and troubleshooting.

**How to use alignment statistics:**

- **Quality control**: Before proceeding to variant calling, check if the alignment meets expected quality thresholds. For human whole genome sequencing, alignment rates above 95% are typically expected.

- **Troubleshooting**: High rates of unmapped reads might indicate sample contamination, wrong reference genome, or poor sequencing quality. 

- **Decision making**: Based on these statistics, you might decide to adjust alignment parameters, try different tools, or even re-sequence samples that show poor alignment metrics.

**Key metrics to look for in the flagstat output:**

- **Mapped reads percentage**: Primary indicator of alignment success. Low values (<80% for DNA-seq) suggest problems.

- **Properly paired percentage**: For paired-end data, this shows read pairs aligning with correct orientation and distance. Values below 85% might indicate structural variations or quality issues.

- **Singletons**: Reads whose mate pair failed to align. High numbers can indicate poor quality in one of the read files.

- **Duplicates**: PCR or optical duplicates. Excessive duplication (>20%) suggests PCR bias or low library complexity.

The alignment statistics serve as a critical checkpoint before investing time in downstream analysis like variant calling.

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

## Next Steps: Analyzing and Visualizing Your Variants

Once you have your VCF file containing variants, there are several ways to explore and interpret this data:

### Visualization Tools

1. **Genome Browsers**
   ```bash
   # Run IGV (Integrative Genomics Viewer) in Docker
   docker pull broadinstitute/igv:latest
   docker run -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $(pwd)/data:/data --rm broadinstitute/igv
   ```
   
   Other popular browsers:
   - UCSC Genome Browser (web-based)
   - JBrowse (can be run locally or on a server)

2. **VCF-specific Viewers**
   - VCF.iobio (web-based tool for VCF exploration)
   - Varsome (online variant interpretation)

### Functional Annotation

```bash
# Annotate variants with SnpEff
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/snpeff:5.0--hdfd78af_1 bash -c \
  "snpEff download -v GRCh38.99 && \
   snpEff ann GRCh38.99 /data/filtered_variants.vcf > /data/annotated_variants.vcf"

# Alternative: Annotate with ANNOVAR
docker pull clinicalgenomics/annovar
docker run -v $(pwd)/data:/data --rm clinicalgenomics/annovar bash -c \
  "cd /annovar && ./annotate_variation.pl -buildver hg38 -downdb refGene humandb/ && \
   ./table_annovar.pl /data/filtered_variants.vcf humandb/ -buildver hg38 -out /data/annotated \
   -remove -protocol refGene -operation g -nastring . -vcfinput"
```

### Statistical Analysis and Reporting

1. **Using R for variant analysis**
   ```bash
   docker pull bioconductor/bioconductor_docker:latest
   docker run -v $(pwd)/data:/data -w /data --rm bioconductor/bioconductor_docker:latest R -e \
     "library(VariantAnnotation); \
      vcf <- readVcf('/data/filtered_variants.vcf', 'hg38'); \
      summary(vcf)"
   ```

2. **Basic VCF statistics**
   ```bash
   docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/bcftools:1.15--h0ea216a_2 \
     bcftools stats /data/filtered_variants.vcf > /data/variant_summary.txt
   ```

### Common Analysis Workflows

1. **For clinical variants**:
   - Filter for rare variants (population frequency < 1%)
   - Focus on protein-altering changes (missense, nonsense, frameshift)
   - Prioritize by disease relevance or gene lists

2. **For research studies**:
   - Compare variant frequencies between groups
   - Identify mutation patterns or signatures
   - Connect variants to phenotypic data

3. **For technical validation**:
   - Examine variant quality and read support in a genome browser
   - Compare results with different variant callers
   - Consider laboratory validation for key findings

### Advanced Visualization

```bash
# Create a VCF summary report with MultiQC
docker pull quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0
docker run -v $(pwd)/data:/data -w /data --rm quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0 \
  multiqc .
```
