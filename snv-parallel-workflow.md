# Scaling SNV Detection: Processing Multiple FASTQ Files with Nextflow and Docker

This guide will help you scale up from processing a single FASTQ file to analyzing thousands of files in parallel using Nextflow and Docker. Nextflow is a workflow manager that's widely used in bioinformatics, relatively easy to learn, and works well with Docker.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Setting Up Nextflow](#setting-up-nextflow)
- [Project Structure](#project-structure)
- [Creating the Workflow](#creating-the-workflow)
- [Running the Pipeline](#running-the-pipeline)
- [Monitoring and Reporting](#monitoring-and-reporting)
- [Advanced Features](#advanced-features)
- [Troubleshooting](#troubleshooting)

## Prerequisites

Before starting, make sure you have:
- Successfully processed one FASTQ file using the previous guide
- Docker installed and working
- Java installed (needed for Nextflow)
- Access to your multiple FASTQ files

## Setting Up Nextflow

### 1. Install Nextflow

```bash
# Download Nextflow
curl -s https://get.nextflow.io | bash

# Make it executable and move it to your PATH
chmod +x nextflow
mv nextflow /usr/local/bin/   # or ~/bin/ if you prefer
```

### 2. Test the Installation

```bash
nextflow -version
```

You should see version information if the installation was successful.

## Project Structure

Create a project structure to organize your workflow:

```bash
mkdir -p ~/snv-practice_batch/{data,workflows,results,configs}
cd ~/snv-practice_batch

# For your input FASTQ files (or create symlinks to save space)
mkdir -p data/fastq_files

# Move reference genome here if you already have it
mkdir -p data/reference
# cp ~/snv-practice/data/hg38.fa data/reference/
```

## Creating the Workflow

Nextflow workflows are defined in files with the `.nf` extension. Let's create a basic workflow for SNV detection.

### 1. Create the Main Workflow File

Create a file named `main.nf` in the workflows directory:

```bash
nano workflows/main.nf
```

Add the following content:

```groovy
#!/usr/bin/env nextflow

// Define parameters with default values
params.fastq_dir = "$baseDir/data/fastq_files"
params.output_dir = "$baseDir/results"
params.reference = "$baseDir/data/reference/hg38.fa"
params.threads = 4
params.help = false

// Show help message
if (params.help) {
    log.info """
    =========================================
    SNV Detection Workflow
    =========================================
    Usage:
    nextflow run workflows/main.nf [options]
    
    Options:
      --fastq_dir     Directory with FASTQ files [default: $params.fastq_dir]
      --output_dir    Output directory [default: $params.output_dir]
      --reference     Reference genome [default: $params.reference]
      --threads       Number of threads [default: $params.threads]
      --help          Show this message
    """
    exit 0
}

// Check if reference genome exists and is indexed
reference_file = file(params.reference)
if (!reference_file.exists()) {
    exit 1, "Reference genome file not found: ${params.reference}"
}

// Create channel for FASTQ files (handles both paired-end and single-end)
// Pattern explanation: Look for files ending with .fastq.gz
//                     Group them by sample name (everything before _R1 or _R2)
fastq_files = Channel
    .fromFilePairs("${params.fastq_dir}/*{_R1,_R2}*.fastq.gz", size: -1)
    .ifEmpty { 
        // If no paired files found, try single-end
        Channel.fromPath("${params.fastq_dir}/*.fastq.gz")
               .map { file -> tuple(file.baseName, file) }
    }

// Index reference genome if needed
process indexReference {
    container 'biocontainers/bwa:v0.7.17_cv1'
    
    input:
    path reference from reference_file
    
    output:
    path "${reference}*" into indexed_reference
    
    script:
    """
    # Check if index exists, if not create it
    if [ ! -f ${reference}.bwt ]; then
        bwa index ${reference}
    fi
    """
}

// Create FAI index if needed
process indexFai {
    container 'quay.io/biocontainers/samtools:1.15--h1170115_1'
    
    input:
    path reference from reference_file
    
    output:
    path "${reference}.fai" into fai_reference
    
    script:
    """
    # Check if index exists, if not create it
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi
    """
}

// Align reads to reference
process alignReads {
    container 'biocontainers/bwa:v0.7.17_cv1'
    publishDir "${params.output_dir}/alignments", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from fastq_files
    path '*' from indexed_reference.collect()
    path reference from reference_file
    
    output:
    tuple val(sample_id), path("${sample_id}.sam") into aligned_sam
    
    script:
    // Determine if single or paired-end from number of files
    if (reads instanceof List && reads.size() > 1) {
        """
        bwa mem -t ${params.threads} ${reference} ${reads[0]} ${reads[1]} > ${sample_id}.sam
        """
    } else {
        // Handle single-end case
        def read = reads instanceof List ? reads[0] : reads
        """
        bwa mem -t ${params.threads} ${reference} ${read} > ${sample_id}.sam
        """
    }
}

// Process SAM files
process processSam {
    container 'quay.io/biocontainers/samtools:1.15--h1170115_1'
    publishDir "${params.output_dir}/bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam_file) from aligned_sam
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai") into sorted_bam
    
    script:
    """
    samtools view -b ${sam_file} > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}

// Call variants
process callVariants {
    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'
    publishDir "${params.output_dir}/vcf", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index) from sorted_bam
    path reference from reference_file
    path fai from fai_reference
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf"), path("${sample_id}.filtered.vcf") into variant_results
    
    script:
    """
    bcftools mpileup -f ${reference} ${bam_file} | bcftools call -mv -o ${sample_id}.vcf
    bcftools filter -i 'QUAL>20 && DP>10' ${sample_id}.vcf > ${sample_id}.filtered.vcf
    bcftools stats ${sample_id}.filtered.vcf > ${sample_id}.stats.txt
    """
}

// Generate report
process generateReport {
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    publishDir "${params.output_dir}/reports", mode: 'copy'
    
    input:
    path '*' from variant_results.map { it[2] }.collect()
    
    output:
    path "multiqc_report.html"
    
    script:
    """
    multiqc . -f
    """
}

workflow.onComplete {
    log.info """
    Pipeline completed successfully!
    Results are in: ${params.output_dir}
    """
}
```

### 2. Create a Configuration File

Create a file named `nextflow.config` in the configs directory:

```bash
nano configs/nextflow.config
```

Add the following content:

```groovy
// Docker configuration
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}

// Process resource allocation
process {
    // Default settings for all processes
    cpus = 1
    memory = 2.GB
    
    // Process-specific settings
    withName: 'alignReads|processSam' {
        cpus = 4
        memory = 8.GB
    }
    
    withName: 'callVariants' {
        memory = 4.GB
    }
}

// Execution reports
report {
    enabled = true
    file = "${params.output_dir}/reports/execution_report.html"
}

timeline {
    enabled = true
    file = "${params.output_dir}/reports/timeline_report.html"
}

trace {
    enabled = true
    file = "${params.output_dir}/reports/trace.txt"
}

// Environment variables
env {
    JAVA_OPTS = '-Xms1g -Xmx4g'
}
```

## Running the Pipeline

### 1. Prepare Your Data

Copy or link your FASTQ files to the input directory:

```bash
# Copy files
cp /path/to/your/fastq/files/*.fastq.gz ~/snv-practice_batch/data/fastq_files/

# Or create symbolic links to save space
ln -s /path/to/your/fastq/files/*.fastq.gz ~/snv-practice_batch/data/fastq_files/

# Copy reference genome if you already have it
cp ~/snv-practice/data/hg38.fa ~/snv-practice_batch/data/reference/
```

### 2. Run the Pipeline

Navigate to your project directory and run the workflow:

```bash
cd ~/snv-practice_batch

# Run the pipeline with default settings
nextflow run workflows/main.nf -c configs/nextflow.config

# To specify custom parameters
nextflow run workflows/main.nf -c configs/nextflow.config --fastq_dir=/path/to/fastq --threads=8
```

### 3. Resume an Interrupted Run

If your pipeline stops for any reason, you can resume from where it left off:

```bash
nextflow run workflows/main.nf -c configs/nextflow.config -resume
```

## Monitoring and Reporting

Nextflow automatically generates several reports:

- **Execution Report**: Shows details about each task execution
- **Timeline Report**: Visualizes the timing of process executions
- **Trace File**: Provides detailed metrics on CPU, memory, and time for each task

To view these reports, open the HTML files in the `results/reports` directory after the pipeline completes.

## Advanced Features

### 1. Running on a Compute Cluster

To run on a cluster, add the appropriate executor configuration to your `nextflow.config`:

```groovy
// For SLURM
executor {
    name = 'slurm'
    queueSize = 100
}

// Process-specific cluster options
process {
    withName: 'alignReads' {
        clusterOptions = '--partition=normal --qos=normal'
    }
}
```

### 2. Using Profiles for Different Environments

Edit your config file to add profiles:

```groovy
profiles {
    standard {
        process.executor = 'local'
    }
    
    cluster {
        process.executor = 'slurm'
        process.queue = 'normal'
    }
    
    cloud {
        // AWS Batch configuration
    }
}
```

Then run with your chosen profile:

```bash
nextflow run workflows/main.nf -c configs/nextflow.config -profile cluster
```

## Troubleshooting

Common issues and solutions:

- **Out of memory errors**: Increase memory allocation in the config file
- **Process fails**: Check the work directory for logs (`work/<process_id>/`)
- **Docker issues**: Ensure Docker is running and you have permission to use it
- **File not found errors**: Verify paths and file permissions

To get detailed debugging information, run Nextflow with the `-with-trace` option:

```bash
nextflow run workflows/main.nf -c configs/nextflow.config -with-trace
```

---

This guide provides a foundation for scaling your SNV detection workflow to process thousands of FASTQ files efficiently using Nextflow and Docker. As you become more comfortable with the workflow, you can extend it with additional steps like variant annotation, quality control, or custom filtering.