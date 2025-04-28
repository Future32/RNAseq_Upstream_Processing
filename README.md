#RNAseq_Upstream_Processing
This project is based on my exploration practice after studying the publicly available course materials and videos of "Intro to Bioinfo & Comp Bio" taught by Professor Shirley Liu from Harvard University.  
During the learning process, I explored RNA-seq upstream processing workflows by applying the software tools, coding practices, and analytical methods introduced in the course.  
The current workflow covers essential steps including quality control, read trimming, alignment, quantification, and quality assessment.  
As I continue learning, more analysis tools and methods will be gradually incorporated to meet diverse research needs.

## 📌 Current Workflow Scripts

### `RNAseq_batch_STAR_RSEM.sh`

- **Workflow Functions**:
  - Batch processing of multiple paired-end Fastq samples
  - Perform initial quality control and trimming using **FastQC** and **Fastp**
  - Conduct post-trimming quality assessment with **FastQC**
  - Perform two-pass alignment to the reference genome using **STAR**
  - Mark duplicate reads using **Picard MarkDuplicates**
  - Assess alignment quality using **RSeQC** and **Qualimap**
  - Quantify gene and transcript expression levels using **RSEM** (TPM/FPKM/count outputs)
  - Summarize quality control results across samples and generate a **MultiQC** report
  - Extract and summarize Mapping Rate and Duplication Rate into dedicated files
  - Automatically generate a comprehensive batch analysis summary and runtime statistics

