.. CNV-SNP Pipeline README

=============================
CNV and SNP Mosquito Pipeline
=============================

Version: 0.1

This bioinformatics pipeline identifies genetic changes (i.e., copy number variations and single nucleotide polymorphisms) in mosquito populations.

Requirements
============

Environment requirements and setup, to be added.

Steps to Run the Pipeline
=========================

1. Clone the pipeline source code. These instructions assume that you are cloning the source code into your Linux home directory:

    ::

        cd ~
        git clone https://github.com/lkothera/cnv-snp-pipeline

2. Create the "samples" directory:

    ::

        cd ~/cnv-snp-pipeline
        mkdir samples

3. Create a directory for the specific group of samples with which you will be working. In this example, the "culex" folder will be used to store samples for Culex pipiens. However, any name can be used.

    ::

        cd ~/cnv-snp-pipeline/samples
        mkdir culex
        cd culex
        mkdir fastq_orig

    Copy all fastq files for the analysis to the "fastq_orig" folder. 

4. Copy a comparison text file to the samples directory. This file describes which phenotypes are to be compared in the analysis. 

5. Create a directory for the reference sequence to which the samples will be aligned. 

    ::

        cd ~/cnv-snp-pipeline
        mkdir reference
        cd reference
        mkdir vectorbase_files

6. Copy VectorBase reference files into the "vectorbase_files" folder. For culex, these can be the following files:

    a. GTF file - e.g., culex-quinquefasciatus-johannesburgbasefeaturescpipj22.gtf
    b. FASTA file - e.g., culex-quinquefasciatus-johannesburgscaffoldscpipj2.fa

7. Copy the list of gene names into a file called "gene_names.txt" within the "reference" folder. This file includes the target list of genes for CNV analysis. 
    
8. Modify the pipeline configuration file:

    ::

        cd ~/cnv-snp-pipeline/scripts
        vi pipeline.sh

Substitute the brackets with appropriate values:

    ::

        BASE=[your base pipeline directory, i.e., ~/cnv-snp-pipeline]
        SAMPLES=[sample folder name, i.e., culex]
        REFERENCE=[reference folder name, i.e., culex]
        VBGTF=[the name of the VectorBase GTF file]
        VBFASTA=[the name of the VectorBase FASTA file]
        RUN_NAME=[custom name of the analysis run, used for organizing output]
        COMPARE=[name of the comparison text file in the sample directory]

9. Build the reference:

    ::

        cd ~/cnv-snp-pipeline
        make reference

10. Process the samples:

    ::

        cd ~/cnv-snp-pipeline
        make samples

11. Call CNVs Using Conifier:

    ::

        cd ~/cnv-snp-pipeline
        make conifer

12. Select Conifer Features:

    ::

        cd ~/cnv-snp-pipeline
        make features

13. Call SNPs Using GATK:

    ::

        cd ~/cnv-snp-pipeline
        make gatk

14. Select SNP Features:

    ::

        cd ~/cnv-snp-pipeline
        make gatk_features

15. Detect Statistically Significant SNPs with LampLink:

    ::

        cd ~/cnv-snp-pipeline
        make lamplink


