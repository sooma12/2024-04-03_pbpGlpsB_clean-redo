# 2024-04-03_pbpGlpsB_clean-redo
Date started: 2024-04-03
Goal: Perform analysis of ∆pbpG and ∆lpsB RNA samples, originally analyzed in Jan 2024

## Sample preparation, library prep, and sequencing
(copied from 2024-01 README)
Followed RNA isolation, reverse transcription, and qPCR protocol to isolate RNA.
Protocol document (followed up to step 22): https://docs.google.com/document/d/1iIDJtlX0jUI4P77D0M6SsOHyaJRPXctn8PhxoLuhKHk/edit
(EG's protocol with MWS's clarifications)
Aliquots of 20 uL of RNA at 100 ng/uL (totaling 2 ug) were made in separate tubes.  These samples were frozen at -80C before shipment to SeqCenter (on ~10 lbs dry ice).

Library preparation and sequencing were performed at SeqCenter.  The following is copy/pasted from their methods pdf (provided with data):
```text
Samples were DNAse treated with Invitrogen DNAse (RNAse free). Library preparation was
performed using Illumina’s Stranded Total RNA Prep Ligation with Ribo-Zero Plus kit and 10bp
unique dual indices (UDI). Sequencing was done on a NovaSeq X Plus, producing paired end
150bp reads. Demultiplexing, quality control, and adapter trimming was performed with bcl-
convert (v4.2.4)1. Sequencing statistics are included in the ‘RNA Sequencing Stats.xlsx’ file.
```

Symlinks to sequencing results (submitted to SeqCenter Jan 11, 2024) were made in `./input/`.

## QC with fastqc

Read-level quality control was performed with FastQC version 0.11.9.

## Trimming

No quality or adapter trimming was performed, as per results of fastQC

## Reference genome prep

A reference genome was prepared using STAR in genomeGenerate mode.
The NZ_CP012004.gff3 file download from NCBI Nucleotide was converted to a gtf file using gffread:
`gffread  REF_GENOMES/17978-mff/NZ_CP012004.gff3 -T -o REF_GENOMES/17978-mff/NZ_CP012004.gtf`



## Alignment

A sample sheet (listing sample names and paired fastq files) was prepared using `prep_sample_sheet.sh`


## Analysis

Functional class overrepresentation was analyzed using the FUNAGE-Pro web server at http://funagepro.molgenrug.nl/.
- Reference genome `Selected RefSeq: Acinetobacter baumannii ATCC 17978-mff ASM107767v1 genomic`

Genes were filtered by padj < 0.1 and either log2FoldChange >1 for upregulated genes or <1 for downregulated genes.
Corresponding lists of DEG locus tags were pasted into the web server text box.  The FUNGAGE-Pro main table for each was downloaded.