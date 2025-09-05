# Maize-Promoter-Motif-Analyzer
A bioinformatics tool that extracts promoter sequences from the maize genome and analyzes the enrichment of transcription factor binding motifs in target genes compared to random sets.

## Overview 

This Python-based tool analyzes promoter sequences of maize genes to identify statistically overrepresented transcription factor binding motifs. It extracts upstream promoter regions from genome FASTA files, processes GFF3 annotations, and compares motif frequencies between user-provided target genes and randomly selected control sets.

## Features

**Genome Processing:** Loads and normalizes chromosome sequences from FASTA files

**GFF3 Parsing:** Extracts gene coordinates and metadata from annotation files

**Promoter Extraction:** Retrieves upstream sequences with customizable length

**Motif Analysis:** Uses regular expressions to identify TF binding sites

**Statistical Comparison:** Compares target genes against multiple random sets

**Comprehensive Output:** Generates tab-separated count files for further analysis

## File Structure

* `promoter_motif_analyzer.py` - Main analysis script 
* `known_motifs.txt` - Database of transcription factor binding motifs 
* `target_genes.txt` - List of gene IDs for analysis 

## Usage

```bash
python promoter_motif_analyzer.py <fasta_directory> <gff3_file> <target_genes_file> <motifs_file> [--promoter_length 500] [--random_sets 5]
```

## Dependencies
Python 3.6+
Biopython
Standard Python libraries: os, logging, argparse, re, random

## Output
The tool generates `motif_counts.txt` containing raw motif counts for both target genes and random control sets, facilitating downstream statistical analysis.
