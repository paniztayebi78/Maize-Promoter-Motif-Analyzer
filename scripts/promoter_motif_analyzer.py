import os
import logging
import argparse
from Bio import SeqIO
import re
import random

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# Function to load genome data from FASTA files
def load_genome(fasta_dir):
    """
    Loads genome data from a directory containing FASTA files and normalizes chromosome keys.
    Args:
        fasta_dir (str): Path to the directory containing genome FASTA files.
    Returns:
        dict: A dictionary where keys are normalized chromosome identifiers and values are Bio.SeqRecord objects.
    """
    genome = {}
    for file in os.listdir(fasta_dir):
        if file.endswith(".fa") or file.endswith(".fasta"):
            filepath = os.path.join(fasta_dir, file)
            for record in SeqIO.parse(filepath, "fasta"):
                # Normalize chromosome names (e.g., remove prefixes and leading zeros)
                normalized_key = record.id.replace("chromosome_", "").replace("chr", "").lstrip("0").upper()
                genome[normalized_key] = record
    logging.info(f"Loaded genome data from {len(genome)} chromosomes.")
    return genome

# Function to parse GFF3 file
def parse_gff3(file_path):
    """
    Parses a GFF3 file to extract gene data.
    Args:
        file_path (str): Path to the GFF3 file.
    Returns:
        dict: A dictionary where keys are gene IDs and values contain gene details.
    """
    genes = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue

            # Extract attributes and normalize gene IDs
            attributes = {
                item.split("=")[0]: item.split("=")[1]
                for item in cols[8].split(";") if "=" in item
            }
            if 'ID' in attributes:
                gene_id = attributes['ID'].split(':')[1].split('.')[0].upper() if ':' in attributes['ID'] else attributes['ID'].split('.')[0].upper()
                seqid = cols[0].replace("chromosome_", "").replace("chr", "").lstrip("0").upper()
                if seqid in [str(i) for i in range(1, 11)]:  # Only consider valid chromosomes
                    genes[gene_id] = {
                        'seqid': seqid,
                        'start': int(cols[3]),
                        'end': int(cols[4]),
                        'strand': cols[6]
                    }
    logging.info(f"Parsed {len(genes)} genes from GFF3 file.")
    return genes

# Function to extract promoter sequences
def extract_promoter(genome, gene_info, length=500):
    """
    Extracts the upstream promoter sequence for a given gene.
    Args:
        genome (dict): A dictionary containing chromosome sequences.
        gene_info (dict): Information about the gene.
        length (int, optional): Length of the upstream promoter region.
    Returns:
        str: The extracted promoter sequence.
    """
    chrom_key = gene_info.get('seqid')
    if chrom_key not in genome:
        logging.warning(f"Chromosome {chrom_key} not found in genome.")
        return None

    # Determine the promoter region based on strand
    if gene_info['strand'] == '+':
        start = max(0, gene_info['start'] - length)
        end = gene_info['start'] - 1
    else:
        start = gene_info['end'] + 1
        end = min(gene_info['end'] + length, len(genome[chrom_key].seq))

    # Extract the sequence from the chromosome and convert to uppercase
    promoter = genome[chrom_key].seq[start:end].upper()

    # Handle missing sequences (Ns) by truncating at the first occurrence of 'N'
    if 'N' in promoter:
        promoter = promoter[:promoter.find('N')]

    # Return reverse complement if the gene is on the negative strand
    return promoter.reverse_complement() if gene_info['strand'] == '-' else promoter

# Function to compile regex patterns for motifs
def compile_motif_patterns(motifs):
    """
    Compiles regex patterns for motifs.
    Args:
        motifs (list): List of motifs.
    Returns:
        dict: Dictionary of compiled regex patterns.
    """
    patterns = {}
    for motif in motifs:
        # Convert motif notation (e.g., [AG]) to regex patterns
        regex = re.sub(r'\[([A-Z]+)\]', lambda x: f"({'|'.join(x.group(1))})", motif.upper())
        patterns[motif] = re.compile(regex)
    return patterns

# Function to read gene IDs
def read_genes(file_path):
    """
    Reads and normalizes gene IDs from the input file.
    Args:
        file_path (str): Path to the file containing gene IDs.
    Returns:
        list: A list of normalized gene IDs.
    """
    if not os.path.exists(file_path):
        logging.error(f"Gene file '{file_path}' not found.")
        return []

    # Open the file and process each line
    with open(file_path, 'r') as file:
        genes = [line.strip().split('.')[0].upper() for line in file.readlines()]  # Normalize gene IDs
    return genes

# Function to read promoter motifs
def read_promoters(file_path):
    """
    Reads promoter motifs from the input file.
    Args:
        file_path (str): Path to the file containing promoter motifs.
    Returns:
        list: A list of motifs.
    """
    if not os.path.exists(file_path):
        logging.error(f"Promoter file '{file_path}' not found.")
        return []

    # Read all motifs from the file, stripping whitespace
    with open(file_path, 'r') as file:
        motifs = [line.strip() for line in file.readlines()]
    return motifs

# Main function
def main():
    parser = argparse.ArgumentParser(description="Analyze maize promoter sequences.")
    parser.add_argument('fasta_dir', help="Directory containing the maize genome FASTA files.")
    parser.add_argument('gff3', help="Path to the maize genome GFF3 file.")
    parser.add_argument('genes', help="Path to the file with target gene IDs.")
    parser.add_argument('promoters', help="Path to the file with known motifs.")
    parser.add_argument('--promoter_length', type=int, default=500, help="Length of upstream promoter region to extract (default: 500 nt).")
    parser.add_argument('--random_sets', type=int, default=5, help="Number of random gene sets to analyze (default: 5).")
    args = parser.parse_args()

    try:
        # Load genome sequences from FASTA files
        genome = load_genome(args.fasta_dir)
        
        # Parse gene information from the GFF3 file
        gene_info = parse_gff3(args.gff3)
        
        # Read the target gene IDs from the input file
        target_genes = read_genes(args.genes)
        
        # Read known transcription factor binding site (TFBS) motifs
        motifs = read_promoters(args.promoters)
        
        # Compile regular expression patterns for each motif
        motif_patterns = compile_motif_patterns(motifs)
    except Exception as e:
        logging.error(e)  # Log any errors encountered during setup
        return

    # Identify unmatched genes (those not found in the GFF3 file)
    unmatched_genes = [gene for gene in target_genes if gene not in gene_info]
    if unmatched_genes:
        logging.warning(f"{len(unmatched_genes)} out of {len(target_genes)} target genes did not match the GFF3 data.")
        logging.warning(f"Examples of unmatched genes: {unmatched_genes[:10]}")

    # Initialize a dictionary to store motif counts for target genes
    target_counts = {}
    for gene in target_genes:
        if gene in gene_info:
            # Extract the promoter sequence for the current gene
            promoter = extract_promoter(genome, gene_info[gene], args.promoter_length)
            if promoter is None:
                continue

            # Count occurrences of each motif in the promoter sequence
            counts = {motif: len(pattern.findall(str(promoter))) for motif, pattern in motif_patterns.items()}
            for motif, count in counts.items():
                target_counts[motif] = target_counts.get(motif, 0) + count

    # Perform analysis on random sets of genes
    random_counts = []
    all_genes = list(gene_info.keys())
    for _ in range(args.random_sets):
        # Randomly select a set of genes equal in size to the target genes
        sample_genes = random.sample(all_genes, len(target_genes))
        sample_counts = {}
        for gene in sample_genes:
            if gene in gene_info:
                # Extract the promoter sequence for the random gene
                promoter = extract_promoter(genome, gene_info[gene], args.promoter_length)
                if promoter:
                    # Count occurrences of motifs in the random promoter sequence
                    counts = {motif: len(pattern.findall(str(promoter))) for motif, pattern in motif_patterns.items()}
                    for motif, count in counts.items():
                        sample_counts[motif] = sample_counts.get(motif, 0) + count
        random_counts.append(sample_counts)

    # Write the results to an output file
    with open("motif_counts.txt", 'w') as output_file:
        # Define header row for the output file
        headers = ["Motif", "Selected_Genes"] + [f"Random_Set_{i+1}" for i in range(args.random_sets)]
        output_file.write("\t".join(headers) + "\n")
        
        # Write motif counts for each motif
        for motif in motifs:
            row = [motif, target_counts.get(motif, 0)] + [counts.get(motif, 0) for counts in random_counts]
            output_file.write("\t".join(map(str, row)) + "\n")

    logging.info("Analysis completed successfully. Results written to 'motif_counts.txt'.")

if __name__ == "__main__":
    main()
