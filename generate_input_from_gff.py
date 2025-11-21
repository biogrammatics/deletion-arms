#!/usr/bin/env python3
"""
Helper script to generate input FASTA files for deletion_arms_designer.py
from genome FASTA and GFF/GTF annotation files.

This extracts 1500 bp upstream + gene + 1500 bp downstream for specified genes.
"""

import sys
import argparse
from pathlib import Path


def parse_fasta(fasta_file):
    """Parse FASTA file and return dict of {chrom: sequence}"""
    sequences = {}
    current_chrom = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom:
                    sequences[current_chrom] = ''.join(current_seq)
                # Take first word as chromosome name
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_chrom:
            sequences[current_chrom] = ''.join(current_seq)

    return sequences


def parse_gff(gff_file, feature_type='gene'):
    """Parse GFF file and return dict of gene annotations"""
    genes = {}

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom = parts[0]
            ftype = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            if ftype.lower() != feature_type.lower():
                continue

            # Parse gene name from attributes
            gene_name = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('Name=') or attr.startswith('gene='):
                    gene_name = attr.split('=')[1]
                    break
                elif attr.startswith('ID='):
                    # Fallback to ID if no Name
                    gene_name = attr.split('=')[1]

            if gene_name:
                genes[gene_name] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                }

    return genes


def extract_gene_region(genome, gene_info, flank=1500):
    """
    Extract gene region with flanking sequences

    Returns: (upstream 1500bp, gene, downstream 1500bp) as single sequence
    """
    chrom = gene_info['chrom']
    start = gene_info['start']
    end = gene_info['end']
    strand = gene_info['strand']

    if chrom not in genome:
        raise ValueError(f"Chromosome {chrom} not found in genome")

    chrom_seq = genome[chrom]

    # Extract with flanks (1-indexed to 0-indexed)
    region_start = max(0, start - 1 - flank)
    region_end = min(len(chrom_seq), end + flank)

    sequence = chrom_seq[region_start:region_end]

    # If on minus strand, reverse complement
    if strand == '-':
        complement = str.maketrans('ATCG', 'TAGC')
        sequence = sequence.translate(complement)[::-1]

    return sequence


def main():
    parser = argparse.ArgumentParser(
        description='Generate input FASTA for deletion_arms_designer.py from genome and GFF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s genome.fasta annotation.gff --genes GENE1,GENE2,GENE3 -o input.fasta
  %(prog)s genome.fasta annotation.gff --all -o all_genes.fasta
        """
    )

    parser.add_argument('genome', help='Genome FASTA file')
    parser.add_argument('gff', help='GFF/GTF annotation file')
    parser.add_argument('--genes', '-g', help='Comma-separated list of gene names to extract')
    parser.add_argument('--all', '-a', action='store_true',
                       help='Extract all genes from GFF')
    parser.add_argument('--output', '-o', required=True,
                       help='Output FASTA file')
    parser.add_argument('--flank', '-f', type=int, default=1500,
                       help='Flanking sequence length (default: 1500)')
    parser.add_argument('--feature-type', '-t', default='gene',
                       help='GFF feature type to extract (default: gene)')

    args = parser.parse_args()

    if not args.genes and not args.all:
        parser.error("Must specify either --genes or --all")

    print(f"Loading genome from {args.genome}...")
    genome = parse_fasta(args.genome)
    print(f"  Loaded {len(genome)} chromosome(s)")

    print(f"Parsing annotations from {args.gff}...")
    genes = parse_gff(args.gff, args.feature_type)
    print(f"  Found {len(genes)} {args.feature_type} features")

    # Determine which genes to extract
    if args.all:
        genes_to_extract = list(genes.keys())
    else:
        genes_to_extract = [g.strip() for g in args.genes.split(',')]

    print(f"\nExtracting {len(genes_to_extract)} gene region(s)...")

    extracted = 0
    with open(args.output, 'w') as out_f:
        for gene_name in genes_to_extract:
            if gene_name not in genes:
                print(f"  Warning: Gene '{gene_name}' not found in GFF, skipping")
                continue

            try:
                gene_info = genes[gene_name]
                sequence = extract_gene_region(genome, gene_info, args.flank)

                # Calculate sizes
                gene_size = gene_info['end'] - gene_info['start'] + 1
                upstream_size = args.flank
                downstream_size = args.flank
                total_size = len(sequence)

                # Write to FASTA
                out_f.write(f">{gene_name} {gene_info['chrom']}:{gene_info['start']}-{gene_info['end']} "
                          f"strand:{gene_info['strand']} "
                          f"[{upstream_size}bp upstream + {gene_size}bp gene + {downstream_size}bp downstream]\n")

                # Write sequence in 60 bp lines
                for i in range(0, len(sequence), 60):
                    out_f.write(sequence[i:i+60] + '\n')

                extracted += 1
                print(f"  ✓ {gene_name}: {total_size} bp "
                     f"({upstream_size} + {gene_size} + {downstream_size})")

            except Exception as e:
                print(f"  Error extracting {gene_name}: {e}")

    print(f"\nSuccessfully extracted {extracted} gene region(s)")
    print(f"Output written to {args.output}")
    print(f"\nNext step:")
    print(f"  python3 deletion_arms_designer.py {args.output}")


if __name__ == '__main__':
    main()
