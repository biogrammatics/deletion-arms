# Gene Knockout Construct Designer for *Pichia pastoris*

A Python tool for designing double homologous recombination gene knockout constructs using blunt-cutting restriction enzymes with half-site strategy.

## Overview

This tool automates the design of knockout constructs by:
1. Analyzing DNA sequences with homology arms flanking a target deletion region
2. Finding optimal "half-sites" of blunt-cutting restriction enzymes
3. Designing constructs where half-sites recreate full restriction sites around a removable stuffer fragment
4. Prioritizing designs by enzyme cost-efficiency

## Strategy

### Construct Architecture

```
[Downstream arm—LeftHalf_RE1] [RightHalf_RE1—Stuffer—LeftHalf_RE2] [RightHalf_RE2—Upstream arm]
```

When the half-sites are joined during cloning:

```
[Downstream arm]—[Full RE1 site]—[Stuffer]—[Full RE2 site]—[Upstream arm]
```

**Benefits:**
- Stuffer (selection marker, reporter) can be removed by digesting with RE1 and RE2
- Clean deletion with precise homology arms
- Flexible stuffer replacement
- Cost-optimized enzyme selection

## Features

- ✅ **Configurable homology arm length** (default: 1500 bp)
- ✅ **Configurable half-site placement range** (default: 900-1300 bp from deletion)
- ✅ **Custom stuffer insert support** (selection markers, reporters, etc.)
- ✅ **Vector compatibility checking** (excludes enzymes that cut your cloning vector)
- ✅ **Cost-based enzyme prioritization** (uses real NEB pricing data)
- ✅ **Multiple design candidates** per gene
- ✅ **Batch processing** of multiple genes

## Installation

### Requirements
- Python 3.6+
- No external dependencies (uses standard library only)

### Setup

```bash
# Clone the repository
git clone <repository-url>
cd deletion-arms-tool

# Make scripts executable
chmod +x deletion_arms_designer.py
chmod +x generate_input_from_gff.py

# Verify installation
./deletion_arms_designer.py --help
```

## Quick Start

### 1. Prepare Your Input Sequences

Input FASTA format:
```
[arm_length bp upstream][Target deletion region][arm_length bp downstream]
```

**Example:**
```fasta
>MY_GENE_1
[1500 bp upstream][gene to delete][1500 bp downstream]
```

### 2. Basic Usage

```bash
# Design knockouts with default settings
./deletion_arms_designer.py input.fasta

# Save output to file
./deletion_arms_designer.py input.fasta --output results.txt

# Process with vector compatibility check
./deletion_arms_designer.py input.fasta --vector pRS426.fasta
```

### 3. Advanced Usage

```bash
# Custom everything
./deletion_arms_designer.py genes.fasta \
  --arm-length 2000 \
  --half-site-range 1000 1500 \
  --stuffer KanMX.fasta \
  --vector pRS426.fasta \
  --max-designs 10 \
  --output results.txt
```

## Configuration Files

### enzymes.tsv

Tab-separated file containing available restriction enzymes:

```tsv
Enzyme_Name	Recognition_Sequence	Cost_Tier	NEB_Price	NEB_Units	Notes
EcoRV-HF	GAT⇅ATC	inexpensive	283	50000	High-fidelity version
SmaI	CCC⇅GGG	inexpensive	308	10000	Time-Saver qualified
```

**Columns:**
- `Enzyme_Name`: NEB enzyme name
- `Recognition_Sequence`: Cut site marked with ⇅
- `Cost_Tier`: inexpensive/moderate/expensive/unknown
- `NEB_Price`: Price in USD
- `NEB_Units`: Units in catalog size
- `Notes`: Optional notes

**To update:** Simply edit the TSV file in any text editor.

### vectors.fa

FASTA file with your cloning vector sequences. The tool will exclude any enzymes that cut these vectors.

```fasta
>pR-B
GCTAGCATGCTCGAGCTTACTTAAGATGCACCGGTGC...

>pR-H
GCTAGCATGCTCGAGCTTACTTAAGATGCACCGGTGC...
```

## Command-Line Options

```
positional arguments:
  fasta                 Input FASTA file with sequences

optional arguments:
  -h, --help            Show help message
  --vector, -v          Vector sequence file (optional)
  --stuffer, -s         Stuffer insert FASTA file (optional)
  --enzymes, -e         Enzyme list TSV file (default: enzymes.tsv)
  --arm-length, -a      Length of homology arms in bp (default: 1500)
  --half-site-range, -r MIN MAX
                        Distance range from deletion for half-sites (default: 900 1300)
  --max-designs, -m     Maximum designs per gene (default: 5)
  --output, -o          Output file (default: stdout)
```

## Examples

### Example 1: Default Settings
```bash
./deletion_arms_designer.py genes.fasta
```

### Example 2: Custom Marker
```bash
./deletion_arms_designer.py genes.fasta \
  --stuffer KanMX_marker.fasta \
  --output designs.txt
```

### Example 3: Shorter Arms
```bash
./deletion_arms_designer.py genes.fasta \
  --arm-length 1000 \
  --half-site-range 600 800
```

### Example 4: From Genome Annotation
```bash
# Generate input from GFF
./generate_input_from_gff.py genome.fasta annotation.gff \
  --genes AOX1,GUT1,PEX5 \
  -o input.fasta

# Design constructs
./deletion_arms_designer.py input.fasta --output results.txt
```

## Output Format

For each valid design, the tool provides:

```
================================================================================
CONSTRUCT DESIGN FOR: Gene_1
================================================================================

Enzyme 1: EcoRV-HF (GAT⇅ATC)
  - Cost tier: inexpensive
  - Price: $283.00 for 50000 units ($0.0057/unit)
  - Half-site in downstream arm: GAT at position 1050

Enzyme 2: SmaI (CCC⇅GGG)
  - Cost tier: inexpensive
  - Price: $308.00 for 10000 units ($0.0308/unit)
  - Half-site in upstream arm: GGG at position 450

--------------------------------------------------------------------------------
CONSTRUCT STRUCTURE:
--------------------------------------------------------------------------------
Downstream arm (1500 bp) ending with: ...AGCTAGCTAGAT
  → Left half RE1: GAT

Stuffer fragment (50 bp): ATC[filler or custom insert]CCC
  → Right half RE1: ATC
  → Left half RE2: CCC

Upstream arm (1500 bp) starting with: GGGATCGATCGAT...
  → Right half RE2: GGG

--------------------------------------------------------------------------------
CLONING FRAGMENT:
--------------------------------------------------------------------------------
...CTAGCTAGCTAGCTAGCTAGCTAGCTAGAT|ATC...CCC|GGGATCGATCGATCGATCGATCGATCGATCGAT...
             ▲ RE1 site ▲          ▲ RE2 site ▲

Deletion size: 500 bp
```

## Helper Tools

### generate_input_from_gff.py

Generates input FASTA files from genome and GFF annotation:

```bash
./generate_input_from_gff.py genome.fasta annotation.gff \
  --genes GENE1,GENE2,GENE3 \
  -o input.fasta
```

**Features:**
- Extracts genes with flanking regions
- Handles both strands (automatically reverse complements minus strand)
- Customizable flank length

## Design Constraints

The tool ensures valid designs by checking:

1. ✅ **Full RE sites absent from arms** - Prevents unwanted cutting
2. ✅ **Full RE sites absent from vector** - Enables cloning
3. ✅ **Half-sites in correct positions** - Ensures proper construct assembly
4. ✅ **Cost optimization** - Prioritizes economical enzymes

## Cost Optimization

Enzymes are automatically prioritized by **cost per unit**:

```
Most Economical:
#1  EcoRV-HF   $0.0057/unit  ⭐ BEST VALUE (50,000 units / $283)
#2  SmaI       $0.0308/unit
#3  PmlI       $0.0331/unit

Least Economical:
#14 SrfI       $0.1664/unit  (29× more expensive)
#15 PsiI-v2    $0.2775/unit  (49× more expensive)
#16 AfeI       $0.3410/unit  (60× more expensive)
```

## Validation

All designs are validated to ensure:
- No restriction sites in construct arms (prevents unwanted cutting)
- No restriction sites in cloning vector (enables successful cloning)
- Half-sites positioned correctly for efficient recombination
- Both enzymes are cost-effective and available

## Project Structure

```
deletion-arms-tool/
├── deletion_arms_designer.py    # Main design tool
├── generate_input_from_gff.py   # Helper: generate input from genome
├── enzymes.tsv                  # Enzyme database (editable)
├── vectors.fa                   # Cloning vectors (4 vectors included)
├── test_example.fasta           # Example input
├── test_stuffer.fasta           # Example stuffer (KanMX)
├── README.md                    # This file
├── FEATURE_SUMMARY.md           # Detailed feature documentation
└── .gitignore                   # Git ignore rules
```

## Troubleshooting

### "Sequence too short" error
- Ensure sequences are at least `2 × arm_length` bp
- Default requires minimum 3000 bp (1500 + 0 + 1500)

### "No valid designs found"
- Try increasing `--max-designs`
- Check if sequences contain too many restriction sites
- Verify vector file doesn't conflict with available enzymes
- Try adjusting `--half-site-range`

### Enzyme pricing outdated
- Edit `enzymes.tsv` with current NEB prices
- Tool automatically recalculates cost per unit

## Contributing

To add new enzymes:
1. Edit `enzymes.tsv`
2. Add line with enzyme details (tab-separated)
3. Tool automatically loads on next run

To add new vectors:
1. Add sequences to `vectors.fa` in FASTA format
2. Rerun enzyme validation if needed

## Citation

If you use this tool in your research, please cite:
- This tool (GitHub repository)
- New England Biolabs (NEB) for restriction enzyme data

## License

MIT License - Free to use, modify, and distribute

## Contact

For issues, questions, or suggestions, please open an issue on GitHub.

## Acknowledgments

- Developed for *Pichia pastoris* gene knockout experiments
- Enzyme data from New England Biolabs (NEB)
- Built with Python standard library only (no dependencies!)
