# Deletion Arms Designer - Usage Guide

A command-line tool for designing gene knockout constructs using half-site restriction enzyme strategy.

## Quick Start

```bash
# Basic usage - design constructs for sequences in input.fasta
python3 deletion_arms_designer.py input.fasta

# Save output to file
python3 deletion_arms_designer.py input.fasta -o results.txt

# Generate FASTA of construct sequences (uses balanced optimization by default)
python3 generate_construct_fasta.py input.fasta -o constructs.fa
```

## Scripts Overview

| Script | Purpose | Default Optimization |
|--------|---------|---------------------|
| `deletion_arms_designer.py` | Detailed design reports | Cost-only (cw=1.0, lw=0.0) |
| `generate_construct_fasta.py` | FASTA output for synthesis | Balanced (cw=0.5, lw=0.5) |

Both scripts accept the same command-line options.

## Input File Format

The input FASTA file should contain sequences in the format:
```
>Gene_Name
[1500bp upstream][deletion region][1500bp downstream]
```

Each sequence must be at least 3000bp (2 × arm_length). The tool will:
1. Split the sequence into upstream arm, deletion region, and downstream arm
2. Search for restriction enzyme half-sites within the specified range
3. Design constructs that remove the deletion region via homologous recombination

### Generating Input from GFF

Use the helper script to extract sequences from a genome + GFF annotation:
```bash
python3 generate_input_from_gff.py genome.fasta annotations.gff -o input.fasta
```

## Command Line Options

### Required Arguments

| Argument | Description |
|----------|-------------|
| `fasta` | Input FASTA file with gene sequences |

### Optional Arguments

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--output` | `-o` | stdout | Output file for results |
| `--vector` | `-v` | none | Additional vector to check for RE site conflicts (see note below) |
| `--stuffer` | `-s` | none | Custom stuffer insert sequence (FASTA) |
| `--enzymes` | `-e` | enzymes.tsv | Enzyme list TSV file |
| `--arm-length` | `-a` | 1500 | Length of homology arms in bp |
| `--half-site-range` | `-r` | 900 1300 | Min/max distance from deletion for half-sites |
| `--max-designs` | `-m` | 5 (designer) / 4 (fasta) | Maximum designs to output per gene |
| `--cost-weight` | `-cw` | 1.0 (designer) / 0.5 (fasta) | Weight for enzyme cost optimization (0-1) |
| `--length-weight` | `-lw` | 0.0 (designer) / 0.5 (fasta) | Weight for arm length optimization (0-1) |

### Note on Vector Checking

The default `enzymes.tsv` file contains enzymes that have been **pre-screened** to NOT cut the following vectors:
- pR-B
- pR-H
- pR-N
- pR-Z

The `--vector` option is only needed if you are using an **additional vector** not in this pre-screened set. It performs a runtime check to exclude any enzymes that would cut your custom vector.

## Optimization Strategies

The tool ranks designs using a combined score based on enzyme cost and arm length. Use `--cost-weight` and `--length-weight` to control priorities.

### Cost-Only (Default for deletion_arms_designer.py)
Minimize enzyme cost, ignore arm length:
```bash
python3 deletion_arms_designer.py input.fasta --cost-weight 1.0 --length-weight 0.0
```

### Length-Only
Maximize arm length for better recombination efficiency:
```bash
python3 deletion_arms_designer.py input.fasta --cost-weight 0.0 --length-weight 1.0
```

### Balanced (Default for generate_construct_fasta.py)
Equal priority to cost and length:
```bash
python3 deletion_arms_designer.py input.fasta --cost-weight 0.5 --length-weight 0.5
```

### Cost-Preferred
Prioritize cost but consider length:
```bash
python3 deletion_arms_designer.py input.fasta --cost-weight 0.7 --length-weight 0.3
```

## Design Types

### Single-Enzyme Designs (Preferred)
When the same enzyme's left half-site is found in the downstream arm AND right half-site in the upstream arm, no stuffer fragment is needed. The construct is simply:

```
[Downstream arm]--[Upstream arm]
              ^^^^
         Full RE site
```

**Advantages:**
- Simpler construct
- Only one enzyme needed for excision
- No stuffer fragment required

### Two-Enzyme Designs
When different enzymes are used for each site, a stuffer fragment connects them:

```
[Downstream arm]--[Stuffer]--[Upstream arm]
              ^^^^        ^^^^
           RE1 site    RE2 site
```

The stuffer contains:
- Right half of RE1
- AT-rich filler (or custom sequence via `--stuffer`)
- Left half of RE2

## Examples

### Basic Design
```bash
python3 deletion_arms_designer.py genes.fasta
```

### Check Against Additional Vector
If using a vector other than pR-B/H/N/Z, check for conflicts:
```bash
python3 deletion_arms_designer.py genes.fasta --vector pUC19.fasta
```

### Custom Arm Length
Use 2000bp arms with half-sites at 1200-1800bp from deletion:
```bash
python3 deletion_arms_designer.py genes.fasta \
    --arm-length 2000 \
    --half-site-range 1200 1800
```

### Custom Stuffer Insert
Include a selectable marker in the stuffer (for two-enzyme designs):
```bash
python3 deletion_arms_designer.py genes.fasta --stuffer marker.fasta
```

### Generate FASTA Output
Create a FASTA file of construct sequences for synthesis:
```bash
# Uses balanced optimization by default
python3 generate_construct_fasta.py genes.fasta --output constructs.fa

# With custom settings
python3 generate_construct_fasta.py genes.fasta \
    --max-designs 4 \
    --arm-length 2000 \
    --cost-weight 0.7 \
    --length-weight 0.3 \
    --output constructs.fa
```

## Output Format

### deletion_arms_designer.py

Outputs detailed design information:

```
================================================================================
CONSTRUCT DESIGN FOR: Gene_1
================================================================================

Enzyme: PmlI (CAC⇅GTG)
  - Cost tier: inexpensive
  - Price: $331.00 for 10000 units ($0.0331/unit)
  - Left half in downstream arm: CAC at position 1136
  - Right half in upstream arm: GTG at position 308

Total enzyme cost: $0.0331/unit
  ★ Single enzyme design - no stuffer fragment needed!

Effective arm lengths:
  - Downstream: 1139 bp
  - Upstream: 1192 bp
  - Total: 2331 bp

--------------------------------------------------------------------------------
CONSTRUCT STRUCTURE:
--------------------------------------------------------------------------------
Downstream arm (1500 bp) ending with: ...GGAAACACTTTTAGAGAAAT
  → Left half: CAC

  [No stuffer - direct junction]

Upstream arm (1500 bp) starting with: TAATACCGATGCAAGAATAA...
  → Right half: GTG
```

### generate_construct_fasta.py

Outputs FASTA format with descriptive headers:

```
>Gene_1_design1_PmlI_down1139bp_up1192bp
ACTATCACTTGGCTAGTCGGTTCCATTTCTGATAGTAATAGATATAAATATAAATTAGCT
ATTTATTGATTGACACTTAT...
```

Header format: `{gene}_{design#}_{enzyme(s)}_{downstream_length}_{upstream_length}`

## Enzyme File Format

The `enzymes.tsv` file defines available restriction enzymes:

```tsv
Enzyme_Name	Recognition_Sequence	Cost_Tier	NEB_Price	NEB_Units	Notes
EcoRV-HF	GAT⇅ATC	inexpensive	283	50000	High-fidelity version, $0.0057/unit
SmaI	CCC⇅GGG	inexpensive	308	10000	Time-Saver qualified, $0.0308/unit
PmlI	CAC⇅GTG	inexpensive	331	10000	$0.0331/unit
```

- `⇅` marks the cut site (blunt cutters only)
- Cost tiers based on actual $/unit:
  - `inexpensive`: < $0.05/unit
  - `moderate`: $0.05 - $0.10/unit
  - `expensive`: > $0.10/unit
- Enzymes are automatically sorted by cost per unit
- Pre-screened to not cut vectors: pR-B, pR-H, pR-N, pR-Z

## Web Interface

A web-based interface is also available for designing constructs interactively.

### Running Locally

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Start the web server
uvicorn web.app:app --host 0.0.0.0 --port 8000

# 3. Open in browser
# http://localhost:8000
```

### Running with Docker

```bash
# Build the Docker image
docker build -t deletion-arms .

# Run the container
docker run -p 8000:8000 deletion-arms

# Open in browser
# http://localhost:8000
```

### Web Interface Features

- **File upload or paste** - Upload FASTA files or paste sequences directly
- **Adjustable parameters** - Arm length, half-site range, max designs
- **Optimization slider** - Balance between enzyme cost and arm length
- **Enzyme filtering** - Select only enzymes you have available in your freezer
- **Download results** - Export as FASTA or detailed text report
- **Interactive results** - Expandable details for each design

### API Endpoints

The web interface also exposes a REST API:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/enzymes` | GET | List all available enzymes |
| `/api/design` | POST | Design constructs (JSON response) |
| `/api/design/fasta` | POST | Design constructs (FASTA download) |
| `/api/design/report` | POST | Design constructs (text report download) |

## Troubleshooting

### No designs found
- Check that input sequences are long enough (≥ 2 × arm_length)
- Expand the half-site search range: `--half-site-range 800 1400`
- Try a different enzyme set

### All designs use expensive enzymes
- The sequence may lack half-sites for cheaper enzymes
- Consider using `--length-weight 0` to strictly prioritize cost

### Arm lengths too short
- Increase `--length-weight` to prioritize longer arms
- Expand the half-site search range toward the deletion boundary

### Need to use a different vector
- Use `--vector your_vector.fasta` to exclude enzymes that cut your vector
- The default enzyme list already excludes enzymes cutting pR-B, pR-H, pR-N, pR-Z
