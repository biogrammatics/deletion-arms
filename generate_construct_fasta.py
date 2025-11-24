#!/usr/bin/env python3
"""
Generate FASTA file of construct DNA sequences from deletion_arms_designer output.
"""

from deletion_arms_designer import DeletionArmsDesigner, ConstructDesign

def get_construct_sequence(design: ConstructDesign) -> str:
    """
    Generate the full construct DNA sequence.

    Structure:
    - Downstream arm (from start to cut site, including left half of RE1)
    - Stuffer fragment (if different enzymes used)
    - Upstream arm (from cut site to end, including right half of RE2)
    """
    # Downstream portion: from start to end of left half-site
    downstream_seq = design.downstream_arm[:design.downstream_arm_length]

    # Stuffer (empty if same enzyme)
    stuffer_seq = design.generate_stuffer()

    # Upstream portion: from right half-site to end
    upstream_seq = design.upstream_arm[design.upstream_half_site.position:]

    return downstream_seq + stuffer_seq + upstream_seq

def get_construct_name(design: ConstructDesign, design_num: int) -> str:
    """Generate a descriptive FASTA header for the construct."""
    if design.needs_stuffer:
        enzymes = f"{design.re1.name}+{design.re2.name}"
    else:
        enzymes = design.re1.name

    return (f"{design.gene_name}_design{design_num}_{enzymes}_"
            f"down{design.downstream_arm_length}bp_up{design.upstream_arm_length}bp")

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Generate FASTA of construct sequences')
    parser.add_argument('fasta', help='Input FASTA file')
    parser.add_argument('--max-designs', '-m', type=int, default=4)
    parser.add_argument('--cost-weight', '-cw', type=float, default=0.5)
    parser.add_argument('--length-weight', '-lw', type=float, default=0.5)
    parser.add_argument('--output', '-o', default='constructs.fa')

    args = parser.parse_args()

    # Run designer
    designer = DeletionArmsDesigner()
    designs = designer.design_constructs(
        args.fasta,
        max_designs=args.max_designs,
        cost_weight=args.cost_weight,
        length_weight=args.length_weight
    )

    # Group designs by gene
    gene_designs = {}
    for design in designs:
        if design.gene_name not in gene_designs:
            gene_designs[design.gene_name] = []
        gene_designs[design.gene_name].append(design)

    # Output FASTA
    with open(args.output, 'w') as f:
        for gene_name in sorted(gene_designs.keys()):
            for i, design in enumerate(gene_designs[gene_name], 1):
                header = get_construct_name(design, i)
                sequence = get_construct_sequence(design)

                f.write(f">{header}\n")
                # Write sequence in 80-character lines
                for j in range(0, len(sequence), 80):
                    f.write(sequence[j:j+80] + "\n")

    print(f"Wrote {len(designs)} constructs to {args.output}")

if __name__ == '__main__':
    main()
