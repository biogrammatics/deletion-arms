#!/usr/bin/env python3
"""
Gene Knockout Construct Designer for Pichia pastoris

This tool designs double homologous recombination knockout constructs by:
1. Taking sequences with 1500bp upstream + deletion + 1500bp downstream
2. Finding half-sites of blunt-cutting restriction enzymes
3. Designing constructs with homology arms and removable stuffer fragments
"""

import re
import sys
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path


@dataclass
class RestrictionEnzyme:
    """Represents a restriction enzyme with its properties"""
    name: str
    recognition_seq: str  # Full sequence with cut site marked
    left_half: str
    right_half: str
    cost_tier: str = "unknown"
    neb_price: float = 0.0  # Price in USD
    neb_units: int = 0  # Number of units

    def __post_init__(self):
        """Parse the recognition sequence to extract half-sites"""
        # Remove the cut site marker to get full sequence
        self.full_site = self.recognition_seq.replace('⇅', '')

    @property
    def has_degenerate_bases(self) -> bool:
        """Check if enzyme has degenerate bases (RYKMSWN)"""
        degenerate = set('RYKMSWN')
        return any(base in degenerate for base in self.full_site)

    @property
    def cost_per_unit(self) -> float:
        """Calculate cost per unit (lower is better)"""
        if self.neb_units > 0:
            return self.neb_price / self.neb_units
        return float('inf')


@dataclass
class HalfSiteMatch:
    """Represents a found half-site in a sequence"""
    enzyme: RestrictionEnzyme
    position: int
    sequence: str
    is_left_half: bool  # True if left half, False if right half
    in_upstream: bool  # True if in upstream arm, False if in downstream arm


@dataclass
class ConstructDesign:
    """Represents a complete knockout construct design"""
    gene_name: str
    re1: RestrictionEnzyme  # First restriction enzyme
    re2: RestrictionEnzyme  # Second restriction enzyme
    downstream_half_site: HalfSiteMatch  # Left half of RE1 in downstream arm
    upstream_half_site: HalfSiteMatch  # Right half of RE2 in upstream arm
    upstream_arm: str
    downstream_arm: str
    deletion_region: str
    custom_stuffer: str = ""  # Optional custom stuffer insert sequence

    @property
    def needs_stuffer(self) -> bool:
        """Check if stuffer fragment is needed (only when enzymes differ)"""
        return self.re1.name != self.re2.name

    @property
    def total_enzyme_cost(self) -> float:
        """
        Calculate total enzyme cost per unit for this design.
        For single-enzyme designs, returns just that enzyme's cost.
        For two-enzyme designs, returns the sum of both costs.
        """
        if self.needs_stuffer:
            return self.re1.cost_per_unit + self.re2.cost_per_unit
        else:
            return self.re1.cost_per_unit

    @property
    def downstream_arm_length(self) -> int:
        """
        Effective downstream arm length (from start to restriction cut site).
        This is the portion that will participate in homologous recombination.
        """
        return self.downstream_half_site.position + len(self.re1.left_half)

    @property
    def upstream_arm_length(self) -> int:
        """
        Effective upstream arm length (from restriction cut site to end).
        This is the portion that will participate in homologous recombination.
        """
        return len(self.upstream_arm) - self.upstream_half_site.position

    @property
    def total_arm_length(self) -> int:
        """Total effective arm length (sum of both arms)."""
        return self.downstream_arm_length + self.upstream_arm_length

    def optimization_score(self, cost_weight: float = 1.0, length_weight: float = 0.0,
                          max_cost: float = 1.0, max_length: int = 3000) -> float:
        """
        Calculate combined optimization score (lower is better).

        Args:
            cost_weight: Weight for enzyme cost (0-1, default 1.0)
            length_weight: Weight for arm length penalty (0-1, default 0.0)
            max_cost: Maximum expected cost for normalization (default $1.00/unit)
            max_length: Maximum possible arm length for normalization (default 3000bp)

        Returns:
            Combined score where lower is better.
            Score = cost_weight * (normalized_cost) + length_weight * (1 - normalized_length)
        """
        # Normalize cost (0 to 1, where 0 is free and 1 is max_cost)
        normalized_cost = min(self.total_enzyme_cost / max_cost, 1.0)

        # Normalize length (0 to 1, where 1 is max_length)
        # We invert this so that longer arms = lower score (better)
        normalized_length = min(self.total_arm_length / max_length, 1.0)
        length_penalty = 1.0 - normalized_length  # Higher penalty for shorter arms

        return cost_weight * normalized_cost + length_weight * length_penalty

    def generate_stuffer(self, length: int = 50) -> str:
        """
        Generate stuffer fragment sequence
        Format: [Right half RE1] + [custom insert OR filler] + [Left half RE2]

        Returns empty string if both restriction sites use the same enzyme,
        since the sites are compatible and no stuffer is needed.
        """
        # No stuffer needed if both sites use the same enzyme
        if not self.needs_stuffer:
            return ""

        if self.custom_stuffer:
            # Use custom stuffer insert (user-provided sequence)
            return self.re1.right_half + self.custom_stuffer + self.re2.left_half
        else:
            # Generate default AT-rich filler
            filler_length = length - len(self.re1.right_half) - len(self.re2.left_half)
            if filler_length < 0:
                filler_length = 10

            # Use simple AT-rich filler to avoid creating unwanted restriction sites
            filler = "AT" * (filler_length // 2)
            if filler_length % 2:
                filler += "A"

            return self.re1.right_half + filler + self.re2.left_half

    def __str__(self) -> str:
        """Generate human-readable output"""
        stuffer = self.generate_stuffer()

        output = []
        output.append(f"\n{'='*80}")
        output.append(f"CONSTRUCT DESIGN FOR: {self.gene_name}")
        output.append(f"{'='*80}")

        if self.needs_stuffer:
            # Two different enzymes
            output.append(f"\nEnzyme 1: {self.re1.name} ({self.re1.recognition_seq})")
            output.append(f"  - Cost tier: {self.re1.cost_tier}")
            if self.re1.neb_price > 0:
                output.append(f"  - Price: ${self.re1.neb_price:.2f} for {self.re1.neb_units} units (${self.re1.cost_per_unit:.4f}/unit)")
            output.append(f"  - Half-site in downstream arm: {self.downstream_half_site.sequence} at position {self.downstream_half_site.position}")
            output.append(f"\nEnzyme 2: {self.re2.name} ({self.re2.recognition_seq})")
            output.append(f"  - Cost tier: {self.re2.cost_tier}")
            if self.re2.neb_price > 0:
                output.append(f"  - Price: ${self.re2.neb_price:.2f} for {self.re2.neb_units} units (${self.re2.cost_per_unit:.4f}/unit)")
            output.append(f"  - Half-site in upstream arm: {self.upstream_half_site.sequence} at position {self.upstream_half_site.position}")
            output.append(f"\nTotal enzyme cost: ${self.total_enzyme_cost:.4f}/unit (combined)")
        else:
            # Same enzyme for both sites - no stuffer needed
            output.append(f"\nEnzyme: {self.re1.name} ({self.re1.recognition_seq})")
            output.append(f"  - Cost tier: {self.re1.cost_tier}")
            if self.re1.neb_price > 0:
                output.append(f"  - Price: ${self.re1.neb_price:.2f} for {self.re1.neb_units} units (${self.re1.cost_per_unit:.4f}/unit)")
            output.append(f"  - Left half in downstream arm: {self.downstream_half_site.sequence} at position {self.downstream_half_site.position}")
            output.append(f"  - Right half in upstream arm: {self.upstream_half_site.sequence} at position {self.upstream_half_site.position}")
            output.append(f"\nTotal enzyme cost: ${self.total_enzyme_cost:.4f}/unit")
            output.append(f"  ★ Single enzyme design - no stuffer fragment needed!")

        # Arm length summary
        output.append(f"\nEffective arm lengths:")
        output.append(f"  - Downstream: {self.downstream_arm_length} bp")
        output.append(f"  - Upstream: {self.upstream_arm_length} bp")
        output.append(f"  - Total: {self.total_arm_length} bp")

        output.append(f"\n{'-'*80}")
        output.append("CONSTRUCT STRUCTURE:")
        output.append(f"{'-'*80}")
        output.append(f"Downstream arm ({len(self.downstream_arm)} bp) ending with: ...{self.downstream_arm[-20:]}")
        output.append(f"  → Left half: {self.re1.left_half}")

        if self.needs_stuffer:
            output.append(f"\nStuffer fragment ({len(stuffer)} bp): {stuffer}")
            output.append(f"  → Right half RE1: {self.re1.right_half}")
            output.append(f"  → Left half RE2: {self.re2.left_half}")
        else:
            output.append(f"\n  [No stuffer - direct junction]")

        output.append(f"\nUpstream arm ({len(self.upstream_arm)} bp) starting with: {self.upstream_arm[:20]}...")
        output.append(f"  → Right half: {self.re2.right_half}")

        output.append(f"\n{'-'*80}")
        output.append("CLONING FRAGMENT:")
        output.append(f"{'-'*80}")

        # Show the junction regions
        if self.needs_stuffer:
            down_end = self.downstream_arm[-30:] + self.re1.left_half
            stuffer_seq = self.re1.right_half + "..." + self.re2.left_half
            up_start = self.re2.right_half + self.upstream_arm[:30]
            output.append(f"...{down_end}|{stuffer_seq}|{up_start}...")
            output.append(f"         {'▲ RE1 site ▲':^20}   {'▲ RE2 site ▲':^20}")
        else:
            down_end = self.downstream_arm[-30:] + self.re1.left_half
            up_start = self.re2.right_half + self.upstream_arm[:30]
            output.append(f"...{down_end}|{up_start}...")
            output.append(f"              {'▲ Single RE site ▲':^30}")

        output.append(f"\nDeletion size: {len(self.deletion_region)} bp")

        return "\n".join(output)


class DeletionArmsDesigner:
    """Main class for designing deletion constructs"""

    def __init__(self, enzyme_file: str = 'enzymes.tsv'):
        """
        Initialize designer with enzyme list from TSV file

        Args:
            enzyme_file: Path to TSV file with enzyme information
        """
        self.enzymes: List[RestrictionEnzyme] = []
        self.enzyme_file = enzyme_file
        self._parse_enzymes()

    def _parse_enzymes(self):
        """
        Parse enzyme list from TSV file

        Expected format:
        Enzyme_Name <tab> Recognition_Sequence <tab> Cost_Tier <tab> NEB_Price <tab> NEB_Units <tab> Notes
        Lines starting with # are comments
        """
        try:
            with open(self.enzyme_file, 'r') as f:
                for line in f:
                    line = line.strip()

                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue

                    # Skip header line
                    if line.startswith('Enzyme_Name'):
                        continue

                    # Parse TSV
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        name = parts[0].strip()
                        rec_seq = parts[1].strip()
                        cost_tier = parts[2].strip()

                        # Try to parse price and units if available
                        neb_price = 0.0
                        neb_units = 0
                        if len(parts) >= 5:
                            try:
                                neb_price = float(parts[3].strip())
                                neb_units = int(parts[4].strip())
                            except (ValueError, IndexError):
                                pass  # Keep defaults if parsing fails

                        # Split at cut site
                        halves = rec_seq.split('⇅')
                        if len(halves) == 2:
                            enzyme = RestrictionEnzyme(
                                name=name,
                                recognition_seq=rec_seq,
                                left_half=halves[0],
                                right_half=halves[1],
                                cost_tier=cost_tier,
                                neb_price=neb_price,
                                neb_units=neb_units
                            )
                            self.enzymes.append(enzyme)

            # Sort by cost per unit (lower is better), then by tier, then by name
            tier_order = {'inexpensive': 0, 'moderate': 1, 'expensive': 2, 'unknown': 3}
            self.enzymes.sort(key=lambda e: (e.cost_per_unit, tier_order.get(e.cost_tier, 3), e.name))

        except FileNotFoundError:
            raise FileNotFoundError(
                f"Enzyme file '{self.enzyme_file}' not found. "
                f"Please provide a TSV file with enzyme information."
            )

    def parse_fasta(self, fasta_file: str) -> Dict[str, str]:
        """Parse FASTA file and return dict of {name: sequence}"""
        sequences = {}
        current_name = None
        current_seq = []

        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_name:
                        sequences[current_name] = ''.join(current_seq)
                    current_name = line[1:].split()[0]  # Take first word after >
                    current_seq = []
                else:
                    current_seq.append(line.upper())

            # Don't forget the last sequence
            if current_name:
                sequences[current_name] = ''.join(current_seq)

        return sequences

    def split_sequence(self, sequence: str, arm_length: int = 1500) -> Tuple[str, str, str]:
        """
        Split sequence into upstream arm, deletion, downstream arm
        Format: [arm_length bp upstream][deletion][arm_length bp downstream]

        Args:
            sequence: Full sequence including arms and deletion
            arm_length: Length of homology arms in bp (default: 1500)
        """
        min_length = arm_length * 2
        if len(sequence) < min_length:
            raise ValueError(f"Sequence too short ({len(sequence)} bp). Need at least {min_length} bp for {arm_length} bp arms.")

        upstream = sequence[:arm_length]
        downstream = sequence[-arm_length:]
        deletion = sequence[arm_length:-arm_length]

        return upstream, deletion, downstream

    def expand_degenerate(self, pattern: str) -> List[str]:
        """Expand degenerate IUPAC codes to all possible sequences"""
        iupac = {
            'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
            'R': ['A', 'G'], 'Y': ['C', 'T'], 'M': ['A', 'C'],
            'K': ['G', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
            'N': ['A', 'C', 'G', 'T'],
            'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T'],
            'V': ['A', 'C', 'G'], 'D': ['A', 'G', 'T']
        }

        def expand_recursive(seq: str) -> List[str]:
            if not seq:
                return ['']

            first = seq[0]
            rest = seq[1:]

            expanded_rest = expand_recursive(rest)
            result = []

            for base in iupac.get(first, [first]):
                for suffix in expanded_rest:
                    result.append(base + suffix)

            return result

        return expand_recursive(pattern)

    def find_half_sites(self, sequence: str, enzyme: RestrictionEnzyme,
                       search_left: bool, start: int, end: int) -> List[HalfSiteMatch]:
        """
        Find half-sites of an enzyme in a sequence region

        Args:
            sequence: The sequence to search in
            enzyme: The restriction enzyme
            search_left: True to search for left half, False for right half
            start: Start position in sequence (0-indexed)
            end: End position in sequence (exclusive)
        """
        matches = []
        half_site = enzyme.left_half if search_left else enzyme.right_half

        # Generate all possible sequences if degenerate
        possible_seqs = self.expand_degenerate(half_site)

        # Search in the specified region
        search_region = sequence[start:end]

        for possible_seq in possible_seqs:
            # Find all occurrences
            for match in re.finditer(possible_seq, search_region):
                pos = start + match.start()
                matches.append(HalfSiteMatch(
                    enzyme=enzyme,
                    position=pos,
                    sequence=possible_seq,
                    is_left_half=search_left,
                    in_upstream=False  # Will be set by caller
                ))

        return matches

    def deduplicate_matches_by_position(self, matches: List[HalfSiteMatch]) -> List[HalfSiteMatch]:
        """
        Deduplicate half-site matches by position, keeping only the cheapest enzyme.

        When degenerate enzyme sites (like BsaAI: YAC⇅GTR) match at the same position
        as more specific sites (like SnaBI: TAC⇅GTA), this keeps only the cheaper enzyme.

        Args:
            matches: List of HalfSiteMatch objects

        Returns:
            Deduplicated list with only the cheapest enzyme per position
        """
        if not matches:
            return matches

        # Group matches by position
        by_position: Dict[int, List[HalfSiteMatch]] = {}
        for m in matches:
            if m.position not in by_position:
                by_position[m.position] = []
            by_position[m.position].append(m)

        # For each position, keep only the cheapest enzyme
        deduplicated = []
        for pos, pos_matches in by_position.items():
            # Sort by cost per unit (cheapest first)
            pos_matches.sort(key=lambda m: m.enzyme.cost_per_unit)
            deduplicated.append(pos_matches[0])

        return deduplicated

    def deduplicate_designs_by_position(self, designs: List[ConstructDesign]) -> List[ConstructDesign]:
        """
        Deduplicate designs that have the same half-site positions but use different enzymes.

        When degenerate sites like BsaAI overlap with specific sites like SnaBI,
        multiple designs may be generated for the same positions. This keeps only
        the cheapest enzyme combination for each position pair.

        Args:
            designs: List of ConstructDesign objects

        Returns:
            Deduplicated list with only the cheapest enzyme combination per position pair
        """
        if not designs:
            return designs

        # Group designs by their half-site positions
        by_positions: Dict[Tuple[int, int], List[ConstructDesign]] = {}
        for d in designs:
            key = (d.downstream_half_site.position, d.upstream_half_site.position)
            if key not in by_positions:
                by_positions[key] = []
            by_positions[key].append(d)

        # For each position pair, keep only the design with lowest total enzyme cost
        deduplicated = []
        for positions, pos_designs in by_positions.items():
            # Sort by total enzyme cost (cheapest first)
            pos_designs.sort(key=lambda d: d.total_enzyme_cost)
            deduplicated.append(pos_designs[0])

        return deduplicated

    def check_full_site_absent(self, sequence: str, enzyme: RestrictionEnzyme,
                              before_position: int) -> bool:
        """
        Check that full restriction site doesn't appear before given position
        """
        search_region = sequence[:before_position]
        possible_sites = self.expand_degenerate(enzyme.full_site)

        for site in possible_sites:
            if site in search_region:
                return False

        return True

    def check_site_absent_in_construct(self, enzyme: RestrictionEnzyme,
                                       upstream: str, downstream: str) -> bool:
        """
        Check that full restriction site doesn't appear anywhere in the construct arms
        This is critical - if the site appears anywhere, digestion will cut the construct
        into pieces!
        """
        possible_sites = self.expand_degenerate(enzyme.full_site)

        for site in possible_sites:
            if site in upstream or site in downstream:
                return False

        return True

    def design_constructs(self, fasta_file: str, vector_file: Optional[str] = None,
                         stuffer_file: Optional[str] = None,
                         arm_length: int = 1500,
                         half_site_min: int = 900,
                         half_site_max: int = 1300,
                         max_designs: int = 5,
                         cost_weight: float = 1.0,
                         length_weight: float = 0.0,
                         available_enzymes: Optional[List[str]] = None) -> List[ConstructDesign]:
        """
        Main function to design knockout constructs

        Args:
            fasta_file: Path to FASTA file with sequences
            vector_file: Optional path to vector sequence file
            stuffer_file: Optional path to FASTA file with custom stuffer insert sequence
            arm_length: Length of homology arms in bp (default: 1500)
            half_site_min: Minimum distance from deletion for half-sites (default: 900)
            half_site_max: Maximum distance from deletion for half-sites (default: 1300)
            max_designs: Maximum number of designs to return per gene
            cost_weight: Weight for enzyme cost optimization (0-1, default 1.0)
            length_weight: Weight for arm length optimization (0-1, default 0.0)
            available_enzymes: Optional list of enzyme names to use (user's freezer stock)
        """
        sequences = self.parse_fasta(fasta_file)
        all_designs = []

        # Filter enzymes if user specified available ones
        if available_enzymes:
            available_set = set(e.lower() for e in available_enzymes)
            enzymes_to_use = [e for e in self.enzymes if e.name.lower() in available_set]
            if not enzymes_to_use:
                raise ValueError(f"None of the specified enzymes found in database: {available_enzymes}")
            print(f"Using {len(enzymes_to_use)} user-specified enzymes: {[e.name for e in enzymes_to_use]}")
        else:
            enzymes_to_use = self.enzymes

        # Parse vector if provided
        vector_seq = None
        if vector_file:
            vector_seqs = self.parse_fasta(vector_file)
            vector_seq = list(vector_seqs.values())[0] if vector_seqs else None

        # Parse stuffer insert if provided
        stuffer_insert = ""
        if stuffer_file:
            stuffer_seqs = self.parse_fasta(stuffer_file)
            stuffer_insert = list(stuffer_seqs.values())[0] if stuffer_seqs else ""
            print(f"Using custom stuffer insert: {len(stuffer_insert)} bp")

        for gene_name, sequence in sequences.items():
            print(f"\nProcessing {gene_name}...")
            print(f"  Total sequence length: {len(sequence)} bp")

            upstream, deletion, downstream = self.split_sequence(sequence, arm_length)
            print(f"  Upstream: {len(upstream)} bp, Deletion: {len(deletion)} bp, Downstream: {len(downstream)} bp")

            designs_for_gene = []

            # For each enzyme pair
            for i, re1 in enumerate(enzymes_to_use):
                # CRITICAL: Check that RE1 doesn't appear anywhere in the construct
                # If it does, digestion will cut the construct into pieces
                if not self.check_site_absent_in_construct(re1, upstream, downstream):
                    continue

                # Find left half of RE1 in downstream arm (positions half_site_min to half_site_max)
                downstream_matches = self.find_half_sites(
                    downstream, re1, search_left=True, start=half_site_min, end=half_site_max
                )

                if not downstream_matches:
                    continue

                # Check vector if provided
                if vector_seq:
                    possible_sites = self.expand_degenerate(re1.full_site)
                    if any(site in vector_seq for site in possible_sites):
                        continue

                for down_match in downstream_matches:
                    down_match.in_upstream = False

                    # Now find enzyme for RE2 (can be same as RE1 for single-enzyme designs)
                    for re2 in enzymes_to_use:
                        # If different enzyme, check that RE2 doesn't appear anywhere in the construct
                        if re2.name != re1.name:
                            if not self.check_site_absent_in_construct(re2, upstream, downstream):
                                continue

                        # Find right half of RE2 in upstream arm
                        # Convert distance from deletion to position from start of upstream arm
                        upstream_start = arm_length - half_site_max
                        upstream_end = arm_length - half_site_min
                        upstream_matches = self.find_half_sites(
                            upstream, re2, search_left=False, start=upstream_start, end=upstream_end
                        )

                        if not upstream_matches:
                            continue

                        # Check vector if provided
                        if vector_seq:
                            possible_sites = self.expand_degenerate(re2.full_site)
                            if any(site in vector_seq for site in possible_sites):
                                continue

                        for up_match in upstream_matches:
                            up_match.in_upstream = True

                            # We have a valid design!
                            design = ConstructDesign(
                                gene_name=gene_name,
                                re1=re1,
                                re2=re2,
                                downstream_half_site=down_match,
                                upstream_half_site=up_match,
                                upstream_arm=upstream,
                                downstream_arm=downstream,
                                deletion_region=deletion,
                                custom_stuffer=stuffer_insert
                            )
                            designs_for_gene.append(design)

            # Deduplicate designs where multiple enzymes match at the same positions
            # (e.g., BsaAI and SnaBI overlapping due to degenerate bases)
            # Keeps only the cheapest enzyme combination for each position pair
            before_dedup = len(designs_for_gene)
            designs_for_gene = self.deduplicate_designs_by_position(designs_for_gene)
            if before_dedup != len(designs_for_gene):
                print(f"  Deduplicated {before_dedup - len(designs_for_gene)} overlapping designs")

            # Sort designs by optimization score (lower is better)
            # Score combines enzyme cost and arm length based on configured weights
            # max_length is 2 * arm_length (theoretical max for both arms combined)
            max_length = 2 * arm_length
            designs_for_gene.sort(key=lambda d: (
                d.optimization_score(cost_weight, length_weight, max_length=max_length),
                d.re1.name,
                d.re2.name
            ))

            # Take only the best designs up to max_designs
            designs_for_gene = designs_for_gene[:max_designs]

            print(f"  Found {len(designs_for_gene)} valid design(s)")
            all_designs.extend(designs_for_gene)

        return all_designs


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Design gene knockout constructs for Pichia pastoris',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s input.fasta
  %(prog)s input.fasta --vector pUC19.fasta --stuffer marker.fasta
  %(prog)s input.fasta --arm-length 2000 --half-site-range 1000 1500
  %(prog)s input.fasta --max-designs 10 --output results.txt

Optimization weights:
  --cost-weight and --length-weight control how designs are ranked.
  Both weights should be between 0 and 1. Higher weight = more importance.

  Examples:
    --cost-weight 1.0 --length-weight 0.0   # Cost only (default)
    --cost-weight 0.0 --length-weight 1.0   # Length only (maximize arm length)
    --cost-weight 0.5 --length-weight 0.5   # Balance cost and length equally
    --cost-weight 0.7 --length-weight 0.3   # Prefer cost, but consider length
        """
    )

    parser.add_argument('fasta', help='Input FASTA file with sequences')
    parser.add_argument('--vector', '-v', help='Vector sequence file (optional)')
    parser.add_argument('--stuffer', '-s', help='Stuffer insert FASTA file (optional, will be flanked by half-sites)')
    parser.add_argument('--enzymes', '-e', default='enzymes.tsv',
                       help='Enzyme list TSV file (default: enzymes.tsv)')
    parser.add_argument('--arm-length', '-a', type=int, default=1500,
                       help='Length of homology arms in bp (default: 1500)')
    parser.add_argument('--half-site-range', '-r', nargs=2, type=int, metavar=('MIN', 'MAX'),
                       default=[900, 1300],
                       help='Distance range from deletion for half-sites in bp (default: 900 1300)')
    parser.add_argument('--max-designs', '-m', type=int, default=5,
                       help='Maximum designs per gene (default: 5)')
    parser.add_argument('--cost-weight', '-cw', type=float, default=1.0,
                       help='Weight for enzyme cost optimization, 0-1 (default: 1.0)')
    parser.add_argument('--length-weight', '-lw', type=float, default=0.0,
                       help='Weight for arm length optimization, 0-1 (default: 0.0)')
    parser.add_argument('--use-enzymes', '-ue', nargs='+', metavar='ENZYME',
                       help='Only use specified enzymes (e.g., --use-enzymes SnaBI EcoRV-HF)')
    parser.add_argument('--output', '-o', help='Output file (default: stdout)')

    args = parser.parse_args()

    # Create designer
    designer = DeletionArmsDesigner(enzyme_file=args.enzymes)

    print(f"Loaded {len(designer.enzymes)} blunt-cutting restriction enzymes from {args.enzymes}")
    print(f"Cost tiers: {sum(1 for e in designer.enzymes if e.cost_tier == 'inexpensive')} inexpensive, "
          f"{sum(1 for e in designer.enzymes if e.cost_tier == 'moderate')} moderate, "
          f"{sum(1 for e in designer.enzymes if e.cost_tier == 'expensive')} expensive")
    print(f"Configuration: {args.arm_length} bp arms, half-sites at {args.half_site_range[0]}-{args.half_site_range[1]} bp from deletion")
    print(f"Optimization: cost_weight={args.cost_weight}, length_weight={args.length_weight}")

    # Design constructs
    designs = designer.design_constructs(
        args.fasta,
        vector_file=args.vector,
        stuffer_file=args.stuffer,
        arm_length=args.arm_length,
        half_site_min=args.half_site_range[0],
        half_site_max=args.half_site_range[1],
        max_designs=args.max_designs,
        cost_weight=args.cost_weight,
        length_weight=args.length_weight,
        available_enzymes=args.use_enzymes
    )

    # Output results
    output_text = []
    for design in designs:
        output_text.append(str(design))

    result = '\n'.join(output_text)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(result)
        print(f"\nResults written to {args.output}")
    else:
        print(result)

    print(f"\n\nTotal designs generated: {len(designs)}")


if __name__ == '__main__':
    main()
