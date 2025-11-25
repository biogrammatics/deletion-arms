"""
API routes for the Deletion Arms Designer
"""

import sys
import os
import tempfile
from typing import Optional
from io import StringIO

from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from fastapi.responses import StreamingResponse

# Add parent directory to path to import deletion_arms_designer
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from deletion_arms_designer import DeletionArmsDesigner, ConstructDesign
from web.models import (
    DesignRequest, DesignResponse, GeneDesigns, DesignResult,
    EnzymeInfo, EnzymeListResponse
)

router = APIRouter()

# Initialize designer with default enzyme file
def get_designer():
    """Get a DeletionArmsDesigner instance"""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    enzyme_file = os.path.join(base_dir, 'enzymes.tsv')
    return DeletionArmsDesigner(enzyme_file=enzyme_file)


def design_to_result(design: ConstructDesign, design_num: int) -> DesignResult:
    """Convert a ConstructDesign to a DesignResult"""
    # Get construct sequence
    downstream_seq = design.downstream_arm[:design.downstream_arm_length]
    stuffer_seq = design.generate_stuffer()
    upstream_seq = design.upstream_arm[design.upstream_half_site.position:]
    construct_sequence = downstream_seq + stuffer_seq + upstream_seq

    return DesignResult(
        design_num=design_num,
        enzyme1_name=design.re1.name,
        enzyme1_recognition=design.re1.recognition_seq,
        enzyme1_cost_per_unit=design.re1.cost_per_unit,
        enzyme2_name=design.re2.name,
        enzyme2_recognition=design.re2.recognition_seq,
        enzyme2_cost_per_unit=design.re2.cost_per_unit,
        single_enzyme=not design.needs_stuffer,
        needs_stuffer=design.needs_stuffer,
        downstream_arm_length=design.downstream_arm_length,
        upstream_arm_length=design.upstream_arm_length,
        total_arm_length=design.total_arm_length,
        total_enzyme_cost=design.total_enzyme_cost,
        downstream_half_site_position=design.downstream_half_site.position,
        upstream_half_site_position=design.upstream_half_site.position,
        construct_sequence=construct_sequence,
        stuffer_sequence=stuffer_seq
    )


@router.get("/api/enzymes", response_model=EnzymeListResponse)
async def list_enzymes():
    """List all available restriction enzymes"""
    designer = get_designer()
    enzymes = [
        EnzymeInfo(
            name=e.name,
            recognition_seq=e.recognition_seq,
            full_site=e.full_site,
            left_half=e.left_half,
            right_half=e.right_half,
            cost_tier=e.cost_tier,
            cost_per_unit=e.cost_per_unit,
            neb_price=e.neb_price,
            neb_units=e.neb_units
        )
        for e in designer.enzymes
    ]
    return EnzymeListResponse(enzymes=enzymes, total=len(enzymes))


@router.post("/api/design", response_model=DesignResponse)
async def design_constructs(
    file: Optional[UploadFile] = File(None),
    sequences: Optional[str] = Form(None),
    arm_length: int = Form(1500),
    half_site_min: int = Form(900),
    half_site_max: int = Form(1300),
    max_designs: int = Form(5),
    cost_weight: float = Form(0.5),
    length_weight: float = Form(0.5),
    available_enzymes: Optional[str] = Form(None)
):
    """
    Design knockout constructs from FASTA sequences.

    Either upload a file or provide sequences as text.
    """
    # Get sequences from file or form data
    fasta_content = None
    if file:
        fasta_content = (await file.read()).decode('utf-8')
    elif sequences:
        fasta_content = sequences
    else:
        raise HTTPException(status_code=400, detail="No sequences provided. Upload a file or provide sequences text.")

    # Write to temp file for processing
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        tmp.write(fasta_content)
        tmp_path = tmp.name

    try:
        designer = get_designer()
        # Parse available enzymes from comma-separated string
        enzyme_list = None
        if available_enzymes:
            enzyme_list = [e.strip() for e in available_enzymes.split(',') if e.strip()]

        designs = designer.design_constructs(
            tmp_path,
            arm_length=arm_length,
            half_site_min=half_site_min,
            half_site_max=half_site_max,
            max_designs=max_designs,
            cost_weight=cost_weight,
            length_weight=length_weight,
            available_enzymes=enzyme_list
        )

        # Group by gene
        gene_designs = {}
        for design in designs:
            if design.gene_name not in gene_designs:
                gene_designs[design.gene_name] = {
                    'deletion_size': len(design.deletion_region),
                    'designs': []
                }
            gene_designs[design.gene_name]['designs'].append(design)

        # Build response
        genes = []
        for gene_name, data in gene_designs.items():
            gene_results = GeneDesigns(
                gene_name=gene_name,
                deletion_size=data['deletion_size'],
                designs=[
                    design_to_result(d, i+1)
                    for i, d in enumerate(data['designs'])
                ]
            )
            genes.append(gene_results)

        return DesignResponse(
            success=True,
            message=f"Generated {len(designs)} designs for {len(genes)} genes",
            genes=genes,
            summary={
                "total_genes": len(genes),
                "total_designs": len(designs),
                "parameters": {
                    "arm_length": arm_length,
                    "half_site_range": [half_site_min, half_site_max],
                    "max_designs": max_designs,
                    "cost_weight": cost_weight,
                    "length_weight": length_weight
                }
            }
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        os.unlink(tmp_path)


@router.post("/api/design/fasta")
async def design_constructs_fasta(
    file: Optional[UploadFile] = File(None),
    sequences: Optional[str] = Form(None),
    arm_length: int = Form(1500),
    half_site_min: int = Form(900),
    half_site_max: int = Form(1300),
    max_designs: int = Form(5),
    cost_weight: float = Form(0.5),
    length_weight: float = Form(0.5),
    available_enzymes: Optional[str] = Form(None)
):
    """
    Design constructs and return as downloadable FASTA file.
    """
    # Get the design response first
    response = await design_constructs(
        file=file, sequences=sequences,
        arm_length=arm_length, half_site_min=half_site_min,
        half_site_max=half_site_max, max_designs=max_designs,
        cost_weight=cost_weight, length_weight=length_weight,
        available_enzymes=available_enzymes
    )

    # Generate FASTA content
    fasta_lines = []
    for gene in response.genes:
        for design in gene.designs:
            if design.single_enzyme:
                enzymes = design.enzyme1_name
            else:
                enzymes = f"{design.enzyme1_name}+{design.enzyme2_name}"

            header = (f">{gene.gene_name}_design{design.design_num}_{enzymes}_"
                     f"down{design.downstream_arm_length}bp_up{design.upstream_arm_length}bp")
            fasta_lines.append(header)

            # Wrap sequence at 80 characters
            seq = design.construct_sequence
            for i in range(0, len(seq), 80):
                fasta_lines.append(seq[i:i+80])

    fasta_content = "\n".join(fasta_lines)

    return StreamingResponse(
        StringIO(fasta_content),
        media_type="text/plain",
        headers={"Content-Disposition": "attachment; filename=constructs.fa"}
    )


@router.post("/api/design/report")
async def design_constructs_report(
    file: Optional[UploadFile] = File(None),
    sequences: Optional[str] = Form(None),
    arm_length: int = Form(1500),
    half_site_min: int = Form(900),
    half_site_max: int = Form(1300),
    max_designs: int = Form(5),
    cost_weight: float = Form(0.5),
    length_weight: float = Form(0.5),
    available_enzymes: Optional[str] = Form(None)
):
    """
    Design constructs and return detailed text report.
    """
    # Get sequences from file or form data
    fasta_content = None
    if file:
        fasta_content = (await file.read()).decode('utf-8')
    elif sequences:
        fasta_content = sequences
    else:
        raise HTTPException(status_code=400, detail="No sequences provided.")

    # Write to temp file for processing
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        tmp.write(fasta_content)
        tmp_path = tmp.name

    try:
        designer = get_designer()
        # Parse available enzymes from comma-separated string
        enzyme_list = None
        if available_enzymes:
            enzyme_list = [e.strip() for e in available_enzymes.split(',') if e.strip()]

        designs = designer.design_constructs(
            tmp_path,
            arm_length=arm_length,
            half_site_min=half_site_min,
            half_site_max=half_site_max,
            max_designs=max_designs,
            cost_weight=cost_weight,
            length_weight=length_weight,
            available_enzymes=enzyme_list
        )

        # Generate report using the __str__ method of ConstructDesign
        report_lines = []
        report_lines.append("DELETION ARMS DESIGNER - CONSTRUCT REPORT")
        report_lines.append("=" * 80)
        report_lines.append(f"\nParameters:")
        report_lines.append(f"  Arm length: {arm_length} bp")
        report_lines.append(f"  Half-site range: {half_site_min}-{half_site_max} bp")
        report_lines.append(f"  Max designs: {max_designs}")
        report_lines.append(f"  Cost weight: {cost_weight}")
        report_lines.append(f"  Length weight: {length_weight}")
        report_lines.append(f"\nTotal designs: {len(designs)}")

        for design in designs:
            report_lines.append(str(design))

        report_content = "\n".join(report_lines)

        return StreamingResponse(
            StringIO(report_content),
            media_type="text/plain",
            headers={"Content-Disposition": "attachment; filename=constructs_report.txt"}
        )
    finally:
        os.unlink(tmp_path)
