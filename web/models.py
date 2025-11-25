"""
Pydantic models for the Deletion Arms Designer API
"""

from typing import List, Optional
from pydantic import BaseModel, Field


class DesignRequest(BaseModel):
    """Request model for construct design"""
    sequences: Optional[str] = Field(None, description="FASTA formatted sequences (alternative to file upload)")
    arm_length: int = Field(1500, ge=500, le=5000, description="Length of homology arms in bp")
    half_site_min: int = Field(900, ge=0, description="Minimum distance from deletion for half-sites")
    half_site_max: int = Field(1300, ge=0, description="Maximum distance from deletion for half-sites")
    max_designs: int = Field(5, ge=1, le=100, description="Maximum designs per gene")
    cost_weight: float = Field(0.5, ge=0, le=1, description="Weight for enzyme cost optimization")
    length_weight: float = Field(0.5, ge=0, le=1, description="Weight for arm length optimization")


class EnzymeInfo(BaseModel):
    """Enzyme information for API response"""
    name: str
    recognition_seq: str
    full_site: str
    left_half: str
    right_half: str
    cost_tier: str
    cost_per_unit: float
    neb_price: float
    neb_units: int


class DesignResult(BaseModel):
    """Single construct design result"""
    design_num: int
    enzyme1_name: str
    enzyme1_recognition: str
    enzyme1_cost_per_unit: float
    enzyme2_name: str
    enzyme2_recognition: str
    enzyme2_cost_per_unit: float
    single_enzyme: bool
    needs_stuffer: bool
    downstream_arm_length: int
    upstream_arm_length: int
    total_arm_length: int
    total_enzyme_cost: float
    downstream_half_site_position: int
    upstream_half_site_position: int
    construct_sequence: str
    stuffer_sequence: str


class GeneDesigns(BaseModel):
    """All designs for a single gene"""
    gene_name: str
    deletion_size: int
    designs: List[DesignResult]


class DesignResponse(BaseModel):
    """Full response for design request"""
    success: bool
    message: str = ""
    genes: List[GeneDesigns] = []
    summary: dict = {}


class EnzymeListResponse(BaseModel):
    """Response for enzyme list endpoint"""
    enzymes: List[EnzymeInfo]
    total: int
