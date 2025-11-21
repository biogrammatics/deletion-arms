# New Configurable Features - Summary

## 1. **Configurable Homology Arm Length**

**Parameter:** `--arm-length` or `-a`  
**Default:** 1500 bp  
**Usage:**
```bash
./deletion_arms_designer.py input.fasta --arm-length 2000
```

**Effect:**
- Input sequences are parsed as: `[arm_length bp upstream][deletion][arm_length bp downstream]`
- Allows flexibility for different experimental requirements
- Minimum sequence length = 2 × arm_length

**Example:**
- With `--arm-length 1000`: Uses 1000 bp arms instead of default 1500 bp
- Gene_1 (3918 bp total) → 1000 bp upstream + 1918 bp deletion + 1000 bp downstream

---

## 2. **Configurable Half-Site Placement Range**

**Parameter:** `--half-site-range` or `-r`  
**Default:** 900 1300 (bp from deletion)  
**Usage:**
```bash
./deletion_arms_designer.py input.fasta --half-site-range 600 800
```

**Effect:**
- Specifies distance from deletion region where half-sites are searched
- Allows optimization of recombination efficiency vs stability
- Range is measured in bp from the deletion boundary

**Technical Details:**
- For downstream arm: searches positions `half_site_min` to `half_site_max` from start
- For upstream arm: searches positions `(arm_length - half_site_max)` to `(arm_length - half_site_min)` from start

**Example:**
- Default (900-1300): Half-sites are 900-1300 bp away from deletion
- Custom (600-800): Half-sites are 600-800 bp away from deletion (closer to deletion)

---

## 3. **Custom Stuffer Insert Sequence**

**Parameter:** `--stuffer` or `-s`  
**Default:** Auto-generated AT-rich filler  
**Usage:**
```bash
./deletion_arms_designer.py input.fasta --stuffer marker.fasta
```

**Effect:**
- Uses your custom sequence (e.g., selection marker, reporter gene) as stuffer insert
- Tool automatically flanks with appropriate restriction enzyme half-sites
- No restriction site checking in stuffer (it's removed in final construct)

**Stuffer Structure:**
```
[Right-half RE1] + [Your Custom Insert] + [Left-half RE2]
```

**Example Output:**
```
Default (no stuffer file):
  Stuffer fragment (50 bp): GGGATATATATATATATATATATATATATATATATATATATATATATATCAC

With custom KanMX marker (795 bp):
  Stuffer fragment (801 bp): GGG + [795 bp KanMX] + CAC
```

---

## Combined Usage Examples

### Example 1: Custom Everything
```bash
./deletion_arms_designer.py genes.fasta \
  --arm-length 2000 \
  --half-site-range 1000 1500 \
  --stuffer KanMX.fasta \
  --vector pRS426.fasta \
  --max-designs 5 \
  --output results.txt
```

### Example 2: Shorter Arms, Closer Half-Sites
```bash
./deletion_arms_designer.py genes.fasta \
  --arm-length 800 \
  --half-site-range 400 600
```

### Example 3: Just Custom Stuffer
```bash
./deletion_arms_designer.py genes.fasta \
  --stuffer GFP_marker.fasta
```

---

## Input File Requirements

### Gene Sequences (required)
- **Format:** FASTA
- **Structure:** `[arm_length bp upstream][deletion][arm_length bp downstream]`
- **Minimum length:** 2 × arm_length bp

### Stuffer Insert (optional)
- **Format:** FASTA
- **Content:** Selection marker, reporter gene, or any insert sequence
- **Note:** Tool adds restriction half-sites automatically

### Vector (optional)
- **Format:** FASTA  
- **Purpose:** Exclude enzymes that cut your cloning vector

---

## Output Changes

When custom stuffer is used, the output shows:

```
Stuffer fragment (801 bp): GGG[custom_795bp_insert]CAC
```

Instead of default:

```
Stuffer fragment (50 bp): GGGATATATATATATATATATATATATATATATATATATATATATATCAC
```

---

## Validation Tests Completed

✅ **Default settings** (1500 bp arms, 900-1300 range, auto stuffer)  
✅ **Custom stuffer** (795 bp KanMX marker)  
✅ **Custom arm length** (1000 bp arms)  
✅ **Custom half-site range** (600-800 bp from deletion)  
✅ **Combined custom settings** (all parameters together)  

All features working correctly!

