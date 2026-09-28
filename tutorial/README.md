# Example VCF preparation
The tutorial VCFs are derived from the [GIAB Ashkenazim and Chinese trio datasets](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/). The commands below select 4 samples from joint GRCh38 VCFs and restrict the variants to the regions in `Twist_Exome_Core_Covered_Targets_hg38.bed`.

Run these commands from the repository root. They require `wget`, `bcftools`. The two source VCFs are large downloads (~ 1.3 GB each).

```bash
cd tutorial

# Download the joint trio VCFs and their indexes
ASH_URL="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/RTG_RTGJointTrio_08122019/GRCh38/family.merged.avr_0.1.vcf.gz"
CHI_URL="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/analysis/RTG_RTGJointTrio_06062019/GRCh38/family.merged.avr_0.1.vcf.gz"

wget -O ashkenazim.family.merged.avr_0.1.vcf.gz "$ASH_URL"
wget -O ashkenazim.family.merged.avr_0.1.vcf.gz.tbi "$ASH_URL.tbi"
wget -O chinese.family.merged.avr_0.1.vcf.gz "$CHI_URL"
wget -O chinese.family.merged.avr_0.1.vcf.gz.tbi "$CHI_URL.tbi"

# Restrict each VCF to the supplied exome regions
BED="Twist_Exome_Core_Covered_Targets_hg38.bed"
bcftools view -R "$BED" -Oz -o ashkenazim.family.merged.avr_0.1_bed.vcf.gz ashkenazim.family.merged.avr_0.1.vcf.gz
bcftools index ashkenazim.family.merged.avr_0.1_bed.vcf.gz

bcftools view -R "$BED" -Oz -o chinese.family.merged.avr_0.1_bed.vcf.gz chinese.family.merged.avr_0.1.vcf.gz
bcftools index chinese.family.merged.avr_0.1_bed.vcf.gz

# Merge the trios and keep the four samples used in the cohort example
bcftools merge -Ou ashkenazim.family.merged.avr_0.1_bed.vcf.gz chinese.family.merged.avr_0.1_bed.vcf.gz \
  | bcftools view -s NA24385,NA24149,NA24631,NA24695 -Oz \
      -o NA24385_NA24149_NA24631_NA24695.vcf.gz

# Rename the samples
bcftools reheader \
  -s <(printf '%s\n' \
    'NA24385 HG002' \
    'NA24149 HG003' \
    'NA24631 HG005' \
    'NA24695 HG007') \
  -o HG002_HG003_HG005_HG007.vcf.gz \
  NA24385_NA24149_NA24631_NA24695.vcf.gz
bcftools index -t HG002_HG003_HG005_HG007.vcf.gz
```

To extract individual VCF (e.g. `HG002`), run this command from `tutorial/` after creating the joint VCF:

```bash
bcftools view --samples HG002 -Oz -o HG002.vcf.gz HG002_HG003_HG005_HG007.vcf.gz
```