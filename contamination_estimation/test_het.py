import pysam
import matplotlib.pyplot as plt
import seaborn as sns

# Filepath to your VCF file
vcf_file = "/gpfs/data/bergstrom/harry/variant_calling/wolf/pseudohaploid/wolf_ph_fix.vcf"

# Name of the sample to process
target_sample = "Neige_2_3"

# Open the VCF file
vcf = pysam.VariantFile(vcf_file)

# Check if the target sample exists in the VCF
if target_sample not in vcf.header.samples:
    raise ValueError(f"Sample '{target_sample}' not found in the VCF file.")

# List to store ref/alt allele ratios
ratios = []

# Counters for skipped sites
zero_depth_count = 0
missing_depth_count = 0
invalid_ad_count = 0

# Process each record in the VCF
for record in vcf:
    # Get the sample data for the target sample
    sample = record.samples[target_sample]

    # Extract allele depths (AD field)
    ad = sample.get("AD")
    if ad and len(ad) == 2:  # Ensure AD field exists and has ref/alt depths
        ref_depth, alt_depth = ad
        if ref_depth is not None and alt_depth is not None:  # Check for valid depths
            if ref_depth + alt_depth > 0:  # Avoid division by zero
                ratio = ref_depth / (ref_depth + alt_depth)
                if ratio >= 0.01:  # Omit ratios less than 0.01
                    ratios.append(ratio)
            else:
                zero_depth_count += 1
        else:
            missing_depth_count += 1
    else:
        invalid_ad_count += 1

# Close the VCF file
vcf.close()

# Print summary of skipped sites
print(f"Summary of skipped sites:")
print(f"  Sites with zero depth: {zero_depth_count}")
print(f"  Sites with missing depth values: {missing_depth_count}")
print(f"  Sites with invalid AD fields: {invalid_ad_count}")

# Check if any ratios were calculated
if not ratios:
    print("No sites with valid allele depth were found for the sample.")
else:
    print(f"Calculated {len(ratios)} ref/alt allele ratios for the sample '{target_sample}'.")

# Plot the distribution of ref/alt allele ratios
sns.histplot(ratios, kde=True, bins=50, color="blue")
plt.axvline(0.5, color="red", linestyle="--", label="Theoretical 50%")
plt.xlabel("Ref/Alt Allele Ratio")
plt.ylabel("Frequency")
plt.title(f"Distribution of Ref/Alt Allele Ratios for {target_sample}")
plt.legend()
plt.savefig(f"het_distribution_{target_sample}.png", dpi=300, bbox_inches="tight")
plt.show()