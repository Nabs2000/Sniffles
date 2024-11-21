import sys

import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd


def parse_vcf(vcf_file):
    """
    Parse a VCF file to extract structural variants.
    """
    vcf = pysam.VariantFile(vcf_file)
    data = []
    for record in vcf:
        chrom = record.chrom
        pos = record.pos
        svtype = record.info.get("SVTYPE", "NA")  # Structural variant type
        svlen = abs(record.info.get("SVLEN", [0])[0])  # Variant length (absolute value)
        data.append({"chrom": chrom, "pos": pos, "svtype": svtype, "svlen": svlen})
    return pd.DataFrame(data)


def plot_ideogram(data, chromosome_lengths):
    """
    Create an ideogram to visualize structural variants.
    """
    # Unique chromosomes
    chromosomes = sorted(data["chrom"].unique())
    chrom_height = 0.6
    chrom_spacing = 1.5

    # Initialize the plot
    fig, ax = plt.subplots(figsize=(10, len(chromosomes) * chrom_spacing + 1))
    ax.set_xlim(0, max(chromosome_lengths.values()))
    ax.set_ylim(-chrom_spacing, len(chromosomes) * chrom_spacing)
    ax.axis("off")

    # Plot chromosomes
    for i, chrom in enumerate(chromosomes):
        y = i * chrom_spacing
        chrom_len = chromosome_lengths.get(chrom, 0)
        ax.add_patch(patches.Rectangle((0, y), chrom_len, chrom_height, color="lightgray"))
        ax.text(-10 ** 6, y + chrom_height / 2, chrom, ha="right", va="center", fontsize=10)

    # Plot structural variants
    colors = {"DEL": "red", "INS": "blue", "DUP": "orange", "INV": "green", "NA": "gray"}
    for _, row in data.iterrows():
        chrom = row["chrom"]
        pos = row["pos"]
        svtype = row["svtype"]
        svlen = row["svlen"]
        y = chromosomes.index(chrom) * chrom_spacing + chrom_height / 2
        ax.scatter(pos, y, s=svlen / 1000, color=colors.get(svtype, "gray"), alpha=0.7, label=svtype)

    # Add title and legend
    ax.set_title("Ideogram of Structural Variants", fontsize=14)
    handles = [patches.Patch(color=color, label=svtype) for svtype, color in colors.items()]
    ax.legend(handles=handles, loc="upper right", title="SV Type")

    plt.tight_layout()
    plt.savefig("ideogram.png")


def main():
    vcf_file = sys.argv[1]  # Path to your VCF file
    chromosome_lengths = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
        "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
        "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468,
        "chrX": 156040895, "chrY": 57227415  # Example lengths; add as needed
    }

    # Parse the VCF and plot
    data = parse_vcf(vcf_file)
    plot_ideogram(data, chromosome_lengths)
