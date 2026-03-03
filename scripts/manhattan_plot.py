from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def make_manhattan(file_name):
    """
    Generate a Manhattan plot from a GWAS TSV file.
    Saves plot as results/manhattan_plot.png.
    """
    # Paths
    script_dir = Path(__file__).resolve().parent
    file_path = script_dir.parent / "data" / file_name
    results_dir = script_dir.parent / "results"
    results_dir.mkdir(exist_ok=True)

    # Load data
    gwas_data = pd.read_csv(file_path, sep="\t")

    # Ensure numeric values
    gwas_data["P-VALUE"] = pd.to_numeric(gwas_data["P-VALUE"], errors="coerce")
    gwas_data = gwas_data.dropna(subset=["P-VALUE"])

    # Compute -log10(P)
    gwas_data["minus_log10_p"] = -np.log10(gwas_data["P-VALUE"])

    # Handle chromosomes: convert to int if possible
    gwas_data = gwas_data[gwas_data["CHR_ID"].notna()]
    gwas_data["CHR_ID"] = pd.to_numeric(gwas_data["CHR_ID"], errors="coerce")
    gwas_data = gwas_data.dropna(subset=["CHR_ID"])
    gwas_data["CHR_ID"] = gwas_data["CHR_ID"].astype(int)

    # Sort by chromosome and position
    gwas_data["CHR_POS"] = pd.to_numeric(gwas_data["CHR_POS"], errors="coerce")
    gwas_data = gwas_data.dropna(subset=["CHR_POS"])
    gwas_data = gwas_data.sort_values(["CHR_ID", "CHR_POS"])

    # Compute cumulative positions
    chroms = gwas_data["CHR_ID"].unique()
    cumulative_pos = []
    cum_offset = 0
    chrom_offsets = {}
    for chrom in chroms:
        chrom_data = gwas_data[gwas_data["CHR_ID"] == chrom]
        chrom_offsets[chrom] = cum_offset
        cumulative_pos.extend(chrom_data["CHR_POS"] + cum_offset)
        cum_offset += chrom_data["CHR_POS"].max() + 1e6  # small gap between chromosomes

    gwas_data["cumulative_pos"] = cumulative_pos

    # Alternate colors for chromosomes
    colors = ["#1f77b4", "#ff7f0e"]
    gwas_data["color"] = [colors[i % 2] for i in gwas_data["CHR_ID"]]

    # Plot
    plt.figure(figsize=(12,6))
    plt.scatter(
        gwas_data["cumulative_pos"],
        gwas_data["minus_log10_p"],
        c=gwas_data["color"],
        s=10
    )
    plt.axhline(y=-np.log10(5e-8), color="red", linestyle="--", label="Genome-wide significance")
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P-value)")
    plt.title("Manhattan Plot")

    # X-ticks at chromosome centers
    chrom_centers = []
    for chrom in chroms:
        chrom_data = gwas_data[gwas_data["CHR_ID"] == chrom]
        center = (chrom_data["cumulative_pos"].min() + chrom_data["cumulative_pos"].max()) / 2
        chrom_centers.append(center)
    plt.xticks(chrom_centers, chroms)

    plt.tight_layout()
    save_path = results_dir / "manhattan_plot.png"
    plt.savefig(save_path, dpi=300)
    plt.show()

    print(f"Manhattan plot saved to: {save_path}")

# Example usage:
make_manhattan("example_gwas_data_prostate_cancer.tsv")