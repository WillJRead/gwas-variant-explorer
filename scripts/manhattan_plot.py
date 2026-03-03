from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def make_manhattan(file_name):

    #ensuring pathing is correct
    script_dir = Path(__file__).resolve().parent
    file_path = script_dir.parent / "data" / file_name
    results_dir = script_dir.parent / "results"
    results_dir.mkdir(exist_ok=True)

    #load data
    gwas_data = pd.read_csv(file_path, sep="\t")

    #ensure numeric values
    gwas_data["P-VALUE"] = pd.to_numeric(gwas_data["P-VALUE"], errors="coerce")
    gwas_data = gwas_data.dropna(subset=["P-VALUE"])

    #compute -log10(P)
    gwas_data["minus_log10_p"] = -np.log10(gwas_data["P-VALUE"])

    #handle chromosomes: convert to int if possible
    gwas_data = gwas_data[gwas_data["CHR_ID"].notna()]
    gwas_data["CHR_ID"] = pd.to_numeric(gwas_data["CHR_ID"], errors="coerce")
    gwas_data = gwas_data.dropna(subset=["CHR_ID"])
    gwas_data["CHR_ID"] = gwas_data["CHR_ID"].astype(int)

    #sort by chromosome and position
    gwas_data["CHR_POS"] = pd.to_numeric(gwas_data["CHR_POS"], errors="coerce")
    gwas_data = gwas_data.dropna(subset=["CHR_POS"])
    gwas_data = gwas_data.sort_values(["CHR_ID", "CHR_POS"])

    #compute cumulative positions
    chroms = gwas_data["CHR_ID"].unique()
    cumulative_pos = []
    cum_offset = 0
    chrom_offsets = {}
    for chrom in chroms:
        chrom_data = gwas_data[gwas_data["CHR_ID"] == chrom]
        chrom_offsets[chrom] = cum_offset
        cumulative_pos.extend(chrom_data["CHR_POS"] + cum_offset)
        cum_offset += chrom_data["CHR_POS"].max() + 1e6  

    gwas_data["cumulative_pos"] = cumulative_pos

    #alternate colors for chromosomes
    colors = ["#1f77b4", "#ff7f0e"]
    gwas_data["color"] = [colors[i % 2] for i in gwas_data["CHR_ID"]]

    #plot
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

    #x-ticks at chromosome centers
    chrom_centers = []
    for chrom in chroms:
        chrom_data = gwas_data[gwas_data["CHR_ID"] == chrom]
        center = (chrom_data["cumulative_pos"].min() + chrom_data["cumulative_pos"].max()) / 2
        chrom_centers.append(center)
    plt.xticks(chrom_centers, chroms)

    plt.tight_layout()
    
    #saving plot
    save_path = results_dir / "manhattan_plot.png"
    plt.savefig(save_path, dpi=300)
    plt.show()

    print(f"Manhattan plot saved to: {save_path}")

#example usage:
#make_manhattan("example_gwas_data_prostate_cancer.tsv")