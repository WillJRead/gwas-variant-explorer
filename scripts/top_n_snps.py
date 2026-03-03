#defining function to take GWAS data and create a csv of top n significantlly signfiicant SNPs 
def top_snps(file_name, n_snps=None):
    from pathlib import Path
    import pandas as pd

    #ensuring pathing is correct for user
    script_dir = Path(__file__).resolve().parent
    file_path = script_dir.parent / "data" / file_name
    results_dir = script_dir.parent / "results"
    results_dir.mkdir(exist_ok=True)
    
    #load data
    gwas_data = pd.read_csv(file_path, sep="\t")
    
    # clean numeric columns
    gwas_data["P-VALUE"] = pd.to_numeric(gwas_data["P-VALUE"], errors="coerce")
    gwas_data = gwas_data.dropna(subset=["P-VALUE"])

    # get top 10 significant SNPs
    top_hits = gwas_data.sort_values("P-VALUE").head(n_snps)
    
    #save CSV
    save_path = results_dir / f"top_{n_snps}_snps.csv"
    top_hits.to_csv(save_path, index=False)
    print(f"Top {n_snps} SNPs saved to: {save_path}")
    
#call the function
top_hits = top_snps("gwas_subset_prostate_cancer.tsv", n_snps=15)
print(top_hits)

