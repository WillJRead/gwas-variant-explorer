# GWAS Variant Explorer 

This is a personal project where I explored **genome-wide association study (GWAS) data** for prostate cancer. The goal was to summarize the most significant SNPs and visualize the results using **Manhattan plots**. This project allowed me to practice Python data analysis, Pandas, and plotting with Matplotlib while working with real-world genomic datasets. I later adapated the project so that the scripts are now custom functions that can take any subset of GWAS associations from the GWAS catalogue and output a .csv of the most signficant SNPs and a basic manhattan plot, therefore allowing for anyone to download and use. 

---

## **Project Structure**
```
Prostate_Cancer_GWAS_Summary/
│
├── data/ <-- download your own data to this folder
│ ├── gwas_subset_prostate_cancer.tsv # Sample GWAS data used as an example in the code
│ └── .gitkeep # placeholder to keep folder in Git
│
├── results/  <-- output files are generated here
│ └── .gitkeep # placeholder 
│
├── scripts/ <-- holds custom scripts
│ ├── top_snps.py # Extract top N SNPs by p-value
│ └── manhattan_plot.py # Generate Manhattan plot
│
└── README.md
```
--- 

## **Usage**

Ensure before starting that you have downlaoded this project as a ZIP file, decompressed it and downloaded the GWAS associations data of interest into /data.

1. **Clone the repository**

```bash
git clone https://github.com/WillJRead/gwas-variant-explorer.git
cd gwas-variant-explorer
```

2. **Install required python packages**
Check requirements.txt and see which packages are required which are not currently installed

3. **Top SNPs Extraction**

```python
from scripts.top_snps import top_snps

# Extract top 10 SNPs by default
top_snps("your.file.name.tsv")

# Or extract a custom number of SNPs
top_snps("your.file.name.tsv", n_snps=20)
```

This will generate top_snps.csv in the results/ folder. Replace "your.file.name.tsv" with your own .tsv file name. 

```
4. **Manhattan Plot**

```
from scripts.manhattan_plot import make_manhattan

make_manhattan("your.file.name.tsv")
```
This will create manhattan_plot.png in the results/ folder.


Aimed as a personal project to practice genomic data analysis and visualization. Make sure the data/ folder contains the GWAS TSV file before running scripts. Results will always be saved in results/ for easy access. The scripts are modular, so you can reuse them with other GWAS datasets. Any suggestions for improvements or collaboration ideas are welcome!