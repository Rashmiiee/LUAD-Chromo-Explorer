
# LUAD Chromo-Explorer

**LUAD Chromo-Explorer** is an interactive web application for exploring **chromatin and gene expression features in lung adenocarcinoma (LUAD)**, focusing on **chromosome 1 genes**. The app integrates **genomic, epigenomic, and transcriptomic data**, enabling users to examine the relationship between chromatin compartments, gene regulation, and expression patterns.  

---

## Features

1. **Gene-centric exploration**  
   - Input any ENSG gene ID from chromosome 1 to retrieve its detailed profile.  
   - View chromatin state (A/B compartment, eigenvalues, domain) and expression summary (tumor vs normal).  

2. **Regulatory scoring**  
   - Computes a **regulatory score** combining chromatin eigenvalues and differential expression.  
   - Highlights genes with unusual regulation patterns, e.g., **open compartment but repressed** in tumors.  

3. **Interactive visualizations**  
   - **Tumor vs Normal Expression**: Boxplot of expression distribution.  
   - **Expression vs Eigenvalue**: Scatterplot linking compartment activity with expression.  
   - **Neighborhood Plot**: Visualizes nearby genes ±200kb, showing expression and compartment context.  

4. **Chromatin-expression anomaly detection**  
   - Detect genes in open compartments that are **unexpectedly downregulated**, highlighting potential regulatory repression events.  

---

## Installation

1. Clone the repository:  
```bash
git clone https://github.com/<your_username>/luad-chromo-explorer.git
cd luad-chromo-explorer
```

2. Create a python environment and install the packages:
   
```bash
conda create -n luad_explorer python=3.10
conda activate luad_explorer
pip install -r requirements.txt
```

3. Run the app
```bash
streamlit run app_full.py
```

Data Sources

Chromosome 1 genes: genes_chr1.bed (Gencode annotation)

Chromatin compartments: compartments.txt (Hi-C eigenvalues/domains)

LUAD expression: TCGA RSEM TPM values from Toil Xena Hub

Phenotype metadata: TCGA PanCanAtlas phenotype table

Methodology

Data preprocessing

Genes are intersected with chromatin compartment data using PyRanges.

LUAD expression data is merged per gene.

Regulatory score calculation
For each gene:

Compute mean tumor vs normal expression → logFC

Combine with compartment eigenvalue:

synergy = eigen * logFC * penalty
norm = tanh(synergy)
regulatory_score = 2 / (1 + exp(-shift_scale * norm)) - 1


Penalty is applied if the gene is in a repressive compartment (B).

High positive scores indicate upregulation in open regions; negative scores can flag repression in open chromatin.

Visualization

Expression boxplots, eigenvalue scatterplots, and neighborhood plots allow users to intuitively assess regulation patterns and genomic context.

Example Use Cases

Identify genes in open compartments that are unexpectedly downregulated in LUAD tumors.

Explore tumor vs normal expression patterns for key oncogenes or tumor suppressors.

Visualize the chromatin landscape of neighboring genes to detect local regulatory trends.

