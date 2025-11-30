# app_full.py
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
from pathlib import Path
import pyranges as pr
import time

sns.set(style="whitegrid")
st.set_page_config(page_title="LUAD Gene Browser (chr1)", layout="wide")

# -------------------------------------------------------------------
# 0️⃣ Preprocessing (create gene_df_chr1.csv if missing)
# -------------------------------------------------------------------
if not Path("data/gene_df_chr1.csv").exists():
    start_time = time.time()
    st.info("Preprocessing gene_df_chr1.csv... This may take ~3-5 minutes depending on internet and CPU speed.")

    # 1️⃣ Load chromosome 1 genes BED
    genes_df = pd.read_csv("data/genes_chr1.bed", sep="\t")
    genes_df = genes_df.rename(columns={
        "chr": "Chromosome",
        "start": "Start",
        "end": "End",
        "gene_name": "gene"
    })
    if "gene_id" not in genes_df.columns:
        genes_df["gene_id"] = genes_df["gene"]  # fallback if BED lacks gene_id
    genes_pr = pr.PyRanges(genes_df[["Chromosome","Start","End","gene","gene_id"]])

    # 2️⃣ Load compartment info (chr1 only)
    comp = pd.read_csv(
        "data/compartments.txt",
        sep="\t",
        header=None,
        skiprows=2,
        names=["chr","start","end","eigen","domain"]
    )
    comp = comp[comp["chr"]=="chr1"].reset_index(drop=True)
    comp["compartment"] = comp["domain"].map({"open":"A","closed":"B"})
    comp_pr = pr.PyRanges(comp.rename(columns={"chr":"Chromosome","start":"Start","end":"End"}))

    # 3️⃣ Join genes with compartments
    joined = genes_pr.join(comp_pr)
    joined_df = joined.df
    collapsed = joined_df.groupby("gene").agg({
        "gene_id":"first",
        "Chromosome":"first",
        "Start":"min",
        "End":"max",
        "compartment": lambda s: s.value_counts().idxmax(),
        "eigen":"mean",
        "domain": lambda s: s.value_counts().idxmax()
    }).reset_index()

    # 4️⃣ Load LUAD expression
    url_expr = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz"
    expr = pd.read_csv(url_expr, sep="\t")
    pheno2 = pd.read_csv(
        "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
        sep="\t"
    )
    luad_pheno = pheno2[pheno2["_primary_disease"]=="lung adenocarcinoma"]
    luad_primary = luad_pheno[luad_pheno["sample_type"]=="Primary Tumor"]["sample"].unique()
    luad_normal = luad_pheno[luad_pheno["sample_type"]=="Solid Tissue Normal"]["sample"].unique()

    expr_cols = list(expr.columns[1:])
    tumor_cols = [s for s in luad_primary if s in expr_cols]
    normal_cols = [s for s in luad_normal if s in expr_cols]

    expr = expr.rename(columns={"sample":"gene_id"})
    expr_subset = expr[["gene_id"] + tumor_cols + normal_cols]
    expr_subset["gene_id_clean"] = expr_subset["gene_id"].str.split(".").str[0]
    collapsed["gene_id_clean"] = collapsed["gene_id"].str.split(".").str[0]

    # 5️⃣ Merge and compute summary columns
    merged = collapsed.merge(expr_subset.drop(columns=["gene_id"]), on="gene_id_clean", how="inner")
    merged["tumor_mean"] = merged[tumor_cols].mean(axis=1)
    merged["normal_mean"] = merged[normal_cols].mean(axis=1)
    merged["logFC"] = merged["tumor_mean"] - merged["normal_mean"]

    # 6️⃣ Save merged gene_df
    Path("data").mkdir(parents=True, exist_ok=True)
    merged.to_csv("data/gene_df_chr1.csv", index=False)
    st.success(f"gene_df_chr1.csv created in {time.time()-start_time:.1f} seconds.")

# -------------------------------------------------------------------
# 1️⃣ Load gene_df
# -------------------------------------------------------------------
@st.cache_data
def load_gene_df():
    df = pd.read_csv("data/gene_df_chr1.csv")
    df = df.set_index("gene_id_clean")
    return df

gene_df = load_gene_df()
tumor_cols = [c for c in gene_df.columns if c.endswith("-01")]
normal_cols = [c for c in gene_df.columns if c.endswith("-11")]

# -------------------------------------------------------------------
# 2️⃣ DESCRIPTION
# -------------------------------------------------------------------
st.title("LUAD Gene Browser (chr1 genes only)")
st.markdown("""
Explore the interplay between chromatin organization and gene expression in lung adenocarcinoma. This app provides a chromosome 1–focused view of how regulatory compartments and tumor-specific expression patterns shape gene activity**.
- Only **chromosome 1 genes** are available.
- You can view:
    - Chromatin state
    - Expression summary in tumor vs normal samples
    - Regulatory score
    - Plots for expression and neighboring genes
- Enter a gene **ENSG ID** below to start.
""")

# -------------------------------------------------------------------
# 3️⃣ UTILITY & PLOT FUNCTIONS
# -------------------------------------------------------------------
def get_chromatin_state(gene_id_clean):
    row = gene_df.loc[gene_id_clean]
    return {
        "gene": row["gene"],
        "gene_id": gene_id_clean,
        "chromosome": row["Chromosome"],
        "start": int(row["Start"]),
        "end": int(row["End"]),
        "compartment": row["compartment"],
        "domain": row["domain"],
        "eigenvalue": float(row["eigen"])
    }

def get_expression_summary(gene_id_clean):
    row = gene_df.loc[gene_id_clean]
    tumor_vals = row[tumor_cols].astype(float)
    normal_vals = row[normal_cols].astype(float)
    return {
        "tumor_mean": float(tumor_vals.mean()),
        "normal_mean": float(normal_vals.mean()),
        "logFC": float(tumor_vals.mean() - normal_vals.mean()),
        "tumor_median": float(tumor_vals.median()),
        "normal_median": float(normal_vals.median())
    }

def compute_regulatory_score(gene_id_clean, pseudocount=0.1, comp_penalty=0.6, shift_scale=3):
    row = gene_df.loc[gene_id_clean]
    chrom_shift = float(row["eigen"])
    compartment = row["compartment"]
    expr = get_expression_summary(gene_id_clean)
    logFC = expr["logFC"]
    penalty = comp_penalty if compartment=="B" else 1.0
    synergy = chrom_shift * logFC * penalty
    norm = np.tanh(synergy)
    final_score = 2 / (1 + np.exp(-shift_scale * norm)) - 1
    return {
        "chrom_shift": chrom_shift,
        "compartment": compartment,
        "logFC": logFC,
        "synergy_raw": synergy,
        "normalized": norm,
        "regulatory_score": final_score
    }

def get_neighborhood(gene_id_clean, window=200000):
    row = gene_df.loc[gene_id_clean]
    chr_ = row["Chromosome"]
    start = row["Start"]
    end = row["End"]
    region = gene_df[(gene_df["Chromosome"]==chr_) & 
                     (gene_df["Start"] >= start - window) &
                     (gene_df["End"] <= end + window)].copy()
    return region

def is_open_but_repressed(gene_id, gene_df, tumor_cols, normal_cols, logfc_threshold=-0.5):
    """
    Check if a single gene is in an open (A) compartment but repressed in tumors.
    
    Parameters
    ----------
    gene_id : str
        The gene_id (ENSG) to check.
    gene_df : pd.DataFrame
        DataFrame containing gene info, compartments, and expression.
    tumor_cols : list
        Columns corresponding to tumor samples.
    normal_cols : list
        Columns corresponding to normal samples.
    logfc_threshold : float
        Log fold change cutoff for considering a gene repressed in tumor (negative values).
        
    Returns
    -------
    bool
        True if the gene is in an open compartment and downregulated, False otherwise.
    dict
        Info about the gene (compartment, tumor mean, normal mean, logFC).
    """
    if gene_id not in gene_df.index:
        return False, None
    
    row = gene_df.loc[gene_id].copy()
    tumor_mean = row[tumor_cols].astype(float).mean()
    normal_mean = row[normal_cols].astype(float).mean()
    logFC = tumor_mean - normal_mean
    
    info = {
        "compartment": row["compartment"],
        "tumor_mean": tumor_mean,
        "normal_mean": normal_mean,
        "logFC": logFC
    }
    
    is_open_repressed = (row["compartment"] == "A") and (logFC < logfc_threshold)
    return is_open_repressed, info


# -------------------------------------------------------------------
# 3️⃣ PLOT FUNCTIONS (improved aesthetics & captions)
# -------------------------------------------------------------------
def plot_tumor_normal_expression(gene_id_clean):
    row = gene_df.loc[gene_id_clean]
    tumor_vals = row[tumor_cols].astype(float).values.flatten()
    normal_vals = row[normal_cols].astype(float).values.flatten()
    
    plt.figure(figsize=(2,2))  # smaller, compact
    sns.boxplot(data=[tumor_vals, normal_vals], palette=["#e41a1c","#377eb8"])
    plt.xticks([0,1], ["Tumor","Normal"])
    plt.ylabel("Expression (TPM)")
    plt.title(f"Tumor vs Normal Expression for {gene_id_clean}")
    st.pyplot(plt.gcf())
    plt.close()
    
    st.markdown("**Interpretation:** This boxplot compares expression levels in tumor (red) vs normal (blue) tissues. "
                "Higher values in tumor indicate potential upregulation in LUAD, while lower values suggest downregulation.")

def plot_eigen_vs_expression(gene_id_clean):
    row = gene_df.loc[gene_id_clean]
    eigen = float(row["eigen"])
    tumor_vals = row[tumor_cols].astype(float).values.flatten()
    normal_vals = row[normal_cols].astype(float).values.flatten()
    
    plt.figure(figsize=(3,3))
    plt.axvspan(-1,0,color="#a6cee3",alpha=0.2,label="B (Closed)")
    plt.axvspan(0,1,color="#fb9a99",alpha=0.2,label="A (Open)")
    sns.scatterplot(x=[eigen]*len(tumor_vals), y=tumor_vals, color="#e41a1c", label="Tumor", s=50)
    sns.scatterplot(x=[eigen]*len(normal_vals), y=normal_vals, color="#377eb8", label="Normal", s=50)
    plt.hlines(tumor_vals.mean(), eigen-0.05, eigen+0.05, colors="#e41a1c", lw=2)
    plt.hlines(normal_vals.mean(), eigen-0.05, eigen+0.05, colors="#377eb8", lw=2)
    plt.xlabel("Eigenvalue")
    plt.ylabel("Expression (TPM)")
    plt.title(f"Expression vs Eigenvalue for {gene_id_clean}")
    st.pyplot(plt.gcf())
    plt.close()
    
    st.markdown("**Interpretation:** Eigenvalue indicates chromatin openness: positive = A (open), negative = B (closed). "
                "Red points (tumor) above blue (normal) show genes upregulated in tumor. The background color indicates compartment type.")

def plot_neighborhood(gene_id_clean, window=200000):
    neigh = get_neighborhood(gene_id_clean, window)
    if neigh.empty:
        st.warning("No neighboring genes found in this window.")
        return
    neigh["tumor_mean"] = neigh[tumor_cols].mean(axis=1)
    neigh["normal_mean"] = neigh[normal_cols].mean(axis=1)
    
    plt.figure(figsize=(8,3))
    sns.scatterplot(x=range(len(neigh)), y=neigh["tumor_mean"], color="#e41a1c", s=50, label="Tumor")
    sns.scatterplot(x=range(len(neigh)), y=neigh["normal_mean"], color="#377eb8", s=50, label="Normal")
    for i, comp in enumerate(neigh["compartment"]):
        plt.axvspan(i-0.4, i+0.4, color="#fb9a99" if comp=="A" else "#a6cee3", alpha=0.2)
    plt.xticks(range(len(neigh)), neigh["gene"], rotation=45, ha="right")
    plt.ylabel("Mean Expression (TPM)")
    plt.title(f"Neighborhood of {gene_id_clean} (±{window//1000} kb)")
    plt.legend()
    st.pyplot(plt.gcf())
    plt.close()
    
    st.markdown("**Interpretation:** This scatterplot shows expression of the queried gene and its neighbors. "
                "Background shading indicates chromatin compartment (pink=A, light blue=B). Differences between red (tumor) and blue (normal) reveal local co-regulation patterns.")


# -------------------------------------------------------------------
# 4️⃣ STREAMLIT INTERFACE
# -------------------------------------------------------------------
gene_input = st.text_input("Enter a gene ENSG ID", "ENSG00000198691")

if gene_input in gene_df.index:
    chrom_state = get_chromatin_state(gene_input)
    expr_summary = get_expression_summary(gene_input)
    reg_score = compute_regulatory_score(gene_input)
    st.subheader("Gene & Chromatin Info")
    st.table(pd.DataFrame([chrom_state]))
    st.subheader("Expression Summary (LUAD)")
    st.table(pd.DataFrame([expr_summary]))
    st.subheader("Regulatory Score")
    st.table(pd.DataFrame([reg_score]))
    st.markdown("""
    The **regulatory score** integrates chromatin organization with tumor-specific expression changes to highlight genes that may be **functionally     dysregulated** in LUAD.  
            
     - **Chromatin context**: Open (A) compartments generally favor transcription, while closed (B) compartments are less accessible.  
     - **Expression change**: Tumor vs normal differential expression highlights genes that are up- or down-regulated in cancer.  
     The score combines these features to prioritize genes that are **both in an active chromatin environment and show tumor-specific   
     dysregulation**, revealing potential drivers or key regulators that would be missed by looking at expression or chromatin alone.
     """)

    if st.button("Check if gene is open but repressed"):
        if gene_input:
            open_repressed, info = is_open_but_repressed(gene_input, gene_df, tumor_cols, normal_cols)
            if open_repressed:
                st.warning(f"{gene_input} is in an open compartment but downregulated in tumors.")
            else:
                st.success(f"{gene_input} is NOT in an open compartment and repressed.")
            st.table(info)
    
    if st.button("Plot Tumor vs Normal Expression"):
        plot_tumor_normal_expression(gene_input)
        
    if st.button("Plot Expression vs Eigenvalue"):
        plot_eigen_vs_expression(gene_input)
        
    if st.button("Plot Neighborhood (±200kb)"):
        plot_neighborhood(gene_input)
else:
    st.error("Gene not found in LUAD chr1 dataset.")
