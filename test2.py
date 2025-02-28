import numpy as np

def compute_gene_contributions(adverse_effects, drug_to_genes, drug_to_ae):
    """
    Computes the contribution of each gene to each adverse effect.
    
    Parameters:
    - adverse_effects: dict mapping AE to list of associated drugs
    - drug_to_genes: dict mapping each drug to {gene: impact_score}
    - drug_to_ae: dict mapping each drug to {AE: probability}
    
    Returns:
    - gene_contributions: dict mapping AE to {gene: contribution_score}
    """
    gene_contributions = {}
    
    for ae, drugs in adverse_effects.items():
        gene_scores = {}
        
        for drug in drugs:
            if drug in drug_to_genes and drug in drug_to_ae and ae in drug_to_ae[drug]:
                for gene, impact in drug_to_genes[drug].items():
                    score = impact * drug_to_ae[drug][ae]  # Compute weighted contribution
                    gene_scores[gene] = gene_scores.get(gene, 0) + score
        
        # Normalize scores
        total_score = sum(gene_scores.values())
        if total_score > 0:
            gene_contributions[ae] = {g: s / total_score for g, s in gene_scores.items()}
        else:
            gene_contributions[ae] = gene_scores  # Avoid division by zero
    
    return gene_contributions

# Example Data
adverse_effects = {
    "Liver Toxicity": ["DrugA", "DrugB"],
    "Cardiac Issues": ["DrugB", "DrugC"]
}

drug_to_genes = {
    "DrugA": {"Gene1": 0.8, "Gene2": 0.6},
    "DrugB": {"Gene2": 0.5, "Gene3": 0.7},
    "DrugC": {"Gene3": 0.9, "Gene4": 0.4}
}

drug_to_ae = {
    "DrugA": {"Liver Toxicity": 0.7},
    "DrugB": {"Liver Toxicity": 0.5, "Cardiac Issues": 0.8},
    "DrugC": {"Cardiac Issues": 0.6}
}

# Compute gene contributions
gene_contributions = compute_gene_contributions(adverse_effects, drug_to_genes, drug_to_ae)

# Print results
for ae, genes in gene_contributions.items():
    print(f"Adverse Effect: {ae}")
    for gene, score in genes.items():
        print(f"  {gene}: {score:.4f}")
