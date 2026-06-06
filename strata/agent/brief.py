def write_brief(drug, assoc, driver_genes, n_samples) -> str:
    return f"""# Commercial-Development Opportunity — {drug}

**Discovered subgroup:** defined by methylation at {driver_genes[0]} — high-methylation
stratum (cluster {assoc.responder_label}, {assoc.grade} responder),
n={n_samples}, AUC effect={assoc.effect_size:.3f}, p={assoc.pvalue:.3g}.

**Methylation signature:** led by {driver_genes[0]}, with supporting markers {", ".join(driver_genes[1:6])}.

**Why it's invisible to current taxonomy:** the stratum is defined by a methylation
state, not by histology or a single mutation — it cuts across the indication.

**Commercial implication:** a methylation companion-marker could define a responder
indication for {drug}, raising trial enrichment and peak-year revenue while reducing
opportunity-identification time and cost.
"""
