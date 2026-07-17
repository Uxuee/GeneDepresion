# MDD WGCNA Biological Annotation Report

This report annotates prioritized non-grey WGCNA modules from the full Step6 MDD/control run. No module eigengene survived BH-FDR < 0.10 for Diagnosis, so all biological interpretations are exploratory and intended for follow-up prioritization.

## Methods and Databases

- Modules annotated: 10 prioritized/secondary modules.
- Background: region-specific WGCNA gene universe from the saved best-solution module assignments.
- GO-like enrichment: custom hypergeometric tests over MSigDB C5 GO:BP, GO:CC, and GO:MF gene sets because clusterProfiler was not installed.
- Pathways: MSigDB Hallmark, C2/CP, and available Reactome/KEGG-style MSigDB sets.
- Cell type: built-in minimal marker sanity-check list because no local marker resource was found.
- AHBA context: integrated from Step 11 published AHBA module annotation tables.

## Top Biological Themes by Module

- Ventral Striatum MElavenderblush3: neuronal/synaptic; confidence: high-priority exploratory; top GO BP: Protein localization to membrane; Localization within membrane; Protein localization to plasma membrane
- Ventral Striatum MElightcyan1: neuronal/synaptic; confidence: high-priority exploratory; top GO BP: Endoplasmic reticulum to golgi vesicle mediated transport; Lipid biosynthetic process; Zymosterol metabolic process
- Ventral Striatum MEdarkred: neuronal/synaptic; confidence: high-priority exploratory; top GO BP: Response to wounding; Neuronal stem cell population maintenance; Hemostasis
- Ventral Striatum MElightgreen: neuronal/synaptic; confidence: high-priority exploratory; top GO BP: Vesicle mediated transport in synapse; Synaptic vesicle maturation; Behavior
- Insula MEorange: neuronal/synaptic; confidence: high-priority exploratory; top GO BP: Commitment of neuronal cell to specific neuron type in forebrain; Forebrain neuron fate commitment; Glial cell derived neurotrophic factor receptor signaling pathway
- Hippocampal Formation MEgreen: mixed/unclear; confidence: high-priority exploratory; top GO BP: Response to lectin; Mesenchymal cell migration; Glial cell differentiation
- Frontal Cortex MEdarkgreen: ciliary/ependymal-like; confidence: high-priority exploratory; top GO BP: Detection of chemical stimulus involved in sensory perception of taste; Sensory perception of bitter taste; Sensory perception of taste
- Frontal Cortex MEblack: mixed/unclear; confidence: high-priority exploratory; top GO BP: Sensory perception of chemical stimulus; Sensory perception of smell; Detection of chemical stimulus
- Frontal Cortex MEsalmon: mixed/unclear; confidence: high-priority exploratory; top GO BP: Mrna 5 splice site recognition; Protein rna complex organization; Mrna splice site recognition
- Insula MEdarkorange: mixed/unclear; confidence: candidate/supportive; top GO BP: Multicellular organismal reproductive process; Regulation of neutrophil extravasation; Acrosomal vesicle exocytosis

## Strong Tier2/ML or AHBA Evidence

- Ventral Striatum MElavenderblush3: best Tier2/ML FDR=0.00159; best AHBA=M6 FDR=9.08e-04
- Ventral Striatum MElightcyan1: best Tier2/ML FDR=0.00159; best AHBA=M18 FDR=0.00177
- Ventral Striatum MEdarkred: best Tier2/ML FDR=0.00638; best AHBA=M31 FDR=1.12e-12
- Ventral Striatum MElightgreen: best Tier2/ML FDR=0.478; best AHBA=M9 FDR=4.43e-12
- Insula MEorange: best Tier2/ML FDR=0.219; best AHBA=M2 FDR=1.08e-04
- Hippocampal Formation MEgreen: best Tier2/ML FDR=1.83e-07; best AHBA=M30 FDR=1.79e-101
- Frontal Cortex MEdarkgreen: best Tier2/ML FDR=1; best AHBA=M26 FDR=4.21e-04
- Frontal Cortex MEblack: best Tier2/ML FDR=1; best AHBA=M17 FDR=9.24e-06
- Frontal Cortex MEsalmon: best Tier2/ML FDR=1.88e-06; best AHBA=M14 FDR=0.251
- Insula MEdarkorange: best Tier2/ML FDR=0.961; best AHBA=M24 FDR=0.121

## Limitations

- These are prioritized exploratory modules, not validated disease modules.
- No module eigengene survived Diagnosis BH-FDR < 0.10.
- org.Hs.eg.db, clusterProfiler, ReactomePA, and enrichplot were unavailable in this environment; the script used MSigDB and custom enrichment fallbacks where possible.
- The cell-type marker analysis is a minimal sanity check and should not be overinterpreted.
- Enrichment terms suggest biological context but do not establish causality or diagnostic utility.
