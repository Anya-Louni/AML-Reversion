# Reversion Landscape: final report

An in-silico hypothesis generator for differentiation therapy research in AML. Course project. All results are computational hypotheses from public data; nothing here is medical advice.

## 1. Introduction

Acute myeloid leukaemia (AML) is a cancer of myeloid progenitors that fail to complete maturation. Differentiation therapy, used clinically with all-trans retinoic acid (ATRA) in acute promyelocytic leukaemia and with IDH inhibitors in IDH-mutant AML, works by letting malignant cells resume maturation. We asked whether public data can be combined to (1) describe where each patient's malignant cells sit on the normal maturation path, (2) rank existing compounds by how strongly their expression signatures point toward maturity for that patient, and (3) check those rankings against real ex vivo drug responses.

## 2. Related work

- **Connectivity Map repurposing.** CMap and LINCS L1000 (Lamb et al. 2006; Subramanian et al. 2017) score compounds by how their expression signatures oppose or mimic a disease signature. Our reversion score follows this idea, using a maturation direction learned from single cells instead of a case-control signature.
- **AML hierarchies.** van Galen et al. 2019 classified malignant AML cells into six normal-like types and showed that monocyte-like cells suppress T cells. Their annotations are the backbone of our reference and validation.
- **Drug response prediction in Beat AML.** Tyner et al. 2018 and Bottomly et al. 2022 linked ex vivo sensitivity to mutations, expression and deconvolved differentiation state. Bottomly et al. already showed that differentiation state shapes drug response, so our deconvolution is a building block, not a new finding.
- **Perturbation modelling.** Methods such as scGen or CellOracle model expression responses to perturbations. We use a far simpler linear shift, chosen for transparency.

## 3. Methods

Full methods and every parameter choice, with justifications, are in the phase reports. In brief:

- **Data (Phase 0).** van Galen GSE116256 (41,090 cells, 43 samples), Beat AML 2.0 open processed data (707 RNA-seq samples, 166 drugs, mutations), CMap LINCS 2020 Level 5 signatures read by HTTP range requests (12,414 signatures), all harmonised to HGNC symbols. Drug names matched by normalisation plus a curated synonym table: 140 of 166 Beat AML drugs matched.
- **Landscape (Phase 1).** 7,956 normal myeloid cells; 2,000 HVGs; PCA; Harmony by sample. Diffusion pseudotime rooted in the healthy HSC with the highest authors' HSC score, checked against Palantir and 10 alternative roots. PHATE embedding; potential surface from pseudotime and density. A fixed PCA + kNN map places any cell on the landscape. Maturation signature: top 200 genes up and down with pseudotime, excluding ribosomal and mitochondrial genes.
- **Patients (Phase 2).** Maturity threshold fixed on normal cells (maximum balanced accuracy). Malignant cells mapped with an off-manifold flag. Beat AML bulk deconvolved by NNLS with van Galen reference profiles, validated on pseudo-bulk samples with the patient held out. Stuck signatures: patient expression minus expected mature expression.
- **Drugs (Phase 3).** Consensus signature per compound from AML/monocytic lines where available. Reversion score: cosine between the drug signature and the patient's push-to-maturity vector; KS connectivity as a second method. Positive controls fixed in the config before scoring.
- **Validation (Phase 4).** Pre-registered plan committed to git before the code. 520 samples, 479 patients, 123 drugs. Patient-level 5-fold CV. Baselines: random, global drug ranking, mutation-only and fraction-only ridge models; our raw and calibrated scores; a combined model.
- **Simulation (Phase 5).** Expression shifted linearly along each drug's signature, rescaled to a common length chosen so that on average at most 20% of cells leave the normal manifold; each frame re-mapped.

## 4. Results

**Landscape.** Pseudotime orders the annotated stages at Spearman 0.920 (Palantir 0.903; agreement 0.900; root sensitivity minimum 0.990). Held-out mapping error is 0.025 on a 0-1 scale, 0.028 with whole samples held out. CD34 and PROM1 fall, CD14, LYZ and CSF1R rise.

**Patients.** The maturity threshold (0.461) separates normal immature from mature cells with 98.7% balanced accuracy. Malignant cells' mapped scores follow the authors' labels at Spearman 0.84. Stuck fractions range from 2% to 100% across 12 patients. Deconvolved mature share is recovered at Pearson 0.895 on held-out pseudo-bulk; fine cell types less well (overall r 0.61). 95% of signature genes change in the same direction in bulk as in single cells.

**Drug ranking.** Tier-1 positive controls rank high among 1,691 compounds in AML/monocytic lines: AM-580 1st, tretinoin 14th, ingenol mebutate and TTNPB above the 95th percentile (AUROC 0.81, p = 0.003). The LSD1 inhibitor tranylcypromine ranks 2nd. Calcitriol scores low (8th percentile), a clear miss. Rankings are very similar across patients (mean between-patient Spearman 0.93 in van Galen, 0.61 in Beat AML).

**Validation.** Reversion scores carry a weak signal about which drugs a sample responds to (P1: mean Spearman 0.145, p = 0.0001), but no patient-specific signal once each drug's average potency is removed (P2: -0.008, p = 0.95). The global drug ranking reaches 0.696. The combined model (0.720) beats it, but removing the reversion score from it changes little (+0.0005 raw, +0.002 drug-centred). In an exploratory analysis, letting each drug learn its own slope gives a drug-centred Spearman of 0.094, because the direction of the association differs between drugs. Mature share correlates with venetoclax resistance (Spearman -0.35), matching published biology. The ex vivo assay shows no extra sensitivity of IDH-mutant samples to the tool IDH inhibitors, so it cannot test IDH-directed differentiation.

**Simulation.** With common-length shifts, 16% of cells end off the normal manifold on average. Outcomes vary widely by patient; tretinoin moves a mean of 17% of stuck cells past the threshold (0% to 100% across patients).

## 5. Limitations

- LINCS signatures come from cell lines after short exposures; most AML/monocytic signatures are from THP-1 alone.
- Ex vivo viability is only a proxy for maturation, and the Beat AML panel contains no clinical differentiation agents.
- The maturation signature and mapping rely on similarity to normal cells; malignant cells can look mature in expression and still behave abnormally.
- The linear shift does not model dose, time, toxicity, resistance mechanisms or cell fate dynamics.
- Only 12 single-cell patients are available, some with few malignant cells.
- Results depend on the pre-specified choices documented in each report; alternative reasonable choices could change rankings.

## 6. Future work

- Wet-lab testing of top-ranked compounds with differentiation readouts (CD11b, CD14, morphology) rather than viability.
- Perturbation response models trained on myeloid cells and longer exposures.
- Ex vivo datasets that include differentiation agents and maturation endpoints.
- Larger single-cell cohorts with genotype information.

## References

- Bottomly D. et al. *Cancer Cell* 2022.
- Lamb J. et al. *Science* 2006.
- Pei S. et al. *Cancer Discovery* 2020.
- Subramanian A. et al. *Cell* 2017.
- Tyner J.W. et al. *Nature* 2018.
- van Galen P. et al. *Cell* 2019.
