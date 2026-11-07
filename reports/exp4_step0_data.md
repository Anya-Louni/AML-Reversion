# Experiment 4, Step F0: perturbation data availability

Written 2026-11-07. Every accession below was fetched from its actual GEO record, and every file marked
*parsed* was downloaded and read. Nothing here is taken from a paper abstract, a search result or the
brief without checking. Where a claim could not be checked, it is marked unverified.

## Why Experiment 4 exists

Experiment 3 closed with four specific dead ends, all of them data problems rather than modelling problems:

1. the simulation has **no time axis**, and three conditionable drugs could not support fitting one (F2);
2. **per-cell resistance is not measurable** from the open data, so F3's pre-registered bar was never met;
3. **one** measured combination exists in the whole open set, so F4 rests on a single example;
4. every result is a **cell line or xenograft**, and the maturation readout is always *expression*, never protein.

The common cause is that LINCS gives us drug signatures, not measured cell states after a defined
intervention. Perturbation data gives the latter. This step asks which perturbation datasets actually
exist, what is in them, and which of the four dead ends each one can and cannot open.

## Headline findings

1. **GSE213957 (THP-1 CaRPool-seq) is the strongest match in the open data, and it was missed by the
   Experiment 3 audit.** It is 30,707 THP-1 cells — the AML line behind most of our own LINCS
   signatures — carrying **187 perturbation conditions over 28 target genes**, of which **28,095 cells are
   dual perturbations**. It also ships a **58-plex surface-protein panel including CD11b, CD14, CD15, CD34,
   CD117, CD33 and HLA-DR**. That is the protein differentiation readout the final report's future-work
   section asks for, already measured.
2. **It turns F4 from one example into a real experiment.** Experiment 3 could test combination prediction
   against exactly one measured combination (ATRA + palbociclib). Here there are **158 dual-perturbation
   conditions whose two single perturbations are also measured**, which is the design a synergy test needs.
   An additive model still cannot express synergy; the difference is that we could now *measure* how often
   it fails, instead of reporting that it cannot in principle.
3. **KDM1A is a perturbation target** (17 conditions, 2,265 cells). Tranylcypromine, a KDM1A/LSD1
   inhibitor, was our **2nd-ranked compound in Phase 3** and the F4 positive control. The genetic and the
   pharmacological arm can be compared directly for the first time. HDAC3 (vorinostat, entinostat) and
   EP300 are likewise targets. **This is a target-level analogy, not the same intervention** — see Limits.
4. **GSE294098 (MOLM-13 Perturb-seq) supplies the menin arm Experiment 3 had to drop.** F0 of Experiment 3
   found that menin inhibitors are absent from LINCS entirely, which killed the menin test. This dataset
   knocks out **MEN1, DOT1L and KAT6A genetically**, so no drug signature is needed. F3's one surviving
   result — maturity rises +0.051 under a menin inhibitor — becomes independently testable.
5. **Replogle 2022 is not on GEO and is not myeloid.** A GEO author search returns only co-authored series,
   none of them the genome-scale dataset. The data is distributed through `gwps.wi.mit.edu` (verified live,
   HTTP 200) and Figshare+. It is K562 and RPE1. Given that `config.yaml` already demotes K562 to a
   lower-priority hematopoietic line because it is CML and erythroid-leaning, Replogle is a **weaker fit than
   the three AML datasets below** and is not recommended as a primary source.
6. **Norman 2019 (GSE133344) is usable but carries the same cell-context problem.** It is K562 CRISPRa.
   Its myeloid TFs are genuinely there (CEBPA, CEBPB, CEBPE, SPI1, KLF1, ETS2, EGR1, JUN, CNN1 confirmed
   in the parsed file, including the pair `CEBPE_RUNX1T1`), but any maturation claim from it inherits an
   erythroid-leaning background that our own config flags.

## Usable datasets

| Accession | Open | Assay | Cells / samples | Design | Targets | Readout | Gene IDs | Verified |
|---|---|---|---|---|---|---|---|---|
| **GSE213957** | yes | CaRPool-seq (Cas13d knockdown) scRNA-seq + CITE-seq | 30,707 THP-1 cells; 4 lanes, 8 HTOs | 187 conditions: 28,095 dual, 2,188 single, 424 NT | 28 genes | 21,052-gene transcriptome **+ 58 surface proteins** | Ensembl + symbol (10x 3-col) | downloaded (810 MB) and parsed |
| **GSE294098** | yes | Perturb-seq (CRISPR KO) | MOLM-13; 2 GSMs (GEX + CRISPR); paper reports 31,015 cells | 69 guides, 3 per target | 22 genes + AAVS1 control | 38,606-gene transcriptome | Ensembl + symbol | feature and barcode files parsed |
| **GSE294096** | yes | bulk RNA-seq (companion to the above) | 8 samples | validation arm | same study | read counts | — | file list verified |
| **GSE305037** | yes | Perturb-seq (Parse Evercode) | 12 samples | 186 sgRNAs, then **A-485 or DMSO**, 24 h | p300/CREBBP acetyl-interactome | raw count matrices | — | record parsed |
| GSE133344 | yes | Perturb-seq (CRISPRa) | 111,445 cells (filtered) | 290 guide-identity strings; paper reports 236 perturbations (105 single, 131 combinatorial) | ~109 genes incl. myeloid TFs | transcriptome | 10x | cell-identity file parsed |

GSE213957 target genes, as parsed from the metadata (not from the paper): ATXN7L3, CD14, DOT1L, EP300,
FLI1, GATA2, GFI1, GSE1, HDAC3, IKZF1, INTS1, KDM1A, MED13L, MED14, MED16, MED19, MED22, MED24, MYB,
RCOR1, RING1, SMARCA5, SMARCD1, SSRP1, SUPT16H, SUPT5H, SUPT6H, VPRBP.

GSE294098 targets, as parsed from the 69 guide features: BCOR, BRD1, BRPF1, BRPF3, DOT1L, EED, EZH2,
ING4, ING5, JADE1, JADE2, JADE3, KAT6A, KAT7, MLLT1, MLLT10, MEN1, PCGF1, PSIP1, RING1, RYBP, SUZ12,
plus AAVS1 safe-harbour controls.

## Datasets examined and excluded

| Accession | Reason |
|---|---|
| GSE264667 | Returned by a Replogle author search, but it is HepG2 and Jurkat (224 samples), not AML and not the genome-scale K562/RPE1 dataset |
| Replogle 2022 genome-scale (gwps.wi.mit.edu, Figshare+) | Not on GEO; K562 and RPE1 only. Non-myeloid, and K562 is already demoted in our config. Kept as a possible ceiling-analysis source only, not a primary one |
| Norman 2019 as a *primary* source | Usable, but K562 CRISPRa; the cell-context mismatch would have to lead any result, as the cell-line limitation led F1 |

## Which datasets support which part

- **Direction, with a protein readout (the F1 question, done properly).** GSE213957 only. It is the one
  dataset where a maturation prediction can be checked against **measured CD11b / CD14 / CD34 / CD117**
  rather than against expression alone. Experiment 3 could never do this.
- **The menin result, independently.** GSE294098 (MEN1, DOT1L, KAT6A knockouts in MOLM-13). This tests
  the single F3 finding that survived its confound, without needing a LINCS signature that does not exist.
- **Combinations and synergy (the F4 question, done properly).** GSE213957: 158 dual conditions with both
  singles measured. 162 conditions have at least 50 cells and 125 have at least 100, so most pairs are
  powered; the weakest are not.
- **Drug-by-genotype interaction.** GSE305037 (MOLM-13 and OCI-AML3, sgRNA library then A-485 or DMSO).
  Both lines are in our `myeloid_cell_lines`. This is the only candidate that crosses a genetic
  perturbation with an actual drug, so it is the natural bridge back to the reversion score.
- **A ceiling analysis for perturbation-response models.** Any of them, and this is where our existing
  ceiling method transfers most directly. It needs replicate structure; GSE213957's four lanes and eight
  hashtags are the most promising, but **whether they constitute genuine technical replicates is unverified**
  and must be established before any ceiling number is quoted.
- **Time.** *No candidate solves this.* Every dataset here is a single timepoint (16 h, 24 h, or one
  harvest). The largest gap Experiment 3 identified is **still not addressable with open data**, and
  Experiment 4 should say so up front rather than quietly dropping it.
- **Patient cells.** Still nothing. Every dataset is a cell line.

## Limits that must lead the Experiment 4 plan, not be buried in it

- **A knockdown is not a drug.** CaRPool-seq uses Cas13d RNA targeting; Perturb-seq uses Cas9 knockout.
  Our reversion score is defined on *drug* signatures. KDM1A knockdown and tranylcypromine share a target,
  not a mechanism, a dose or a time course. Any comparison between the two arms is an analogy and has to be
  labelled as one, the way F2's cell-line-matched result was labelled a lead rather than a gain.
- **Most of these targets are chromatin and transcription machinery, not differentiation TFs.** INTS1,
  the Mediator subunits, SSRP1 and SUPT5H/SUPT16H are core transcriptional components whose knockdown may
  cause general collapse rather than directed maturation. A maturity score will move; that movement is not
  evidence of differentiation. Pre-specify which targets are differentiation hypotheses (GFI1, MYB, GATA2,
  FLI1, IKZF1, KDM1A, DOT1L, EP300) and which are expected-collapse controls, **before** scoring anything.
- **INTS1 is an anchor in the design.** It is the most frequent partner across pairs (its conditions are the
  largest in the dataset). The pair set is therefore not a uniform 28x28 grid, and any synergy analysis has
  to account for the design rather than assume balance.
- **Single-perturbation power is thin.** Only 2,188 cells are single perturbations across 29 conditions,
  from 290 cells (INTS1) down to 15 (MYB). Several single-gene effects will not be estimable, which
  directly limits how many pairs have usable singles.
- **The maturation-axis overlap check is deferred.** Our Phase 1 signature lives in `data/processed`, which
  a fresh clone does not have. What F0 can state is that all candidates are **whole-transcriptome**
  (21,052 and 38,606 genes), not LINCS landmark space, so signature coverage is not a risk. The exact
  overlap is an F1 step after `make phase1`.
- **Experiment 4 inherits Experiment 2's reuse problem if it is careless.** These are new datasets, so a
  genuinely untouched test set is possible here. It is worth taking: pre-register the split before looking
  at any outcome, as Experiment 2 did.

## Two corrections to the framing this step started from

- **"Perturb-seq (Norman 2019, Replogle 2022)" is the wrong pair of sources for this project.** Both are
  K562/RPE1. The AML perturbation datasets (GSE213957, GSE294098, GSE305037) are a better match on cell
  context, and two of them did not exist publicly when the Experiment 3 audit was written. Norman stays as
  a secondary cross-cell-type check, not a primary source.
- **GSE294098's own records disagree about its size.** The series design says "Perturb seq against 16 AML
  target genes" and the summary says 16 epigenetic regulators, but the guide feature list contains
  **22 distinct targets plus AAVS1 controls** (69 guides, exactly 3 per label). The discrepancy is
  unexplained; F1 should use the parsed guide list, not the stated count.

## Recommendation

Proceed with **GSE213957 as the primary dataset** and **GSE294098 as the menin/MOLM-13 replication**, with
GSE305037 held for the drug-by-genotype bridge. Drop Replogle as a primary source and keep Norman as a
secondary check only.

The honest framing for Experiment 4 is narrower than "model perturbation responses". It is: *we now have
measured cell states after defined interventions in the right cell lines, with a protein maturation
readout and real combinations — so the F1 and F4 questions can be asked properly, the menin result can be
replicated, and the time-axis gap remains unsolved.*

## Unverified

- Whether GSE213957's lanes and hashtags give genuine technical replicates, which a ceiling analysis would
  need. Only the metadata columns were read, not the experimental design in the paper.
- The Figshare+ record for Replogle 2022 returned HTTP 202 (a challenge page), so its contents were not
  inspected; only `gwps.wi.mit.edu` was confirmed to be live.
- GSE294098's per-cell guide assignment was not computed. The CRISPR Guide Capture features are present in
  the matrix, so it is recoverable, but the matrix itself (3.2M raw barcodes) was not downloaded here.
- Whether DOT1L, EP300 or MEN1 inhibitors exist in our Phase 3 compound set. Tranylcypromine, vorinostat
  and entinostat are listed as present by the Experiment 3 audit; the others were not checked.
