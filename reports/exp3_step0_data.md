# Experiment 3, Step F0: data availability

Written 2026-10-15. Every accession below was checked against its actual GEO record and, where marked verified, its actual supplementary files were downloaded and parsed. Nothing here is taken from the task description without checking.

## Headline findings

1. **GSE216730 is not the dataset described in the brief.** It is bulk RNA-seq of **patient-derived xenografts (PDX) treated with SNDX-5613 for one week** (12 samples), not paired before/after patient trial samples. There are no responder/non-responder labels in the GEO record. **The F1 responder test has no data behind it** and should be dropped unless another source appears.
2. **Menin inhibitors are absent from our LINCS signature set.** VTP50469, SNDX-5613 and revumenib have no consensus signature, so our linear simulation cannot be run for them at all. (`VTP-27999` in our set is a renin inhibitor, an unrelated compound.) This removes the menin datasets from the core F1 direction test.
3. **Two ATRA datasets are ideal and fully usable**, both in HL-60, which is one of the cell lines behind our own tretinoin signature.
4. **A real combination dataset exists** (ATRA + palbociclib), with both drugs in our LINCS set. F4 can be validated against measured data rather than reported as untested hypotheses.
5. The ATRA time course timepoints are **0 h, 4 h, 2 d, 4 d**, not days 1, 3 and 6 as the brief assumed.
6. The published ATRA + LSD1-inhibitor positive control for F4 is verified, and the specific drug used (tranylcypromine) is in our LINCS set.

## Usable datasets

| Accession | Open | Type | Samples | Design | Drugs | In our LINCS set | Gene IDs | Verified |
|---|---|---|---|---|---|---|---|---|
| **GSE227542** | yes | bulk RNA-seq | 12 | HL-60 ATRA time course: 0 h, 4 h, 2 d, 4 d, 3 replicates each | tretinoin | **yes** (AML/monocytic: HL60, SKM1, THP1, U937) | per-sample count files | file list verified |
| **GSE276457** | yes | bulk RNA-seq | 12 | HL-60, 48 h: untreated / ATRA / palbociclib / **combination**, 3 replicates each | tretinoin, Palbociclib | **yes, both** (Palbociclib: AML/monocytic, THP1) | gene symbols | downloaded and parsed |
| **GSE228321** | yes | bulk RNA-seq | 18 | MOLM13, MV4-11, OCI-AML3, 16 h: control / VTP50469 / OTX015 / **VTP+OTX combination** (combination in MOLM13 only) | OTX-015 **yes** (non-myeloid tier); VTP50469 **no** | partial | Ensembl with version | downloaded and parsed |
| **GSE228322** | yes | scRNA-seq | 2 | PDX AML cells, control vs 500 nM VTP50469, 16 h | VTP50469 | **no** | 10x barcodes/features/matrix | file list verified |
| van Galen GSE116256 | already in project | scRNA-seq | 19 treated samples from 8 patients with paired day 0 | diagnosis plus follow-up days 14 to 171 | treatment not annotated in GEO | n/a | HGNC (harmonized) | in repo |

## Datasets examined and excluded

| Accession | Reason |
|---|---|
| GSE216730 | PDX, not patient trial; drug (SNDX-5613) absent from our LINCS set; no responder labels |
| GSE228323, GSE228325 | ChIP-seq and scATAC-seq, not expression |
| GSE201664 | Murine cells (I1DN, cKit+/Cd11b−); AG-120 absent from our set; cross-species transfer adds an untestable confound |
| enasidenib trial single-cell | No open GEO dataset found by direct search. Consistent with the brief's statement that it is controlled-access (EGA); I verified only its absence from GEO, not the EGA listing itself |

## Which datasets support which part

- **F1 (does the linear simulation predict real change?)**
  - Direction and maturity tests: **GSE227542** (ATRA, 3 treated timepoints) and **GSE276457** (ATRA and palbociclib single-agent arms). Both have the drug in our set.
  - Secondary: **GSE228321** OTX015 arms in MOLM13 (non-myeloid-tier signature, so a weaker test of transfer).
  - Trajectory test: **GSE227542** only, with the caveat that its timepoints are 4 h, 2 d and 4 d.
  - Responder test: **no data**, see finding 1.
  - Cell-line versus patient split: all usable F1 data is cell line or PDX. **There is no open patient before/after dataset with a drug in our LINCS set.** The cell-line transfer question therefore cannot be answered by F1 as specified, and F1 should report cell-line results only, saying plainly that the patient case is untested.
- **F2 (learned dynamics):** the same datasets. Total usable treated/control arms is small (about 10 drug-cell-line conditions, plus 3 ATRA timepoints). Holding out entire drugs is possible but the held-out set will be tiny; this constrains how complex a model can honestly be fitted and should be stated in the F2 plan.
- **F3 (resistance prediction):** **GSE228322** is the only open single-cell before/after AML dataset found (2 samples, 1 PDX, 1 drug). The van Galen treated samples add 8 patients but under chemotherapy, where malignant cells mostly disappear rather than mature (for example AML556: 2,062 malignant cells at day 0, 0 at day 15; AML707B: 1,370 to 61 at day 18). Both are thin; F3 should be labelled exploratory from the start.
- **F4 (combinations):** **GSE276457** gives a measured ATRA + palbociclib combination arm with both drugs in our set, so a genuine validation is possible. **GSE228321** adds VTP+OTX in MOLM13, but only one of the two drugs is in our set. The published positive control (ATRA + LSD1 inhibitor) is verified below.

## Positive control for F4, verified

According to PubMed, Schenk et al. 2012 (*Nature Medicine* 18:605-611, [DOI](https://doi.org/10.1038/nm.2661)) report that LSD1 inhibitors, **specifically tranylcypromine**, unlock an ATRA differentiation response in non-APL AML, with the combination superior to either drug alone in vivo. Kahl et al. 2019 (*Leukemia* 33:2628-2639, [DOI](https://doi.org/10.1038/s41375-019-0581-y)) report the same direction, with LSD1 inhibition potentiating a GCN5-inhibition-plus-ATRA response. Both tretinoin and tranylcypromine have AML/monocytic-tier signatures in our set, so the pair can be scored. Note this evidence is for differentiation and engraftment endpoints, not for a transcriptome combination measurement.

## Drug coverage check against our Phase 3 set

In our set (AML/monocytic tier unless noted): tretinoin, AM-580, TTNPB, calcitriol, tranylcypromine, decitabine, AGI-5198, AGI-6780, JQ1, I-BET-151, I-BET-762, vorinostat, entinostat, brequinar, Palbociclib. Non-myeloid tier only: OTX-015, alitretinoin, Panobinostat. **Absent:** all menin inhibitors, enasidenib, ivosidenib, azacitidine, GSK-2879552, ORY-1001/iadademstat.

## Two corrections to the brief

- **There is no "Experiment 3" in this repository.** The project contains Phase 4 (first pre-registered validation), the Experiment 2 ceiling analysis, and Experiment 2. Step F5 should say so rather than referring to a third experiment.
- **There is no `paper/` directory.** The write-ups are `README.md`, `reports/final_report.md` and the per-step reports. F5 should either extend `reports/final_report.md` or create `paper/` explicitly as new work.

## Recommendation

Proceed with F1 using GSE227542 and GSE276457 as the primary tests and GSE228321 (OTX015) as a secondary transfer test, and drop the responder test. Expect F1 to be a cell-line-only validation; this is a real limitation of the open data, not a choice, and the F1 report should lead with it.

Unverified: the per-patient treatment in van Galen is not annotated in the GEO files, so the brief's statement that all but one patient received standard chemotherapy could not be confirmed from the open data. What is verifiable is that malignant cells largely disappear at follow-up, which is consistent with cytotoxic therapy rather than differentiation.
