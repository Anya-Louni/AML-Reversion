# Experiment 2, Step 3: development-set results

Generated 2026-10-15 by `python -m pipeline.run --phase exp2-3`. Plan: [exp2_plan.md](exp2_plan.md). Only development data were used; the locked test set was not loaded into any model or metric.

- Development: 342 samples from 310 patients; 104 drugs (82 with AML/monocytic-line signatures, used by HM).
- Mutation/clinical features (chosen on development data): 25.
- Penalty per model chosen from {0.01, 0.1, 1, 10, 100, 1000} by 5-fold patient-level CV on the primary metric. Because the same folds choose the penalty and report the score, development numbers are slightly optimistic; only the Step 4 test result counts.
- HM note: the Phase 3 consensus already builds a drug's signature only from AML/monocytic lines whenever such lines exist, so recomputing with only those lines gives identical scores. HM is therefore HR restricted to those drugs.

## Cross-validated results (development)

| model                                                       |    alpha |   primary_mean |   ci_lo |   ci_hi |   n_patients |   pooled_pearson |   top5_hit_rate |
|:------------------------------------------------------------|---------:|---------------:|--------:|--------:|-------------:|-----------------:|----------------:|
| A global ranking                                            |  nan     |          0.000 |   0.000 |   0.000 |          310 |            0.000 |           0.200 |
| B mutation + clinical                                       | 1000.000 |          0.195 |   0.170 |   0.221 |          284 |            0.258 |           0.400 |
| B mutation + clinical [low-signal excluded]                 | 1000.000 |          0.187 |   0.158 |   0.215 |          238 |            0.246 |           0.388 |
| H hierarchy (PRIMARY)                                       | 1000.000 |          0.192 |   0.167 |   0.217 |          284 |            0.172 |           0.374 |
| H hierarchy (PRIMARY) [low-signal excluded]                 | 1000.000 |          0.194 |   0.166 |   0.222 |          238 |            0.165 |           0.356 |
| HR hierarchy + reversion                                    | 1000.000 |          0.198 |   0.173 |   0.222 |          284 |            0.189 |           0.364 |
| HR hierarchy + reversion [low-signal excluded]              | 1000.000 |          0.203 |   0.175 |   0.230 |          238 |            0.183 |           0.359 |
| HM hierarchy + myeloid-line reversion                       | 1000.000 |          0.205 |   0.179 |   0.232 |          283 |            0.199 |           0.364 |
| HM hierarchy + myeloid-line reversion [low-signal excluded] | 1000.000 |          0.213 |   0.186 |   0.242 |          237 |            0.193 |           0.359 |
| Phase 4 raw reversion score (descriptive)                   |  nan     |         -0.012 |  -0.025 |   0.002 |          284 |           -0.010 |           0.237 |

`primary_mean`: mean over patients of per-sample Spearman between predicted and observed residuals, with patient-level bootstrap 95% CI. `top5_hit_rate`: share of a sample's 5 highest-predicted drugs in its observed top 20% (chance 0.20).

## Paired differences vs the primary model (development CV, descriptive)

|                                               |   mean_diff |   p_two_sided |
|:----------------------------------------------|------------:|--------------:|
| H minus B mutation + clinical                 |     -0.0028 |        0.8773 |
| H minus HR hierarchy + reversion              |     -0.0054 |        0.1155 |
| H minus HM hierarchy + myeloid-line reversion |     -0.0141 |        0.0019 |

## Penalty grid (CV primary metric)

|           |   B mutation + clinical |   H hierarchy (PRIMARY) |   HR hierarchy + reversion |   HM hierarchy + myeloid-line reversion |
|----------:|------------------------:|------------------------:|---------------------------:|----------------------------------------:|
|    0.0100 |                  0.1823 |                  0.1872 |                     0.1875 |                                  0.1960 |
|    0.1000 |                  0.1823 |                  0.1872 |                     0.1875 |                                  0.1960 |
|    1.0000 |                  0.1827 |                  0.1872 |                     0.1878 |                                  0.1962 |
|   10.0000 |                  0.1848 |                  0.1882 |                     0.1890 |                                  0.1974 |
|  100.0000 |                  0.1933 |                  0.1917 |                     0.1957 |                                  0.2041 |
| 1000.0000 |                  0.1949 |                  0.1921 |                     0.1975 |                                  0.2054 |

All four models chose the largest penalty in the pre-registered grid. The curves have flattened there (H gains 0.0004 from 100 to 1000), so extending the grid would change little; the grid is kept as registered.

## Exploratory diagnostic (development only, not pre-registered)

The target removes each drug's mean but not each patient's overall sensitivity. A generally sensitive patient therefore tends to have larger residuals on drugs whose response varies most between patients, and a model could score well just by detecting general sensitivity, without matching specific drugs to patients. The checks below estimate how much each model relies on that pattern. The last rows remove the pattern from the target using each patient's fitted overall effect (an oracle quantity, used only for this diagnostic).

|                                                                          |   value |
|:-------------------------------------------------------------------------|--------:|
| observed residuals vs sign(patient effect) x drug SD                     |   0.175 |
| B mutation + clinical: predictions vs that pattern                       |   0.167 |
| B mutation + clinical: CV score with the pattern removed from the target |   0.185 |
| H hierarchy (PRIMARY): predictions vs that pattern                       |   0.069 |
| H hierarchy (PRIMARY): CV score with the pattern removed from the target |   0.191 |

If a model's score barely drops when the pattern is removed, its signal is mostly drug-specific.

The primary model (H) was fixed in the plan and is not changed by these results.
