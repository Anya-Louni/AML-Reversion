// Data layer: typed loaders for the precomputed exports in /public/data.

export interface Manifest {
  generated: string;
  threshold: number;
  frames: number;
  alpha_max: number;
  like_types: string[];
  quantization: { x: [number, number]; y: [number, number]; score: [number, number] };
  status_codes: Record<string, string>;
  patients: PatientSummary[];
  disclaimer: string;
}

export interface PatientSummary {
  id: string;
  n_malignant_d0: number;
  n_shown: number;
  mutations: string[];
  stuck_fraction: number;
  median_score: number;
  file: string;
  top_drug: string;
}

export interface RankedDrug {
  drug: string;
  rank: number;
  cosine: number;
  connectivity: number;
  tier: string;
  cells: string;
  n_sigs: number;
  in_beataml_panel: boolean;
  beat_median_auc?: number | null;
  beat_n?: number | null;
  validation_rho?: number | null;
}

export interface SimDrug {
  drug: string;
  file: string;
  status: number[];
  n: number;
  already_mature: number;
  stuck_start: number;
  matured: number;
  resistant: number;
  transient: number;
  pct_matured_of_stuck: number;
}

export interface PatientDetail extends PatientSummary {
  cell_types: number[];
  simulated_drugs: SimDrug[];
  ranking: RankedDrug[];
}

export interface Landscape {
  n: number;
  x: [number, number];
  y: [number, number];
  height: number[];
  score?: number[];
  support?: number[];
  height_vis?: number[];
  threshold: number;
  labels: { name: string; x: number; y: number }[];
  reference: { x: number[]; y: number[]; score: number[]; type: string[] };
}

export interface Trajectory {
  frames: number;
  cells: number;
  // frame-major: [frame][cell][x, y, score], dequantized
  data: Float32Array;
}

const BASE = `${import.meta.env.BASE_URL}data/`;
const cache = new Map<string, Promise<unknown>>();

function getJSON<T>(path: string): Promise<T> {
  if (!cache.has(path)) cache.set(path, fetch(BASE + path).then((r) => {
    if (!r.ok) throw new Error(`${path}: ${r.status}`);
    return r.json();
  }));
  return cache.get(path) as Promise<T>;
}

export const loadManifest = () => getJSON<Manifest>("manifest.json");
export const loadLandscape = () => getJSON<Landscape>("landscape.json");
export const loadPatient = (file: string) => getJSON<PatientDetail>(file);
export const loadValidation = () => getJSON<Validation>("validation.json");
export const loadHierarchy = () => getJSON<Hierarchy>("hierarchy.json");
export const loadExp3 = () => getJSON<Exp3>("exp3.json");

export interface Exp3 {
  f1: {
    direction: { id: string; drug: string; role: string; cosine: number; rank: number; n_null: number; passes: boolean; variance_explained: number }[];
    n_pass: number; n_total: number;
    trajectory: { timepoint: string; off_axis_fraction: number; parallel: number; cosine: number }[];
    maturity: { id: string; drug: string; real_shift: number; sim_shift: number; same_direction: boolean }[];
    n_sim_genes: number; n_null_drugs: number;
  };
  f2: {
    models: Record<string, string | number | null>[];
    verdict: Record<string, { mean_improvement: number; beats_linear: boolean; n_conditions: number }>;
    m3: Record<string, { cosine: number | null }>;
    n_map_training_drugs: number;
  };
  f3: { overall: boolean; verdict: Record<string, { full: number; over_S1: number; over_S2: number }>;
        direction: Record<string, Record<string, number>>; cells: Record<string, { untreated: number; treated: number }> };
  f4: {
    test0: Record<string, number>;
    test1: { table: { prediction: string; cosine: number }[]; beats_both_singles: boolean;
             additive_rank_among_single_drugs: number; best_single_drug: string; best_single_cosine: number };
    test2: { pair: string; score: number; rank: number; n_pairs: number; percentile: number };
    top_pairs: [string, string, number][];
    test3: Record<string, number>;
  };
  validation_label: { headline: string; detail: string; level: string };
}

export interface Hierarchy {
  cell_types: string[];
  sensitivity: Record<string, number[]>;   // drug -> sensitivity per cell type (minus AUC; higher = more sensitive)
  drug_means: Record<string, number>;
  patients: Record<string, { fractions: number[]; ranking: { drug: string; pred_residual: number }[] }>;
}

export interface Exp2Model { model: string; primary_mean: number; ci_lo: number; ci_hi: number; top5_hit_rate: number; pooled_pearson: number }
export interface Exp2 {
  ceiling: { central: number; bias_corrected: number; test_retest: number; variance_model: number; shares: Record<string, number> };
  primary: { mean: number; ci: [number, number]; p_vs_A: number; supported: boolean; fraction_of_ceiling: Record<string, number> };
  models: Exp2Model[];
  secondary: Record<string, { mean_diff: number; p: number }>;
  n_test_patients: number; n_test_samples: number;
  phase4_descriptive: { primary_mean: number; ci_lo: number; ci_hi: number };
  venetoclax_sensitivity: Record<string, number> | null;
}

export function loadTrajectory(file: string, m: Manifest, cells: number): Promise<Trajectory> {
  const key = `traj:${file}`;
  if (!cache.has(key)) {
    cache.set(key, fetch(BASE + file).then((r) => r.arrayBuffer()).then((buf) => {
      const u = new Uint16Array(buf);
      const out = new Float32Array(u.length);
      const [x0, x1] = m.quantization.x, [y0, y1] = m.quantization.y;
      for (let i = 0; i < u.length; i += 3) {
        out[i] = x0 + (u[i] / 65535) * (x1 - x0);
        out[i + 1] = y0 + (u[i + 1] / 65535) * (y1 - y0);
        out[i + 2] = u[i + 2] / 65535;
      }
      return { frames: m.frames, cells, data: out };
    }));
  }
  return cache.get(key) as Promise<Trajectory>;
}

export interface Validation {
  summary: {
    n_samples: number; n_patients: number; n_drugs: number; n_pairs: number;
    P1: { statistic: number; p: number; significant: boolean };
    P2: { statistic: number; p: number; significant: boolean };
    claims: Record<string, boolean>;
  };
  models: { model: string; raw_mean: number; raw_lo: number; raw_hi: number; centred_mean: number; centred_lo: number; centred_hi: number; per_drug_mean: number; per_drug_lo: number; per_drug_hi: number }[];
  per_drug: { drug: string; spearman: number; q: number; tier: string; n: number }[];
  scatter: { s: number; y: number; d: string }[];
  controls: { drug: string; tier: string; percentile: number | null; cosine: number | null }[];
  control_tests: Record<string, { auroc: number; perm_p: number; n_controls: number; n_background: number }>;
  idh: { drug: string; mutation: string; n_mut: number; n_wt: number; median_sens_mut: number; median_sens_wt: number; p_mut_more_sensitive: number }[];
  venetoclax?: { spearman: number; p: number; n: number };
  generic_hist: { bin: number; count: number }[];
  exp2?: Exp2;
}

// ---------------------------------------------------------------- landscape helpers
/** Bilinear sample of a row-major n x n grid at landscape coordinates. */
export function sampleGrid(L: Landscape, field: number[], x: number, y: number): number {
  const n = L.n;
  const gx = ((x - L.x[0]) / (L.x[1] - L.x[0])) * (n - 1);
  const gy = ((y - L.y[0]) / (L.y[1] - L.y[0])) * (n - 1);
  const i0 = Math.max(0, Math.min(n - 2, Math.floor(gx))), j0 = Math.max(0, Math.min(n - 2, Math.floor(gy)));
  const fx = Math.min(1, Math.max(0, gx - i0)), fy = Math.min(1, Math.max(0, gy - j0));
  const a = field[j0 * n + i0], b = field[j0 * n + i0 + 1], c = field[(j0 + 1) * n + i0], d = field[(j0 + 1) * n + i0 + 1];
  return a * (1 - fx) * (1 - fy) + b * fx * (1 - fy) + c * (1 - fx) * fy + d * fx * fy;
}

/** World height of the terrain surface at landscape (x, y); matches the mesh in Terrain.tsx. */
export function terrainY(L: Landscape, x: number, y: number, HEIGHT: number): number {
  return sampleGrid(L, L.height_vis ?? L.height, x, y) * HEIGHT;
}

/** Landscape coordinates -> world (x, z), keeping aspect, spanning about 20 units. */
export function makeWorld(L: Landscape) {
  const sx = L.x[1] - L.x[0], sy = L.y[1] - L.y[0];
  const s = 20 / Math.max(sx, sy);
  const cx = (L.x[0] + L.x[1]) / 2, cy = (L.y[0] + L.y[1]) / 2;
  return {
    W: sx * s, D: sy * s, HEIGHT: 9.0,   // vertical exaggeration, stated in the UI
    toX: (x: number) => (x - cx) * s,
    toZ: (y: number) => -(y - cy) * s,
    // inverse, for turning a pointer hit on the terrain back into landscape coordinates
    fromX: (wx: number) => wx / s + cx,
    fromY: (wz: number) => -wz / s + cy,
  };
}
