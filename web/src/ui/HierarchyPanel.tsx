import { useEffect, useMemo, useState } from "react";
import { loadHierarchy, type Hierarchy } from "../data";

// Experiment 2 primary model, applied to this patient's cell-type composition.
// Predicted relative ex vivo response = sum over cell types of fraction x learned sensitivity,
// minus the drug's average. Trained on Beat AML development patients only.

const LABEL: Record<string, string> = {
  "HSC-like": "HSC-like", "Prog-like": "Prog-like", "GMP-like": "GMP-like", "ProMono-like": "ProMono-like",
  "Mono-like": "Mono-like", "cDC-like": "cDC-like", T_NK: "T/NK", B: "B", Ery: "Erythroid", pDC: "pDC",
};
const COLOR: Record<string, string> = {
  "HSC-like": "#3b4cc0", "Prog-like": "#6a8bef", "GMP-like": "#9abbff", "ProMono-like": "#f2a481", "Mono-like": "#d8553f",
  "cDC-like": "#8e44ad", T_NK: "#5a5a72", B: "#7a7a92", Ery: "#9a7a72", pDC: "#6a8a7a",
};

export function HierarchyPanel({ patient }: { patient: string }) {
  const [H, setH] = useState<Hierarchy | null>(null);
  const [err, setErr] = useState(false);
  const [sel, setSel] = useState<string | null>(null);
  useEffect(() => { loadHierarchy().then(setH).catch(() => setErr(true)); }, []);
  const P = H?.patients[patient];
  useEffect(() => { setSel(P?.ranking[0]?.drug ?? null); }, [P]);
  const max = useMemo(() => Math.max(...(P?.ranking.map((r) => Math.abs(r.pred_residual)) ?? [1]), 1e-6), [P]);
  if (err) return <div className="empty">Hierarchy model data not available.</div>;
  if (!H || !P) return <div className="empty">Loading…</div>;
  const sens = sel ? H.sensitivity[sel] : null;
  const sMean = sens ? sens.reduce((a, b) => a + b, 0) / sens.length : 0;
  const sMax = sens ? Math.max(...sens.map((v) => Math.abs(v - sMean)), 1) : 1;
  return (
    <>
      <div style={{ padding: "6px 10px 10px" }}>
        <div className="sub" style={{ marginBottom: 6 }}>Cell-type composition (single-cell labels)</div>
        <div style={{ display: "flex", height: 10, borderRadius: 6, overflow: "hidden", border: "1px solid var(--line)" }} aria-label="Cell-type composition">
          {H.cell_types.map((c, i) => <i key={c} title={`${LABEL[c]} ${(P.fractions[i] * 100).toFixed(1)}%`} style={{ width: `${P.fractions[i] * 100}%`, background: COLOR[c] }} />)}
        </div>
        <div className="tags" style={{ marginTop: 6 }}>
          {H.cell_types.map((c, i) => P.fractions[i] >= 0.03 && (
            <span key={c} className="tag"><span style={{ display: "inline-block", width: 7, height: 7, borderRadius: 2, background: COLOR[c], marginRight: 4 }} />{LABEL[c]} {Math.round(P.fractions[i] * 100)}%</span>
          ))}
        </div>
      </div>
      {sens && (
        <div style={{ padding: "4px 10px 10px", borderTop: "1px solid var(--line)" }}>
          <div className="sub" style={{ margin: "6px 0" }}>{sel}: learned sensitivity by cell type (relative to its average)</div>
          {H.cell_types.map((c, i) => {
            const v = sens[i] - sMean, w = (Math.abs(v) / sMax) * 50;
            return (
              <div key={c} style={{ display: "grid", gridTemplateColumns: "78px 1fr", alignItems: "center", gap: 6, fontSize: 11, color: "var(--text-2)" }}>
                <span>{LABEL[c]}</span>
                <span style={{ position: "relative", height: 6, background: "rgba(255,255,255,.05)", borderRadius: 4 }}>
                  <i style={{ position: "absolute", top: 0, bottom: 0, left: v >= 0 ? "50%" : `${50 - w}%`, width: `${w}%`, borderRadius: 4, background: v >= 0 ? "#ffc46b" : "#6f62e8" }} />
                </span>
              </div>
            );
          })}
        </div>
      )}
      <div className="scroll" style={{ borderTop: "1px solid var(--line)" }}>
        {P.ranking.slice(0, 30).map((r, i) => {
          const w = Math.min(50, (Math.abs(r.pred_residual) / max) * 50);
          return (
            <button key={r.drug} className="drug" aria-pressed={r.drug === sel} onClick={() => setSel(r.drug)}>
              <span className="rk">{i + 1}</span>
              <span className="dn">{r.drug}</span>
              <span className="sc">{r.pred_residual >= 0 ? "+" : ""}{r.pred_residual.toFixed(1)}</span>
              <span className="bar"><i style={r.pred_residual >= 0
                ? { left: "50%", width: `${w}%`, background: "linear-gradient(90deg, rgba(143,240,176,.4), #8ff0b0)" }
                : { left: `${50 - w}%`, width: `${w}%`, background: "#6f62e8" }} /></span>
            </button>
          );
        })}
      </div>
      <div className="hint">
        Experiment 2 primary model: predicted ex vivo response relative to each drug's average (minus AUC units; higher = more sensitive than usual).
        On held-out Beat AML patients it ranks drugs modestly better than chance (Spearman 0.20, about 30% of the estimated ceiling). It was trained on bulk
        deconvolution and is applied here to single-cell counts, so treat it as a ranking hypothesis, not a prediction of clinical benefit.
      </div>
    </>
  );
}
