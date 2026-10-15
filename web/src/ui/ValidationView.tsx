import { motion } from "framer-motion";
import { useEffect, useState } from "react";
import {
  Bar, BarChart, CartesianGrid, Cell, ErrorBar, ReferenceLine, ResponsiveContainer, Scatter, ScatterChart, Tooltip, XAxis, YAxis, ZAxis,
} from "recharts";
import { loadValidation, type Validation } from "../data";

const f3 = (x: unknown) => (typeof x === "number" ? x.toFixed(3) : String(x));
const pct = (x: unknown) => (typeof x === "number" ? `${x.toFixed(1)} percentile` : String(x));
const fmtP = (p: number) => (p < 0.001 ? p.toExponential(1) : p.toFixed(3));

function short(m: string) {
  return m.replace(" (raw)", "").replace(" (calibrated)", " cal.").replace("reversion score", "reversion").replace("differentiation fractions", "diff. fractions").replace("global drug ranking", "global ranking");
}

export function ValidationView() {
  const [v, setV] = useState<Validation | null>(null);
  const [err, setErr] = useState<string | null>(null);
  useEffect(() => { loadValidation().then(setV).catch((e) => setErr(String(e))); }, []);
  if (err) return <div className="validation"><div className="empty">Validation data not found ({err}).</div></div>;
  if (!v) return <div className="validation"><div className="empty">Loading validation…</div></div>;
  const s = v.summary;
  const models = v.models.map((m) => ({ name: short(m.model), raw: m.raw_mean, rawErr: [m.raw_mean - m.raw_lo, m.raw_hi - m.raw_mean],
    cen: m.centred_mean, cenErr: [m.centred_mean - m.centred_lo, m.centred_hi - m.centred_mean], ours: m.model.startsWith("M") }));
  const pd = [...v.per_drug].sort((a, b) => b.spearman - a.spearman);
  const t1 = v.control_tests["tier1_cosine"];
  const tip = { contentStyle: { background: "#120f26", border: "1px solid rgba(190,180,255,.24)", borderRadius: 8, fontSize: 12 }, itemStyle: { color: "#ece9fb" }, labelStyle: { color: "#b2acd3" } };
  return (
    <motion.div className="validation" initial={{ opacity: 0, y: 24 }} animate={{ opacity: 1, y: 0 }} exit={{ opacity: 0, y: 24 }} transition={{ duration: 0.5 }}>
      <div className="v-grid">
        <div className="card span-12">
          <h3>Does the reversion score predict real patients' drug response?</h3>
          <div className="sub">Beat AML ex vivo sensitivity, pre-registered analysis, patient-level 5-fold cross-validation.</div>
          <div className="kpi">
            <div><div className="v">{s.n_samples}</div><div className="k">samples</div></div>
            <div><div className="v">{s.n_drugs}</div><div className="k">drugs</div></div>
            <div><div className="v">{s.n_pairs.toLocaleString()}</div><div className="k">patient-drug pairs</div></div>
            <div><div className="v">{s.P1.statistic.toFixed(3)}</div><div className="k">P1 mean Spearman (p {fmtP(s.P1.p)})</div></div>
            <div><div className="v">{s.P2.statistic.toFixed(3)}</div><div className="k">P2 patient-specific (p {fmtP(s.P2.p)})</div></div>
            <div><div className="v">{t1 ? t1.auroc.toFixed(2) : "–"}</div><div className="k">positive-control AUROC</div></div>
          </div>
          <div className="verdict">
            <b>Pre-registered verdict.</b> Information about response: <b>{s.claims.information_about_response ? "yes" : "no"}</b>. Patient-specific value beyond the global drug ranking: <b>{s.claims.patient_specific_value ? "yes" : "no"}</b>.
            Ex vivo assays measure 3-day viability, not maturation, so this is a proxy test, and the panel contains no clinical differentiation agents.
          </div>
        </div>

        <div className="card span-6">
          <h3>Ranking drugs within a patient</h3>
          <div className="sub">Mean per-patient Spearman (raw sensitivity), 95% bootstrap CI. Global ranking is the bar to beat.</div>
          <ResponsiveContainer width="100%" height={280}>
            <BarChart data={models} layout="vertical" margin={{ left: 30, right: 20 }}>
              <CartesianGrid horizontal={false} />
              <XAxis type="number" domain={["auto", "auto"]} tickFormatter={(x) => x.toFixed(2)} />
              <YAxis type="category" dataKey="name" width={150} />
              <ReferenceLine x={0} stroke="#8a84ae" />
              <Tooltip {...tip} formatter={f3} />
              <Bar dataKey="raw" radius={[0, 4, 4, 0]}>
                {models.map((m) => <Cell key={m.name} fill={m.ours ? "#ffc46b" : "#6f62e8"} />)}
                <ErrorBar dataKey="rawErr" width={4} stroke="#ece9fb" direction="x" />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="card span-6">
          <h3>Patient-specific signal</h3>
          <div className="sub">Same, after removing each drug's average potency. Global ranking scores 0 here by construction.</div>
          <ResponsiveContainer width="100%" height={280}>
            <BarChart data={models} layout="vertical" margin={{ left: 30, right: 20 }}>
              <CartesianGrid horizontal={false} />
              <XAxis type="number" domain={["auto", "auto"]} tickFormatter={(x) => x.toFixed(2)} />
              <YAxis type="category" dataKey="name" width={150} />
              <ReferenceLine x={0} stroke="#8a84ae" />
              <Tooltip {...tip} formatter={f3} />
              <Bar dataKey="cen" radius={[0, 4, 4, 0]}>
                {models.map((m) => <Cell key={m.name} fill={m.ours ? "#ffc46b" : "#6f62e8"} />)}
                <ErrorBar dataKey="cenErr" width={4} stroke="#ece9fb" direction="x" />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="card span-7">
          <h3>Predicted vs observed</h3>
          <div className="sub">Reversion score vs ex vivo sensitivity (−AUC) for a random sample of patient-drug pairs.</div>
          <ResponsiveContainer width="100%" height={300}>
            <ScatterChart margin={{ left: 10, right: 20, bottom: 10 }}>
              <CartesianGrid />
              <XAxis type="number" dataKey="s" name="reversion score" tickFormatter={(x) => x.toFixed(2)} label={{ value: "reversion score", position: "insideBottom", offset: -4, fill: "#8a84ae", fontSize: 11 }} />
              <YAxis type="number" dataKey="y" name="sensitivity" label={{ value: "sensitivity (−AUC)", angle: -90, position: "insideLeft", fill: "#8a84ae", fontSize: 11 }} />
              <ZAxis range={[10, 10]} />
              <Tooltip {...tip} formatter={f3} />
              <Scatter data={v.scatter} fill="#6f62e8" fillOpacity={0.35} />
            </ScatterChart>
          </ResponsiveContainer>
        </div>

        <div className="card span-5">
          <h3>Positive controls</h3>
          <div className="sub">Generic maturation score of known differentiation agents among all compounds profiled in AML/monocytic lines (percentile).</div>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={v.controls.filter((c) => c.percentile != null)} layout="vertical" margin={{ left: 20, right: 20 }}>
              <CartesianGrid horizontal={false} />
              <XAxis type="number" domain={[0, 100]} />
              <YAxis type="category" dataKey="drug" width={120} interval={0} tick={{ fontSize: 10 }} />
              <ReferenceLine x={50} stroke="#8a84ae" strokeDasharray="3 3" />
              <Tooltip {...tip} formatter={pct} />
              <Bar dataKey="percentile" radius={[0, 4, 4, 0]}>
                {v.controls.filter((c) => c.percentile != null).map((c) => <Cell key={c.drug} fill={c.tier === "tier1" ? "#ffc46b" : "#7fe7ff"} />)}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="card span-12">
          <h3>Per drug, across patients</h3>
          <div className="sub">Spearman between reversion score and sensitivity for each drug. Most panel drugs kill cells rather than mature them, so weak correlations are expected.</div>
          <ResponsiveContainer width="100%" height={260}>
            <BarChart data={pd} margin={{ left: 10, right: 10, bottom: 70 }}>
              <CartesianGrid vertical={false} />
              <XAxis dataKey="drug" interval={0} angle={-70} textAnchor="end" height={80} tick={{ fontSize: 8 }} />
              <YAxis tickFormatter={(x) => x.toFixed(2)} />
              <ReferenceLine y={0} stroke="#8a84ae" />
              <Tooltip {...tip} formatter={f3} />
              <Bar dataKey="spearman">
                {pd.map((d) => <Cell key={d.drug} fill={d.q < 0.05 ? (d.spearman > 0 ? "#ffc46b" : "#ff7a90") : "rgba(111,98,232,.55)"} />)}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="card span-6">
          <h3>IDH subgroup check</h3>
          <div className="sub">Tool IDH inhibitors: are IDH-mutant samples more sensitive ex vivo?</div>
          <table style={{ width: "100%", fontSize: 13, borderCollapse: "collapse" }}>
            <thead><tr style={{ color: "#8a84ae", textAlign: "left" }}><th>drug</th><th>mutation</th><th>n mut / wt</th><th>median sens. mut / wt</th><th>p</th></tr></thead>
            <tbody>{v.idh.map((r) => (
              <tr key={r.drug} style={{ borderTop: "1px solid rgba(190,180,255,.12)" }}>
                <td>{r.drug}</td><td>{r.mutation}</td><td>{r.n_mut} / {r.n_wt}</td>
                <td>{r.median_sens_mut.toFixed(1)} / {r.median_sens_wt.toFixed(1)}</td><td>{fmtP(r.p_mut_more_sensitive)}</td>
              </tr>))}</tbody>
          </table>
        </div>
        {v.venetoclax && (
          <div className="card span-6">
            <h3>Sanity check: venetoclax</h3>
            <div className="sub">Monocytic AML is known to resist venetoclax. Our deconvolved mature share should correlate negatively with sensitivity.</div>
            <div className="kpi"><div><div className="v">{v.venetoclax.spearman.toFixed(2)}</div><div className="k">Spearman (n = {v.venetoclax.n}, p {fmtP(v.venetoclax.p)})</div></div></div>
          </div>
        )}
      </div>
    </motion.div>
  );
}
