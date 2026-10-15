import { motion } from "framer-motion";
import { useEffect, useState } from "react";
import {
  Bar, BarChart, CartesianGrid, Cell, LabelList, Legend, ReferenceLine,
  ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import { loadExp3, type Exp3 } from "../data";

// Chart colours are the validated categorical steps declared in index.css.
const S1 = "#3987e5", S3 = "#199e70", S4 = "#c98500";
const PASS = "#199e70", FAIL = "#d55181", MUTED = "rgba(190,180,255,.22)";
const f3 = (x: unknown) => (typeof x === "number" ? x.toFixed(3) : String(x));
const tipStyle = {
  contentStyle: { background: "#120f26", border: "1px solid rgba(190,180,255,.24)", borderRadius: 8, fontSize: 12 },
  itemStyle: { color: "#ece9fb" }, labelStyle: { color: "#b2acd3" },
};

export function Exp3View() {
  const [e, setE] = useState<Exp3 | null>(null);
  const [err, setErr] = useState<string | null>(null);
  useEffect(() => { loadExp3().then(setE).catch((x) => setErr(String(x))); }, []);
  if (err) return <div className="validation"><div className="empty">Experiment 3 data not found ({err}).</div></div>;
  if (!e) return <div className="validation"><div className="empty">Loading…</div></div>;

  const dir = [...e.f1.direction].sort((a, b) => b.cosine - a.cosine)
    .map((d) => ({ ...d, label: d.id.replace(/^D\d /, ""), captured: d.variance_explained * 100,
                   rankLabel: d.passes ? `#${d.rank} of ${d.n_null.toLocaleString()}` : "" }));
  // Share of the measured response captured = cosine squared (variance explained), the same
  // quantity as the headline figure. Not 1 - off_axis_fraction, which is a different number.
  const traj = e.f1.trajectory.map((t) => ({
    timepoint: t.timepoint,
    along: t.cosine * t.cosine * 100,
  }));
  const f2rows = e.f2.models.map((m) => ({
    id: String(m.id ?? "").replace(/^D\d /, ""),
    linear: m["B1 linear shift"] as number,
    matched: m["M1 matched cell line"] as number | null,
    learned: m["M2 learned map"] as number | null,
  }));
  const combo = [
    { rule: "ATRA alone", v: e.f4.test1.table.find((t) => t.prediction === "ATRA only")?.cosine ?? 0, single: true },
    { rule: "palbociclib alone", v: e.f4.test1.table.find((t) => t.prediction === "palbociclib only")?.cosine ?? 0, single: true },
    { rule: "both, combined", v: e.f4.test1.table.find((t) => t.prediction === "additive")?.cosine ?? 0, single: false },
  ];
  const m3 = e.f2.m3;

  return (
    <motion.div className="validation" initial={{ opacity: 0, y: 20 }} animate={{ opacity: 1, y: 0 }}
      exit={{ opacity: 0, y: 20 }} transition={{ duration: 0.45 }}>
      <div className="v-grid">
        <div className="card span-12">
          <h3>Experiment 3: is the simulation real?</h3>
          <div className="sub">
            Everything else on this site assumes the animation means something. Here we test it against measured
            transcriptomes of AML cells actually treated with these drugs.
          </div>
          <div className="kpi">
            <div><div className="v">{e.f1.n_pass}/{e.f1.n_total}</div><div className="k">direction tests passed</div></div>
            <div><div className="v">{Math.round(Math.min(...dir.filter((d) => d.passes).map((d) => d.captured)))}–{Math.round(Math.max(...dir.map((d) => d.captured)))}%</div>
              <div className="k">of real movement captured</div></div>
            <div><div className="v">{Math.round(e.f1.trajectory[2].off_axis_fraction * 100)}%</div><div className="k">of the ATRA response is off our axis</div></div>
            <div><div className="v">0</div><div className="k">patient datasets available</div></div>
          </div>
          <div className="verdict">
            <b>{e.validation_label.headline}.</b> {e.validation_label.detail}
          </div>
        </div>

        <div className="card span-7">
          <h3>Does the predicted direction match reality?</h3>
          <div className="sub">
            Cosine between the simulated and the measured change, per experiment. Each drug was also ranked against
            all {e.f1.n_null_drugs.toLocaleString()} other compounds scored the same way; that rank is printed on the bar.
          </div>
          <ResponsiveContainer width="100%" height={260}>
            <BarChart data={dir} layout="vertical" margin={{ left: 10, right: 96 }}>
              <CartesianGrid horizontal={false} />
              <XAxis type="number" domain={[-0.1, 0.7]} tickFormatter={f3} />
              <YAxis type="category" dataKey="label" width={118} interval={0} tick={{ fontSize: 11 }} />
              <ReferenceLine x={0} stroke="#8a84ae" />
              <Tooltip {...tipStyle} formatter={f3} />
              <Bar dataKey="cosine" radius={[0, 4, 4, 0]} barSize={16}>
                {dir.map((d) => <Cell key={d.id} fill={d.passes ? PASS : FAIL} />)}
                <LabelList dataKey="rankLabel" position="right" fontSize={10} fill="#b2acd3" />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
          <div className="legend-row">
            <span><i style={{ background: PASS }} />beat chance</span>
            <span><i style={{ background: FAIL }} />failed (4 h, too early)</span>
          </div>
        </div>

        <div className="card span-5">
          <h3>But how much does it capture?</h3>
          <div className="sub">
            The honest number. Of everything ATRA really does to these cells, this is the share lying along our
            straight-line prediction, at each measured timepoint.
          </div>
          <div className="capture-rows">
            {traj.map((t) => (
              <div key={t.timepoint}>
                <span>{t.timepoint} after ATRA</span>
                <div className="meter wide"><i style={{ width: `${t.along}%`, background: S1 }} /></div>
                <b>{t.along.toFixed(0)}%</b>
              </div>
            ))}
          </div>
          <div className="verdict" style={{ marginTop: 16 }}>
            Even at its best, {(100 - traj[2].along).toFixed(0)}% of what the cells actually do at 4 days is movement
            our model does not represent. At 4 hours the prediction points the wrong way entirely.
          </div>
        </div>

        <div className="card span-7">
          <h3>Can anything do better?</h3>
          <div className="sub">
            Two alternatives tested against the same measurements. Using the drug's signature from the matching cell
            line helps; a map learned from {e.f2.n_map_training_drugs} other drugs does not. Gaps show where a model could not be applied.
          </div>
          <ResponsiveContainer width="100%" height={250}>
            <BarChart data={f2rows} margin={{ left: 0, right: 10, bottom: 4 }}>
              <CartesianGrid vertical={false} />
              <XAxis dataKey="id" tick={{ fontSize: 10 }} interval={0} />
              <YAxis domain={[-0.1, 0.7]} tickFormatter={f3} />
              <ReferenceLine y={0} stroke="#8a84ae" />
              <Tooltip {...tipStyle} formatter={f3} />
              <Legend wrapperStyle={{ fontSize: 11 }} />
              <Bar dataKey="linear" name="our linear shift" fill={S1} barSize={14} radius={[3, 3, 0, 0]} />
              <Bar dataKey="matched" name="matched cell line" fill={S3} barSize={14} radius={[3, 3, 0, 0]} />
              <Bar dataKey="learned" name="learned map" fill={S4} barSize={14} radius={[3, 3, 0, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="card span-5">
          <h3>The missing ingredient is time</h3>
          <div className="sub">
            Our model predicts the same change at 4 hours and 4 days. Interpolating between two real timepoints
            predicts the held-out middle one far better, on identical genes.
          </div>
          <div className="compare-rows">
            <div>
              <span>our static prediction</span>
              <div className="meter"><i style={{ width: `${(m3["B1 linear shift"].cosine ?? 0) * 100}%`, background: S1 }} /></div>
              <b>{(m3["B1 linear shift"].cosine ?? 0).toFixed(2)}</b>
            </div>
            <div>
              <span>knowing the time course</span>
              <div className="meter"><i style={{ width: `${(m3["M3 time interpolation"].cosine ?? 0) * 100}%`, background: S3 }} /></div>
              <b>{(m3["M3 time interpolation"].cosine ?? 0).toFixed(2)}</b>
            </div>
          </div>
          <div className="verdict" style={{ marginTop: 14 }}>
            This is not a fair contest: interpolation gets to see real measurements of the same drug. It measures how
            much a time axis would be worth, and the answer is <b>a great deal</b>.
          </div>
        </div>

        <div className="card span-6">
          <h3>Two drugs at once</h3>
          <div className="sub">
            One real combination exists in open data: ATRA plus palbociclib in HL-60. Predicting from both drugs beats
            predicting from either alone.
          </div>
          <ResponsiveContainer width="100%" height={190}>
            <BarChart data={combo} layout="vertical" margin={{ left: 10, right: 30 }}>
              <CartesianGrid horizontal={false} />
              <XAxis type="number" domain={[0, 0.6]} tickFormatter={f3} />
              <YAxis type="category" dataKey="rule" width={128} tick={{ fontSize: 11 }} />
              <Tooltip {...tipStyle} formatter={f3} />
              <Bar dataKey="v" radius={[0, 4, 4, 0]} barSize={18}>
                {combo.map((c) => <Cell key={c.rule} fill={c.single ? MUTED : S3} />)}
                <LabelList dataKey="v" position="right" fontSize={10} fill="#b2acd3" formatter={f3} />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
          <div className="verdict">
            One combination, one cell line, one timepoint. And an unrelated single drug still predicts it better than
            we do, so this is a weak result, not a validated capability.
          </div>
        </div>

        <div className="card span-6">
          <h3>What we could not test</h3>
          <div className="sub">Limits of the open data, not choices.</div>
          <ul className="limits">
            <li><b>Patient cells.</b> No public dataset treats patient AML cells with a drug whose signature we hold. Every
              result here is a cell line.</li>
            <li><b>Which cells resist.</b> Sequencing destroys the cell, so no cell is seen before and after. In the one
              usable dataset, treatment and sequencing run are perfectly confounded.</li>
            <li><b>Synergy.</b> Our model adds two drug effects together, so it cannot express the interaction that makes
              combinations worth studying.</li>
            <li><b>Menin inhibitors, IDH inhibitors.</b> Absent from the signature library, so the drugs with the best
              clinical differentiation evidence cannot be simulated at all.</li>
          </ul>
        </div>
      </div>
    </motion.div>
  );
}
