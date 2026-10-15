import { motion } from "framer-motion";
import type { PatientDetail, PatientSummary, RankedDrug } from "../data";

const HOT = /^(IDH1|IDH2|FLT3|NPM1|TP53)/;

export function PatientPanel({ patients, selected, onSelect }: { patients: PatientSummary[]; selected: string | null; onSelect: (id: string) => void }) {
  return (
    <motion.aside className="panel left" initial={{ x: -30, opacity: 0 }} animate={{ x: 0, opacity: 1 }} transition={{ duration: 0.6, delay: 0.3 }} aria-label="Patients">
      <header>
        <h2>Patients</h2>
        <div className="sub">van Galen et al. 2019, malignant cells at diagnosis</div>
      </header>
      <div className="scroll" role="list">
        {patients.map((p) => (
          <button key={p.id} className="patient" aria-pressed={p.id === selected} onClick={() => onSelect(p.id)} role="listitem">
            <div className="row"><span className="name">{p.id}</span><span className="n">{p.n_malignant_d0.toLocaleString()} cells</span></div>
            <div className="tags">
              {p.mutations.length === 0 && <span className="tag">no genotyped mutation detected</span>}
              {p.mutations.slice(0, 5).map((m) => <span key={m} className={`tag ${HOT.test(m) ? "hot" : ""}`}>{m}</span>)}
              {p.mutations.length > 5 && <span className="tag">+{p.mutations.length - 5}</span>}
            </div>
            <div className="stuck-label"><span>stuck fraction</span><b>{Math.round(p.stuck_fraction * 100)}%</b></div>
            <div className="stuckbar"><i style={{ width: `${p.stuck_fraction * 100}%` }} /></div>
          </button>
        ))}
      </div>
      <div className="hint">Stuck = score below the maturity threshold learned on normal marrow. Mutation tags: variants seen in the authors' single-cell genotyping.</div>
    </motion.aside>
  );
}

function scoreBar(v: number, max: number) {
  const w = Math.min(50, (Math.abs(v) / max) * 50);
  return v >= 0
    ? { left: "50%", width: `${w}%`, background: "linear-gradient(90deg, rgba(255,196,107,.5), #ffc46b)" }
    : { left: `${50 - w}%`, width: `${w}%`, background: "linear-gradient(90deg, #6f62e8, rgba(111,98,232,.5))" };
}

export function DrugPanel({ patient, active, compare, compareMode, onPick }: {
  patient: PatientDetail | null; active: string | null; compare: string | null; compareMode: boolean; onPick: (d: string) => void;
}) {
  if (!patient) return <aside className="panel right"><div className="empty">Select a patient.</div></aside>;
  const simulated = new Set(patient.simulated_drugs.map((d) => d.drug));
  const list: RankedDrug[] = [...patient.ranking];
  for (const s of patient.simulated_drugs) if (!list.some((r) => r.drug === s.drug)) {
    list.push({ drug: s.drug, rank: NaN, cosine: NaN, connectivity: NaN, tier: "AML/monocytic", cells: "", n_sigs: 0, in_beataml_panel: false });
  }
  const max = Math.max(...list.filter((d) => isFinite(d.cosine)).map((d) => Math.abs(d.cosine)), 0.01);
  return (
    <motion.aside className="panel right" initial={{ x: 30, opacity: 0 }} animate={{ x: 0, opacity: 1 }} transition={{ duration: 0.6, delay: 0.4 }} aria-label="Ranked drugs">
      <header>
        <h2>Ranked drugs</h2>
        <div className="sub">{compareMode ? "Compare: pick two simulated drugs" : `Reversion score for ${patient.id}. Highlighted drugs can be simulated.`}</div>
      </header>
      <div className="scroll">
        {list.map((d) => {
          const can = simulated.has(d.drug);
          const pressed = d.drug === active || d.drug === compare;
          return (
            <button key={d.drug} className={`drug ${d.drug === compare ? "compare" : ""}`} aria-pressed={pressed} disabled={!can}
              style={{ opacity: can ? 1 : 0.62, cursor: can ? "pointer" : "default" }} onClick={() => can && onPick(d.drug)}
              title={can ? "Play simulation" : "Not simulated (only each patient's top drugs and tretinoin are exported)"}>
              <span className="rk">{isFinite(d.rank) ? d.rank : "ref"}</span>
              <span className="dn">{d.drug}</span>
              <span className="sc">{isFinite(d.cosine) ? d.cosine.toFixed(3) : "control"}</span>
              <span className="bar">{isFinite(d.cosine) && <i style={scoreBar(d.cosine, max)} />}</span>
              <span className="meta">
                {can && <span className="badge sim">simulated</span>}
                <span>{d.cells ? d.cells.replace(/,/g, " · ") : d.tier}</span>
                {d.in_beataml_panel && (
                  <span className="badge beat" title="Drug is in the Beat AML ex vivo panel: cohort median AUC and the across-patient correlation between reversion score and sensitivity (Phase 4)">
                    Beat AML{d.beat_median_auc != null ? ` AUC ${Math.round(d.beat_median_auc)}` : ""}{d.validation_rho != null ? ` · ρ ${d.validation_rho.toFixed(2)}` : ""}
                  </span>
                )}
              </span>
            </button>
          );
        })}
      </div>
      <div className="hint">Score = cosine between the drug's LINCS signature and this patient's push-to-maturity vector. Lower Beat AML AUC = more sensitive ex vivo (cohort data; these single-cell patients were not in Beat AML).</div>
    </motion.aside>
  );
}
