import { AnimatePresence, motion } from "framer-motion";
import { lazy, Suspense, useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  loadLandscape, loadManifest, loadPatient, loadTrajectory,
  type Landscape, type Manifest, type PatientDetail, type Trajectory,
} from "./data";
import { FlatView } from "./scene/FlatView";
import { Scene } from "./scene/Scene";
import { HowItWorks } from "./ui/HowItWorks";
import { DrugPanel, PatientPanel } from "./ui/Panels";
import { Counters, Timeline } from "./ui/Timeline";
const ValidationView = lazy(() => import("./ui/ValidationView").then((m) => ({ default: m.ValidationView })));

type View = "explore" | "compare" | "validation";
const SPEEDS = [0.5, 1, 2, 4];
const BASE_FPS = 6; // simulation frames per second at 1x

function webglOK() {
  try { const c = document.createElement("canvas"); return !!(c.getContext("webgl2") || c.getContext("webgl")); } catch { return false; }
}

export default function App() {
  const [manifest, setManifest] = useState<Manifest | null>(null);
  const [L, setL] = useState<(Landscape & { support?: number[] }) | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [view, setView] = useState<View>("explore");
  const [how, setHow] = useState(false);
  const [pid, setPid] = useState<string | null>(null);
  const [patient, setPatient] = useState<PatientDetail | null>(null);
  const [drugA, setDrugA] = useState<string | null>(null);
  const [drugB, setDrugB] = useState<string | null>(null);
  const [trajA, setTrajA] = useState<Trajectory | null>(null);
  const [trajB, setTrajB] = useState<Trajectory | null>(null);
  const [frame, setFrame] = useState(0);
  const [playing, setPlaying] = useState(false);
  const [speed, setSpeed] = useState(1);
  const prefersReduced = useMemo(() => window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ?? false, []);
  const [reduced, setReduced] = useState(prefersReduced);
  const lowPower = useMemo(() => !webglOK() || window.innerWidth < 760 || (navigator.hardwareConcurrency ?? 8) <= 2, []);
  const [flat, setFlat] = useState(lowPower);
  const [effects, setEffects] = useState(!lowPower);
  const [intro, setIntro] = useState(!prefersReduced);
  const frameRef = useRef(0);

  useEffect(() => {
    Promise.all([loadManifest(), loadLandscape()]).then(([m, l]) => { setManifest(m); setL(l); setPid(m.patients[0]?.id ?? null); })
      .catch((e) => setError(String(e)));
  }, []);

  useEffect(() => {
    if (!manifest || !pid) return;
    const s = manifest.patients.find((p) => p.id === pid)!;
    loadPatient(s.file).then((p) => {
      setPatient(p);
      const first = p.simulated_drugs[0]?.drug ?? null;
      setDrugA(first);
      setDrugB(p.simulated_drugs.find((d) => d.drug === "tretinoin" && d.drug !== first)?.drug ?? p.simulated_drugs[1]?.drug ?? null);
      frameRef.current = 0; setFrame(0);
    });
  }, [manifest, pid]);

  const simA = patient?.simulated_drugs.find((d) => d.drug === drugA) ?? null;
  const simB = patient?.simulated_drugs.find((d) => d.drug === drugB) ?? null;
  useEffect(() => { if (manifest && patient && simA) loadTrajectory(simA.file, manifest, patient.n_shown).then(setTrajA); else setTrajA(null); }, [manifest, patient, simA]);
  useEffect(() => { if (manifest && patient && simB) loadTrajectory(simB.file, manifest, patient.n_shown).then(setTrajB); else setTrajB(null); }, [manifest, patient, simB]);

  // playback loop
  useEffect(() => {
    if (!playing || !manifest) return;
    let raf = 0, last = performance.now();
    const tick = (t: number) => {
      const dt = (t - last) / 1000; last = t;
      let f = frameRef.current + dt * BASE_FPS * speed;
      if (f >= manifest.frames - 1) { f = manifest.frames - 1; setPlaying(false); }
      frameRef.current = f; setFrame(f);
      raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(raf);
  }, [playing, speed, manifest]);

  const scrub = useCallback((f: number) => { frameRef.current = f; setFrame(f); }, []);
  const toggle = useCallback(() => {
    if (!manifest) return;
    if (!playing && frameRef.current >= manifest.frames - 1) scrub(0);
    setPlaying((p) => !p);
  }, [playing, manifest, scrub]);
  const pickDrug = useCallback((d: string) => {
    if (view === "compare") {
      if (d === drugA) return;
      setDrugB(d);
    } else setDrugA(d);
    scrub(0); setPlaying(true);
  }, [view, drugA, scrub]);

  // keyboard
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.target as HTMLElement)?.tagName === "INPUT" && e.key !== " ") return;
      if (!manifest) return;
      if (e.key === " ") { e.preventDefault(); toggle(); }
      else if (e.key === "ArrowRight") scrub(Math.min(manifest.frames - 1, Math.round(frameRef.current) + 1));
      else if (e.key === "ArrowLeft") scrub(Math.max(0, Math.round(frameRef.current) - 1));
      else if (e.key === "]") setSpeed((s) => SPEEDS[Math.min(SPEEDS.length - 1, SPEEDS.indexOf(s) + 1)]);
      else if (e.key === "[") setSpeed((s) => SPEEDS[Math.max(0, SPEEDS.indexOf(s) - 1)]);
      else if (e.key === "Escape") setHow(false);
      else if (e.key.toLowerCase() === "h") setHow((h) => !h);
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [manifest, toggle, scrub]);

  // start playing once the intro fly-over ends
  const onIntroDone = useCallback(() => { setIntro(false); setTimeout(() => setPlaying(true), 400); }, []);
  useEffect(() => { if (!intro && reduced) setPlaying(false); }, [intro, reduced]);

  const support = useMemo(() => new Float32Array(L?.support ?? (L ? L.height.map(() => 1) : [])), [L]);
  const counters = (traj: Trajectory | null, status: number[] | undefined) => {
    if (!traj || !L || !status) return { matured: 0, resistant: 0, total: 0, mature0: 0 };
    const f = Math.min(traj.frames - 1, Math.round(frame)), n = traj.cells, D = traj.data;
    let matured = 0, resistant = 0, mature0 = 0;
    for (let c = 0; c < n; c++) {
      const s0 = D[c * 3 + 2], s = D[(f * n + c) * 3 + 2];
      if (s0 >= L.threshold) mature0++;
      else if (s >= L.threshold) matured++;
      if (status[c] === 2) resistant++;
    }
    return { matured, resistant, total: n, mature0 };
  };

  if (error) return <div className="loading"><div><p>Could not load data.<br /><small style={{ fontSize: 13 }}>{error}</small></p></div></div>;
  const ready = manifest && L && patient;

  const renderView = (traj: Trajectory | null, status: number[] | null, tag?: string, showIntro = false) => (
    <div className="view" style={{ position: view === "compare" ? "relative" : "absolute", inset: 0 }}>
      {tag && <div className="view-tag">{tag}</div>}
      {flat
        ? <FlatView L={L!} traj={traj} status={status} frameRef={frameRef} />
        : <Scene L={L!} support={support} traj={traj} status={status} frameRef={frameRef} reduced={reduced} effects={effects}
            intro={showIntro && intro} onIntroDone={onIntroDone} showLabels={view !== "compare"} />}
    </div>
  );

  const cA = counters(trajA, simA?.status);
  return (
    <div className="app">
      <AnimatePresence>{!ready && (
        <motion.div className="loading" exit={{ opacity: 0 }}><div style={{ display: "grid", justifyItems: "center" }}><div className="ring" /><p>Growing the landscape…</p></div></motion.div>
      )}</AnimatePresence>

      {ready && view !== "validation" && (
        <div className={`stage ${view === "compare" ? "split" : ""}`}>
          {view === "compare"
            ? <>{renderView(trajA, simA?.status ?? null, drugA ?? "")}{renderView(trajB, simB?.status ?? null, drugB ?? "")}</>
            : renderView(trajA, simA?.status ?? null, undefined, true)}
        </div>
      )}

      <header className="topbar">
        <div className="brand"><h1>Reversion <em>Landscape</em></h1><span>pushing leukaemic cells to mature, in silico</span></div>
        <nav className="nav" aria-label="Views">
          {(["explore", "compare", "validation"] as View[]).map((v) => (
            <button key={v} aria-pressed={view === v} onClick={() => { setView(v); if (v === "compare") setPlaying(false); }}>{v[0].toUpperCase() + v.slice(1)}</button>
          ))}
          <button onClick={() => setHow(true)} aria-label="How it works (H)">How it works</button>
        </nav>
        <div style={{ display: "flex", gap: 14 }}>
          <label className="toggle"><input type="checkbox" checked={reduced} onChange={(e) => setReduced(e.target.checked)} /> reduced motion</label>
          <label className="toggle"><input type="checkbox" checked={flat} onChange={(e) => setFlat(e.target.checked)} /> 2D</label>
          {!flat && <label className="toggle"><input type="checkbox" checked={effects} onChange={(e) => setEffects(e.target.checked)} /> effects</label>}
        </div>
      </header>

      {ready && view !== "validation" && (
        <>
          <PatientPanel patients={manifest!.patients} selected={pid} onSelect={(id) => { setPlaying(false); setPid(id); }} />
          <DrugPanel patient={patient} active={drugA} compare={view === "compare" ? drugB : null} compareMode={view === "compare"} onPick={pickDrug} />
          {view === "explore" && <Counters {...cA} />}
          <Timeline frame={frame} frames={manifest!.frames} playing={playing} speed={speed} onToggle={toggle} onScrub={scrub}
            label={view === "compare" ? `${drugA} vs ${drugB}` : drugA ?? ""}
            onSpeed={() => setSpeed((s) => SPEEDS[(SPEEDS.indexOf(s) + 1) % SPEEDS.length])} />
        </>
      )}

      <Suspense fallback={<div className="validation"><div className="empty">Loading validation…</div></div>}>
        <AnimatePresence>{view === "validation" && <ValidationView key="val" />}</AnimatePresence>
      </Suspense>
      <AnimatePresence>{how && <HowItWorks onClose={() => setHow(false)} />}</AnimatePresence>

      <footer className="disclaimer" role="note">
        <b>Research simulation</b>
        <span>In-silico hypotheses from public data and a linear expression-shift model. Not medical advice. No clinical use.</span>
      </footer>
    </div>
  );
}
