import { AnimatePresence, motion } from "framer-motion";

export function Timeline({ frame, frames, playing, speed, onToggle, onScrub, onSpeed, label }: {
  frame: number; frames: number; playing: boolean; speed: number; label: string;
  onToggle: () => void; onScrub: (f: number) => void; onSpeed: () => void;
}) {
  return (
    <div className="timeline" role="group" aria-label="Simulation playback">
      <button className="play" onClick={onToggle} aria-label={playing ? "Pause (space)" : "Play (space)"}>
        {playing
          ? <svg width="14" height="14" viewBox="0 0 14 14"><rect x="2" y="1" width="3.5" height="12" rx="1" fill="currentColor" /><rect x="8.5" y="1" width="3.5" height="12" rx="1" fill="currentColor" /></svg>
          : <svg width="14" height="14" viewBox="0 0 14 14"><path d="M3 1.5v11l9-5.5z" fill="currentColor" /></svg>}
      </button>
      <input className="scrub" type="range" min={0} max={frames - 1} step={0.01} value={frame}
        onChange={(e) => onScrub(parseFloat(e.target.value))} aria-label="Simulation step (left/right arrows)" />
      <div className="tl-meta">
        <span>{label}</span>
        <span>α {((frame / (frames - 1))).toFixed(2)}</span>
        <button className="speed" onClick={onSpeed} aria-label="Playback speed">{speed}×</button>
      </div>
    </div>
  );
}

export function Counters({ matured, resistant, total, mature0 }: { matured: number; resistant: number; total: number; mature0: number }) {
  const pct = (v: number) => (total ? Math.round((v / total) * 100) : 0);
  return (
    <div className="counters" aria-live="polite">
      <AnimatePresence>
        {[
          { k: "matured", v: pct(matured), cls: "m" },
          { k: "resistant", v: pct(resistant), cls: "r" },
          { k: "mature at start", v: pct(mature0), cls: "" },
        ].map((c) => (
          <motion.div key={c.k} className={`counter ${c.cls}`} initial={{ y: 10, opacity: 0 }} animate={{ y: 0, opacity: 1 }}>
            <div className="v">{c.v}%</div>
            <div className="k">{c.k}</div>
          </motion.div>
        ))}
      </AnimatePresence>
    </div>
  );
}
