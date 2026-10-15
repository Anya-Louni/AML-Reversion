import { motion } from "framer-motion";

const steps = [
  {
    t: "Learn the normal path",
    p: "From 7,956 healthy-marrow myeloid cells we learn how stem cells mature into monocytes: a pseudotime score and a 2D map, turned into a landscape where stem cells sit high and mature cells low.",
    svg: (
      <svg viewBox="0 0 200 90"><defs><linearGradient id="g1" x1="0" x2="1"><stop offset="0" stopColor="#6f62e8" /><stop offset="1" stopColor="#ffc46b" /></linearGradient></defs>
        <path d="M5 20 C60 20 80 75 195 80" stroke="url(#g1)" strokeWidth="3" fill="none" />
        {[10, 45, 90, 140, 185].map((x, i) => <circle key={x} cx={x} cy={[20, 26, 52, 72, 79][i]} r="4" fill={i < 2 ? "#6f62e8" : i < 3 ? "#c2569b" : "#ffc46b"} />)}</svg>
    ),
  },
  {
    t: "Place each patient",
    p: "Leukaemic cells are placed on the same landscape by expression similarity. Cells below a maturity threshold, fixed on normal cells, count as stuck.",
    svg: (
      <svg viewBox="0 0 200 90"><path d="M5 20 C60 20 80 75 195 80" stroke="rgba(190,180,255,.25)" strokeWidth="2" fill="none" />
        <line x1="105" y1="5" x2="105" y2="88" stroke="#ff7a90" strokeDasharray="3 3" />
        {Array.from({ length: 22 }, (_, i) => <circle key={i} cx={10 + (i * 37) % 90} cy={18 + (i * 13) % 30} r="2.4" fill="#d6d1ff" />)}</svg>
    ),
  },
  {
    t: "Score drugs",
    p: "Each drug's LINCS L1000 expression signature is compared with the change this patient needs to mature. Drugs pointing the same way score high.",
    svg: (
      <svg viewBox="0 0 200 90"><line x1="40" y1="70" x2="150" y2="25" stroke="#ffc46b" strokeWidth="3" markerEnd="url(#a)" />
        <line x1="40" y1="70" x2="165" y2="55" stroke="#6f62e8" strokeWidth="3" />
        <path d="M70 58 A35 35 0 0 0 75 50" stroke="#b2acd3" fill="none" />
        <text x="150" y="20" fill="#ffc46b" fontSize="10">drug</text><text x="160" y="68" fill="#6f62e8" fontSize="10">needed</text></svg>
    ),
  },
  {
    t: "Simulate and test",
    p: "Cells are shifted along the drug signature step by step and re-placed on the map. Rankings are tested against real patients' ex vivo responses (Beat AML).",
    svg: (
      <svg viewBox="0 0 200 90"><path d="M5 20 C60 20 80 75 195 80" stroke="rgba(190,180,255,.25)" strokeWidth="2" fill="none" />
        {[0, 1, 2, 3, 4].map((i) => <circle key={i} cx={30 + i * 32} cy={22 + i * 13} r={3.2} fill="#ffc46b" opacity={0.3 + i * 0.17} />)}
        <circle cx="40" cy="30" r="3.2" fill="#7fe7ff" /></svg>
    ),
  },
];

export function HowItWorks({ onClose }: { onClose: () => void }) {
  return (
    <motion.div className="overlay" role="dialog" aria-modal="true" aria-label="How it works" initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }} onClick={onClose}>
      <motion.div className="sheet" initial={{ y: 30, scale: 0.98 }} animate={{ y: 0, scale: 1 }} exit={{ y: 30 }} onClick={(e) => e.stopPropagation()}>
        <button className="close" onClick={onClose} aria-label="Close (Esc)">Close</button>
        <h2>How it works</h2>
        <p className="lede">Acute myeloid leukaemia cells are stuck in an immature state. Differentiation therapy tries to push them to finish maturing. This project asks, entirely in silico, which existing drugs push each patient's cells that way.</p>
        <div className="steps">
          {steps.map((s, i) => (
            <div className="step" key={s.t}><div className="num">0{i + 1}</div><h3>{s.t}</h3>{s.svg}<p>{s.p}</p></div>
          ))}
        </div>
        <p className="lede" style={{ marginTop: 22, marginBottom: 0 }}>
          The motion is a linear expression-shift model: a visualization of a hypothesis, not a prediction of how real cells behave.
          LINCS signatures come from cell lines after short exposure. Nothing here is a treatment recommendation.
        </p>
      </motion.div>
    </motion.div>
  );
}
