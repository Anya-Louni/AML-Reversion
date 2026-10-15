import { OrbitControls } from "@react-three/drei";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { Bloom, DepthOfField, EffectComposer, Vignette } from "@react-three/postprocessing";
import { Suspense, useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import type { Landscape, Trajectory } from "../data";
import { makeWorld, terrainY } from "../data";
import { Cells, ReferenceDust } from "./Cells";
import { Terrain, type Probe } from "./Terrain";

// Clearance for the side panels, so floating labels never sit underneath them.
const GUTTER_LEFT = 340, GUTTER_RIGHT = 380;

const LABEL_NAMES: Record<string, [string, string]> = {
  HSC: ["HSC", "stem"], Prog: ["Progenitor", "early"], GMP: ["GMP", "committed"],
  ProMono: ["Promonocyte", "maturing"], Mono: ["Monocyte", "mature"], cDC: ["Dendritic", "mature branch"],
};

/** Floating stage labels: plain DOM elements positioned each frame by projecting world points.
 *  (Avoids creating a React root per label inside the canvas.) */
function LabelProjector({ L, host }: { L: Landscape; host: React.RefObject<HTMLDivElement | null> }) {
  const w = useMemo(() => makeWorld(L), [L]);
  const pts = useMemo(() => L.labels.map((l) => new THREE.Vector3(w.toX(l.x), terrainY(L, l.x, l.y, w.HEIGHT) + 0.9, w.toZ(l.y))), [L, w]);
  const v = useMemo(() => new THREE.Vector3(), []);
  useFrame(({ camera, size }) => {
    const el = host.current;
    if (!el) return;
    const taken: [number, number][] = [];
    pts.forEach((p, i) => {
      const node = el.children[i] as HTMLElement | undefined;
      if (!node) return;
      v.copy(p).project(camera);
      const px = ((v.x + 1) / 2) * size.width;
      const py = ((1 - v.y) / 2) * size.height;
      // fade labels that would slide under the side panels rather than letting them clip
      const clear = px > GUTTER_LEFT && px < size.width - GUTTER_RIGHT;
      // at a raking angle stages project on top of each other; drop the later one
      const free = !taken.some(([tx, ty]) => Math.abs(tx - px) < 96 && Math.abs(ty - py) < 34);
      const visible = v.z < 1 && Math.abs(v.x) < 1.1 && Math.abs(v.y) < 1.1 && clear && free;
      if (visible) taken.push([px, py]);
      node.style.opacity = visible ? "1" : "0";
      node.style.transform = `translate(${px}px, ${py}px) translate(-50%, -100%)`;
    });
  });
  return null;
}

function LabelLayer({ L, host }: { L: Landscape; host: React.RefObject<HTMLDivElement | null> }) {
  return (
    <div ref={host} style={{ position: "absolute", inset: 0, pointerEvents: "none", zIndex: 2 }} aria-hidden>
      {L.labels.map((l) => {
        const [name, stage] = LABEL_NAMES[l.name] ?? [l.name, ""];
        return <div key={l.name} className="label3d" style={{ position: "absolute", left: 0, top: 0, transition: "opacity .4s" }}>{name}<small>{stage}</small></div>;
      })}
    </div>
  );
}

/** Nearest labelled stage to a maturity score, for the hover readout. */
function stageOf(score: number, _threshold: number): string {
  const stops: [number, string][] = [[0.03, "HSC-like"], [0.05, "progenitor"], [0.25, "GMP-like"],
                                     [0.65, "promonocyte"], [0.9, "monocyte"], [1.01, "mature monocyte"]];
  return stops.find(([v]) => score <= v)?.[1] ?? "mature";
}

function Intro({ enabled, onDone }: { enabled: boolean; onDone: () => void }) {
  const { camera } = useThree();
  const t0 = useRef<number | null>(null);
  const done = useRef(false);
  const from = useMemo(() => new THREE.Vector3(24, 13, -16), []);
  const to = useMemo(() => new THREE.Vector3(2.5, 5.2, 16.5), []);
  useEffect(() => {
    if (!enabled) { camera.position.copy(to); camera.lookAt(0, 2.1, 0); done.current = true; onDone(); }
  }, [enabled, camera, to, onDone]);
  useFrame((state) => {
    if (done.current) return;
    if (t0.current === null) t0.current = state.clock.elapsedTime;
    const t = Math.min(1, (state.clock.elapsedTime - t0.current) / 4.2);
    const e = 1 - Math.pow(1 - t, 3);
    // sweep along an arc over the stem plateau toward the mature valleys
    const p = from.clone().lerp(to, e);
    p.y += Math.sin(e * Math.PI) * 4;
    camera.position.copy(p);
    camera.lookAt(0, 2.1, 0);
    if (t >= 1) { done.current = true; onDone(); }
  });
  return null;
}

export interface SceneProps {
  L: Landscape;
  support: Float32Array;
  traj: Trajectory | null;
  status: number[] | null;
  frameRef: React.MutableRefObject<number>;
  reduced: boolean;
  effects: boolean;
  intro: boolean;
  onIntroDone: () => void;
  showLabels?: boolean;
}

export function Scene(p: SceneProps) {
  const labelHost = useRef<HTMLDivElement>(null);
  const [probe, setProbe] = useState<Probe | null>(null);
  const stage = probe ? stageOf(probe.score, p.L.threshold) : null;
  return (
    <>
    {p.showLabels !== false && (
      <div className={`probe ${probe ? "on" : ""}`} aria-live="polite">
        {probe ? (
          <>
            <b>{probe.score.toFixed(2)}</b>
            <span>maturity score</span>
            <i className={probe.score >= p.L.threshold ? "mature" : "stuck"}>
              {probe.score >= p.L.threshold ? "past threshold" : "below threshold"}
            </i>
            <em>{stage}</em>
          </>
        ) : <span className="hint-text">hover the landscape to read it</span>}
      </div>
    )}
    {p.showLabels !== false && <div className="scale-note">height = potential, vertically exaggerated for relief</div>}
    {p.showLabels !== false && <LabelLayer L={p.L} host={labelHost} />}
    <Canvas
      dpr={[1, 2]}
      gl={{ antialias: true, powerPreference: "high-performance" }}
      camera={{ position: [2.5, 5.2, 16.5], fov: 48, near: 0.1, far: 240 }}
      onCreated={({ gl }) => { gl.setClearColor("#07060f"); gl.toneMapping = THREE.ACESFilmicToneMapping; }}
    >
      <fog attach="fog" args={["#07060f", 16, 72]} />
      <Suspense fallback={null}>
        <Terrain L={p.L} support={p.support} onProbe={setProbe} />
        <ReferenceDust L={p.L} />
        <Cells L={p.L} traj={p.traj} status={p.status} frameRef={p.frameRef} reduced={p.reduced} />
        {p.showLabels !== false && <LabelProjector L={p.L} host={labelHost} />}
      </Suspense>
      <Intro enabled={p.intro && !p.reduced} onDone={p.onIntroDone} />
      <OrbitControls
        enabled={!p.intro}
        enableDamping dampingFactor={0.06}
        autoRotate={!p.reduced && !p.intro} autoRotateSpeed={0.42}
        minDistance={7} maxDistance={42} minPolarAngle={Math.PI * 0.12} maxPolarAngle={Math.PI * 0.495}
        target={[0, 2.1, 0]}
      />
      {p.effects && (
        <EffectComposer multisampling={0}>
          <Bloom intensity={0.55} luminanceThreshold={0.62} luminanceSmoothing={0.3} mipmapBlur />
          <DepthOfField focusDistance={0.02} focalLength={0.05} bokehScale={1.6} />
          <Vignette eskil={false} offset={0.22} darkness={0.75} />
        </EffectComposer>
      )}
    </Canvas>
    </>
  );
}
