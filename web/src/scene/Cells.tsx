import { useFrame } from "@react-three/fiber";
import { useMemo, useRef } from "react";
import * as THREE from "three";
import type { Landscape, Trajectory } from "../data";
import { makeWorld, terrainY } from "../data";

// Glowing particles resting on the surface. Colour encodes state:
//   stuck (score below threshold)  -> cool lilac-white
//   above threshold               -> warm gold (brighter the more mature)
//   resistant (never crosses)     -> ice cyan with a slow pulse
// A short trail of ghost particles follows each cell while the simulation plays.

const vert = /* glsl */ `
  attribute float aScore;
  attribute float aStatus;
  attribute float aAlpha;
  uniform float uTime;
  uniform float uSize;
  uniform float uPixelRatio;
  uniform float uThreshold;
  uniform float uReduced;
  varying vec3 vColor;
  varying float vAlpha;
  void main() {
    vec4 mv = modelViewMatrix * vec4(position, 1.0);
    gl_Position = projectionMatrix * mv;
    float pulse = 1.0;
    vec3 stuck = vec3(0.72, 0.70, 1.00);
    vec3 gold = mix(vec3(1.0, 0.70, 0.38), vec3(1.0, 0.92, 0.66), clamp((aScore - uThreshold) / (1.0 - uThreshold), 0.0, 1.0));
    vec3 c = aScore >= uThreshold ? gold : stuck;
    if (aStatus > 1.5 && aStatus < 2.5) {
      c = vec3(0.50, 0.91, 1.0);
      pulse = 1.0 + (1.0 - uReduced) * 0.35 * sin(uTime * 3.0 + position.x * 2.0);
    }
    vColor = c;
    vAlpha = aAlpha;
    gl_PointSize = uSize * pulse * uPixelRatio * (24.0 / -mv.z);
  }
`;

const frag = /* glsl */ `
  varying vec3 vColor;
  varying float vAlpha;
  void main() {
    vec2 p = gl_PointCoord - 0.5;
    float d = length(p);
    if (d > 0.5) discard;
    float core = smoothstep(0.5, 0.0, d);
    float glow = pow(core, 2.2);
    gl_FragColor = vec4(pow(vColor * (0.5 + 0.7 * glow), vec3(2.2)), vAlpha * core * 0.7);
  }
`;

export const TRAIL = 6;

export function Cells({ L, traj, status, frameRef, reduced, size = 1 }: {
  L: Landscape; traj: Trajectory | null; status: number[] | null;
  frameRef: React.MutableRefObject<number>; reduced: boolean; size?: number;
}) {
  const world = useMemo(() => makeWorld(L), [L]);
  const n = traj?.cells ?? 0;
  const total = n * (TRAIL + 1);
  const { geometry, material } = useMemo(() => {
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(new Float32Array(Math.max(total, 1) * 3), 3));
    g.setAttribute("aScore", new THREE.BufferAttribute(new Float32Array(Math.max(total, 1)), 1));
    g.setAttribute("aStatus", new THREE.BufferAttribute(new Float32Array(Math.max(total, 1)), 1));
    const alpha = new Float32Array(Math.max(total, 1));
    for (let k = 0; k <= TRAIL; k++) for (let c = 0; c < n; c++) alpha[k * n + c] = k === 0 ? 0.9 : 0.22 * (1 - k / (TRAIL + 1));
    g.setAttribute("aAlpha", new THREE.BufferAttribute(alpha, 1));
    const m = new THREE.ShaderMaterial({
      vertexShader: vert, fragmentShader: frag, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending,
      uniforms: { uTime: { value: 0 }, uSize: { value: 3.4 * size }, uPixelRatio: { value: Math.min(window.devicePixelRatio, 2) },
        uThreshold: { value: L.threshold }, uReduced: { value: reduced ? 1 : 0 } },
    });
    return { geometry: g, material: m };
  }, [total, n, L, size, reduced]);

  const statusArr = useMemo(() => {
    const st = geometry.attributes.aStatus as THREE.BufferAttribute;
    if (status) for (let k = 0; k <= TRAIL; k++) for (let c = 0; c < n; c++) st.setX(k * n + c, status[c]);
    st.needsUpdate = true;
    return st;
  }, [status, geometry, n]);

  const last = useRef(-1);
  const lastTraj = useRef<Trajectory | null>(null);
  useFrame((state) => {
    material.uniforms.uTime.value = state.clock.elapsedTime;
    if (!traj || n === 0) return;
    const f = frameRef.current;
    if (Math.abs(f - last.current) < 1e-4 && lastTraj.current === traj) return;
    last.current = f;
    lastTraj.current = traj;
    const pos = geometry.attributes.position as THREE.BufferAttribute;
    const sc = geometry.attributes.aScore as THREE.BufferAttribute;
    const T = traj.frames, D = traj.data;
    for (let k = 0; k <= TRAIL; k++) {
      const fk = Math.max(0, f - k * 0.55);
      const f0 = Math.min(T - 1, Math.floor(fk)), f1 = Math.min(T - 1, f0 + 1), t = fk - f0;
      for (let c = 0; c < n; c++) {
        const a = (f0 * n + c) * 3, b = (f1 * n + c) * 3;
        const x = D[a] + (D[b] - D[a]) * t, y = D[a + 1] + (D[b + 1] - D[a + 1]) * t, s = D[a + 2] + (D[b + 2] - D[a + 2]) * t;
        const v = k * n + c;
        pos.setXYZ(v, world.toX(x), terrainY(L, x, y, world.HEIGHT) + 0.14, world.toZ(y));
        sc.setX(v, s);
      }
    }
    pos.needsUpdate = true; sc.needsUpdate = true;
    void statusArr;
  });
  if (!traj) return null;
  return <points geometry={geometry} material={material} frustumCulled={false} />;
}

/** Faint dust of normal reference cells, coloured by stage. */
export function ReferenceDust({ L }: { L: Landscape }) {
  const world = useMemo(() => makeWorld(L), [L]);
  const geometry = useMemo(() => {
    const r = L.reference, n = r.x.length;
    const pos = new Float32Array(n * 3), col = new Float32Array(n * 3);
    const tint: Record<string, [number, number, number]> = {
      HSC: [0.45, 0.42, 1.0], Prog: [0.55, 0.55, 1.0], GMP: [0.7, 0.72, 1.0],
      ProMono: [1.0, 0.7, 0.5], Mono: [1.0, 0.82, 0.55], cDC: [0.85, 0.55, 1.0],
    };
    for (let i = 0; i < n; i++) {
      pos.set([world.toX(r.x[i]), terrainY(L, r.x[i], r.y[i], world.HEIGHT) + 0.05, world.toZ(r.y[i])], i * 3);
      col.set(tint[r.type[i]] ?? [0.7, 0.7, 0.9], i * 3);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    g.setAttribute("color", new THREE.BufferAttribute(col, 3));
    return g;
  }, [L, world]);
  return (
    <points geometry={geometry} frustumCulled={false}>
      <pointsMaterial size={0.06} vertexColors transparent opacity={0.16} depthWrite={false} blending={THREE.AdditiveBlending} sizeAttenuation />
    </points>
  );
}
