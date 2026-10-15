import { useMemo } from "react";
import * as THREE from "three";
import type { Landscape } from "../data";
import { makeWorld, sampleGrid } from "../data";

// Height-coloured terrain with contour lines. High = stem-like plateau (deep indigo/violet),
// low = mature valleys (warm, luminous gold). Lighting is a soft key light plus rim.
const vert = /* glsl */ `
  attribute float aHeight;
  attribute float aSupport;
  attribute float aScore;
  varying float vH;
  varying float vSupport;
  varying float vScore;
  varying vec3 vNormal;
  varying vec3 vWorld;
  void main() {
    vH = aHeight;
    vSupport = aSupport;
    vScore = aScore;
    vNormal = normalize(normalMatrix * normal);
    vec4 w = modelMatrix * vec4(position, 1.0);
    vWorld = w.xyz;
    gl_Position = projectionMatrix * viewMatrix * w;
  }
`;

const frag = /* glsl */ `
  uniform float uThreshold;
  uniform vec3 uFogColor;
  uniform float uFogNear;
  uniform float uFogFar;
  uniform float uMaturity;
  varying float vH;
  varying float vSupport;
  varying float vScore;
  varying vec3 vNormal;
  varying vec3 vWorld;

  vec3 ramp(float h) {
    vec3 c0 = vec3(0.95, 0.66, 0.30); // valley gold
    vec3 c1 = vec3(0.80, 0.36, 0.22); // amber
    vec3 c2 = vec3(0.52, 0.20, 0.42); // magenta
    vec3 c3 = vec3(0.30, 0.22, 0.62); // violet
    vec3 c4 = vec3(0.14, 0.12, 0.38); // indigo
    vec3 c5 = vec3(0.07, 0.06, 0.20); // rim
    if (h < 0.18) return mix(c0, c1, h / 0.18);
    if (h < 0.38) return mix(c1, c2, (h - 0.18) / 0.20);
    if (h < 0.58) return mix(c2, c3, (h - 0.38) / 0.20);
    if (h < 0.80) return mix(c3, c4, (h - 0.58) / 0.22);
    return mix(c4, c5, (h - 0.80) / 0.20);
  }

  void main() {
    vec3 n = normalize(vNormal);
    vec3 L = normalize(vec3(0.75, 0.55, -0.35));
    float diff = clamp(dot(n, L), 0.0, 1.0);
    float diff2 = clamp(dot(n, normalize(vec3(-0.5, 0.35, 0.6))), 0.0, 1.0);
    float rim = pow(1.0 - clamp(n.y, 0.0, 1.0), 2.0);
    vec3 base = ramp(vH);
    float slope = 1.0 - clamp(n.y, 0.0, 1.0);
    vec3 col = base * (0.13 + 1.08 * diff + 0.20 * diff2) * (1.0 - 0.52 * slope)
             + rim * vec3(0.26, 0.22, 0.62) * 0.55;

    // contour lines every 1/18 of height, antialiased with screen-space derivatives
    float lv = vH * 18.0;
    float w = fwidth(lv);
    float line = 1.0 - smoothstep(0.0, w * 1.4, abs(fract(lv) - 0.5) - 0.5 + w);
    col = mix(col, col * 1.8 + 0.07, line * 0.6 * vSupport);

    // the maturity threshold, learned on normal marrow: the line cells must cross to count as matured
    float m = vScore - uMaturity;
    float mw = fwidth(vScore) * 1.5 + 1e-5;
    float mline = 1.0 - smoothstep(0.0, mw, abs(m));
    col = mix(col, vec3(1.0, 0.86, 0.55), mline * 0.55 * vSupport);

    // glow in the deepest valleys (bloom picks this up)
    col += base * smoothstep(0.16, 0.0, vH) * 0.22 * vSupport;
    // fade outside data support toward the night rim
    col = mix(vec3(0.05, 0.04, 0.13), col, 0.25 + 0.75 * vSupport);

    float d = length(vWorld - cameraPosition);
    float f = smoothstep(uFogNear, uFogFar, d);
    // colours above are authored in sRGB; the renderer converts linear -> sRGB on output
    gl_FragColor = vec4(pow(mix(col, uFogColor, f), vec3(2.2)), 1.0);
  }
`;

export interface Probe { score: number; height: number; x: number; y: number; }

export function Terrain({ L, support, onProbe }: {
  L: Landscape; support: Float32Array; onProbe?: (p: Probe | null) => void;
}) {
  const { geometry, material } = useMemo(() => {
    const world = makeWorld(L);
    const n = L.n;
    const g = new THREE.PlaneGeometry(world.W, world.D, n - 1, n - 1);
    g.rotateX(-Math.PI / 2);
    const pos = g.attributes.position as THREE.BufferAttribute;
    const hAttr = new Float32Array(n * n);
    const sAttr = new Float32Array(n * n);
    const scAttr = new Float32Array(n * n);
    const scoreGrid = L.score ?? [];
    // After rotateX(-pi/2), PlaneGeometry row 0 lies at world z = -D/2. World z = -(y - cy) * s, so
    // row 0 is landscape y max: grid row j = n - 1 - r. Columns run x min -> x max.
    for (let r = 0; r < n; r++) {
      for (let c = 0; c < n; c++) {
        const v = r * n + c;
        const j = n - 1 - r, i = c;
        const h = L.height[j * n + i];
        // outside the data the surface falls away, so the landscape is an island in the dark
        const sup = support[j * n + i];
        pos.setY(v, (L.height_vis ?? L.height)[j * n + i] * world.HEIGHT);
        void sup;
        hAttr[v] = h;
        sAttr[v] = support[j * n + i];
        scAttr[v] = scoreGrid[j * n + i] ?? 0;
      }
    }
    g.setAttribute("aHeight", new THREE.BufferAttribute(hAttr, 1));
    g.setAttribute("aSupport", new THREE.BufferAttribute(sAttr, 1));
    g.setAttribute("aScore", new THREE.BufferAttribute(scAttr, 1));
    g.computeVertexNormals();
    const m = new THREE.ShaderMaterial({
      vertexShader: vert, fragmentShader: frag,
      uniforms: {
        uThreshold: { value: L.threshold },
        uMaturity: { value: L.threshold },
        uFogColor: { value: new THREE.Color("#07060f") },
        uFogNear: { value: 16 }, uFogFar: { value: 72 },
      },
    });
    return { geometry: g, material: m };
  }, [L, support]);
  const world = useMemo(() => makeWorld(L), [L]);
  return (
    <mesh
      geometry={geometry}
      material={material}
      onPointerMove={(e) => {
        if (!onProbe) return;
        e.stopPropagation();
        const x = world.fromX(e.point.x), y = world.fromY(e.point.z);
        const inside = x >= L.x[0] && x <= L.x[1] && y >= L.y[0] && y <= L.y[1];
        onProbe(inside ? { x, y, score: sampleGrid(L, L.score ?? L.height, x, y), height: e.point.y } : null);
      }}
      onPointerOut={() => onProbe?.(null)}
    />
  );
}
