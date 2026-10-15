import { useEffect, useRef } from "react";
import type { Landscape, Trajectory } from "../data";

// 2D fallback for phones, low-power devices and browsers without WebGL:
// the height map drawn as an image with contours, cells as dots.

function ramp(h: number): [number, number, number] {
  const stops: [number, [number, number, number]][] = [
    [0, [255, 219, 140]], [0.18, [247, 143, 82]], [0.38, [199, 84, 148]], [0.58, [115, 84, 209]], [0.8, [51, 46, 133]], [1, [18, 15, 51]],
  ];
  for (let i = 1; i < stops.length; i++) {
    if (h <= stops[i][0]) {
      const [a, ca] = stops[i - 1], [b, cb] = stops[i];
      const t = (h - a) / (b - a);
      return [0, 1, 2].map((k) => ca[k] + (cb[k] - ca[k]) * t) as [number, number, number];
    }
  }
  return stops[stops.length - 1][1];
}

export function FlatView({ L, traj, status, frameRef }: { L: Landscape; traj: Trajectory | null; status: number[] | null; frameRef: React.MutableRefObject<number> }) {
  const ref = useRef<HTMLCanvasElement>(null);
  const bg = useRef<HTMLCanvasElement | null>(null);

  useEffect(() => {
    const n = L.n, c = document.createElement("canvas");
    c.width = n; c.height = n;
    const ctx = c.getContext("2d")!, img = ctx.createImageData(n, n);
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
      const h = L.height[j * n + i], lines = Math.abs((h * 18) % 1 - 0.5) < 0.06 ? 1.25 : 1;
      const [r, g, b] = ramp(h);
      const k = ((n - 1 - j) * n + i) * 4;
      img.data.set([Math.min(255, r * lines), Math.min(255, g * lines), Math.min(255, b * lines), 255], k);
    }
    ctx.putImageData(img, 0, 0);
    bg.current = c;
  }, [L]);

  useEffect(() => {
    let raf = 0;
    const draw = () => {
      const cv = ref.current, b = bg.current;
      if (cv && b) {
        const size = Math.min(cv.parentElement!.clientWidth, cv.parentElement!.clientHeight) * devicePixelRatio;
        if (cv.width !== size) { cv.width = size; cv.height = size; }
        const ctx = cv.getContext("2d")!;
        ctx.imageSmoothingEnabled = true;
        ctx.drawImage(b, 0, 0, size, size);
        if (traj) {
          const f = Math.min(traj.frames - 1, Math.max(0, Math.round(frameRef.current)));
          const D = traj.data, n = traj.cells;
          for (let c = 0; c < n; c++) {
            const a = (f * n + c) * 3;
            const x = (D[a] - L.x[0]) / (L.x[1] - L.x[0]) * size, y = (1 - (D[a + 1] - L.y[0]) / (L.y[1] - L.y[0])) * size;
            const st = status?.[c] ?? 0;
            ctx.fillStyle = st === 2 ? "#7fe7ff" : D[a + 2] >= L.threshold ? "#ffd27a" : "#d6d1ff";
            ctx.beginPath(); ctx.arc(x, y, 2.2 * devicePixelRatio, 0, Math.PI * 2); ctx.fill();
          }
        }
      }
      raf = requestAnimationFrame(draw);
    };
    raf = requestAnimationFrame(draw);
    return () => cancelAnimationFrame(raf);
  }, [traj, status, L, frameRef]);

  return <div className="flat"><canvas ref={ref} style={{ width: "min(100%, 100vh)", aspectRatio: "1" }} aria-label="2D landscape view" /></div>;
}
