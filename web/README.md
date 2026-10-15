# Reversion Landscape web app

Static single-page app (React, Vite, TypeScript, react-three-fiber) that visualises the pipeline outputs in `public/data/`.

```bash
npm install
npm run dev      # local preview
npm run build    # static build in dist/
```

Data files are produced by `python -m pipeline.run --phase 5` in the repository root. Keyboard: space play/pause, left/right step, [ and ] speed, H how it works.
