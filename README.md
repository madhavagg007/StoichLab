# StoichLab

Explainable chemical equation balancer — offline-first PWA, no dependencies.

**Features:** exact integer balancing (rational nullspace, BigInt) with step-by-step explanation · atom/charge verification · ions, hydrates, brackets · answer checker · molar masses & limiting-reagent calculator · oxidation-state/redox analysis · element picker · 19-reaction library · history · shareable links · installable & offline.

**Run:** open `index.html`, or serve over HTTPS/localhost for install + offline (GitHub Pages works).
**Test:** `node tests/tests.mjs`

```
index.html  manifest.webmanifest  sw.js
css/style.css   js/engine.js (chemistry core)   js/app.js (UI)
data/library.js   icons/icon.svg   tests/tests.mjs
```
