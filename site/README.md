# MethylGKB · Strata — demo site

The static front end (no build step) for **Strata**: an AI agent that scans a methylome against a drug's response, surfaces the marker-defined responder subgroup, recovers known biology as a positive control, and writes an opportunity brief.

**Live:** https://methylkb.vercel.app/strata.html · **System reference:** `../docs/SYSTEM_MAP.md`

> **Scope (read this):** the canonical current demo is **`strata.html`** — GDSC tumor cell-line, drug-response, **research-use proof-of-method** (not clinical decision-making). The other pages (`index.html`, `science.html`, `explore.html`, `report.html`) are **earlier consumer-MethylGKB pages** from before the Strata pivot; their whole-blood / "what will work, what won't" framing is **legacy** and does not describe the current system.

## Pages

| File | Status |
|------|--------|
| `strata.html` + `strata.js` + `strata.css` | **Current** — the Strata discovery-agent demo (replays pre-generated runs from `demo/`) |
| `demo/` | pre-generated run artifacts (JSON + figures) for Palbociclib + Luminespib, and `scale.json` |
| `index.html`, `science.html`, `explore.html`, `report.html` | **Legacy** — pre-pivot consumer-MethylGKB pages; kept, not canonical |
| `styles.css` | shared design tokens (IBM Carbon × Bloomberg D2: `#F4F4F5` canvas, electric-blue `#0F62FE`, Inter + JetBrains Mono) |
| `app.js` | legacy hero methylation-β track (used by the older pages) |

## Run locally

```bash
python3 -m http.server 8000 --directory site   # then open http://localhost:8000/strata.html
```

## Deploy (static, no build)

```bash
vercel deploy --prod --yes --cwd site
```

---

*Research use only. Not for clinical or diagnostic use. Illustrative prototype.*
