# Strata Report API

Thin FastAPI server over precomputed standardized JSON reports. The auto-generated
OpenAPI schema (`/openapi.json`, `/docs`) is the standardized contract an agent calls.

## Reports

Reports are produced by the batch `scripts/build_reports.py` (run from repo root)
and written under `api/reports/`:

```
api/reports/drug/<slug>.json          # one per drug          -> DrugReport
api/reports/cancer_type/<slug>.json   # one per cancer type   -> CancerTypeReport
api/reports/manifest.json             # index                 -> Manifest
```

`api/reports/` is **.gitignored** (generated data, not source). The server reads
these files from disk at request time, so regenerating the reports updates the API
with no code change. If `manifest.json` is missing, the service stays healthy and
returns empty lists.

> **Deploy note:** because reports are gitignored, the Docker build context must
> contain `api/reports/`. Regenerate before deploying:
> `PYTHONPATH=. uv run python scripts/build_reports.py`

## Endpoints

| Method | Path | Returns |
|--------|------|---------|
| GET | `/` | name, docs link, counts |
| GET | `/health` | `{"status":"ok"}` |
| GET | `/manifest` | full manifest |
| GET | `/drugs` | `manifest.drugs` |
| GET | `/cancer_types` | `manifest.cancer_types` |
| GET | `/report/drug/{slug}` | `DrugReport` (404 if missing) |
| GET | `/report/cancer_type/{slug}` | `CancerTypeReport` (404 if missing) |
| GET | `/docs`, `/openapi.json` | interactive docs + schema |

## Run locally

Isolated env (does not touch the main project's deps):

```bash
cd api
uv run --with fastapi --with 'uvicorn[standard]' uvicorn main:app --port 8099
# -> http://127.0.0.1:8099/docs
```

Or with a venv:

```bash
cd api
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn main:app --port 8099
```

## Deploy (Fly.io)

The human runs this — reports must be regenerated first (they are gitignored):

```bash
# from repo root, refresh the served data
PYTHONPATH=. uv run python scripts/build_reports.py

# then deploy from api/
cd api
fly deploy
```

App: `strata-report-api`, internal port 8080, force HTTPS, auto stop/start machines.
