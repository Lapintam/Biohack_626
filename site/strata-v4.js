/* ──────────────────────────────────────────────────────────────────
   Strata v4 — interactive report explorer + agent-transcript replay.
   Vanilla JS, no build, no external libs.

   Data flow:
     1. Live API at API_BASE (the Fly URL) — tried first.
     2. On any fetch error, fall back to the static copies bundled under
        demo-v4/reports/ so the demo is bulletproof even if the API is down.
     3. Transcripts are read from demo-v4/transcripts.json (always static).
   ────────────────────────────────────────────────────────────────── */

const API_BASE = 'https://strata-report-api.fly.dev';
const STATIC_BASE = 'demo-v4/reports';

/* Reflects which source answered the most recent report fetch. */
let lastSource = 'live'; // 'live' | 'cached'

const $ = (sel) => document.querySelector(sel);

/* ── fetch helper: try live, fall back to static ─────────────────── */
async function fetchReportJSON(livePath, staticPath) {
  // livePath relative to API_BASE; staticPath relative to STATIC_BASE
  try {
    const res = await fetch(`${API_BASE}${livePath}`, { cache: 'no-store' });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    lastSource = 'live';
    return data;
  } catch (err) {
    // Any network/HTTP error → static fallback.
    const res = await fetch(`${STATIC_BASE}${staticPath}`, { cache: 'no-store' });
    if (!res.ok) throw new Error(`static fallback failed: HTTP ${res.status}`);
    const data = await res.json();
    lastSource = 'cached';
    return data;
  }
}

function getManifest() {
  return fetchReportJSON('/manifest', '/manifest.json');
}
function getDrugReport(slug) {
  return fetchReportJSON(`/report/drug/${slug}`, `/drug/${slug}.json`);
}
function getCancerTypeReport(slug) {
  return fetchReportJSON(`/report/cancer_type/${slug}`, `/cancer_type/${slug}.json`);
}

/* ── small DOM/format utilities ──────────────────────────────────── */
function esc(s) {
  return String(s == null ? '' : s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}
function fmt(x, digits = 3) {
  if (x === null || x === undefined || Number.isNaN(x)) return '—';
  return Number(x).toFixed(digits);
}
function fmtP(p) {
  if (p === null || p === undefined || Number.isNaN(p)) return '';
  if (p === 0) return 'p<1e-4';
  if (p < 0.001) return `p=${p.toExponential(1)}`;
  return `p=${Number(p).toFixed(3)}`;
}

/* ──────────────────────────────────────────────────────────────────
   Interpret this result — POST /interpret, applies the fixed key to the
   actual on-screen result. Shared by §03 (report) and §04 (annotation).
   ────────────────────────────────────────────────────────────────── */

/* What is currently rendered, so "Interpret this result" knows the ref. */
let currentReportKind = 'drug'; // 'drug' | 'cancer_type' (mirrors currentKind)
let currentReportSlug = null;
let currentDemoId = null;
const interpretBusy = { report: false, annotate: false };

async function fetchInterpret(kind, ref) {
  const res = await fetch(`${API_BASE}/interpret`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ kind, ref }),
    cache: 'no-store',
  });
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  return res.json();
}

function interpretBadge(usedLlm) {
  return usedLlm
    ? '<span class="v4-interpret-badge is-live" title="written by the live LLM">● live</span>'
    : '<span class="v4-interpret-badge is-scripted" title="deterministic templated interpretation">○ scripted</span>';
}

/* Render the interpretation into the host element (replaces its contents). */
function renderInterpretInto(host, data) {
  host.hidden = false;
  host.innerHTML =
    `<div class="v4-interpret-head">` +
    `<span class="v4-interpret-eyebrow">INTERPRETATION</span>` +
    `${interpretBadge(!!data.used_llm)}` +
    `</div>` +
    `<div class="v4-answer-body v4-md">${mdLite(data.interpretation || '')}</div>`;
}

/* Wire one "Interpret this result ▸" button: btnSel triggers a /interpret
   call for refFn()'s {kind, ref} and renders into hostSel. */
function runInterpret(slot, btn, host, kind, ref) {
  if (!btn || !host) return;
  if (!ref) {
    host.hidden = false;
    host.innerHTML = `<div class="v4-interpret-head"><span class="v4-interpret-eyebrow">INTERPRETATION</span></div>` +
      `<div class="v4-answer-body">Nothing to interpret yet.</div>`;
    return;
  }
  if (interpretBusy[slot]) return;
  interpretBusy[slot] = true;
  btn.disabled = true;
  const prev = btn.textContent;
  btn.textContent = 'Interpreting…';
  host.hidden = false;
  host.innerHTML = `<div class="v4-interpret-head"><span class="v4-interpret-eyebrow">INTERPRETATION</span></div>` +
    `<div class="v4-answer-body is-loading">Applying the interpretation key to this result…</div>`;
  fetchInterpret(kind, ref)
    .then((data) => renderInterpretInto(host, data))
    .catch((err) => {
      host.hidden = false;
      host.innerHTML = `<div class="v4-interpret-head"><span class="v4-interpret-eyebrow">INTERPRETATION</span></div>` +
        `<div class="v4-answer-body">Could not interpret (live API required): ${esc(err.message)}</div>`;
    })
    .finally(() => {
      interpretBusy[slot] = false;
      btn.disabled = false;
      btn.textContent = prev;
    });
}

/* ── source indicator ────────────────────────────────────────────── */
function updateSourceIndicator() {
  const el = $('#sourceIndicator');
  if (!el) return;
  el.classList.remove('is-live', 'is-cached');
  const txt = el.querySelector('.v4-source-txt');
  if (lastSource === 'live') {
    el.classList.add('is-live');
    txt.textContent = '● live API';
  } else {
    el.classList.add('is-cached');
    txt.textContent = '○ cached';
  }
}

/* ── honesty rails (shared by reports) ───────────────────────────── */
const FALLBACK_HONESTY = [
  'Single-gene methylation→response is a weak, diffuse signal; the claim is directional concordance across independent screens + functional silencing, not a large effect.',
  'Cell-line results: defensible novel hypotheses (evidence ladder ~L3–L5), not clinical claims.',
  'The agent expands the hypothesis space the expert adjudicates; it removes the structural gate, it does not replace expert judgment.',
];
function renderHonesty(rails) {
  const list = (rails && rails.length) ? rails : FALLBACK_HONESTY;
  const host = $('#honestyRails');
  host.innerHTML = list.map((r, i) => `
    <div class="honesty-rail">
      <span class="hr-n">${String(i + 1).padStart(2, '0')}</span>
      <span>${esc(r)}</span>
    </div>`).join('');
}

/* ──────────────────────────────────────────────────────────────────
   §02 — Transcript replay
   ────────────────────────────────────────────────────────────────── */
let TRANSCRIPTS = [];
let activeTranscript = 0;
let replayTimers = [];

function clearReplayTimers() {
  replayTimers.forEach((t) => clearTimeout(t));
  replayTimers = [];
}

function renderTranscriptTabs() {
  const host = $('#transcriptTabs');
  host.innerHTML = TRANSCRIPTS.map((t, i) => `
    <button class="v4-tab ${i === activeTranscript ? 'is-active' : ''}" data-idx="${i}" type="button">
      <span class="v4-tab-idx">Run ${String(i + 1).padStart(2, '0')}</span>
      <span class="v4-tab-q">${esc(t.query)}</span>
    </button>`).join('');
  host.querySelectorAll('.v4-tab').forEach((b) => {
    b.addEventListener('click', () => {
      activeTranscript = Number(b.dataset.idx);
      renderTranscriptTabs();
      selectTranscript();
    });
  });
}

/* Build a colorized terminal-style tool-call line, e.g.
   → get_drug_report("olaparib") */
function toolCallHTML(step) {
  const args = step.args || {};
  const argKeys = Object.keys(args);
  let argStr;
  if (argKeys.length === 0) {
    argStr = '';
  } else {
    argStr = argKeys.map((k) => {
      const v = args[k];
      const vStr = typeof v === 'string' ? `"${esc(v)}"` : esc(JSON.stringify(v));
      return `<span class="t-arg-v">${vStr}</span>`;
    }).join('<span class="t-paren">, </span>');
  }
  return `<span class="t-paren">→ </span><span class="t-fn">${esc(step.tool)}</span>` +
    `<span class="t-paren">(</span>${argStr}<span class="t-paren">)</span>`;
}

function selectTranscript() {
  clearReplayTimers();
  const t = TRANSCRIPTS[activeTranscript];
  $('#transcriptQuery').textContent = t ? t.query : '';
  $('#transcript').innerHTML = '';
  const ans = $('#transcriptAnswer');
  ans.hidden = true;
  ans.innerHTML = '';
}

function replayTranscript() {
  clearReplayTimers();
  const t = TRANSCRIPTS[activeTranscript];
  if (!t) return;
  const host = $('#transcript');
  host.innerHTML = '';
  const ans = $('#transcriptAnswer');
  ans.hidden = true;
  ans.innerHTML = '';

  const steps = t.steps || [];
  const stagger = 650; // ms between reveals

  steps.forEach((step, i) => {
    // tool call line
    replayTimers.push(setTimeout(() => {
      const call = document.createElement('div');
      call.className = 't-step t-call';
      call.innerHTML = `<div class="t-label">Tool call</div>` +
        `<div class="t-cmd">${toolCallHTML(step)}</div>`;
      host.appendChild(call);
    }, i * 2 * stagger));

    // result summary
    replayTimers.push(setTimeout(() => {
      const result = document.createElement('div');
      result.className = 't-step t-result';
      result.innerHTML = `<div class="t-label">Result</div>` +
        `<div class="t-body">${esc(step.result_summary || '')}</div>`;
      host.appendChild(result);
    }, (i * 2 + 1) * stagger));
  });

  // final answer narration
  replayTimers.push(setTimeout(() => {
    ans.innerHTML = `<span class="v4-answer-lbl">Agent answer</span>` +
      `<div class="v4-answer-body">${esc(t.answer || '')}</div>`;
    ans.hidden = false;
  }, steps.length * 2 * stagger + 200));
}

async function initTranscripts() {
  try {
    const res = await fetch('demo-v4/transcripts.json', { cache: 'no-store' });
    TRANSCRIPTS = await res.json();
  } catch (err) {
    TRANSCRIPTS = [];
  }
  if (!TRANSCRIPTS.length) {
    $('#transcriptTabs').innerHTML =
      '<div class="v4-state">No transcripts available.</div>';
    return;
  }
  activeTranscript = 0;
  renderTranscriptTabs();
  selectTranscript();
  $('#replayBtn').addEventListener('click', replayTranscript);
  // auto-play the first run once on load
  replayTranscript();
}

/* ──────────────────────────────────────────────────────────────────
   §02b — Ask the agent (live /ask, with canned-transcript fallback)
   ────────────────────────────────────────────────────────────────── */
const SUGGESTED_QUERIES = [
  "What's actionable for olaparib?",
  "What's actionable in lung cancer?",
  "Annotate patient nci-h209",
];
let askBusy = false;

/* Light markdown-ish formatter: ## headers, **bold**, line breaks.
   Escapes first, then applies a tiny safe subset of inline markup. */
function mdLite(src) {
  const lines = String(src == null ? '' : src).split(/\r?\n/);
  const out = [];
  for (const raw of lines) {
    const line = raw.trimEnd();
    let m;
    if ((m = line.match(/^\s*(#{1,4})\s+(.*)$/))) {
      const lvl = Math.min(m[1].length, 4);
      out.push(`<div class="v4-md-h v4-md-h${lvl}">${inlineMd(m[2])}</div>`);
    } else if (line.trim() === '') {
      out.push('<div class="v4-md-gap"></div>');
    } else if ((m = line.match(/^\s*[-*]\s+(.*)$/))) {
      out.push(`<div class="v4-md-li">${inlineMd(m[1])}</div>`);
    } else {
      out.push(`<div class="v4-md-p">${inlineMd(line)}</div>`);
    }
  }
  return out.join('');
}
function inlineMd(s) {
  // escape, then re-introduce **bold** only
  let t = esc(s);
  t = t.replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>');
  return t;
}

function renderAskChips() {
  const host = $('#askChips');
  if (!host) return;
  host.innerHTML = SUGGESTED_QUERIES.map((q) =>
    `<button class="v4-ask-chip" type="button" data-q="${esc(q)}">${esc(q)}</button>`).join('');
  host.querySelectorAll('.v4-ask-chip').forEach((b) => {
    b.addEventListener('click', () => {
      $('#askInput').value = b.dataset.q;
      runAsk(b.dataset.q);
    });
  });
}

function setAskMode(mode) {
  // mode: 'live' | 'scripted' | null
  const el = $('#askMode');
  if (!el) return;
  if (!mode) { el.hidden = true; return; }
  el.hidden = false;
  el.classList.remove('is-live', 'is-scripted');
  if (mode === 'live') {
    el.classList.add('is-live');
    el.textContent = '● live';
    el.title = 'answered by the live LLM agent';
  } else {
    el.classList.add('is-scripted');
    el.textContent = '○ scripted';
    el.title = 'served from a bundled example transcript';
  }
}

function renderAskSteps(host, steps) {
  host.innerHTML = '';
  (steps || []).forEach((step) => {
    const call = document.createElement('div');
    call.className = 't-step t-call';
    call.innerHTML = `<div class="t-label">Tool call</div>` +
      `<div class="t-cmd">${toolCallHTML(step)}</div>`;
    host.appendChild(call);
    const result = document.createElement('div');
    result.className = 't-step t-result';
    result.innerHTML = `<div class="t-label">Result</div>` +
      `<div class="t-body">${esc(step.result_summary || '')}</div>`;
    host.appendChild(result);
  });
}

/* Best-effort match of a query to a canned transcript for the offline path. */
function fallbackTranscriptFor(query) {
  const q = (query || '').toLowerCase();
  let best = null, bestScore = -1;
  for (const t of TRANSCRIPTS) {
    const tq = (t.query || '').toLowerCase();
    let score = 0;
    for (const w of q.split(/\W+/)) {
      if (w.length > 2 && tq.includes(w)) score += 1;
    }
    if (score > bestScore) { bestScore = score; best = t; }
  }
  return best || TRANSCRIPTS[0] || null;
}

async function runAsk(query) {
  const q = (query || '').trim();
  if (!q || askBusy) return;
  askBusy = true;
  const shell = $('#askShell');
  const stepsHost = $('#askTranscript');
  const ans = $('#askAnswer');
  shell.hidden = false;
  $('#askQuery').textContent = q;
  setAskMode(null);
  stepsHost.innerHTML = '<div class="t-step t-call"><div class="t-label">Working</div>' +
    '<div class="t-cmd"><span class="t-paren">→ </span><span class="t-fn">agent.run</span>' +
    `<span class="t-paren">(</span><span class="t-arg-v">"${esc(q)}"</span><span class="t-paren">)</span>…</div></div>`;
  ans.hidden = true;
  ans.innerHTML = '';
  $('#askBtn').disabled = true;

  try {
    const res = await fetch(`${API_BASE}/ask`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query: q }),
      cache: 'no-store',
    });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    renderAskSteps(stepsHost, data.steps);
    setAskMode(data.used_llm ? 'live' : 'scripted');
    ans.innerHTML = `<span class="v4-answer-lbl">Agent answer</span>` +
      `<div class="v4-answer-body v4-md">${mdLite(data.answer || '')}</div>`;
    ans.hidden = false;
  } catch (err) {
    // Offline fallback: nearest canned transcript.
    const t = fallbackTranscriptFor(q);
    if (t) {
      renderAskSteps(stepsHost, t.steps);
      setAskMode('scripted');
      ans.innerHTML = `<span class="v4-answer-lbl">Agent answer</span>` +
        `<div class="v4-answer-body v4-md">${mdLite(t.answer || '')}</div>` +
        `<div class="v4-ask-fallback">Live agent unreachable — showing a bundled example run.</div>`;
      ans.hidden = false;
    } else {
      stepsHost.innerHTML = `<div class="t-step t-result"><div class="t-label">Error</div>` +
        `<div class="t-body">Could not reach the agent: ${esc(err.message)}</div></div>`;
    }
  } finally {
    askBusy = false;
    $('#askBtn').disabled = false;
  }
}

function wireAsk() {
  renderAskChips();
  const form = $('#askForm');
  if (form) {
    form.addEventListener('submit', (e) => {
      e.preventDefault();
      runAsk($('#askInput').value);
    });
  }
  const prefill = $('#askPrefillBtn');
  if (prefill) {
    prefill.addEventListener('click', () => {
      const t = TRANSCRIPTS[activeTranscript];
      if (!t) return;
      $('#askInput').value = t.query;
      runAsk(t.query);
      $('#transcriptSection').scrollIntoView({ behavior: 'smooth', block: 'start' });
    });
  }
}

/* ──────────────────────────────────────────────────────────────────
   §04 — Patient methylation annotation
   ────────────────────────────────────────────────────────────────── */
let DEMO_PROFILES = [];
let lastAnnotateSource = 'live';
const ANNOTATE_STATIC_BASE = 'demo-v4/annotations';
const ANNOTATE_LEAKAGE_NOTE =
  'Illustrative: this line was in the discovery cohort (training leakage). A real patient tumor/PBMC sample is leakage-free.';

function updateAnnotateSource() {
  const el = $('#annotateSource');
  if (!el) return;
  el.classList.remove('is-live', 'is-cached');
  const txt = el.querySelector('.v4-source-txt');
  if (lastAnnotateSource === 'live') {
    el.classList.add('is-live');
    txt.textContent = '● live API';
  } else {
    el.classList.add('is-cached');
    txt.textContent = '○ cached';
  }
}

async function getDemoProfiles() {
  try {
    const res = await fetch(`${API_BASE}/demo_profiles`, { cache: 'no-store' });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    // handle both a bare list and {profiles:[...]}
    return Array.isArray(data) ? data : (data.profiles || []);
  } catch (err) {
    // fall back to the bundled annotation ids (labels are best-effort)
    const ids = ['nci-h209', 'dohh-2', 'l-363', 'nci-h660'];
    return ids.map((id) => ({ id, label: id, tissue: null, source: 'cell_line_in_cohort' }));
  }
}

function populateProfileSelect() {
  const sel = $('#profileSelect');
  if (!sel) return;
  sel.innerHTML = DEMO_PROFILES.map((p) => {
    const t = p.tissue ? ` · ${esc(p.tissue)}` : '';
    return `<option value="${esc(p.id)}">${esc(p.label || p.id)}${t}</option>`;
  }).join('');
}

/* annotate via demo_id, live-first with bundled JSON fallback */
async function annotateDemo(demoId) {
  try {
    const res = await fetch(`${API_BASE}/annotate`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ demo_id: demoId }),
      cache: 'no-store',
    });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    lastAnnotateSource = 'live';
    return data;
  } catch (err) {
    const res = await fetch(`${ANNOTATE_STATIC_BASE}/${demoId}.json`, { cache: 'no-store' });
    if (!res.ok) throw new Error(`static fallback failed: HTTP ${res.status}`);
    const data = await res.json();
    lastAnnotateSource = 'cached';
    return data;
  }
}

/* annotate a pasted profile (live only — no static analogue) */
async function annotateProfile(profile) {
  const res = await fetch(`${API_BASE}/annotate`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ profile }),
    cache: 'no-store',
  });
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  lastAnnotateSource = 'live';
  return res.json();
}

const CALL_ORDER = { sensitive: 0, neutral: 1, resistant: 2 };

function annotationPredictionHTML(pred, gtAuc) {
  const call = (pred.call || 'neutral').toLowerCase();
  const drivers = (pred.drivers || []).slice(0, 4);
  // contribution bars are scaled to the max |contribution| in this row
  const maxC = drivers.reduce((m, d) => Math.max(m, Math.abs(d.contribution || 0)), 0) || 1;
  const driversHTML = drivers.map((d) => {
    const c = d.contribution || 0;
    const pct = Math.min(100, (Math.abs(c) / maxC) * 100);
    const sign = c < 0 ? 'neg' : 'pos';
    return `
      <div class="v4-driver">
        <span class="v4-driver-gene">${esc(d.gene)}</span>
        <span class="v4-driver-lvl">${esc(d.evidence_level || '')}</span>
        <span class="v4-driver-bar-wrap">
          <span class="v4-driver-bar ${sign}" style="width:${pct.toFixed(0)}%;"></span>
        </span>
        <span class="v4-driver-pct" title="patient percentile within cohort">${d.patient_pct != null ? Math.round(d.patient_pct * 100) + '%' : '—'}</span>
      </div>`;
  }).join('');

  // measured AUC (illustrative, in-cohort) if present
  let gt = '';
  if (gtAuc && Object.prototype.hasOwnProperty.call(gtAuc, pred.drug) && gtAuc[pred.drug] != null) {
    gt = `<span class="v4-pred-gt" title="measured drug-response AUC for this line — illustrative, in-cohort (training leakage)">
        measured AUC ${fmt(gtAuc[pred.drug], 2)} <em>illustrative, in-cohort</em></span>`;
  }

  return `
    <div class="opp-card v4-pred ${call}">
      <div class="opp-card-head">
        <div class="opp-card-title">
          <span class="opp-card-marker">${esc(pred.drug)}</span>
          <span class="opp-meth-tag ${call === 'resistant' ? 'resistant' : (call === 'sensitive' ? 'sensitive' : '')} v4-call ${call}">${esc(call)}</span>
        </div>
        <div class="opp-card-pathway">score ${fmt(pred.score, 2)} · ${pred.n_leads_scored != null ? pred.n_leads_scored : (pred.drivers || []).length} leads</div>
      </div>
      ${gt ? `<div class="v4-pred-gt-row">${gt}</div>` : ''}
      <div class="v4-drivers">
        <div class="v4-drivers-lbl">Driver genes <span class="v4-drivers-sub">gene · level · contribution · patient %ile</span></div>
        ${driversHTML || '<div class="v4-drivers-empty">No driver genes for this drug.</div>'}
      </div>
    </div>`;
}

function renderAnnotation(r, profileMeta) {
  const out = $('#annotateOut');
  const preds = (r.predictions || []).slice().sort((a, b) => {
    const ca = CALL_ORDER[(a.call || 'neutral').toLowerCase()] ?? 1;
    const cb = CALL_ORDER[(b.call || 'neutral').toLowerCase()] ?? 1;
    if (ca !== cb) return ca - cb;
    return (b.score || 0) - (a.score || 0); // within group, higher score first
  });
  const gt = r.ground_truth_auc || null;
  const note = (profileMeta && profileMeta.source === 'cell_line_in_cohort') ? ANNOTATE_LEAKAGE_NOTE : null;

  const label = r.source_label || (profileMeta && profileMeta.label) || 'Pasted profile';
  const tissue = r.tissue || (profileMeta && profileMeta.tissue) || null;

  out.innerHTML = `
    <div class="v4-rep-head">
      <div class="v4-rep-title">
        <span class="v4-rep-drug">${esc(label)}</span>
        <span class="v4-rep-type">Annotation</span>
      </div>
      <div class="v4-rep-meta">
        ${tissue ? `<span><span class="v4-rm-k">tissue</span> <span class="v4-rm-v">${esc(tissue)}</span></span>` : ''}
        <span><span class="v4-rm-k">drugs scored</span> <span class="v4-rm-v">${r.n_drugs_scored != null ? r.n_drugs_scored : preds.length}</span></span>
      </div>
    </div>

    ${note ? `<div class="v4-leakage-note"><span class="v4-leakage-icon">⚠</span> ${esc(note)}</div>` : ''}

    <p class="v4-sub-marker">Predicted sensitivities · ${preds.length} drugs · most-sensitive first</p>
    <div class="opp-cards v4-pred-cards">
      ${preds.map((p) => annotationPredictionHTML(p, gt)).join('') || '<div class="v4-state">No predictions for this profile.</div>'}
    </div>
  `;

  renderHonesty(r.honesty);
}

function hideAnnotateInterpret() {
  const host = $('#annotateInterpretOut');
  if (host) { host.hidden = true; host.innerHTML = ''; }
  const btn = $('#annotateInterpretBtn');
  if (btn) btn.hidden = true;
}

async function loadSelectedProfile() {
  const sel = $('#profileSelect');
  if (!sel || !sel.value) return;
  const demoId = sel.value;
  const meta = DEMO_PROFILES.find((p) => p.id === demoId) || { id: demoId };
  const out = $('#annotateOut');
  hideAnnotateInterpret();
  out.innerHTML = `<div class="v4-state is-loading">Annotating ${esc(meta.label || demoId)}…</div>`;
  try {
    const r = await annotateDemo(demoId);
    updateAnnotateSource();
    renderAnnotation(r, meta);
    // track the demo id on screen and offer interpretation
    currentDemoId = demoId;
    const btn = $('#annotateInterpretBtn');
    if (btn) btn.hidden = false;
  } catch (err) {
    out.innerHTML = `<div class="v4-state is-error">Could not annotate: ${esc(err.message)}</div>`;
  }
}

async function runPasteAnnotate() {
  const ta = $('#pasteArea');
  const out = $('#annotateOut');
  let profile;
  try {
    profile = JSON.parse(ta.value);
    if (!profile || typeof profile !== 'object' || Array.isArray(profile)) {
      throw new Error('expected a JSON object of {gene: beta}');
    }
  } catch (err) {
    out.innerHTML = `<div class="v4-state is-error">Invalid JSON: ${esc(err.message)}</div>`;
    return;
  }
  hideAnnotateInterpret();
  out.innerHTML = `<div class="v4-state is-loading">Annotating pasted profile…</div>`;
  try {
    const r = await annotateProfile(profile);
    updateAnnotateSource();
    renderAnnotation(r, { source: 'user', label: 'Pasted profile' });
    // pasted profiles have no ref the /interpret endpoint can re-fetch; keep it hidden
    currentDemoId = null;
  } catch (err) {
    out.innerHTML = `<div class="v4-state is-error">Could not annotate (live API required for custom profiles): ${esc(err.message)}</div>`;
  }
}

async function initAnnotate() {
  DEMO_PROFILES = await getDemoProfiles();
  populateProfileSelect();
  const sel = $('#profileSelect');
  if (sel) sel.addEventListener('change', loadSelectedProfile);
  const pasteBtn = $('#pasteBtn');
  if (pasteBtn) pasteBtn.addEventListener('click', runPasteAnnotate);
  const aibtn = $('#annotateInterpretBtn');
  if (aibtn) {
    aibtn.addEventListener('click', () =>
      runInterpret('annotate', aibtn, $('#annotateInterpretOut'), 'annotation', currentDemoId));
  }
  if (sel && DEMO_PROFILES.length) {
    // prefer nci-h209 as the starter
    const pref = DEMO_PROFILES.find((p) => p.id === 'nci-h209');
    sel.value = pref ? pref.id : DEMO_PROFILES[0].id;
    loadSelectedProfile();
  } else {
    $('#annotateOut').innerHTML = '<div class="v4-state is-error">No demo profiles available from API or cache.</div>';
  }
}

/* ──────────────────────────────────────────────────────────────────
   §03 — Live report explorer
   ────────────────────────────────────────────────────────────────── */
let MANIFEST = { drugs: [], cancer_types: [], gates: [] };
let currentKind = 'drug'; // 'drug' | 'cancer_type'

function populateSelect() {
  const sel = $('#reportSelect');
  const entries = currentKind === 'drug' ? MANIFEST.drugs : MANIFEST.cancer_types;
  sel.innerHTML = entries.map((e) => {
    const meta = currentKind === 'drug'
      ? `${e.n_leads} leads · ${e.replicated} replicated`
      : `${e.n_leads} leads`;
    return `<option value="${esc(e.slug)}">${esc(e.name)} — ${esc(meta)}</option>`;
  }).join('');
}

/* default-select a sensible starter for each kind */
function defaultSlug() {
  const entries = currentKind === 'drug' ? MANIFEST.drugs : MANIFEST.cancer_types;
  const preferred = currentKind === 'drug' ? 'olaparib' : 'lung';
  const hit = entries.find((e) => e.slug === preferred);
  return hit ? hit.slug : (entries[0] && entries[0].slug);
}

function hideReportInterpret() {
  const host = $('#reportInterpretOut');
  if (host) { host.hidden = true; host.innerHTML = ''; }
  const btn = $('#reportInterpretBtn');
  if (btn) btn.hidden = true;
}

async function loadSelectedReport() {
  const slug = $('#reportSelect').value;
  if (!slug) return;
  const out = $('#reportOut');
  hideReportInterpret();
  out.innerHTML = `<div class="v4-state is-loading">Fetching report for ${esc(slug)}…</div>`;
  try {
    let report;
    if (currentKind === 'drug') {
      report = await getDrugReport(slug);
      updateSourceIndicator();
      renderDrugReport(report);
    } else {
      report = await getCancerTypeReport(slug);
      updateSourceIndicator();
      renderCancerTypeReport(report);
    }
    renderHonesty(report.honesty);
    // track what's on screen and offer interpretation
    currentReportKind = currentKind;
    currentReportSlug = slug;
    const btn = $('#reportInterpretBtn');
    if (btn) btn.hidden = false;
  } catch (err) {
    out.innerHTML = `<div class="v4-state is-error">Could not load report: ${esc(err.message)}</div>`;
  }
}

/* ── drug report render ──────────────────────────────────────────── */
function renderDrugReport(r) {
  const f = r.funnel || {};
  const funnelStages = [
    { label: 'Genes scanned', note: 'every promoter gene enters', n: f.genes_scanned },
    { label: 'Response candidates', note: 'association survives tissue + methylome-PC adjustment', n: f.response_candidates },
    { label: 'Functionally silenced', note: 'methylation tracks expression loss', n: f.functionally_silenced },
    { label: 'Externally replicated', note: 'same direction in independent screens', n: f.externally_replicated },
  ];

  const funnelHTML = funnelStages.map((s, i) => `
    ${i > 0 ? '<div class="funnel-arrow">↓</div>' : ''}
    <div class="funnel-stage">
      <span class="funnel-idx">${String(i + 1).padStart(2, '0')}</span>
      <div class="funnel-main">
        <span class="funnel-label">${esc(s.label)}</span>
        <span class="funnel-note">${esc(s.note)}</span>
      </div>
      <span class="funnel-n">${s.n != null ? Number(s.n).toLocaleString() : '—'}</span>
    </div>`).join('');

  const heroGene = r.hero ? r.hero.gene : null;
  const cardsHTML = (r.leads || []).map((lead) =>
    drugLeadCardHTML(lead, lead.gene === heroGene)
  ).join('');

  $('#reportOut').innerHTML = `
    <div class="v4-rep-head">
      <div class="v4-rep-title">
        <span class="v4-rep-drug">${esc(r.drug)}</span>
        <span class="v4-rep-type">Drug report</span>
      </div>
      <div class="v4-rep-meta">
        ${r.target ? `<span><span class="v4-rm-k">target</span> <span class="v4-rm-v">${esc(r.target)}</span></span>` : ''}
        ${r.pathway ? `<span><span class="v4-rm-k">pathway</span> <span class="v4-rm-v">${esc(r.pathway)}</span></span>` : ''}
        ${f.lines != null ? `<span><span class="v4-rm-k">cell lines</span> <span class="v4-rm-v">${Number(f.lines).toLocaleString()}</span></span>` : ''}
      </div>
    </div>

    <div class="v4-funnel-wrap">
      <p class="v4-sub-marker">Discovery funnel</p>
      <div class="funnel">${funnelHTML}</div>
    </div>

    <p class="v4-sub-marker">Leads · ${(r.leads || []).length}${heroGene ? ` · hero ${esc(heroGene)}` : ''}</p>
    <div class="opp-cards">${cardsHTML || '<div class="v4-state">No leads in this report.</div>'}</div>
  `;
}

function drugLeadCardHTML(lead, isHero) {
  const dir = lead.direction === 'resistant' ? 'resistant' : 'sensitive';
  const rep = lead.recurrence || {};

  // per-screen replication strip
  const screens = [];
  if (lead.gdsc_r != null) screens.push({ name: 'GDSC', r: lead.gdsc_r, p: null });
  if (lead.prism_r != null) screens.push({ name: 'PRISM', r: lead.prism_r, p: lead.prism_p });
  if (lead.ctrp_r != null) screens.push({ name: 'CTRP', r: lead.ctrp_r, p: lead.ctrp_p });

  const screensHTML = screens.map((s) => {
    const sign = s.r < 0 ? 'neg' : 'pos';
    return `
      <div class="opp-screen">
        <span class="opp-screen-name">${esc(s.name)}</span>
        <span class="opp-screen-r ${sign}">${fmt(s.r)}</span>
        ${s.p != null ? `<span class="opp-screen-p">${esc(fmtP(s.p))}</span>` : ''}
      </div>`;
  }).join('');

  const repStrength = lead.replicated ? 'strong' : 'partial';
  const repLabel = lead.replicated
    ? `replicated · ${lead.n_screens || screens.length} screens`
    : `partial · ${lead.n_screens || screens.length} screens`;

  // recurrence mini-bars
  const tested = rep.tissues_tested || 0;
  const same = rep.same_direction || 0;
  const sig = rep.significant || 0;
  let bars = '';
  for (let i = 0; i < tested; i++) {
    const isSame = i < same;
    const isSig = i < sig;
    const cls = (dir === 'resistant' && isSame) ? 'pos' : '';
    bars += `<span class="opp-recur-bar ${isSame ? cls : ''} ${isSig ? 'sig' : ''}"
      style="${isSame ? '' : 'opacity:0.25;'}"></span>`;
  }

  return `
    <div class="opp-card ${isHero ? 'is-lead' : ''}">
      <div class="opp-card-head">
        <div class="opp-card-title">
          <span class="opp-card-marker">${esc(lead.gene)}</span>
          <span class="opp-meth-tag ${dir}">${esc(dir)}</span>
          <span class="opp-level-badge">${esc(lead.evidence_level)}</span>
          ${isHero ? '<span class="v4-hero-flag">★ hero</span>' : ''}
        </div>
        <div class="opp-card-pathway">GDSC r = ${fmt(lead.gdsc_r)}</div>
      </div>

      ${lead.gene_note ? `
      <div class="opp-gene-note">
        <span class="opp-gene-note-lbl">Why an expert would miss it</span>
        ${esc(lead.gene_note)}
      </div>` : ''}

      <div class="opp-rep-strip">
        <div class="opp-rep-strip-head">
          <span class="opp-rep-strip-lbl">External replication</span>
          <span class="opp-rep-strength ${repStrength}">${esc(repLabel)}</span>
        </div>
        <div class="opp-rep-screens">${screensHTML}</div>
      </div>

      <div class="opp-silencing">
        <span class="opp-silencing-lbl">methylation <span class="arrow">→</span> expression loss</span>
        <span class="opp-silencing-val">${fmt(lead.silencing_r)}</span>
      </div>

      <div class="opp-recurrence">
        <div class="opp-recurrence-lbl">
          same direction in <strong>${same}/${tested}</strong> cancer types
          ${sig ? ` · ${sig} significant` : ''}
        </div>
        <div class="opp-recur-bars">${bars}</div>
      </div>
    </div>`;
}

/* ── cancer-type report render ───────────────────────────────────── */
function renderCancerTypeReport(r) {
  const leads = r.leads || [];
  const rows = leads.map((l) => {
    const dir = l.direction === 'resistant' ? 'resistant' : 'sensitive';
    const rhoSign = (l.tissue_rho != null && l.tissue_rho < 0) ? 'neg' : 'pos';
    const also = (l.also_in && l.also_in.length)
      ? l.also_in.join(', ')
      : '—';
    return `
      <tr>
        <td><span class="v4-ct-drug">${esc(l.drug)}</span></td>
        <td><span class="v4-ct-gene">${esc(l.gene)}</span></td>
        <td><span class="v4-dir-tag ${dir}">${esc(dir)}</span></td>
        <td class="num ${rhoSign}">${fmt(l.tissue_rho)}</td>
        <td class="num">${esc(fmtP(l.tissue_p))}</td>
        <td><span class="v4-lvl-badge">${esc(l.evidence_level)}</span></td>
        <td>${l.replicated ? '<span class="v4-rep-yes">✓ yes</span>' : '<span class="v4-rep-no">no</span>'}</td>
        <td><span class="v4-ct-also" title="${esc(also)}">${esc(also)}</span></td>
      </tr>`;
  }).join('');

  const connects = (r.connects_to || []).map((c) =>
    `<span class="v4-chip">${esc(c)}</span>`).join('');

  $('#reportOut').innerHTML = `
    <div class="v4-rep-head">
      <div class="v4-rep-title">
        <span class="v4-rep-drug">${esc(r.cancer_type)}</span>
        <span class="v4-rep-type">Cancer-type report</span>
      </div>
      <div class="v4-rep-meta">
        <span><span class="v4-rm-k">leads</span> <span class="v4-rm-v">${r.n_leads != null ? r.n_leads : leads.length}</span></span>
        <span><span class="v4-rm-k">connects to</span> <span class="v4-rm-v">${(r.connects_to || []).length} cancer types</span></span>
      </div>
    </div>

    <p class="v4-sub-marker">Leads · drug → gene</p>
    <div class="v4-ct-table-wrap">
      <table class="v4-ct-table">
        <thead>
          <tr>
            <th>Drug</th><th>Gene</th><th>Direction</th>
            <th class="num">tissue ρ</th><th class="num">p</th>
            <th>Evidence</th><th>Replicated</th><th>Also in</th>
          </tr>
        </thead>
        <tbody>${rows || '<tr><td colspan="8">No leads in this report.</td></tr>'}</tbody>
      </table>
    </div>

    ${connects ? `
    <p class="v4-sub-marker">Connects to</p>
    <div class="v4-connects">
      <div class="v4-connects-lbl">
        The same molecular states surface across <strong>${(r.connects_to || []).length}</strong> other cancer types:
      </div>
      <div class="v4-chips">${connects}</div>
    </div>` : ''}
  `;
}

/* ── explorer wiring ─────────────────────────────────────────────── */
function wireExplorer() {
  // type toggle
  $('#reportToggle').querySelectorAll('.v4-toggle-btn').forEach((b) => {
    b.addEventListener('click', () => {
      if (b.classList.contains('is-active')) return;
      $('#reportToggle').querySelectorAll('.v4-toggle-btn').forEach((x) => x.classList.remove('is-active'));
      b.classList.add('is-active');
      currentKind = b.dataset.kind;
      populateSelect();
      const d = defaultSlug();
      if (d) $('#reportSelect').value = d;
      loadSelectedReport();
    });
  });
  // select change
  $('#reportSelect').addEventListener('change', loadSelectedReport);
  // interpret this result
  const ibtn = $('#reportInterpretBtn');
  if (ibtn) {
    ibtn.addEventListener('click', () =>
      runInterpret('report', ibtn, $('#reportInterpretOut'), currentReportKind, currentReportSlug));
  }
}

/* ──────────────────────────────────────────────────────────────────
   Init
   ────────────────────────────────────────────────────────────────── */
async function init() {
  await initTranscripts();
  wireAsk();
  wireExplorer();
  initAnnotate();

  try {
    MANIFEST = await getManifest();
    updateSourceIndicator();
  } catch (err) {
    MANIFEST = { drugs: [], cancer_types: [], gates: [] };
  }

  // render honesty immediately from fallback so §04 is never blank;
  // it gets refreshed from each report's own honesty array.
  renderHonesty(null);

  populateSelect();
  const d = defaultSlug();
  if (d) {
    $('#reportSelect').value = d;
    loadSelectedReport();
  } else {
    $('#reportOut').innerHTML = '<div class="v4-state is-error">Manifest empty — no reports available from API or cache.</div>';
  }
}

document.addEventListener('DOMContentLoaded', init);
