/* MethylGKB — minimal interactions.
   The only scripted thing: render the CpG/beta methylation track in the
   hero. Each mark is a CpG; its fill encodes beta (teal=unmethylated →
   rose=methylated), the same data colors as the polymerbio.org viewer. */

(function () {
  const track = document.getElementById('betaTrack');
  if (!track) return;

  // A plausible MGMT-promoter-window beta profile: mostly low, with a
  // hypermethylated stretch in the middle (the actionable signal).
  const N = 44;
  const betas = [];
  for (let i = 0; i < N; i++) {
    const center = N * 0.58;
    const bump = Math.exp(-Math.pow((i - center) / (N * 0.12), 2)); // gaussian hotspot
    const base = 0.06 + 0.10 * Math.abs(Math.sin(i * 1.7));         // low, textured baseline
    let b = base + 0.9 * bump;
    b = Math.max(0.02, Math.min(0.98, b));
    betas.push(b);
  }

  const TEAL = [8, 160, 151];     // --cpg-teal  #08A097
  const AMBER = [180, 83, 9];     // --accent-amber #B45309 (midpoint)
  const ROSE = [190, 18, 60];     // --methyl-rose #BE123C

  function lerp(a, b, t) { return Math.round(a + (b - a) * t); }
  function betaColor(b) {
    // two-stop ramp through amber, matching the .beta-scale gradient
    if (b < 0.5) {
      const t = b / 0.5;
      return `rgb(${lerp(TEAL[0], AMBER[0], t)},${lerp(TEAL[1], AMBER[1], t)},${lerp(TEAL[2], AMBER[2], t)})`;
    }
    const t = (b - 0.5) / 0.5;
    return `rgb(${lerp(AMBER[0], ROSE[0], t)},${lerp(AMBER[1], ROSE[1], t)},${lerp(AMBER[2], ROSE[2], t)})`;
  }

  const reduce = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  betas.forEach((b, i) => {
    const cpg = document.createElement('div');
    cpg.className = 'cpg';

    const stem = document.createElement('div');
    stem.className = 'stem';
    // stem height tracks beta — taller = more methylated
    stem.style.flex = `0 0 ${10 + b * 32}px`;

    const dot = document.createElement('div');
    dot.className = 'dot';
    const color = betaColor(b);
    dot.style.borderColor = color;
    // fill proportional to beta: open at low beta, solid at high
    dot.style.background = b > 0.55 ? color : 'var(--bg-elevated)';
    if (b > 0.55) dot.style.borderColor = color;
    dot.style.animationDelay = reduce ? '0s' : (i * 22) + 'ms';
    dot.title = `CpG ${i + 1} · β = ${b.toFixed(2)}`;

    cpg.appendChild(stem);
    cpg.appendChild(dot);
    track.appendChild(cpg);
  });
})();


/* ──────────────────────────────────────────────────────────────────
   Patient report (report.html). Hit "Analyze sample" → the methylome
   scans → an evidence-graded PGx report renders, grouped by call.
   Each row reuses the formal-claim shape under an expandable detail.
   ────────────────────────────────────────────────────────────────── */
(function () {
  const runBtn = document.getElementById('runBtn');
  if (!runBtn) return;

  const reduce = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  // Seed associations. Grades are honest: MGMT→TMZ is clinical (1A);
  // the rest are candidate/research and labelled as such.
  const ASSOC = [
    {
      drug: 'Temozolomide', cls: 'Alkylating agent · oncology',
      locus: 'MGMT promoter', cpg: 'cg12434587', beta: 0.86, signal: 'hypermethylated',
      call: 'favorable', grade: '1A',
      note: 'Promoter hypermethylation silences MGMT, the repair enzyme that reverses temozolomide damage — predicting favorable response. The canonical, clinically deployed pharmaco-epigenomic biomarker (glioblastoma).',
      prov: 'Hegi/Stupp lineage · PharmGKB 1A',
      region: 'chr10:129,466,683–129,467,448',
    },
    {
      drug: 'Carboplatin', cls: 'Platinum agent · oncology',
      locus: 'Methylation sensitivity signature', cpg: 'multi-CpG (n=214)', beta: 0.63, signal: 'intermediate',
      call: 'favorable', grade: '2B',
      note: 'A methylation-only model across 987 GDSC cell lines associates this profile with platinum sensitivity. Research-grade; awaiting prospective confirmation.',
      prov: 'PLOS ONE 2021 · candidate',
      region: 'genome-wide signature',
    },
    {
      drug: 'Fluorouracil (5-FU)', cls: 'Antimetabolite · oncology',
      locus: 'TYMS regulatory region', cpg: 'cg23145901', beta: 0.21, signal: 'hypomethylated',
      call: 'caution', grade: 'candidate',
      note: 'Hypomethylation here tracks with higher thymidylate-synthase expression, a candidate marker of reduced 5-FU efficacy. Exploratory signal — confirm with orthogonal testing.',
      prov: 'MethylGKB candidate · unverified',
      region: 'chr18:657,604–657,790',
    },
    {
      drug: 'Irinotecan', cls: 'Topoisomerase-I inhibitor · oncology',
      locus: 'UGT1A1 promoter', cpg: 'cg08176849', beta: 0.74, signal: 'hypermethylated',
      call: 'avoid', grade: 'candidate',
      note: 'Promoter hypermethylation is a candidate epigenetic correlate of reduced UGT1A1 activity — the same axis as the *28 toxicity allele. Elevated neutropenia / diarrhea risk; candidate-grade only.',
      prov: 'MethylGKB candidate · unverified',
      region: 'chr2:233,757,013–233,757,300',
    },
    {
      drug: 'Cytarabine', cls: 'Antimetabolite · hematology',
      locus: 'Whole-blood methylation clock', cpg: 'clock (n=353)', beta: 0.51, signal: 'age-adjusted',
      call: 'monitor', grade: 'research',
      note: 'A whole-blood methylation-age offset is being evaluated as a modifier of clearance and marrow tolerance. Research signal — included to show longitudinal direction.',
      prov: 'MethylGKB research lane',
      region: 'genome-wide clock',
    },
  ];

  const CALLS = {
    favorable: { label: 'Favorable', cls: 'favorable', blurb: 'predicted response / standard dosing' },
    caution:   { label: 'Caution',   cls: 'caution',   blurb: 'possible reduced efficacy' },
    avoid:     { label: 'Avoid / adjust', cls: 'avoid', blurb: 'elevated toxicity risk' },
    monitor:   { label: 'Monitor',   cls: 'monitor',   blurb: 'longitudinal signal' },
  };
  const ORDER = ['avoid', 'caution', 'monitor', 'favorable'];

  // Beta color ramp (shared with the hero track)
  const TEAL = [8, 160, 151], AMBER = [180, 83, 9], ROSE = [190, 18, 60];
  const lerp = (a, b, t) => Math.round(a + (b - a) * t);
  function betaColor(b) {
    if (b < 0.5) { const t = b / 0.5; return `rgb(${lerp(TEAL[0],AMBER[0],t)},${lerp(TEAL[1],AMBER[1],t)},${lerp(TEAL[2],AMBER[2],t)})`; }
    const t = (b - 0.5) / 0.5; return `rgb(${lerp(AMBER[0],ROSE[0],t)},${lerp(AMBER[1],ROSE[1],t)},${lerp(AMBER[2],ROSE[2],t)})`;
  }

  // Build the scan track
  const scanTrack = document.getElementById('scanTrack');
  function buildScan() {
    const N = 60;
    for (let i = 0; i < N; i++) {
      const b = Math.max(0.03, Math.min(0.97, 0.1 + 0.85 * Math.exp(-Math.pow((i - N * 0.42) / (N * 0.1), 2)) + 0.12 * Math.abs(Math.sin(i * 2.1))));
      const cpg = document.createElement('div'); cpg.className = 'cpg';
      const stem = document.createElement('div'); stem.className = 'stem'; stem.style.flex = `0 0 ${8 + b * 30}px`;
      const dot = document.createElement('div'); dot.className = 'dot';
      const c = betaColor(b); dot.style.borderColor = c; dot.style.background = b > 0.55 ? c : 'var(--bg-elevated)';
      dot.style.animationDelay = '0s'; dot.style.transform = 'scale(1)';
      cpg.appendChild(stem); cpg.appendChild(dot); scanTrack.appendChild(cpg);
    }
  }

  function recRow(a) {
    const c = CALLS[a.call];
    const color = betaColor(a.beta);
    const row = document.createElement('div');
    row.className = 'rec';
    row.innerHTML = `
      <button class="rec-head" aria-expanded="false">
        <span class="rec-call ${c.cls}">${c.label}</span>
        <span class="rec-main">
          <span class="rec-drug">${a.drug}</span>
          <span class="rec-cls">${a.cls}</span>
        </span>
        <span class="rec-signal">
          <span class="readout-dot" style="border-color:${color};background:${a.beta > 0.55 ? color : 'transparent'}"></span>
          <span class="rec-locus">${a.locus}</span>
          <span class="rec-beta mono">β ${a.beta.toFixed(2)}</span>
        </span>
        <span class="rec-grade"><span class="badge ${a.grade === '1A' ? 'licensed' : 'grade'}">${a.grade === '1A' ? 'GRADE 1A' : a.grade.toUpperCase()}</span></span>
        <span class="rec-chev">▾</span>
      </button>
      <div class="rec-detail" hidden>
        <p class="rec-note">${a.note}</p>
        <div class="claim" style="margin-top:14px;">
          <div class="claim-head"><span>claim · methylation/pgx</span><span>${c.blurb}</span></div>
          <div class="claim-body">
            <div class="kv"><span class="k">subject</span><span class="v">genomic_region · <span class="hl">${a.locus}</span> · ${a.region}</span></div>
            <div class="kv"><span class="k">probe</span><span class="v">${a.cpg} · β=${a.beta.toFixed(2)} (${a.signal})</span></div>
            <div class="kv"><span class="k">predicate</span><span class="v">predicts-response-to · <span class="hl">${a.drug}</span></span></div>
            <div class="kv"><span class="k">evidence</span><span class="v">${a.prov}</span></div>
          </div>
        </div>
      </div>`;
    const head = row.querySelector('.rec-head');
    const detail = row.querySelector('.rec-detail');
    head.addEventListener('click', () => {
      const open = head.getAttribute('aria-expanded') === 'true';
      head.setAttribute('aria-expanded', String(!open));
      detail.hidden = open;
      row.classList.toggle('is-open', !open);
    });
    return row;
  }

  function renderReport() {
    // summary band
    const counts = { favorable: 0, caution: 0, avoid: 0, monitor: 0 };
    ASSOC.forEach(a => counts[a.call]++);
    const summary = document.getElementById('reportSummary');
    summary.className = 'metrics report-summary';
    summary.innerHTML = `
      <div class="metric"><div class="val" style="color:#15803D">${counts.favorable}</div><div class="label">Favorable — predicted response / standard dosing.</div><div class="src">recommend</div></div>
      <div class="metric"><div class="val" style="color:var(--accent-amber)">${counts.caution + counts.monitor}</div><div class="label">Caution / monitor — reduced efficacy or longitudinal signal.</div><div class="src">review</div></div>
      <div class="metric"><div class="val rose">${counts.avoid}</div><div class="label">Avoid / adjust — elevated toxicity risk.</div><div class="src">flag</div></div>`;

    const mount = document.getElementById('reportMount');
    mount.innerHTML = '<p class="section-marker" style="margin-top:40px;">§02 — Recommendations</p>';
    const list = document.createElement('div'); list.className = 'rec-list';
    const sorted = [...ASSOC].sort((x, y) => ORDER.indexOf(x.call) - ORDER.indexOf(y.call));
    sorted.forEach((a, i) => {
      const row = recRow(a);
      if (!reduce) { row.style.opacity = '0'; row.style.transform = 'translateY(8px)'; }
      list.appendChild(row);
      if (!reduce) setTimeout(() => { row.style.transition = 'opacity .4s, transform .4s'; row.style.opacity = '1'; row.style.transform = 'none'; }, 90 * i);
    });
    mount.appendChild(list);

    const section = document.getElementById('reportSection');
    section.hidden = false;
    section.scrollIntoView({ behavior: reduce ? 'auto' : 'smooth', block: 'start' });
  }

  function runScan(done) {
    const scan = document.getElementById('scan');
    const pct = document.getElementById('scanPct');
    const label = document.getElementById('scanLabel');
    const sweep = document.getElementById('scanSweep');
    scan.hidden = false;
    if (reduce) { done(); return; }
    let p = 0;
    sweep.style.animation = 'sweep 1.6s linear forwards';
    const labels = ['Reading methylome…', 'Calling β-values…', 'Matching MethylGKB…', 'Grading evidence…'];
    const t = setInterval(() => {
      p += Math.random() * 14 + 6;
      if (p >= 100) { p = 100; clearInterval(t); pct.textContent = '100%'; label.textContent = 'Report ready'; setTimeout(done, 350); }
      else { pct.textContent = Math.floor(p) + '%'; label.textContent = labels[Math.min(labels.length - 1, Math.floor(p / 26))]; }
    }, 180);
  }

  buildScan();
  runBtn.addEventListener('click', () => {
    runBtn.disabled = true;
    runBtn.textContent = 'Analyzing…';
    runScan(() => { runBtn.textContent = 'Re-run analysis'; runBtn.disabled = false; renderReport(); });
  });
})();
