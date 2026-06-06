/* Strata v3 — three structural gates.
   Loads demo-v3/cards.json and renders the static pitch sections:
   thesis, three gates, discovery funnel, hero cards, honesty rails.
   Vanilla JS, no dependencies. Needs http (uses fetch). */

(function () {
  'use strict';

  var DATA_URL = 'demo-v3/cards.json';

  var els = {
    thesisLede:  document.getElementById('thesisLede'),
    gatesGrid:   document.getElementById('gatesGrid'),
    funnel:      document.getElementById('funnel'),
    heroCards:   document.getElementById('heroCards'),
    honestyRails:document.getElementById('honestyRails')
  };

  // ── helpers (same shape as strata-v2.js) ───────────────────────────
  function esc(s) {
    return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
    });
  }
  function fmtR(x) {
    if (x == null || isNaN(Number(x))) return '—';
    var n = Number(x);
    return (n >= 0 ? '+' : '') + n.toFixed(3);
  }
  function fmtP(p) {
    if (p == null || isNaN(Number(p))) return '';
    var n = Number(p);
    if (n === 0) return 'p<0.001';
    if (n < 0.001) return 'p<0.001';
    return 'p=' + n.toFixed(3);
  }

  // ── §02 three gates ─────────────────────────────────────────────────
  function renderGates(gates) {
    if (!els.gatesGrid) return;
    els.gatesGrid.innerHTML = (gates || []).map(function (g, i) {
      return '<div class="gate">' +
        '<div class="gate-head">' +
          '<span class="gate-num">Gate ' + (i + 1) + ' &middot; ' + esc(g.id) + '</span>' +
          '<div class="gate-title">' + esc(g.title) + '</div>' +
        '</div>' +
        '<div class="gate-body">' +
          '<div class="gate-side gate-side-expert">' +
            '<span class="gate-role">An expert</span>' +
            '<div class="gate-text">' + esc(g.expert) + '</div>' +
          '</div>' +
          '<div class="gate-side gate-side-agent">' +
            '<span class="gate-role">The agent</span>' +
            '<div class="gate-text">' + esc(g.agent) + '</div>' +
          '</div>' +
        '</div>' +
        '<div class="gate-evidence">' +
          '<span class="gate-ev-lbl">Evidence</span>' + esc(g.evidence) +
        '</div>' +
      '</div>';
    }).join('');
  }

  // ── §03 discovery funnel ────────────────────────────────────────────
  function renderFunnel(funnel) {
    if (!els.funnel) return;
    var stages = funnel || [];
    var n = stages.length;
    // funnel narrows visually: width steps down per stage
    var html = stages.map(function (s, i) {
      var w = Math.round(100 - (i * (45 / Math.max(1, n - 1))));
      var arrow = i < n - 1
        ? '<div class="funnel-arrow">&#9660;</div>'
        : '';
      return '<div class="funnel-stage" style="width:' + w + '%;">' +
          '<span class="funnel-idx">' + (i + 1) + '</span>' +
          '<div class="funnel-main">' +
            '<span class="funnel-label">' + esc(s.stage) + '</span>' +
            '<span class="funnel-note">' + esc(s.note) + '</span>' +
          '</div>' +
          '<span class="funnel-n">' + esc(s.n) + '</span>' +
        '</div>' + arrow;
    }).join('');
    els.funnel.innerHTML = html;
  }

  // ── replication strip for one card ─────────────────────────────────
  function repStrip(card) {
    var screens = card.screens || [];
    var strong = card.n_screens >= 3;
    // map screen name → {r, p}
    var byScreen = {
      GDSC: { r: card.gdsc_r, p: null },
      PRISM: { r: card.prism_r, p: card.prism_p },
      CTRP: { r: card.ctrp_r, p: card.ctrp_p }
    };
    var cells = screens.map(function (name) {
      var d = byScreen[name] || { r: null, p: null };
      if (d.r == null) return ''; // skip null screens (e.g. missing CTRP)
      var sign = Number(d.r) < 0 ? 'neg' : 'pos';
      var pStr = name === 'GDSC' ? 'discovery' : fmtP(d.p);
      return '<div class="opp-screen">' +
          '<span class="opp-screen-name">' + esc(name) + '</span>' +
          '<span class="opp-screen-r ' + sign + '">' + fmtR(d.r) + '</span>' +
          '<span class="opp-screen-p">' + esc(pStr) + '</span>' +
        '</div>';
    }).join('');

    var badge = strong
      ? '<span class="opp-rep-strength strong">3-screen · GDSC+PRISM+CTRP</span>'
      : '<span class="opp-rep-strength partial">' + esc(card.n_screens) + '-screen</span>';

    return '<div class="opp-rep-strip">' +
        '<div class="opp-rep-strip-head">' +
          '<span class="opp-rep-strip-lbl">Replication · Spearman r per screen</span>' +
          badge +
        '</div>' +
        '<div class="opp-rep-screens">' + cells + '</div>' +
      '</div>';
  }

  // ── recurrence bar strip ───────────────────────────────────────────
  function recurBars(rec) {
    var bt = (rec && rec.by_tissue) || [];
    if (!bt.length) return '';
    var maxAbs = 0;
    bt.forEach(function (t) { var a = Math.abs(Number(t.rho)); if (a > maxAbs) maxAbs = a; });
    if (maxAbs === 0) maxAbs = 1;
    var bars = bt.map(function (t) {
      var rho = Number(t.rho);
      var h = Math.max(8, Math.round((Math.abs(rho) / maxAbs) * 100));
      var sign = rho < 0 ? '' : ' pos';
      var sig = (t.p != null && Number(t.p) < 0.05) ? ' sig' : '';
      var title = esc(t.tissue) + ' · n=' + esc(t.n) + ' · ρ=' + Number(rho).toFixed(3) +
        (t.p != null ? ' · ' + fmtP(t.p) : '');
      return '<span class="opp-recur-bar' + sign + sig + '" style="height:' + h + '%;" title="' + title + '"></span>';
    }).join('');
    return '<div class="opp-recur-bars">' + bars + '</div>';
  }

  // ── one hero card (reuses .opp-card chrome) ─────────────────────────
  function buildHeroCard(card, isLead) {
    var dir = card.direction || '';
    var rec = card.recurrence || {};
    var note = card.gene_note && card.gene_note.trim()
      ? card.gene_note
      : 'No prior link between this gene and ' + card.drug + ' response — outside any domain prior.';

    return '<div class="opp-card' + (isLead ? ' is-lead' : '') + '">' +

      // head: drug → gene, evidence badge, meth tag
      '<div class="opp-card-head">' +
        '<div class="opp-card-title">' +
          '<span class="opp-card-drug">' + esc(card.drug) + '</span>' +
          '<span class="opp-card-arrow">&rarr;</span>' +
          '<span class="opp-card-marker">' + esc(card.gene) + '</span>' +
          '<span class="opp-level-badge">' + esc(card.evidence_level) + '</span>' +
        '</div>' +
        '<span class="opp-meth-tag ' + esc(dir) + '">meth &rarr; ' + esc(dir) + '</span>' +
      '</div>' +

      // gene note — the punchline
      '<div class="opp-gene-note">' +
        '<span class="opp-gene-note-lbl">Why no expert would prioritize it</span>' +
        esc(note) +
      '</div>' +

      // replication strip
      repStrip(card) +

      // functional silencing
      '<div class="opp-silencing">' +
        '<span class="opp-silencing-lbl">promoter methylation ' +
          '<span class="arrow">&darr;</span> expression</span>' +
        '<span class="opp-silencing-val">' + fmtR(card.silencing_r) + '</span>' +
      '</div>' +

      // cross-indication recurrence
      '<div class="opp-recurrence">' +
        '<div class="opp-recurrence-lbl">' +
          'Same direction in <strong>' + esc(rec.same_direction) + '/' + esc(rec.tissues_tested) +
          '</strong> cancer types (' + esc(rec.significant) + ' significant)' +
        '</div>' +
        recurBars(rec) +
      '</div>' +

    '</div>';
  }

  function renderHeroCards(cards) {
    if (!els.heroCards) return;
    els.heroCards.innerHTML = (cards || []).map(function (c, i) {
      return buildHeroCard(c, i === 0);
    }).join('');
  }

  // ── §05 honesty rails ───────────────────────────────────────────────
  function renderHonesty(honesty) {
    if (!els.honestyRails) return;
    els.honestyRails.innerHTML = (honesty || []).map(function (h, i) {
      return '<div class="honesty-rail">' +
          '<span class="hr-n">0' + (i + 1) + '</span>' +
          '<div>' + esc(h) + '</div>' +
        '</div>';
    }).join('');
  }

  // ── data loading ─────────────────────────────────────────────────────
  function render(data) {
    if (data.thesis && els.thesisLede) els.thesisLede.textContent = data.thesis;
    renderGates(data.gates);
    renderFunnel(data.funnel);
    renderHeroCards(data.hero_cards);
    renderHonesty(data.honesty);
  }

  function renderError(msg) {
    var target = els.gatesGrid;
    if (!target) return;
    target.innerHTML =
      '<div class="neg-control" style="grid-column:1/-1;margin-top:0;">' +
      '<span class="neg-control-label">error</span>' +
      '<p>Could not load ' + esc(DATA_URL) + ' — serve the site over http (it uses fetch). ' +
      esc(msg) + '</p></div>';
  }

  fetch(DATA_URL)
    .then(function (r) {
      if (!r.ok) throw new Error('HTTP ' + r.status);
      return r.json();
    })
    .then(render)
    .catch(function (err) { renderError(err.message); });

})();
