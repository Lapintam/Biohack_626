/* Strata — the discovery-agent demo player.
   Loads site/demo/{drug}.json, plays the agent transcript step by step,
   then reveals the figures, top-markers table, and the opportunity brief.
   Vanilla JS, no dependencies. Needs the site served over http (fetch). */

(function () {
  'use strict';

  var KNOWN_CONTROL = { CDKN2A: 1, MTAP: 1 };
  var STEP_DELAY = 750; // ms between revealed transcript steps
  var reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  var els = {
    runBar: document.getElementById('runBar'),
    runBtn: document.getElementById('runBtn'),
    goalLine: document.getElementById('goalLine'),
    transcript: document.getElementById('transcript'),
    results: document.getElementById('resultsSection'),
    figGrid: document.getElementById('figGrid'),
    markerBody: document.querySelector('#markerTable tbody'),
    brief: document.getElementById('brief'),
    scaleGrid: document.getElementById('scaleGrid'),
    scaleProse: document.getElementById('scaleProse')
  };

  var state = { run: 'Palbociclib', data: {}, playing: false, timers: [] };

  // ── helpers ────────────────────────────────────────────────────────
  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
    });
  }
  function fmtRho(x) { return (x >= 0 ? '+' : '') + Number(x).toFixed(3); }
  function fmtQ(x) {
    x = Number(x);
    if (x === 0) return '0';
    if (x >= 0.001) return x.toFixed(3);
    return x.toExponential(1);
  }
  function clearTimers() { state.timers.forEach(clearTimeout); state.timers = []; }

  // ── tiny markdown -> html (headings, bold, code, paragraphs) ─────────
  function renderMarkdown(md) {
    var lines = String(md).replace(/\r\n/g, '\n').split('\n');
    var html = '', para = [];
    function flush() {
      if (para.length) { html += '<p>' + inline(para.join(' ')) + '</p>'; para = []; }
    }
    function inline(t) {
      t = esc(t);
      t = t.replace(/`([^`]+)`/g, '<code>$1</code>');
      t = t.replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>');
      return t;
    }
    for (var i = 0; i < lines.length; i++) {
      var ln = lines[i];
      var h = /^(#{1,3})\s+(.*)$/.exec(ln);
      if (h) { flush(); var lvl = Math.min(h[1].length, 2); html += '<h' + lvl + '>' + inline(h[2]) + '</h' + lvl + '>'; continue; }
      if (!ln.trim()) { flush(); continue; }
      para.push(ln.trim());
    }
    flush();
    return html;
  }

  // ── transcript step renderers ────────────────────────────────────────
  function renderToolCall(step) {
    var args = step.args || {};
    var parts = Object.keys(args).map(function (k) {
      var v = args[k];
      var vs = typeof v === 'string' ? '"' + esc(v) + '"' : esc(v);
      return '<span class="t-arg-k">' + esc(k) + '</span>=<span class="t-arg-v">' + vs + '</span>';
    });
    return '<div class="t-label">tool call</div>' +
      '<div class="t-cmd">$ <span class="t-fn">' + esc(step.name) + '</span>' +
      '<span class="t-paren">(</span>' + parts.join('<span class="t-paren">, </span>') +
      '<span class="t-paren">)</span></div>';
  }

  function renderToolResult(step) {
    var res = step.result || {};
    var body;
    if (step.name === 'rank_markers' && Array.isArray(res.markers)) {
      var chips = res.markers.slice(0, 8).map(function (m) {
        var ctrl = KNOWN_CONTROL[String(m.gene).toUpperCase()] ? ' is-control' : '';
        return '<span class="res-chip' + ctrl + '"><span class="g">' + esc(m.gene) +
          '</span> <span class="m">rho=' + fmtRho(m.rho) + ' q=' + fmtQ(m.qvalue) + '</span></span>';
      }).join('');
      var nSig = res.n_significant != null ? res.n_significant : res.markers.length;
      var nLines = res.n != null ? res.n : null;
      body = '<div class="res-markers">' + chips + '</div>' +
        '<div class="m mono" style="margin-top:10px;font-size:11px;color:var(--text-muted)">' +
          nSig + ' FDR-significant markers' +
          (nLines != null ? ' · n=' + esc(nLines) + ' cell lines' : '') +
        '</div>';
    } else if (step.name === 'evaluate_subgroup') {
      var kv = [
        ['marker', res.marker], ['grade', res.grade], ['p-value', res.pvalue != null ? fmtQ(res.pvalue) : null],
        ['effect (ΔAUC)', res.effect_size != null ? Number(res.effect_size).toFixed(3) : null],
        ['n', res.n], ['responder cluster', res.responder_label]
      ].filter(function (r) { return r[1] != null; });
      body = '<div class="res-kv">' + kv.map(function (r) {
        return '<span class="rk">' + esc(r[0]) + '</span><span class="rv">' + esc(r[1]) + '</span>';
      }).join('') + '</div>';
    } else {
      body = '<div class="mono" style="font-size:12px;color:var(--text-muted)">' + esc(JSON.stringify(res)) + '</div>';
    }
    return '<div class="t-label">tool result · ' + esc(step.name) + '</div>' +
      '<div class="t-body">' + body + '</div>';
  }

  function renderThought(step, withCursor) {
    return '<div class="t-label">agent thought</div>' +
      '<div class="t-body">' + esc(step.text) +
      (withCursor ? '<span class="t-cursor"></span>' : '') + '</div>';
  }

  function stepEl(step, isLast) {
    var div = document.createElement('div');
    var cls = step.type === 'tool_call' ? 't-call' : step.type === 'tool_result' ? 't-result' : 't-thought';
    div.className = 't-step ' + cls;
    if (step.type === 'tool_call') div.innerHTML = renderToolCall(step);
    else if (step.type === 'tool_result') div.innerHTML = renderToolResult(step);
    else div.innerHTML = renderThought(step, isLast);
    return div;
  }

  // ── results (figures, markers, brief) ────────────────────────────────
  var FIG_CAPS = {
    umap: 'UMAP · methylome embedding (responder stratum)',
    auc: 'Drug response · AUC by companion-marker split',
    heatmap: 'Top markers · methylation heatmap by stratum'
  };
  function renderResults(data) {
    // figures
    els.figGrid.innerHTML = '';
    ['umap', 'auc', 'heatmap'].forEach(function (k) {
      var fn = data.figures && data.figures[k];
      if (!fn) return;
      var f = document.createElement('div');
      f.className = 'fig';
      f.innerHTML = '<div class="fig-cap">' + esc(FIG_CAPS[k]) + '</div>' +
        '<img src="demo/' + esc(fn) + '" alt="' + esc(FIG_CAPS[k]) + '" loading="lazy" />';
      els.figGrid.appendChild(f);
    });
    // markers
    els.markerBody.innerHTML = '';
    (data.top_markers || []).forEach(function (m) {
      var ctrl = KNOWN_CONTROL[String(m.gene).toUpperCase()];
      var tr = document.createElement('tr');
      if (ctrl) tr.className = 'is-control';
      var negCls = m.rho < 0 ? ' neg' : '';
      tr.innerHTML = '<td class="gene">' + esc(m.gene) +
        (ctrl ? ' <span class="ctrl-flag">◆ control</span>' : '') + '</td>' +
        '<td class="num' + negCls + '">' + fmtRho(m.rho) + '</td>' +
        '<td class="num">' + fmtQ(m.qvalue) + '</td>';
      els.markerBody.appendChild(tr);
    });
    // control_check rigor line (Palbociclib only)
    var rigorEl = document.getElementById('controlRigor');
    if (!rigorEl) {
      rigorEl = document.createElement('div');
      rigorEl.id = 'controlRigor';
      rigorEl.className = 'control-rigor';
      var markerWrap = document.querySelector('.marker-table-wrap');
      if (markerWrap) markerWrap.parentNode.insertBefore(rigorEl, markerWrap.nextSibling);
    }
    var cc = data.control_check;
    if (cc) {
      var rhoStr = (cc.rho >= 0 ? '+' : '') + Number(cc.rho).toFixed(3);
      var qStr = Number(cc.qvalue).toExponential(1);
      rigorEl.innerHTML = '<span class="rigor-label">Canonical-biomarker check</span> ' +
        'CDKN2A itself is recovered at rank ' + esc(cc.rank) + '/' + esc(cc.total) +
        ' (ρ=' + esc(rhoStr) + ', q=' + esc(qStr) + ') — ' +
        'the #1 marker MTAP is its 9p21 co-deletion partner.';
      rigorEl.hidden = false;
    } else {
      rigorEl.hidden = true;
    }
    // brief
    els.brief.innerHTML = renderMarkdown(data.brief || '');
    els.results.hidden = false;
  }

  // ── scale / NVIDIA section (loads demo/scale.json) ───────────────────
  function commas(n) { return Math.round(Number(n)).toLocaleString('en-US'); }
  function compact(n) {
    n = Number(n);
    if (n >= 1e9) return (n / 1e9).toFixed(1).replace(/\.0$/, '') + 'B';
    if (n >= 1e6) return (n / 1e6).toFixed(1).replace(/\.0$/, '') + 'M';
    if (n >= 1e3) return (n / 1e3).toFixed(0) + 'k';
    return String(Math.round(n));
  }
  function renderScale() {
    if (!els.scaleGrid) return;
    fetch('demo/scale.json')
      .then(function (r) { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); })
      .then(function (s) {
        var cur = s.current, bio = s.biobank_projection;
        var stats = [
          [compact(cur.full_panel_correlations), 'correlations · current GDSC panel'],
          [cur.full_panel_seconds_cpu + ' s', 'one full pass on CPU'],
          ['~' + commas(bio.cpu_days) + ' CPU-days', 'one pass at biobank scale'],
          ['≈' + compact(bio.scale_factor_vs_gdsc) + '×', 'more compute than today — fully parallel → GPU']
        ];
        els.scaleGrid.innerHTML = stats.map(function (st) {
          return '<div class="scale-stat"><span class="ss-n">' + esc(st[0]) +
            '</span><span class="ss-l">' + esc(st[1]) + '</span></div>';
        }).join('');
        var a = bio.assumption;
        els.scaleProse.innerHTML =
          'At GDSC scale the screen is <strong>' + commas(cur.full_panel_correlations) +
          ' correlations in ' + cur.full_panel_seconds_cpu + ' s</strong> on CPU (' +
          esc(cur.microseconds_per_correlation) + ' µs each). Project to a biobank — ' +
          commas(a.genes) + ' genes × ' + commas(a.drugs_or_outcomes) + ' drugs/outcomes × ' +
          compact(a.cohort_n) + ' patients — and a single pass is <strong>~' +
          commas(bio.cpu_days) + ' CPU-days</strong>, about ' + compact(bio.scale_factor_vs_gdsc) +
          '× today. The work is independent per (gene × drug × cohort): the exact shape ' +
          'cuDF / cuML accelerate. We don’t <em>want</em> GPUs — the unit economics require them.';
      })
      .catch(function () {
        els.scaleGrid.innerHTML = '<div class="scale-stat"><span class="ss-l">Scale figures load when the site is served over http.</span></div>';
      });
  }

  // ── playback ─────────────────────────────────────────────────────────
  function play(data) {
    clearTimers();
    state.playing = true;
    els.runBtn.disabled = true;
    els.runBtn.textContent = '… running';
    els.transcript.innerHTML = '';
    els.results.hidden = true;

    var steps = data.transcript || [];
    var delay = reduceMotion ? 0 : STEP_DELAY;

    steps.forEach(function (step, i) {
      var t = setTimeout(function () {
        // drop cursor from previous thought
        var prevCur = els.transcript.querySelector('.t-cursor');
        if (prevCur) prevCur.remove();
        var isLastThought = step.type === 'thought' && i === steps.length - 1;
        els.transcript.appendChild(stepEl(step, isLastThought && !reduceMotion));
        if (i === steps.length - 1) {
          renderResults(data);
          state.playing = false;
          els.runBtn.disabled = false;
          els.runBtn.textContent = '▶ Replay discovery';
          els.results.scrollIntoView({ behavior: reduceMotion ? 'auto' : 'smooth', block: 'start' });
        }
      }, delay * i);
      state.timers.push(t);
    });

    if (steps.length === 0) {
      state.playing = false;
      els.runBtn.disabled = false;
      els.runBtn.textContent = '▶ Run discovery';
    }
  }

  // ── data loading ─────────────────────────────────────────────────────
  function loadRun(drug) {
    if (state.data[drug]) return Promise.resolve(state.data[drug]);
    return fetch('demo/' + encodeURIComponent(drug) + '.json')
      .then(function (r) { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); })
      .then(function (d) { state.data[drug] = d; return d; });
  }

  function setGoal(data) {
    els.goalLine.innerHTML = '<span class="gk">Goal</span>' + esc(data.goal || '');
  }

  function selectRun(drug, autoplay) {
    state.run = drug;
    clearTimers();
    Array.prototype.forEach.call(els.runBar.querySelectorAll('.strata-runtab'), function (b) {
      var on = b.getAttribute('data-run') === drug;
      b.classList.toggle('is-active', on);
      b.setAttribute('aria-selected', on ? 'true' : 'false');
    });
    els.transcript.innerHTML = '';
    els.results.hidden = true;
    els.runBtn.disabled = false;
    els.runBtn.textContent = '▶ Run discovery';

    loadRun(drug).then(function (data) {
      setGoal(data);
      if (autoplay) play(data);
    }).catch(function (err) {
      els.goalLine.innerHTML = '<span class="gk">Error</span>Could not load demo/' +
        esc(drug) + '.json — serve the site over http (it uses fetch). ' + esc(err.message);
    });
  }

  // ── wiring ───────────────────────────────────────────────────────────
  els.runBar.addEventListener('click', function (e) {
    var btn = e.target.closest('.strata-runtab');
    if (!btn) return;
    selectRun(btn.getAttribute('data-run'), false);
  });
  els.runBtn.addEventListener('click', function () {
    if (state.playing) return;
    loadRun(state.run).then(play);
  });

  selectRun('Palbociclib', false);
  renderScale();
})();
