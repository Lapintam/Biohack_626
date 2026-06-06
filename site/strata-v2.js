/* Strata v2 — cross-indication opportunity engine demo player.
   Loads demo-v2/opportunities.json, plays the agent transcript step by
   step, then reveals §03-§06. Vanilla JS, no dependencies. Needs http. */

(function () {
  'use strict';

  var DATA_URL = 'demo-v2/opportunities.json';
  var STEP_DELAY = 750; // ms between transcript steps
  var reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  var els = {
    runBtn:       document.getElementById('runBtn'),
    thesisLede:   document.getElementById('thesisLede'),
    noveltyBox:   document.getElementById('noveltyBox'),
    noveltyText:  document.getElementById('noveltyText'),
    transcript:   document.getElementById('transcript'),
    featured:     document.getElementById('featuredSection'),
    negSection:   document.getElementById('negSection'),
    tableSection: document.getElementById('tableSection'),
    ladderSection:document.getElementById('ladderSection'),
    oppCards:     document.getElementById('oppCards'),
    negControl:   document.getElementById('negControl'),
    oppTable:     document.getElementById('oppTable'),
    ladder:       document.getElementById('evidenceLadder')
  };

  var state = { data: null, playing: false, timers: [] };

  // ── helpers ────────────────────────────────────────────────────────
  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
    });
  }
  function fmtNum(x, places) {
    var n = Number(x);
    var s = (n >= 0 ? '+' : '') + n.toFixed(places != null ? places : 3);
    return s;
  }
  function clearTimers() { state.timers.forEach(clearTimeout); state.timers = []; }

  // ── tiny markdown → html (headings, bold, code, paragraphs) ───────
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
      if (h) {
        flush();
        var lvl = Math.min(h[1].length, 2);
        html += '<h' + lvl + '>' + inline(h[2]) + '</h' + lvl + '>';
        continue;
      }
      if (!ln.trim()) { flush(); continue; }
      para.push(ln.trim());
    }
    flush();
    return html;
  }

  // ── transcript step renderers (same pattern as v1 strata.js) ──────
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
    var body = '<div class="mono" style="font-size:12px;color:var(--text-muted)">' +
      esc(JSON.stringify(res)) + '</div>';
    return '<div class="t-label">tool result &middot; ' + esc(step.name) + '</div>' +
      '<div class="t-body">' + body + '</div>';
  }

  function renderThought(step, withCursor) {
    var isMarkdown = step.text && /^#/.test(step.text.trimStart());
    var body = isMarkdown
      ? '<div class="brief" style="border:none;padding:0;background:transparent;">' +
          renderMarkdown(step.text) + '</div>'
      : esc(step.text) + (withCursor ? '<span class="t-cursor"></span>' : '');
    return '<div class="t-label">agent thought</div>' +
      '<div class="t-body">' + body + '</div>';
  }

  function stepEl(step, isLast) {
    var div = document.createElement('div');
    var cls = step.type === 'tool_call' ? 't-call'
            : step.type === 'tool_result' ? 't-result'
            : 't-thought';
    div.className = 't-step ' + cls;
    if (step.type === 'tool_call')       div.innerHTML = renderToolCall(step);
    else if (step.type === 'tool_result') div.innerHTML = renderToolResult(step);
    else                                  div.innerHTML = renderThought(step, isLast && !reduceMotion);
    return div;
  }

  // ── tissue recurrence bars ─────────────────────────────────────────
  function buildTissuePanel(perTissue) {
    if (!perTissue || !perTissue.length) return '';
    var maxAbs = 0;
    perTissue.forEach(function (t) { var a = Math.abs(t.rho); if (a > maxAbs) maxAbs = a; });
    if (maxAbs === 0) maxAbs = 1;
    var rows = perTissue.map(function (t) {
      var pct = Math.round((Math.abs(t.rho) / maxAbs) * 100);
      var barCls = t.rho < 0 ? 'neg-bar' : 'pos-bar';
      var sigCls = t.sig ? ' is-sig' : '';
      return '<div class="opp-tissue-row">' +
        '<span class="opp-tissue-name" title="' + esc(t.tissue) + '">' + esc(t.tissue) + '</span>' +
        '<span class="opp-tissue-sig' + sigCls + '">' + (t.sig ? '&#10003;' : '') + '</span>' +
        '<div class="opp-tissue-bar-wrap" style="flex:1;background:var(--bg-deep);border-radius:1px;" title="' +
          'n=' + esc(t.n) + ' &rho;=' + Number(t.rho).toFixed(3) + '">' +
          '<span class="opp-tissue-bar ' + barCls + '" style="width:' + pct + '%;display:block;"></span>' +
        '</div>' +
      '</div>';
    }).join('');
    return '<div class="opp-tissue-head">Tissue recurrence &mdash; &rho; by cancer type (neg = sensitive, &#10003; = sig)</div>' +
      '<div class="opp-tissue-list">' + rows + '</div>';
  }

  // ── single opportunity card ─────────────────────────────────────────
  function buildCard(card, tagHtml) {
    var rAdj = Number(card.r_adj);
    var rhoP = Number(card.rho_pooled);
    var mc = card.mechanism_card || {};
    var reps = card.replicated_by || [];

    var repChips = reps.length
      ? '<div class="opp-replicated">' +
          '<div class="opp-replicated-lbl">also recovered in</div>' +
          '<div class="opp-rep-chips">' +
          reps.map(function (r) {
            return '<span class="opp-rep-chip"><strong>' + esc(r.drug) + '</strong> &middot; ' + esc(r.note) + '</span>';
          }).join('') +
          '</div></div>'
      : '';

    return '<div class="opp-card">' +

      // head
      '<div class="opp-card-head">' +
        '<div class="opp-card-title">' +
          '<span class="opp-card-drug">' + esc(card.drug) + '</span>' +
          '<span class="opp-card-arrow">&rarr;</span>' +
          '<span class="opp-card-marker">' + esc(card.marker) + '</span>' +
          '<span class="opp-level-badge">' + esc(card.level) + '</span>' +
        '</div>' +
        tagHtml +
        '<div class="opp-card-pathway">' + esc(card.pathway) + '</div>' +
      '</div>' +

      // stats
      '<div class="opp-stats">' +
        '<div class="opp-stat">' +
          '<span class="opp-stat-val' + (rAdj < 0 ? ' neg' : '') + '">' + fmtNum(rAdj, 3) + '</span>' +
          '<span class="opp-stat-lbl">tissue-adjusted r</span>' +
        '</div>' +
        '<div class="opp-stat">' +
          '<span class="opp-stat-val' + (rhoP < 0 ? ' neg' : '') + '">' + fmtNum(rhoP, 3) + '</span>' +
          '<span class="opp-stat-lbl">pooled &rho;</span>' +
        '</div>' +
        '<div class="opp-stat">' +
          '<span class="opp-stat-val">' + esc(card.within_sig) + ' / ' + esc(card.n_tissues) + '</span>' +
          '<span class="opp-stat-lbl">sig within tissues</span>' +
        '</div>' +
        '<div class="opp-stat">' +
          '<span class="opp-stat-val">' + esc(card.n) + '</span>' +
          '<span class="opp-stat-lbl">cell lines</span>' +
        '</div>' +
      '</div>' +

      // tissue panel
      buildTissuePanel(card.per_tissue) +

      // mechanism card
      '<div class="opp-mech">' +
        '<div class="opp-mech-state">' + esc(mc.state || '') + '</div>' +
        '<div class="opp-mech-kv">' +
          '<div><span class="opp-mech-k">Known biology</span>' + esc(mc.known_biology || '') + '</div>' +
          '<div><span class="opp-mech-k">Implication</span>' + esc(mc.implication || '') + '</div>' +
        '</div>' +
        '<div class="opp-level-note">' + esc(card.level_note || '') + '</div>' +
      '</div>' +

      repChips +

    '</div>';
  }

  // ── render §03 featured cards ──────────────────────────────────────
  function renderFeatured(data) {
    var f = data.featured || {};
    var html = '';

    if (f.positive_control) {
      html += buildCard(f.positive_control,
        '<span class="opp-card-tag opp-tag-control">Positive control &middot; confirmation</span>');
    }
    (f.novel_leads || []).forEach(function (lead) {
      html += buildCard(lead,
        '<span class="opp-card-tag opp-tag-novel">Novel lead &middot; hypothesis</span>');
    });

    els.oppCards.innerHTML = html;
    showReveal(els.featured);
  }

  // ── render §04 negative control ───────────────────────────────────
  function renderNeg(data) {
    var nc = (data.featured || {}).negative_control || {};
    els.negControl.innerHTML =
      '<span class="neg-control-label">Negative control &mdash; the engine isn\'t forcing results</span>' +
      '<p>' + esc(nc.note || '') + '</p>';
    showReveal(els.negSection);
  }

  // ── render §05 opportunity table ──────────────────────────────────
  function renderTable(data) {
    var tbody = els.oppTable.querySelector('tbody');
    tbody.innerHTML = '';
    (data.opportunity_table || []).forEach(function (row) {
      var rAdj = Number(row.r_adj);
      var rhoP = Number(row.rho_pooled);
      var tr = document.createElement('tr');
      tr.innerHTML =
        '<td class="td-drug">' + esc(row.drug) + '</td>' +
        '<td style="color:var(--text-tertiary);font-size:12px;">' + esc(row.pathway) + '</td>' +
        '<td class="td-marker">' + esc(row.marker) + '</td>' +
        '<td><span class="opp-level-badge" style="font-size:9px;">' + esc(row.level) + '</span></td>' +
        '<td class="num' + (rAdj < 0 ? ' neg' : ' pos') + '">' + fmtNum(rAdj, 3) + '</td>' +
        '<td class="num' + (rhoP < 0 ? ' neg' : ' pos') + '">' + fmtNum(rhoP, 3) + '</td>' +
        '<td class="num">' + esc(row.within_sig) + '</td>';
      tbody.appendChild(tr);
    });
    showReveal(els.tableSection);
  }

  // ── render §06 evidence ladder ─────────────────────────────────────
  function renderLadder(data) {
    var rungs = data.evidence_ladder || [];
    // v2 targets L2-L3; mark those as active
    var active = { L2: true, L3: true };
    els.ladder.innerHTML = rungs.map(function (r) {
      var cls = active[r.level] ? ' is-active' : '';
      return '<div class="ev-rung' + cls + '">' +
        '<span class="ev-rung-level">' + esc(r.level) + '</span>' +
        '<span class="ev-rung-label">' + esc(r.label) + '</span>' +
      '</div>';
    }).join('');
    showReveal(els.ladderSection);
  }

  // ── show a hidden reveal section ───────────────────────────────────
  function showReveal(el) {
    if (!el) return;
    el.hidden = false;
    el.style.opacity = '0';
    requestAnimationFrame(function () {
      requestAnimationFrame(function () {
        el.style.opacity = '1';
      });
    });
  }

  // ── playback ───────────────────────────────────────────────────────
  function play(data) {
    clearTimers();
    state.playing = true;
    els.runBtn.disabled = true;
    els.runBtn.textContent = '… running';
    els.transcript.innerHTML = '';

    // hide reveal sections
    [els.featured, els.negSection, els.tableSection, els.ladderSection].forEach(function (s) {
      if (s) { s.hidden = true; s.style.opacity = '0'; }
    });

    var steps = data.transcript || [];
    var delay = reduceMotion ? 0 : STEP_DELAY;

    steps.forEach(function (step, i) {
      var t = setTimeout(function () {
        // remove cursor from previous thought
        var prevCur = els.transcript.querySelector('.t-cursor');
        if (prevCur) prevCur.remove();

        var isLastThought = step.type === 'thought' && i === steps.length - 1;
        els.transcript.appendChild(stepEl(step, isLastThought));

        if (i === steps.length - 1) {
          // reveal all result sections
          renderFeatured(data);
          renderNeg(data);
          renderTable(data);
          renderLadder(data);

          state.playing = false;
          els.runBtn.disabled = false;
          els.runBtn.textContent = '&#9654; Replay';

          if (els.featured) {
            els.featured.scrollIntoView({
              behavior: reduceMotion ? 'auto' : 'smooth',
              block: 'start'
            });
          }
        }
      }, delay * i);
      state.timers.push(t);
    });

    if (steps.length === 0) {
      state.playing = false;
      els.runBtn.disabled = false;
      els.runBtn.textContent = '&#9654; Run discovery';
    }
  }

  // ── data loading ───────────────────────────────────────────────────
  function loadData() {
    if (state.data) return Promise.resolve(state.data);
    return fetch(DATA_URL)
      .then(function (r) {
        if (!r.ok) throw new Error('HTTP ' + r.status);
        return r.json();
      })
      .then(function (d) {
        state.data = d;
        return d;
      });
  }

  function applyStaticContent(data) {
    // thesis lede
    if (data.thesis && els.thesisLede) {
      els.thesisLede.textContent = data.thesis;
    }
    // novelty box
    if (data.novelty && els.noveltyText && els.noveltyBox) {
      els.noveltyText.textContent = data.novelty;
      els.noveltyBox.hidden = false;
    }
  }

  // ── wiring ─────────────────────────────────────────────────────────
  els.runBtn.addEventListener('click', function () {
    if (state.playing) return;
    loadData()
      .then(function (data) {
        applyStaticContent(data);
        play(data);
      })
      .catch(function (err) {
        els.transcript.innerHTML =
          '<div class="t-step t-thought" style="opacity:1;transform:none;">' +
          '<div class="t-label">error</div>' +
          '<div class="t-body" style="color:var(--methyl-rose);">Could not load ' +
          esc(DATA_URL) + ' — serve the site over http (it uses fetch). ' +
          esc(err.message) + '</div></div>';
      });
  });

  // pre-load static content on page open (graceful — no autoplay)
  loadData()
    .then(applyStaticContent)
    .catch(function () { /* fail silently on static open */ });

})();
