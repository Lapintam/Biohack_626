/* ──────────────────────────────────────────────────────────────────
   Polymer tutorial — interaction layer.
   - lightweight Python syntax highlighting for .nb-code blocks
   - per-cell "hide code" + "copy" toggles
   - sticky stage rail: scrollspy active state + click-to-scroll
   - "Run tutorial" guided walkthrough (sequential scroll + highlight)
   - on-scroll reveal of notebook cells
   ────────────────────────────────────────────────────────────────── */
(function () {
  "use strict";
  var reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  /* ── Python syntax highlighter ────────────────────────────────── */
  var PY_KW = new RegExp(
    "\\b(?:import|from|as|def|class|return|if|elif|else|for|while|in|not|and|or|" +
    "is|None|True|False|with|try|except|finally|lambda|global|nonlocal|pass|" +
    "break|continue|yield|raise|assert|del|self)\\b", "g");
  var PY_BI = new RegExp(
    "\\b(?:print|range|len|set|sorted|map|int|float|str|list|dict|tuple|enumerate|" +
    "zip|abs|min|max|sum|super|isinstance|np|pd|plt|torch|nn)\\b", "g");

  function esc(s) {
    return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

  // Tokenise a single physical line; comments/strings take priority.
  function highlightLine(line) {
    var out = "";
    var i = 0;
    while (i < line.length) {
      var ch = line[i];
      // comment to end of line
      if (ch === "#") {
        out += '<span class="tok-cmt">' + esc(line.slice(i)) + "</span>";
        break;
      }
      // string literal (single or double, naive — no triple/escapes needed here)
      if (ch === '"' || ch === "'") {
        var q = ch, j = i + 1;
        while (j < line.length && line[j] !== q) j++;
        var lit = line.slice(i, Math.min(j + 1, line.length));
        out += '<span class="tok-str">' + esc(lit) + "</span>";
        i = j + 1;
        continue;
      }
      // run of "other" chars until next special starter
      var k = i;
      while (k < line.length && line[k] !== "#" && line[k] !== '"' && line[k] !== "'") k++;
      var seg = line.slice(i, k);
      seg = esc(seg)
        .replace(PY_KW, '<span class="tok-kw">$&</span>')
        .replace(PY_BI, '<span class="tok-bi">$&</span>')
        .replace(/\b(\d+\.?\d*(?:e-?\d+)?)\b/g, '<span class="tok-num">$&</span>')
        // function call name immediately before "("
        .replace(/\b([A-Za-z_]\w*)(?=\s*\()/g, function (m, name) {
          if (/tok-/.test(name)) return m;
          return '<span class="tok-fn">' + name + "</span>";
        });
      out += seg;
      i = k;
    }
    return out;
  }

  function highlight(codeEl) {
    var raw = codeEl.textContent;
    var lines = raw.split("\n");
    codeEl.innerHTML = lines.map(highlightLine).join("\n");
  }

  /* ── Wire up notebook cells ───────────────────────────────────── */
  var cells = Array.prototype.slice.call(document.querySelectorAll(".nb-cell"));
  cells.forEach(function (cell) {
    var codeEl = cell.querySelector(".nb-code code");
    if (codeEl) highlight(codeEl);

    var toggle = cell.querySelector('[data-act="toggle"]');
    if (toggle) {
      toggle.addEventListener("click", function () {
        var hidden = cell.classList.toggle("code-hidden");
        toggle.textContent = hidden ? "show code" : "hide code";
      });
    }
    var copy = cell.querySelector('[data-act="copy"]');
    if (copy && codeEl) {
      copy.addEventListener("click", function () {
        navigator.clipboard && navigator.clipboard.writeText(codeEl.textContent).then(function () {
          var prev = copy.textContent;
          copy.textContent = "copied";
          setTimeout(function () { copy.textContent = prev; }, 1200);
        });
      });
    }
  });

  /* ── On-scroll reveal ─────────────────────────────────────────── */
  if (!reduceMotion && "IntersectionObserver" in window) {
    var revObs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) { e.target.classList.add("is-in"); revObs.unobserve(e.target); }
      });
    }, { threshold: 0.08 });
    cells.forEach(function (c) { c.classList.add("reveal"); revObs.observe(c); });
  }

  /* ── Stage rail: scrollspy + click-to-scroll ──────────────────── */
  var chips = Array.prototype.slice.call(document.querySelectorAll(".sr-chip"));
  var stages = chips.map(function (c) { return document.getElementById(c.dataset.target); }).filter(Boolean);
  var progressBar = document.querySelector(".sr-progress");

  function setActive(idx) {
    chips.forEach(function (c, i) {
      c.classList.toggle("is-active", i === idx);
      c.classList.toggle("is-done", i < idx);
    });
  }

  chips.forEach(function (chip) {
    chip.addEventListener("click", function () {
      var el = document.getElementById(chip.dataset.target);
      if (el) el.scrollIntoView({ behavior: reduceMotion ? "auto" : "smooth", block: "start" });
    });
  });

  function onScroll() {
    var mid = window.scrollY + window.innerHeight * 0.32;
    var active = 0;
    stages.forEach(function (s, i) { if (s.offsetTop <= mid) active = i; });
    setActive(active);
    var doc = document.documentElement;
    var scrollable = doc.scrollHeight - window.innerHeight;
    if (progressBar) progressBar.style.width = (scrollable > 0 ? (window.scrollY / scrollable) * 100 : 0) + "%";
  }
  window.addEventListener("scroll", onScroll, { passive: true });
  window.addEventListener("resize", onScroll, { passive: true });
  onScroll();

  /* ── "Run tutorial" guided walkthrough ────────────────────────── */
  var runBtn = document.getElementById("runTutorial");
  var running = false;
  if (runBtn) {
    runBtn.addEventListener("click", function () {
      if (running || !stages.length) return;
      running = true;
      runBtn.disabled = true;
      var origLabel = runBtn.textContent;
      var i = 0;
      var step = function () {
        if (i >= stages.length) {
          running = false;
          runBtn.disabled = false;
          runBtn.textContent = origLabel;
          return;
        }
        setActive(i);
        runBtn.textContent = "▶ Stage " + (i + 1) + " / " + stages.length;
        stages[i].scrollIntoView({ behavior: reduceMotion ? "auto" : "smooth", block: "start" });
        i++;
        setTimeout(step, reduceMotion ? 350 : 2200);
      };
      step();
    });
  }
})();
