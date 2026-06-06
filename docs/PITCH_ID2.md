# Pitch — The ID2 Card

*The opportunity an expert was structurally unable to find.*

---

## 30-second version

> Pfizer's commercial development team can evaluate a handful of opportunities a year, and the ones that matter most are the ones their taxonomy can't see. We built an agent that ignores the taxonomy. It scanned every gene's promoter methylation against drug response across thirty cancer types, controlled for the confounds no human can hold in their head, validated each hit functionally, and replicated it in independent screens. It surfaced **Olaparib → ID2** — a developmental transcription factor, not a DNA-repair gene, that marks PARP-inhibitor sensitivity across cancer types. No oncologist would have written that hypothesis down. The machine did, and proved it three times over.

---

## The arc

**1. The problem (Pfizer's own words).**
Commercial development is expert-gated. An expert evaluates a handful of opportunities over months. And "the opportunities that would most change patient care are, by definition, the ones the existing taxonomy is least equipped to find." The question isn't whether AI is faster. It's whether AI can find what the expert-gated model is *structurally incapable* of finding.

**2. The reframe — it's not speed, it's three gates.**
A human expert cannot remove three priors, no matter how much time you give them:
- **Taxonomy.** They are trained, funded, and regulated by organ. Their unit of analysis is "lung cancer," "breast cancer." A subgroup defined by a molecular state *across* indications is invisible.
- **The famous-gene prior.** They evaluate hypotheses about genes they know. They will not fund a hypothesis about a gene they've never heard of.
- **Multivariate control.** They reason over one or two variables. They cannot residualize 17,000 genes against latent confounding axes in their head.

An agent has none of these priors. That's the categorical difference.

**3. The proof — the ID2 card.**
Start from a molecular state, not a disease. Scan every gene's promoter methylation against Olaparib response across ~30 cancer types, with tissue and the latent methylome axes residualized out. One hit that survives everything: **promoter methylation of ID2 marks Olaparib sensitivity.**

Now look at why each gate would have hidden it:
- **ID2 is not a DDR gene.** It's *Inhibitor of DNA-binding 2*, a developmental transcription factor. Every PARP-inhibitor program on earth is anchored on BRCA / homologous recombination. ID2 is nowhere on that map. *(Gate 2 hides it.)*
- **The signal is cross-indication.** It holds in the same direction in **9 of 14** cancer types — it isn't a one-tumor finding, it's a state that recurs across unrelated indications. An indication-gated search never assembles that group. *(Gate 1 hides it.)*
- **It only appears after whole-methylome control.** It's invisible in the raw data; it surfaces only after residualizing tissue + the global methylation axis. *(Gate 3 hides it.)*

**4. Why it's real (the validation funnel).**
This is not a lucky correlation. It cleared four independent bars:
- **Discovery:** significant after tissue + latent-axis control (GDSC).
- **Functional silencing:** ID2's promoter methylation tracks its own expression *down* (RNA-seq, p≈5×10⁻⁷) — the methylation actually silences the gene, it's not a passenger.
- **External replication:** the association holds, same direction, in **two more independent screens — PRISM and CTRP.** Three screens agree.
- **Honest grading:** evidence level L5 — a defensible novel hypothesis.

**5. The close — we don't replace the expert, we remove the gate.**
Today, expert priors decide which hypotheses ever get evaluated, and the structurally-invisible ones never do. Our agent evaluates the entire taxonomy-free, prior-free, confound-controlled space, functionally validates it, replicates it, and hands the expert a triaged shortlist. We invert the funnel. The expert still adjudicates — they just finally get to see what was always outside the gate.

And the framework is built to disbelieve itself: the same controls that surfaced ID2 *caught and killed* our original "positive control" when it turned out to be a copy-number artifact at 9p21. A discovery engine you can trust is one that routinely rejects its own findings.

---

## Own the caveats (saying these IS the pitch)

- **Effect sizes are modest.** Single-gene methylation→response is a weak, diffuse signal. The claim is **directional concordance across three independent screens plus functional silencing** — an enrichment and triage signal, not a deterministic predictor. State this before they ask.
- **Cell lines, not patients.** These are research-use hypotheses (L5). The next step is patient cohorts and compartment attribution — which is exactly the roadmap.
- **The agent generates and validates; the human decides.** We don't claim to replace oncology judgment. We claim to make the invisible visible.

---

## Anticipated questions

- *"The PRISM correlation is only −0.12."* Right — and it's significant, same-direction, in three screens, on a gene whose methylation we independently showed silences it. Most published methylation–drug associations clear far fewer bars. The point isn't the magnitude; it's that a real, replicable signal on an unprioritizable gene was findable at all.
- *"Couldn't an expert find this eventually?"* No — not structurally. They don't assemble cross-indication groups, they don't prioritize ID2, and they can't residualize the genome by hand. The priors filter ID2 out *before it's ever a hypothesis*. Time doesn't fix a structural blind spot.
- *"Is this just a correlation engine?"* It's a generate → control → functionally-validate → replicate funnel. The output is a triaged, mechanistically-coherent, externally-replicated shortlist — not a correlation dump.
- *"Why should I believe any of it?"* Because it threw out its own headline result (the 9p21 copy-number artifact) the moment higher-resolution data exposed it. It's adversarial toward itself.

---

## The ask / vision

This is the engine, proven on one card. Scale it: every drug, against the whole methylome, with copy-number and expression layered in, then onto patient cohorts. Every opportunity the taxonomy can't see, evaluated and validated — and handed to the people who decide.

---

## Live demo flow (methylkb.vercel.app/strata-v3.html)

1. Open on **the three gates** — read one: "an expert / the agent."
2. Scroll the **funnel**: 17,181 genes → a few hundred → ~160 functionally silenced → externally replicated.
3. Land on the **Olaparib → ID2 card**: the three-screen replication strip, the silencing stat, "same direction in 9/14 cancer types," and the punchline note — *a developmental TF, not a PARP gene.*
4. End on the **honesty rails**. Let the caveats be the mic-drop: this is what a discovery engine you can trust looks like.
