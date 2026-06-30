export const meta = {
  name: 'adversarial-paper-review',
  description: 'Deep multi-persona adversarial review of the markov-chain-sampler manuscript. Reader-base personas (ecologist, mathematician, student, interested public, computer-scientist) plus specialists (reproducibility, implementation, abstract-vs-open-questions, Reviewer-2 skeptic) read the RENDERED submission PDF first; every finding is adversarially refuted before an editor synthesis. Heavy; not for routine use.',
  phases: [
    { title: 'Review', detail: 'each persona reads the rendered PDF pages (some also the code) and files findings' },
    { title: 'Verify', detail: 'adversarially refute every finding against the rendered PDF / source / code' },
    { title: 'Synthesize', detail: 'editor dedups, ranks by severity x consensus, writes the review report' },
  ],
}

// ───────────────────────── inputs (args overrides; sane defaults) ─────────────────────────
const A          = args || {}
const PAPER      = A.paperDir   || '../markov-chain-sampler-paper'
const CODE       = A.codeDir    || '.'                     // GaussNiche repo root (cwd)
const USEMCMC    = A.useMcmcDir  || '../USE.MCMC'           // the sampler implementation
const PAGES      = A.pagesDir   || null                    // dir of page-NN.png renders of main.pdf (PRIMARY artifact)
const PDF        = A.pdf        || (PAPER + '/main.pdf')    // the source PDF (cannot be rendered on a host w/o poppler)
const FOCUS      = A.focus      || null                    // optional array of persona keys to run a subset
const MAX_FIND   = A.maxFindings || 12                     // findings verified per persona (top-severity slice)
const SECTION_PAGES = A.sectionPages || null               // optional "abstract pp1-2, methods pp4-12, ..." human note
const EXPLORE_WEB = A.exploreWeb !== false                 // expert personas may WebSearch/WebFetch further literature (default true)

if (!PAGES) {
  log('⚠ No pagesDir given. This skill is PDF-FIRST: render ' + PDF + ' to page PNGs and pass {pagesDir}.')
  log('  Portable:  pdftoppm -png -r 150 ' + PDF + ' <dir>/page   (or pdftocairo / magick -density 150)')
  log('  On Euler:  render with R pdftools inside the apptainer container (libpoppler, no ghostscript) — see SKILL.md.')
  log('  Proceeding in DEGRADED mode: personas review the .tex source only (final-submission/typesetting issues will be MISSED).')
}

// ───────────────────────── shared context every reviewer gets ─────────────────────────
const PDF_PRIMARY = PAGES
  ? `PRIMARY ARTIFACT — the rendered submission PDF, as page images in \`${PAGES}/\` (page-01.png, page-02.png, …).
**Read those PNGs — you are reviewing the typeset paper a reader/typesetter actually sees, NOT the LaTeX.** Use the Read tool on each
page image. ${SECTION_PAGES ? 'Section→page guide: ' + SECTION_PAGES + '. ' : ''}Read every page in your remit; skim the rest.
You may open the LaTeX source under \`${PAPER}/text/\` and \`${PAPER}/main.tex\` ONLY to pin a precise locus (file:line) for a fix —
never as a substitute for looking at the rendered page.`
  : `DEGRADED MODE — no rendered pages were supplied, so read the LaTeX source under \`${PAPER}/text/\` and \`${PAPER}/main.tex\`.
Explicitly note that you could NOT inspect the typeset PDF, so figure-legibility / overfull-box / float-placement / page-break /
print-size issues are OUT of your reach this run.`

const PAPER_MAP = `THE MANUSCRIPT (sibling repo \`${PAPER}\`, ~26 pages, Wiley NJD two-column template, target journal: Methods in Ecology and Evolution):
- Title + author block + ABSTRACT (one long paragraph) + keywords: \`${PAPER}/main.tex\` (the \\abstract{...} macro).
- Body, \\include'd in this order: \`text/introduction.tex\`, \`text/methods.tex\` (longest — the MCMC/GMM method, the equations,
  the virtual-species experimental design, the benchmark samplers), \`text/results.tex\` (sampler-comparison boxplots 2-D & 5-D +
  a "Computational performance" subsection of profiling numbers), \`text/discussion.tex\`, \`text/conclusion.tex\`.
- Appendix: \`text/appendix/appendix.tex\` (MCMC diagnostics — autocorrelation/Gelman-Rubin/trace; pseudocode; thinning justification;
  NN sampling) with \`text/appendix/ablation.tex\` \\input at its end (the uniform+ settings ablation).
- The method is "uniform+": a Metropolis-Hastings MCMC + Gaussian-Mixture-Model target that extends the grid-based \`uniform\` (USE)
  pseudo-absence sampler of dare_use_2023 to >2 environmental dimensions. Implemented elsewhere (this repo only writes it up).`

const RENDER_CHECKLIST = PAGES ? `FINAL-SUBMISSION / TYPESETTING CHECKS the rendered PDF makes possible (look for these as you read the pages):
- Figures: legible at PRINT size? resolution (any raster/pixelated figure — e.g. a .png among vector PDFs)? colour-blind-safe & greyscale-robust?
  caption matches what the panel shows? float placed near its first \\ref (not pages away / not orphaned at a page top with no text)?
- Overfull/underfull boxes: text or math running into the margin/gutter; equations overflowing the column; ugly inter-word "rivers".
- Tables: overflow past the column; mis-aligned columns; caption above/below per template.
- Page economy: orphan/widow lines; a heading at the very bottom of a page; a near-empty page; awkward column breaks; page count vs journal limit.
- Front matter: abstract block length/formatting in the template; author/affiliation/corresponding-author/keywords formatting; running heads.
- Stray artifacts actually RENDERED into the PDF: a visible "TODO", placeholder "XX", "Author Type", a Lato font fallback, a broken/blank figure box.
- List of figures / cross-references resolve to the right float.` : ''

const CITES_OK = !!A.citationsComplete
const CITE_CAVEAT = CITES_OK
  ? `- Citations ARE resolved in this build (biber ran) — review them NORMALLY: a "[?]" that still appears is a genuinely missing key (a real
  finding); also flag a duplicated / mis-attributed / unsupported citation. The bibliography is in scope.`
  : `- Unresolved citations rendering as "[?]" and a single "undefined references" warning are a no-biber local-build artifact (the citation-complete
  PDF comes from a biber build / Overleaf). **Ignore every [?] / "missing citation" / "empty bibliography" issue** — they are NOT manuscript defects here.`
const KNOWN_CAVEATS = `KNOWN BUILD CAVEATS — do NOT report these as findings unless noted otherwise:
${CITE_CAVEAT}
- "Article Type" in the header and the empty received/revised/accepted dates are journal-template placeholders, not content errors.
- The "% TODO" notes in Author contributions / Acknowledgments / Financial disclosure are LaTeX comments — flag them ONLY if they actually
  render as visible text in the PDF (they should not).`

const REVIEW_RULES = `HOW TO FILE A FINDING (quality bar — the verification phase will DELETE weak or wrong findings, so be precise, not voluminous):
- Ground every finding in something concrete: a quoted phrase + the PDF PAGE it is on (and, if helpful for the fix, the \`text/*.tex\` file:line).
- Assign honest severity: blocker (would sink the paper / is factually wrong / fatal flaw) · major (a reviewer would demand a change) ·
  minor (should fix) · nitpick (polish). Don't inflate.
- Prefer a few high-confidence, specific findings over many vague ones. "Section X could be clearer" with no locus is not a finding.
- Stay IN your persona's lens (below). The editor will merge overlaps across personas; you don't need to cover everything.
- A finding may be a STRENGTH worth preserving if it bears on your verdict, but your job is adversarial: find what is wrong, weak, unclear,
  unsupported, over-claimed, or unreproducible.`

const CODE_MAP = `THE IMPLEMENTATION (for code-reading specialists only):
- GaussNiche driver/engine: \`${CODE}\` — \`virtualSpecies_nd_fn.R\` (the k-D engine + the \`pa_random\`/\`pa_mcmc\`/\`pa_nn\` samplers),
  \`virtualSpecies_fn.R\` (2-D), the \`run_*_experiment.R\` / \`run_5d_hsm.R\` / \`tune_5d_hsm.R\` drivers, \`profile_samplers.R\`
  (the source of the "Computational performance" numbers → \`results/profile/profile_samplers.csv\`), \`hsm_eval.R\` (downstream HSMs).
  Read \`${CODE}/.claude/CLAUDE.md\` for the authoritative map + the settled findings.
- The sampler itself: \`${USEMCMC}/R/paSamplingMcmc.R\`, \`mclustDensityFunction.R\`, \`mcmcSampling.R\`, \`addHighDimGaussian.R\`,
  \`precomputeMcmcEnvironment.R\`, \`paSampling.R\` (the grid USE baseline).
KNOWN SUBTLETIES to check the prose against (each is a real, easy-to-misstate detail):
- env.rast MUST be the ORIGINAL rasters, never the PC scores (re-PCA-ing rotates E-space and makes uniform indistinguishable from random).
- \`species.cutoff.threshold\`: HIGHER quantile = WEAKER presence exclusion; =1 skips the presence GMM → exactly-uniform target (continuous limit).
- Methods says the presence cutoff is q_0.95 and the env cutoff q_0.001; check these against the code defaults and the ablation grid.
- Methods says the presence cap is N_max = 1000; the engine's documented default \`max_pres\` is 500 — check which the run scripts actually use.
- The downstream GLM is LINEAR (a quadratic GLM is near-oracle on a Gaussian niche → sampler-blind); Maxent = maxnet (no Java).
- SETTLED empirical finding (consistency check against the prose): uniform+ BEATS the naive samplers on discrimination (AUC/TSS) but LOSES
  on truth-recovery (correlation/RMSE to the known suitability surface, CBI) — random pseudo-absences reconstruct the true surface best.
  The ablation shows tuning the cutoffs narrows but never closes that gap. Watch for the main text over-claiming past this.`

// ───────────────────────── literature access (check claims against sources) ─────────────────────────
const LIT_PRIMARY = `LITERATURE — check the manuscript's claims against the ACTUAL cited sources (read the primary papers themselves; do NOT rely on any
second-hand digest or summary of them — those can inherit the manuscript's own framing):
- PRIMARY, on disk: the cited source PDFs at \`${PAPER}/sources/<citekey>.pdf\` — one per cited work, named by its \`references.bib\` key. Map a claim's
  \\cite{key} to \`sources/key.pdf\`, open it (Read renders/extracts the PDF), and read the relevant section of the SOURCE to judge whether it genuinely
  supports the claim, is over-attributed, or is mis-stated. Quote the source page in your finding. If a source PDF is absent, fall back to an open-access
  copy (below) or mark the claim "unverified — source not on hand"; never assert from memory of what a paper "probably" says.
- DOIs / URLs of every cited work: \`${PAPER}/references.bib\`.`
const LIT_WEB = EXPLORE_WEB
  ? `
- FURTHER literature (use it to check a claim against the broader field, find a MISSING comparison / counter-result, or confirm whether a value or method is
  standard): WebSearch + WebFetch. If they are not directly callable, load them first with ToolSearch (\`select:WebSearch,WebFetch\`). WebSearch returns
  results + OA/preprint links; WebFetch reads OA pages, abstracts and metadata. **Paywalled full text is NOT reachable** (WebFetch returns HTTP 402 on
  paywalled publishers — there is no institutional login here). For a cited paper, an open-access copy (ecoevorxiv / bioRxiv / ResearchGate / the author's
  GitHub) or the abstract + the on-disk source PDF is usually enough; if paywalled full text is genuinely needed to settle a point, SAY SO in the finding
  (the user can pull it). \`curl\` (Bash) from this host also reaches OA pages + Crossref/OpenAlex metadata. NEVER fabricate a source, quote, or result — "unverified" beats a guess.`
  : `
- Web exploration is DISABLED for this run (exploreWeb=false): use the cited source PDFs (\`sources/\`) + references.bib only; do not WebSearch/WebFetch.`
const LITERATURE = LIT_PRIMARY + LIT_WEB

// per-lens literature encouragement (woven into each persona's prompt; scales by how much that reader should dig)
const LIT_GUIDANCE = {
  ecologist: 'Read the load-bearing ecological sources themselves (bedia_dangers_2013 = buffer-out, meynard_testing_2019 = virtual species, dare_use_2023 = USE) and check the paper states them faithfully. WebSearch whether a relevant newer method/comparison is missing (e.g. the 2024 "pseudo-absences in ecological space" follow-up).',
  mathematician: 'Read the MCMC / optimal-scaling sources (roberts_weak_1997, roberts_examples_2009, metropolis_equation_1953, gelfand_sampling-based_1990, norris_markov_1998) and confirm the paper states their results correctly. WebSearch whether the 0.234 result\'s regularity conditions actually apply to this multimodal GMM target.',
  student: 'Light touch: if a cited definition or result is unclear in the text, glance at the cited source PDF to see what it actually says.',
  public: 'Light touch: sanity-check that an abstract claim is not stronger than the cited source supports.',
  cs: 'WebSearch the obvious algorithmic alternatives and prior art — rejection sampling / quasi-Monte-Carlo (Sobol) / HMC-NUTS on a differentiable GMM target / direct sampling of a bounded density — is a standard, cheaper method left un-cited and un-compared?',
  reproducibility: 'Read the methods sources to confirm the paper describes them faithfully; confirm the manuscript\'s data/code links (the GitHub repo) actually resolve (WebFetch/curl).',
  implementation: 'Use the cited source PDFs for any algorithm the code implements; otherwise keep the focus on code↔text consistency.',
  abstract: 'For each abstract sentence, confirm BOTH the body and the cited source support it; flag any claim the literature does not back.',
  skeptic: 'STRONGLY use the literature: read the load-bearing cited sources and check the paper does not over-attribute or mis-state them. WebSearch for missing baselines, a counter-result, or newer work that supersedes the contribution.',
}

// ───────────────────────── the persona roster ─────────────────────────
// readsCode = also inspects the implementation. Everyone reads the rendered PDF first.
const ALL_PERSONAS = [
  {
    key: 'ecologist', label: 'reader:ecologist', readsCode: false, effort: 'high',
    who: 'a working species-distribution / habitat-suitability modeller — the core Methods in Ecology & Evolution reader who might actually adopt this sampler',
    mandate: `Judge the paper as a practitioner deciding whether to use uniform+ on their own data.
- Is the problem real and honestly framed? Pseudo-absence choice genuinely matters for SDMs — is that motivated without hand-waving?
- Are the virtual species ecologically defensible caricatures? A bivariate-Gaussian niche on PCA axes; the 2x2 generalist/specialist x
  common/rare design — does it span the cases a practitioner cares about, or is it a toy that flatters the method?
- THE adoption question: the paper sells "highest environmental range coverage" and "drop-in replacement for USE", but would an ecologist
  actually switch? Does sampling ENVIRONMENTAL space (vs geographic) deliver what they want? Is the objective-dependent trade-off (uniform+
  helps discrimination but random recovers the true surface better) surfaced in the main text, or buried so an adopter would be misled?
- The real-data gap: everything is virtual species on a clean grid — no real occurrences, sampling bias, spatial autocorrelation, or detection
  error. Will a practitioner believe it transfers? Is interpretability lost by working on PC axes rather than named gradients (the paper admits this)?`,
  },
  {
    key: 'mathematician', label: 'reader:mathematician', readsCode: false, effort: 'xhigh',
    who: 'a statistician / applied probabilist who scrutinises the MCMC and mixture-model machinery for correctness',
    mandate: `Audit the mathematical claims for rigor and correctness.
- Is the target f(p) a proper (integrable) density on its support? Is the Metropolis-Hastings step valid (symmetric Gaussian proposal → the
  proposal ratio must cancel; acceptance = f(p')/f(p)) — and is detailed balance / the stationary-∝-f claim stated correctly?
- Does "uniform in environmental space" actually hold? The chain targets f≡1 on the support, but the kept points are nearest-neighbour-REMAPPED
  to discrete background cells, then thinned and de-duplicated. Nearest-neighbour remap of a continuous-uniform sample is Voronoi-area weighted,
  NOT uniform over the cells — is that approximation acknowledged and bounded, or quietly equated with "uniform"?
- The 0.234 optimal-acceptance result and Robbins-Monro adaptation: are the cited regularity conditions appropriate for this multimodal,
  GMM-shaped, bounded target? Does adapting only during burn-in preserve ergodicity (diminishing adaptation)?
- Convergence evidence: Gelman-Rubin + autocorrelation → thinning 30. Is effective sample size reported? Is a diagnostic on ONE setting
  generalised to all species/dimensions? GMM/EM local optima & component-count cap (≤9) — does fit instability perturb the target?
- Notation/precision: e.g. the "Monte-Carlo integration retrieves a density ∝ f" sentence; the env-threshold equation \\tau_env = q_0.001;
  the within-species paired comparisons over R=50 non-independent realisations (pseudo-replication?). Flag anything imprecise or unjustified.`,
  },
  {
    key: 'student', label: 'reader:student', readsCode: false, effort: 'medium',
    who: 'a first-year ecology/stats graduate student encountering pseudo-absence sampling and MCMC for the first time',
    mandate: `Judge clarity, pedagogy, and whether a newcomer can follow and conceptually reproduce the method.
- Is every term defined before first use (HSM, VS, GMM, PCA, pseudo-absence, suitability, prevalence, hypervolume, thinning, burn-in)?
- Is the notation consistent (f(p) vs f(x); m_env / \\tau_env; the jump from the f≡1 base target to the τ_env/1000-floored full equation)?
- Does the Methods narrative repeat itself or over-explain in places while under-explaining others? Where exactly would a student get lost?
- Do the figures stand alone (self-contained captions)? Are dense figures (the methods-overview multi-panel) legible and labelled?
- Does the appendix pseudocode match the prose method? Could the student redraw the pipeline from the paper alone?`,
  },
  {
    key: 'public', label: 'reader:public', readsCode: false, effort: 'medium',
    who: 'a scientifically-literate non-specialist / science communicator skimming the title, abstract, and figures',
    mandate: `Judge accessibility and whether the framing oversells.
- Is the ABSTRACT followable by a non-specialist? It is one long paragraph mixing motivation, method, and a three-part limitation list —
  is there a clear "we did X, found Y, it matters because Z", or just dense jargon?
- Does the TITLE deliver? It promises "high-dimensional environmental spaces", but the experiments stop at five dimensions — would a critical
  lay reader feel the claim is bigger than the evidence?
- Is the significance communicated — why anyone outside SDM should care (biodiversity / conservation impact)? Or is impact assumed?
- Jargon density in the opening: are MCMC / GMM / PCA given any intuition before they are leaned on?`,
  },
  {
    key: 'cs', label: 'reader:cs', readsCode: false, effort: 'high',
    who: 'a computer scientist / ML practitioner evaluating the algorithm, its complexity, and its design choices',
    mandate: `Judge the algorithmic and software claims.
- Scalability is the paper's headline, but is it QUANTIFIED? Grids cost ~res^d (stated); what is uniform+'s cost as a function of d, and is
  there any scaling curve, or only isolated 2-D & 5-D timings? "High-dimensional" with a 5-D ceiling — does the evidence match the claim?
- Algorithm choice: the target is f≡1 on a bounded support, lowered near presences — i.e. near-uniform on a region you can evaluate cheaply.
  Why random-walk Metropolis-Hastings at all, rather than direct REJECTION sampling on that bounded target, or low-discrepancy/quasi-Monte-Carlo
  (Sobol) coverage, or gradient-based HMC/NUTS (the GMM target is differentiable)? Is the MCMC choice justified against the obvious alternatives?
- "Drop-in replacement" is a software claim — is API/behavioural compatibility actually evidenced? Determinism/seeding/versioning for reproducibility?
- Minor tension: "kept in pure R so ecologists can read it" vs the compiled C++ inner loop. Complexity of the NN remap (KD-tree? O(n) scan?).`,
  },
  {
    key: 'reproducibility', label: 'spec:reproducibility', readsCode: true, effort: 'high',
    who: 'a reproducibility editor cross-checking the manuscript against the actual code and the numbers it reports',
    mandate: `Can a third party reproduce the results, and does the text match the implementation?
- Cross-check every reported NUMBER against its producer: the "Computational performance" figures (~7.7 s / 12.6 s env-GMM fit, 2.2 / 2.6 s
  presence fit, ~10 / 15 s full call, chain <1%, ~380 s USE grid search, 0.234, thinning 30) vs \`${CODE}/results/profile/profile_samplers.csv\`;
  the 82.5% variance and the 5-D species-parameter table vs the env5d build. Flag any number with no on-disk source or a mismatch.
- Methods↔implementation consistency: do the equations (env threshold q_0.001, presence threshold q_0.95, the τ_env/1000 floor, the
  half-largest-NN-distance remap rejection) match \`${USEMCMC}/R/paSamplingMcmc.R\` & \`mclustDensityFunction.R\`? Does N_max=1000 match the run scripts?
- Data & code availability: the GitHub link (danddr/GaussNiche) — is the env data reproducible (geodata downloads), are seeds documented,
  is the "reproducible vignette covering the four VS" real and runnable? Are package/software versions (USE.MCMC, mclust, hypervolume) pinned anywhere?
- Figure regenerability: is each paper figure regenerable from a documented producer? (Cross-ref the GaussNiche publication-figures map.)`,
  },
  {
    key: 'implementation', label: 'spec:implementation', readsCode: true, effort: 'high',
    who: 'a code reviewer checking that the sampler does what the paper says it does',
    mandate: `Read \`${USEMCMC}/R/paSamplingMcmc.R\` + \`mclustDensityFunction.R\` (+ the GaussNiche \`pa_mcmc\` wrapper) and confront them with the Methods text.
- Does the implemented target match equation "full density"? Is the presence term a SMOOTH roll-off (1 − m_pres/τ_pres) with a hard floor,
  and does "5% of points lie in an area not sampled at all" (q_0.95) actually follow from max(floor, 1 − m_pres/τ_pres)? Any off-by-one in the quantile sense?
- The species-quantile DIRECTION: paper must be consistent with HIGHER = WEAKER exclusion and =1 = skip-presence-model (exactly uniform). Does the prose imply the opposite anywhere?
- The remap-rejection rule ("half the largest nearest-neighbour distance, excluding outliers") — does the code do exactly that? The thinning ("every n-th, keep 2x desired") and de-duplication — does the prose match?
- Default parameter values quoted in the paper (proposal scale c=0.075, ≤9 GMM components, burn-in adaptation γ_t=(t+1)^{-0.6}) — verify each against the code.
- Report any place the manuscript describes behaviour the code does not implement, or vice-versa.`,
  },
  {
    key: 'abstract', label: 'spec:abstract-openq', readsCode: false, effort: 'high',
    who: 'an analyst assessing how faithfully the abstract (and conclusion) cover the paper\'s real contributions AND its open questions/limitations',
    mandate: `This is the dedicated "abstract vs. the rest of the paper" lens. Read the abstract (main.tex), then the results/discussion/conclusion/ablation, and adjudicate coverage.
- COVERAGE: does the abstract name every genuine contribution (MCMC+GMM target, drop-in USE replacement, modular custom-chain components,
  PCA for tractability, adaptive tuning, removal of the three USE limitations)? Anything claimed in the abstract the body does not deliver?
- OVER-CLAIM: "flexible and scalable framework … in complex environmental settings" — scalability shown only to 5-D; "complex settings" are
  virtual species, not real data. Does the abstract promise more than the experiments support?
- HIDDEN CAVEAT: the central downstream result is objective-dependent — uniform+ does NOT improve truth-recovery (random is best); the ablation
  shows tuning never closes the gap. The abstract presents only upside. Is omitting the main limitation honest for a methods paper, or should the
  abstract state it? Does the abstract's "for species distribution modeling" imply a downstream-accuracy benefit the appendix contradicts?
- OPEN QUESTIONS the paper raises but leaves unresolved (sensitivity to d — "left for future work"; real-data validation; MCMC vs cheaper
  alternatives; the remap/uniformity approximation error; GMM component-count sensitivity). Are these acknowledged in the conclusion, or silently dropped?
- Produce an explicit "what a handling editor would demand the abstract add/soften before acceptance" list.`,
  },
  {
    key: 'skeptic', label: 'spec:reviewer-2', readsCode: true, effort: 'xhigh',
    who: 'the harshest defensible peer reviewer ("Reviewer 2") hunting for fatal or near-fatal flaws',
    mandate: `Try to SINK the paper with the strongest fair criticisms a top-journal reviewer could raise. Be adversarial but honest (a flaw you invent will be refuted).
- CIRCULARITY / rigged evaluation: the virtual species are bivariate Gaussian and the sampler's target is built from Gaussian MIXTURES fitted to
  those presences — the method is tested on data matching its own modelling family. Would it survive non-Gaussian niches? Is this a stacked deck?
- OVERSOLD SCOPE: "high-dimensional" demonstrated at d=5, where the grid-explosion argument barely bites. No d=10/20/50 demonstration in the
  regime that motivates the whole paper. Is the central premise actually tested?
- FAIRNESS of the comparison: uniform+ gets two tunable cutoffs (and a whole appendix ablation); do RND / buffer-out (fixed 50 km) get any
  comparable tuning? Is the headline metric ("range coverage") almost TAUTOLOGICALLY won by an environmental sampler — i.e. are the diagnostics chosen to favour the method?
- SELECTIVE EMPHASIS: the main text foregrounds the metric uniform+ wins (coverage) and backgrounds the one it loses (downstream truth-recovery,
  where random wins — appendix). For a method whose POINT is better SDMs, is omitting the downstream result from the main Results a problem?
- STATS: R=50 Bernoulli realisations on ONE fixed background = pseudo-replication; boxplots but where are effect sizes / CIs in the main text?
- Any reported number resting on a known-contaminated cell (the env=0.25 ablation cell is ~24% random-fallback-inflated per the code's own notes)?
- Name the 1–3 issues that, unaddressed, would justify rejection, and what minimal change would defuse each.`,
  },
]

const PERSONAS = FOCUS ? ALL_PERSONAS.filter(p => FOCUS.includes(p.key)) : ALL_PERSONAS

// ───────────────────────── schemas ─────────────────────────
const SEVERITY = { type: 'string', enum: ['blocker', 'major', 'minor', 'nitpick'] }
const REVIEW_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    persona: { type: 'string' },
    overall: { type: 'string', description: '2–4 sentence verdict from this lens: would this reader accept/reject, and why' },
    findings: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        properties: {
          title:      { type: 'string', description: 'one-line summary of the issue' },
          severity:   SEVERITY,
          page:       { type: 'string', description: 'PDF page number(s) the issue appears on (or "n/a — source only")' },
          tex_locus:  { type: 'string', description: 'optional text/*.tex file:line for the fix' },
          render_only:{ type: 'boolean', description: 'true if this is a typeset/visual issue only visible in the rendered PDF' },
          claim:      { type: 'string', description: 'what the paper says / shows' },
          problem:    { type: 'string', description: 'why it is wrong / weak / unclear / unsupported, in this persona\'s voice' },
          evidence:   { type: 'string', description: 'the quoted phrase or cross-reference that grounds the finding' },
          suggestion: { type: 'string', description: 'a concrete, specific fix' },
        },
        required: ['title', 'severity', 'page', 'problem', 'suggestion'],
      },
    },
  },
  required: ['persona', 'overall', 'findings'],
}
const VERDICT_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    verdict:    { type: 'string', enum: ['confirmed', 'partly', 'refuted'] },
    confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
    corrected_severity: { type: 'string', enum: ['blocker', 'major', 'minor', 'nitpick', 'none'] },
    reasoning:  { type: 'string', description: 'what the manuscript/code actually says, and why the finding stands, is overstated, or is wrong' },
  },
  required: ['verdict', 'reasoning'],
}

// ───────────────────────── prompt builders ─────────────────────────
function reviewPrompt(p) {
  return `You are reviewing a scientific manuscript as ${p.who}.

${PDF_PRIMARY}

${PAPER_MAP}

${p.readsCode ? CODE_MAP + '\n' : ''}${RENDER_CHECKLIST ? RENDER_CHECKLIST + '\n' : ''}
${KNOWN_CAVEATS}

${LITERATURE}

YOUR LENS — ${p.key}:
${p.mandate}${LIT_GUIDANCE[p.key] ? '\n\nLiterature for your lens: ' + LIT_GUIDANCE[p.key] : ''}

${REVIEW_RULES}

Return your structured review: an \`overall\` verdict in this persona's voice, and a list of \`findings\` (your most important first;
the synthesis keeps roughly the top ${MAX_FIND}). Read the rendered pages before you write anything.`
}

function refutePrompt(f, persona) {
  return `You are an ADVERSARIAL fact-checker auditing one criticism raised against a scientific manuscript. Your default is to REFUTE.

The manuscript:
${PDF_PRIMARY}
${PAPER_MAP}
${KNOWN_CAVEATS}

${LIT_PRIMARY}
If the criticism concerns a CITATION or whether a claim is supported, open the cited source itself at \`${PAPER}/sources/<citekey>.pdf\` (or an open-access
copy): REFUTE if the source supports the claim; CONFIRM if it does not. Do not invent source content.

A "${persona}" reviewer raised this criticism:
- title:    ${f.title}
- severity: ${f.severity}
- page:     ${f.page}${f.tex_locus ? '\n- tex:      ' + f.tex_locus : ''}
- claim:    ${f.claim || '(see problem)'}
- problem:  ${f.problem}
- evidence: ${f.evidence || '(none given)'}
- proposed fix: ${f.suggestion}

Go to the cited page/locus in the rendered PDF (and the .tex source or code if needed) and check it FOR YOURSELF. Decide:
- refuted  — the criticism misreads the paper, is factually wrong, contradicts what the page actually shows, is already addressed elsewhere,
             or is one of the KNOWN BUILD CAVEATS above. Refute unless the evidence on the page clearly supports the criticism.
- partly   — there is a real but smaller/narrower issue than stated (give the corrected severity).
- confirmed — the criticism holds as stated, verifiable on the page.
Be specific about what the manuscript actually says. A plausible-sounding but unverifiable criticism is 'refuted', not 'confirmed'.`
}

// ───────────────────────── run: review → adversarial verify (pipelined) → synthesize ─────────────────────────
phase('Review')
log(`Adversarial review: ${PERSONAS.length} personas over ${PAGES ? 'the rendered PDF (' + PAGES + ')' : 'the LaTeX source (DEGRADED — no PDF render)'}.`)

const perPersona = await pipeline(
  PERSONAS,
  // stage 1 — the persona reads the rendered paper and files findings
  (p) => agent(reviewPrompt(p), { label: p.label, phase: 'Review', schema: REVIEW_SCHEMA, effort: p.effort }),
  // stage 2 — every finding from this persona is independently, adversarially refuted (no barrier vs other personas)
  (review, p) => {
    if (!review || !Array.isArray(review.findings) || review.findings.length === 0) return []
    const rank = { blocker: 0, major: 1, minor: 2, nitpick: 3 }
    const top = review.findings
      .slice()
      .sort((a, b) => (rank[a.severity] ?? 9) - (rank[b.severity] ?? 9))
      .slice(0, MAX_FIND)
    return parallel(top.map((f) => () =>
      agent(refutePrompt(f, p.key), { label: `verify:${p.key}`, phase: 'Verify', schema: VERDICT_SCHEMA, effort: 'low' })
        .then((v) => ({ ...f, persona: p.key, verdict: v }))
        .catch(() => ({ ...f, persona: p.key, verdict: { verdict: 'partly', reasoning: 'verifier error — left for the editor' } }))
    ))
  }
)

// collect what survived (pipeline returns stage-2 output: one array of verified findings per persona)
const verified = (perPersona || []).filter(Boolean).flat().filter(Boolean)
const confirmed = verified.filter((f) => f.verdict && f.verdict.verdict !== 'refuted')
const refutedCount = verified.length - confirmed.length

const tally = (arr, sev) => arr.filter((f) => (f.verdict?.corrected_severity && f.verdict.corrected_severity !== 'none' ? f.verdict.corrected_severity : f.severity) === sev).length
log(`Findings: ${verified.length} raised · ${confirmed.length} survived verification · ${refutedCount} refuted.`)
log(`Surviving by severity — blocker ${tally(confirmed, 'blocker')} · major ${tally(confirmed, 'major')} · minor ${tally(confirmed, 'minor')} · nitpick ${tally(confirmed, 'nitpick')}.`)

phase('Synthesize')
// compact payload for the editor — keep it lean
const payload = confirmed.map((f) => ({
  persona: f.persona,
  title: f.title,
  severity: (f.verdict?.corrected_severity && f.verdict.corrected_severity !== 'none') ? f.verdict.corrected_severity : f.severity,
  raised_severity: f.severity,
  verdict: f.verdict?.verdict,
  confidence: f.verdict?.confidence,
  page: f.page,
  tex_locus: f.tex_locus || '',
  render_only: !!f.render_only,
  problem: f.problem,
  suggestion: f.suggestion,
  verifier_note: f.verdict?.reasoning || '',
}))

const editorPrompt = `You are the HANDLING EDITOR assembling a single adversarial review of the manuscript "${'Pseudo-absences generation through a Markov Chain sampler for uniform sampling in high-dimensional environmental spaces'}" for Methods in Ecology and Evolution.

You are given ${payload.length} findings that already SURVIVED an adversarial verification pass (weak/wrong ones were culled), each tagged with the
reader-persona that raised it, a verifier-corrected severity, and the PDF page. Reader personas: ${PERSONAS.map(p => p.key).join(', ')}.

FINDINGS (JSON):
${JSON.stringify(payload, null, 1)}

Write a Markdown review report with these sections:
1. **Editor's summary & recommendation** — 4–8 sentences: the paper's core claim, its strongest issues, and a recommendation
   (accept / minor revision / major revision / reject), justified. Note explicitly if the review was PDF-based or source-only.
2. **Blocking & major issues** — a ranked list. MERGE findings that multiple personas raised (note the consensus, e.g. "raised by skeptic +
   mathematician + abstract") — cross-persona agreement is the strongest signal. For each: the issue, why it matters, the page/locus, and the concrete fix.
3. **Final-submission / typesetting** — the render-only issues (figure legibility/resolution, overfull boxes, float placement, page breaks, front matter). Omit this section if there were none.
4. **Abstract & open-questions assessment** — fold in the abstract-coverage lens: what the abstract over-claims, what limitation it hides, and the open questions the paper leaves unresolved.
5. **By reader audience** — a short paragraph each for ecologist, mathematician, student, interested public, computer-scientist: would this reader be satisfied, and the one thing to fix for them.
6. **Reproducibility & implementation** — does the text match the code and the reported numbers; what a reproducer would stumble on.
7. **Citation integrity** — any claim unsupported by, or over-attributed to, its cited source (checked against the source PDF), a mis-rendered/duplicated
   citation, or a missing comparison the literature demands. Omit if there were none. (For a methods paper these are serious — surface them prominently.)
8. **Prioritised fix list** — a numbered, do-this-first checklist the authors can work through, each item one line.

Be specific and cite pages. Do not invent findings beyond those given (you may merge, rank, and contextualise them). Return ONLY the Markdown.`

const report = await agent(editorPrompt, { label: 'editor:synthesis', phase: 'Synthesize', effort: 'high' })

return {
  report: report || '(editor produced no report)',
  stats: {
    personas: PERSONAS.map((p) => p.key),
    pdf_based: !!PAGES,
    raised: verified.length,
    confirmed: confirmed.length,
    refuted: refutedCount,
    by_severity: {
      blocker: tally(confirmed, 'blocker'),
      major: tally(confirmed, 'major'),
      minor: tally(confirmed, 'minor'),
      nitpick: tally(confirmed, 'nitpick'),
    },
  },
  findings: payload,
}
