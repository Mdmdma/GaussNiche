---
name: adversarial-paper-review
description: Deep, heavy, multi-persona ADVERSARIAL review of the markov-chain-sampler manuscript (../markov-chain-sampler-paper). Spawns reader-base personas (ecologist, mathematician, student, interested public, computer-scientist) plus specialists (reproducibility, implementation, abstract-vs-open-questions, a Reviewer-2 skeptic); they review the RENDERED submission PDF first (typeset figures, overfull boxes, float placement, page breaks — the things LaTeX-reading misses), check claims against the primary literature (the cited source PDFs in the paper's sources/ dir, NOT second-hand digests; experts also WebSearch/WebFetch open literature — no paywalled access), then every finding is adversarially refuted before an editor synthesis. Use when the user wants a thorough pre-submission critique of the paper from its target audiences, a referee-style report, or to stress-test the abstract/claims/reproducibility/citations. NOT routine — it fans out many agents and reads page images, so it is token-heavy by design. Runs via the Workflow tool. Pairs with the publication-figures skill (figure↔producer map) and euler-* skills (the PDF render needs the apptainer container).
---

# Adversarial paper review (multi-persona, PDF-first)

A heavy, non-routine harness that reviews the manuscript the way a journal does: from the
**rendered submission PDF**, through the eyes of its **target readers**, with every criticism
**adversarially fact-checked** before it reaches the report. It is deliberately expensive — run
it on demand (pre-submission, after a major revision, when a reviewer's angle needs simulating),
not in a loop.

The manuscript is the sibling repo `../markov-chain-sampler-paper/` (a single-paper Wiley NJD
LaTeX project, ~26 pages, target journal *Methods in Ecology and Evolution*). The method under
review is **uniform+**: a Metropolis–Hastings MCMC + Gaussian-Mixture-Model target that extends
the grid-based `uniform` (USE) pseudo-absence sampler of `dare_use_2023` to higher-dimensional
environmental spaces. The implementation lives in `../USE.MCMC` and is driven from this repo.

## Why PDF-first (and what that costs)

**Reviewers read the typeset PDF, not the LaTeX.** Reading the rendered pages is the whole point:
it catches the *final-submission* problems source-reading is blind to — figures illegible or
pixelated at print size (e.g. a raster `.png` among vector PDFs), overfull/underfull boxes spilling
into the margin, equations or tables overflowing the column, floats stranded pages from their
`\ref`, orphan/widow lines, a heading at the foot of a page, abstract-block formatting, stray
rendered placeholders. So this skill's **primary artifact is `main.pdf`**, and the source `.tex`
and the code are consulted only to pin a precise locus for a fix.

The catch on this cluster: **the login node has no PDF rasterizer** (no `pdftoppm`/`gs`/`pdftocairo`;
the spack poppler/imagemagick modules are library-only; `pdftotext` extracts nothing from the Lato
XeLaTeX PDF). The Read tool's own PDF mode therefore fails here. **So you must pre-render `main.pdf`
to page PNGs** (PNGs need no poppler — the Read tool shows them as images), then point the workflow
at that directory. This is the one required preflight.

> Token cost: each persona *looks at* the rendered pages (vision tokens), ×9 personas, plus
> per-finding verification and an editor pass. That is the intended weight — the user opts into it.
> To trim, pass `focus` (a subset of personas) or a `sectionPages` guide so personas read only their pages.

## Preflight A — build a citation-complete PDF (so the citation lens works)

Build with **biber** so `\cite`s resolve (otherwise ~33 render as `[?]` and the citation review is wasted). On Euler the
texlive module ships **no biber binary**, so a matching standalone biber is installed at **`~/.local/bin/biber`** (version **2.19**,
to match this texlive's biblatex **3.19** — version match is mandatory; a mismatched biber refuses the `.bcf`). Build from the repo
root (cwd = root for the Lato fonts + `graphics/` paths):

```bash
cd ../markov-chain-sampler-paper
module load stack/2025-06 gcc/12.2.0 texlive/20240312
export PATH="$HOME/.local/bin:$PATH"          # put biber on PATH
xelatex -interaction=nonstopmode main.tex     # generates main.bcf
biber main                                    # resolves citations → main.bbl
xelatex -interaction=nonstopmode main.tex     # pulls in the bibliography
xelatex -interaction=nonstopmode main.tex     # settles \ref/\cite cross-refs
# verify: grep -c "undefined references" main.log  → 0 ;  grep -c '\\entry{' main.bbl  → ~30
```

Then pass **`citationsComplete: true`** to the workflow so the personas review citations normally (a remaining `[?]` is then a genuinely
missing key, a real finding). If you cannot build with biber, skip this and leave `citationsComplete` false — the workflow then tells
personas to ignore `[?]`. (To review citations exactly as the journal sees them, build the Overleaf PDF instead.) If biber is ever
missing again, re-fetch the version-matched binary: the SourceForge `biblatex-biber/<ver>/binaries/Linux/biber-linux_x86_64.tar.gz`,
where `<ver>` matches the texlive's biblatex (`grep bltxversion main.bcf`).

## Preflight B — render `main.pdf` → page PNGs (REQUIRED)

Produce `page-01.png … page-NN.png` in a directory, then pass it as `pagesDir`.

**Recommended everywhere — PyMuPDF** (a self-contained binary wheel that bundles mupdf: no poppler, no
ghostscript, no container, no system libs). This is the **verified Euler login-node path** (the texlive module
ships no rasterizer and the rocker SIF has no `pdftools`):

```bash
# one-time: bootstrap pip for the system python, then install pymupdf (both are light wheel installs)
python3 -m pip --version || { curl -sL https://bootstrap.pypa.io/get-pip.py | python3 - --user; }
python3 -m pip install --user pymupdf
# render the citation-complete main.pdf → page PNGs on scratch
python3 - <<'PY'
import fitz, os, glob
OUT = "/cluster/scratch/%s/GaussNiche/review_pages" % os.environ["USER"]
os.makedirs(OUT, exist_ok=True)
for f in glob.glob(os.path.join(OUT, "page-*.png")): os.remove(f)   # clear stale
d = fitz.open(os.path.expanduser("~/markov-chain-sampler-paper/main.pdf"))
for i, p in enumerate(d):
    p.get_pixmap(dpi=150).save(os.path.join(OUT, "page-%02d.png" % (i + 1)))
print("rendered", d.page_count, "pages")
PY
```

Rendering ~27 pages at 150 dpi is a few seconds, single-core — light, fine on the login node (not "real compute").

**Alternatives if a rasterizer is already on PATH** (workstations):

```bash
pdftoppm  -png -r 150 main.pdf OUT/page      # poppler  → page-01.png …
pdftocairo -png -r 150 main.pdf OUT/page      # poppler (cairo)
magick -density 150 main.pdf OUT/page-%02d.png   # ImageMagick (needs a ghostscript delegate)
```

`pdftools` inside the apptainer rocker SIF (R `pdftools::pdf_convert`, libpoppler) also works **once installed** —
but it is NOT installed in the SIF by default, so that route needs the **euler-r-spack-setup** skill first; PyMuPDF
above avoids all of that. **Whichever path: verify before launching** that `page-01.png` exists and the Read tool
can open it as an image (and that citations show as author-year, not `[?]`).

### Two build caveats the personas are told to ignore

1. **Unresolved citations `[?]`.** Only if you build **without** biber. With the Preflight-A biber build (and
   `citationsComplete: true`) the citations resolve and are in scope. If you skip biber, the workflow tells every
   persona and verifier to **ignore `[?]` / missing-bibliography issues** (a no-biber artifact, not a defect).
2. **Template placeholders.** "Article Type" in the header, empty received/accepted dates, and the `% TODO`
   comments in the back-matter are journal-template/comment artifacts — flagged only if they actually render.

## The persona roster

Nine reviewers. The five **reader-base** personas reflect the audience the paper must satisfy and read the
PDF only; the four **specialists** go deeper and three of them also read the code. (The roster lives in the
workflow script — edit it there to add/tune personas.)

| key | persona | reads | lens (what they hunt) |
| --- | --- | --- | --- |
| `ecologist` | SDM/HSM practitioner (core MEE reader) | PDF | Is the problem real & honestly framed? Are the virtual species defensible? Would they actually adopt uniform+ given it loses on truth-recovery? The real-data gap. |
| `mathematician` | statistician / applied probabilist | PDF | Is the target a proper density & the M–H step valid? Is "uniform" actually achieved after NN-remap/thinning (Voronoi weighting)? 0.234/Robbins–Monro conditions, convergence evidence, pseudo-replication. |
| `student` | first-year grad newcomer | PDF | Terms defined before use? Notation consistent? Where would a newcomer get lost? Do figures & pseudocode stand alone? |
| `public` | literate non-specialist / communicator | PDF | Abstract followable? Does the title ("high-dimensional") oversell a 5-D result? Is the significance communicated? Jargon. |
| `cs` | computer scientist / ML | PDF | Is "scalable" quantified (cost vs d)? Why M–H rather than rejection / quasi-MC / HMC on this near-uniform target? "Drop-in replacement" evidenced? |
| `reproducibility` | reproducibility editor | PDF + **code** | Do the reported numbers match `results/profile/*.csv` & the env build? Methods↔implementation consistency. Data/seed/version availability. Figure regenerability. |
| `implementation` | code reviewer | PDF + **code** | Does `paSamplingMcmc.R` actually implement the equations (thresholds, the τ_env/1000 floor, remap rejection, thinning)? Cutoff *direction*. Default values (c=0.075, ≤9 comps, N_max). |
| `abstract` | abstract-vs-open-questions analyst | PDF | The dedicated lens the user asked for: does the abstract cover the real contributions AND the open questions? What does it over-claim or hide (the objective-dependent downstream caveat)? Editor demands. |
| `skeptic` | "Reviewer 2" | PDF + **code** | The fatal-flaw hunt: circular Gaussian-on-Gaussian evaluation, oversold "high-D" at d=5, unfair tuning, selective emphasis, stats, contaminated cells. The 1–3 reject-worthy issues. |

### The pipeline

`Review` → `Verify` → `Synthesize`, pipelined (a persona's findings are verified while other personas
still read):

1. **Review** — each persona reads the rendered pages (code personas also read `../USE.MCMC` + this repo)
   and files structured findings (title, severity, **PDF page**, claim, problem, evidence, fix), staying in lens.
2. **Verify (adversarial)** — every finding (top ~12 per persona by severity) gets an independent skeptic
   that **defaults to refuted**: it goes to the cited page and refutes the criticism unless the page clearly
   supports it. This culls plausible-but-wrong findings and the known build caveats. The editor never sees a refuted finding.
3. **Synthesize** — a handling-editor agent dedups across personas (cross-persona consensus = strongest
   signal), ranks by severity × consensus, and writes the report: editor recommendation, blocking/major issues,
   a **final-submission/typesetting** section, the **abstract & open-questions** assessment, a by-audience
   readout, reproducibility, and a prioritised fix list.

## Literature access — reviewers check claims against sources

Reviewers do not take the manuscript's citations on faith; they verify load-bearing claims against the cited
literature, and the experts explore further. What is actually available here (all wired into the workflow):

- **Primary, on disk — the cited source PDFs.** `../markov-chain-sampler-paper/sources/<citekey>.pdf` (one per cited
  work, named by its `references.bib` key). Reviewers map a claim's `\cite{key}` → `sources/key.pdf`, **read the
  actual source paper**, and judge whether it genuinely supports the claim. Every reviewer and every verifier is told
  to read the source itself before challenging a citation — and **explicitly NOT to rely on any second-hand digest**
  (e.g. the repo's `sidecars/`), which can inherit the manuscript's own framing and contaminate an adversarial read.
  `references.bib` supplies the DOIs/URLs. (NOTE: the `sources/` PDFs live on the *author's* machine, not the Euler
  box — so run the review where they are present, e.g. the laptop; otherwise reviewers fall back to OA copies + mark
  unverifiable claims as such.)
- **Further/open literature — WebSearch + WebFetch** (encouraged for the expert lenses: mathematician, cs,
  ecologist, reproducibility, skeptic). WebSearch works and surfaces OA/preprint links; WebFetch reads OA pages,
  abstracts, and metadata. The expert mandates point them at specific look-ups (e.g. cs → rejection / quasi-MC /
  HMC prior art; skeptic → missing baselines & newer work; mathematician → the 0.234 regularity conditions).
- **Honest limit — no institutional/paywall access.** WebFetch runs Anthropic-side (not the ETH IP) and returns
  **HTTP 402 on paywalled publishers**; there is no headless browser with ETH Shibboleth/EZproxy login, and the
  skill does **not** scrape paywalled journals through the login node's ETH IP (that needs your personal
  credentials and is brittle). For cited papers the on-disk source PDF (or an OA copy) is almost always enough; when a
  reviewer truly needs a paywalled full text it doesn't have, the report flags it **for you to pull**.

Set **`exploreWeb: false`** to restrict reviewers to the source PDFs + `references.bib` (no web) — cheaper/faster, and
fully deterministic. Default is `true`.

## How to run

The skill ships a turnkey workflow script. After the preflight, run it via the **Workflow** tool with
`scriptPath` (do not paste the script inline — point at the file so edits/iteration use the same path):

```
Workflow({
  scriptPath: ".claude/skills/adversarial-paper-review/adversarial_review.workflow.js",
  args: {
    pagesDir: "/cluster/scratch/$USER/GaussNiche/review_pages",  // REQUIRED for PDF-first review
    citationsComplete: true,                     // set true after the Preflight-A biber build (else [?] is ignored)
    paperDir: "../markov-chain-sampler-paper",   // default
    codeDir:  ".",                                // default (GaussNiche root, the cwd)
    useMcmcDir: "../USE.MCMC",                    // default
    // optional:
    // exploreWeb: true,                            // experts may WebSearch/WebFetch further literature (default true; source PDFs always used)
    // focus: ["skeptic","mathematician","abstract"],  // run a subset
    // sectionPages: "abstract pp1, intro pp2-3, methods pp3-10, results pp10-13, discussion pp13-15, appendix pp16-26",
    // maxFindings: 12
  }
})
```

The workflow runs in the background and returns `{ report, stats, findings }`. **Write the `report`
Markdown to a file and surface the headline** to the user — suggested destination (heavy artifacts on
scratch, and the review is regenerable):

```
/cluster/scratch/$USER/GaussNiche/reviews/adversarial-review-<UTC-timestamp>.md
```

(stamp the timestamp yourself — the workflow can't). Also dump `findings` as JSON beside it. If the user
wants the review inside the paper repo, copy it by hand — but **do not commit it to the paper repo without
asking** (that repo syncs to Overleaf and has strict git rules; see its `CLAUDE.md`).

### Degraded mode

If you launch without `pagesDir`, the workflow still runs but personas review the **LaTeX source only** and
explicitly flag that typeset/figure/page issues were out of reach. That defeats the PDF-first purpose — only
use it when rendering is genuinely impossible, and tell the user what was skipped.

## Extending / tuning

- **Add or sharpen a persona**: edit the `ALL_PERSONAS` array in the workflow script (key, `who`, `readsCode`,
  `effort`, `mandate`). Keep mandates *specific to this paper* — the existing ones cite the actual equations,
  numbers, and the settled findings so reviewers are sharp, not generic. Then re-run with the same `scriptPath`.
- **Subset run**: `args.focus = ["skeptic","abstract"]` for a fast, targeted pass (e.g. just the harshest
  critic + the abstract analyst after a revision).
- **Citation review**: render the Overleaf (citation-complete) PDF, pass it, and tell the user citations are in
  scope (otherwise `[?]` is ignored by design).
- **Resume after an edit**: the Workflow tool result includes a `runId`; relaunch with `resumeFromRunId` to reuse
  cached persona reviews and only re-run what changed.

## Pitfalls

- **No render → no submission review.** The single most important step is the preflight. If the Read tool can't
  open `page-01.png`, fix that before launching, or the whole point is lost.
- **Don't over-trust a single verifier.** Verification defaults to *refuted* to kill hallucinated criticism, but a
  genuinely subtle finding can be wrongly refuted. The editor sees the verifier's reasoning; if a known-true issue
  is missing from the report, re-run that persona with `focus` and a sharper mandate.
- **The numbers are checkable.** The reproducibility/implementation personas cross-read `../USE.MCMC` and this
  repo's `results/` — keep `codeDir`/`useMcmcDir` pointed correctly or those lenses go blind.
- **Heavy by design.** This is not `/code-review`; it spawns ~9 reviewers + dozens of verifiers + an editor and
  reads page images. Run it deliberately.
