// Build the wet-lab handoff package (docs/HANDOFF.md) for issue #11.
//
//   node docs/build_handoff.mjs           # writes docs/HANDOFF.md
//
// Joins the Stage 3 sweep ranking (results/sweep/*.csv) with candidate
// delivery metadata (data/vicinity_molecules.js) and the plain-English
// brain-access layer (data/brain_access.js), then emits a self-contained
// markdown handoff for an external wet-lab partner. Every score in the
// output is pulled live from the sweep CSV so it cannot drift from the
// pipeline; the per-molecule mechanism / assay / tier notes are editorial
// and live in HANDOFF_NOTES below.
//
// Run from the repo root. Output lives in docs/ alongside the results page.

import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '..');

const CSV_PATH = join(root, 'results', 'sweep', 'fusco_parallel_3mer_core70-88_relaxed_sweep.csv');
const MOL_PATH = join(root, 'data', 'vicinity_molecules.js');
const BRAIN_PATH = join(root, 'data', 'brain_access.js');
const OUT_PATH = join(here, 'HANDOFF.md');
const REPO = 'https://github.com/xag/asyn-oligomer-screen';

function loadConst(path, name) {
  const body = readFileSync(path, 'utf8').replace(/^export\s+/gm, '');
  // eslint-disable-next-line no-new-func
  return new Function(`${body}\nreturn ${name};`)();
}

function parseCsv(text) {
  const rows = [];
  let field = '', row = [], inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (inQuotes) {
      if (c === '"') {
        if (text[i + 1] === '"') { field += '"'; i++; } else inQuotes = false;
      } else field += c;
    } else if (c === '"') inQuotes = true;
    else if (c === ',') { row.push(field); field = ''; }
    else if (c === '\n' || c === '\r') {
      if (c === '\r' && text[i + 1] === '\n') i++;
      if (field !== '' || row.length) { row.push(field); rows.push(row); row = []; field = ''; }
    } else field += c;
  }
  if (field !== '' || row.length) { row.push(field); rows.push(row); }
  const header = rows.shift();
  return rows.map((r) => Object.fromEntries(header.map((h, i) => [h, r[i]])));
}

const metaById = new Map(loadConst(MOL_PATH, 'VICINITY_MOLECULES').map((m) => [m.id, m]));
const BRAIN = loadConst(BRAIN_PATH, 'BRAIN_ACCESS');
const byId = new Map(parseCsv(readFileSync(CSV_PATH, 'utf8')).map((r) => [r.mol_id, r]));

// --- editorial layer: which molecules go in the handoff, and why ---------
// `assay`   the orthogonal-mechanism readout this molecule needs (caveat #29)
// `mech`    one-line predicted mechanism behind the score
// `supply`  one-line supply / regulatory note feeding the priority tier
// `starter` true → part of the recommended first-pass trio
const PROTECTIVE = ['dhea', 'retinoic-acid', 'allopregnanolone', 'thc', 'piperine', 'trehalose', 'urolithin-a'];
const HARMFUL = ['malondialdehyde', 'acrolein', '4-hne', 'methylglyoxal'];

const NOTES = {
  dhea: {
    starter: true,
    mech: 'Steroid scaffold docks into the exposed hydrophobic β-core groove of the oligomer; predicted to raise the toxic shape’s free energy.',
    assay: 'ThT aggregation kinetics + DLS oligomer-size distribution on recombinant α-syn; dose–response from reported brain concentrations.',
    supply: 'OTC supplement, cheap, well-characterised pharmacology — lowest-friction candidate.',
  },
  'retinoic-acid': {
    mech: 'Nuclear-receptor ligand; here scored as a direct binder of the toxic shape independent of its transcriptional role.',
    assay: 'ThT + DLS under amber light (photolabile); pair with all-trans/13-cis comparison.',
    supply: 'Active form is light/oxygen-sensitive; therapeutic isotretinoin is regulated. Use dietary-retinol framing for the actionable readout.',
  },
  allopregnanolone: {
    mech: 'Neurosteroid; predicted direct binder of the hydrophobic groove, same scaffold class as DHEA.',
    assay: 'ThT + DLS; precursor (pregnenolone) arm for the actionable lever.',
    supply: 'Injectable in trials; the actionable handle is oral pregnenolone, not allopregnanolone itself.',
  },
  thc: {
    mech: 'Lipophilic cannabinoid; scored as a direct binder. CB1/CB2 signalling is a separate, untested axis.',
    assay: 'ThT + DLS in vitro (receptor-independent), so the binding hypothesis is testable without cannabinoid signalling.',
    supply: 'Controlled substance in most jurisdictions — highest regulatory friction; deprioritise unless a licensed partner is available.',
  },
  piperine: {
    starter: true,
    mech: 'Black-pepper alkaloid; scored as a direct binder. Note the separate, well-known bioavailability-enhancer role.',
    assay: 'ThT + DLS; clean food-compound chemistry, no special handling.',
    supply: 'Dietary, cheap, OTC — alongside DHEA the easiest to source.',
  },
  trehalose: {
    mech: 'Chemical chaperone / autophagy inducer. The reported in-vivo benefit is peripheral/autophagic — it does not reach the brain (see brain-access note).',
    assay: 'In-vitro ThT + DLS still tests the direct-binding hypothesis; do NOT read a positive in-vitro result as a brain-delivery claim.',
    supply: 'Food/supplement, cheap, but oral CNS entry is the limiting step — mechanism is the open question, not supply.',
  },
  'urolithin-a': {
    starter: true,
    mech: 'Gut metabolite of ellagitannins; mitophagy inducer with reported α-syn reduction in worm/cell models. Scored as a direct binder here.',
    assay: 'ThT + DLS for the binding hypothesis; cell-based mitophagy readout for the orthogonal mechanism.',
    supply: 'Direct supplement (Urolithin A / Mitopure) sidesteps the ~40% microbiome-converter limitation — well-supplied.',
  },
  malondialdehyde: {
    mech: 'Lipid-peroxidation dialdehyde; predicted Lys-Schiff cross-linking on the disordered N-tail.',
    assay: 'ThT (does the adduct accelerate aggregation?) + LC-MS adduct mapping to localise the modified residues.',
    supply: 'Lever: lower oxidative load; avoid oxidised / high-PUFA and ultra-processed fats.',
  },
  acrolein: {
    mech: 'Most reactive common α,β-unsaturated aldehyde; Michael adducts cross-link α-syn and accelerate oligomerisation (published anchor).',
    assay: 'ThT + LC-MS adduct mapping; positive control for the covalent channel given existing literature.',
    supply: 'Lever: avoid tobacco smoke and overheated cooking oils.',
  },
  '4-hne': {
    mech: 'ω-6 lipid-peroxidation aldehyde; Michael adducts on α-syn nucleophiles, elevated in PD substantia nigra (published anchor).',
    assay: 'ThT + LC-MS adduct mapping; benchmark against published 4-HNE–α-syn adduct sites.',
    supply: 'Lever: antioxidant status and dietary fat quality.',
  },
  methylglyoxal: {
    mech: 'Reactive dicarbonyl; glycates and cross-links α-syn directly (published anchor; hyperglycaemia→aggregation bridge).',
    assay: 'ThT + LC-MS adduct mapping (CEL/Schiff sites); glucose-load arm for the actionable lever.',
    supply: 'Lever: glycaemic control; cut high-sugar / high-AGE intake.',
  },
};

const ROUTE_CATEGORY = {
  supplement: 'Supplement',
  precursor: 'Diet / precursor',
  diet: 'Diet',
  endogenous: 'Endogenous / exposure',
};

function provenance(meta) {
  if (meta?.validation_holdout) return 'known (blind hold-out)';
  if (meta?.role === 'anchor' || meta?.role === 'both') return 'documented';
  return 'novel — no prior α-syn evidence';
}

function brainLevel(meta) {
  const c = meta?.cns_conc;
  if (!c) return '—';
  if (c.low && c.high) return `${c.low}–${c.high}`;
  return c.low || c.high || '—';
}

function row(cells) { return `| ${cells.join(' | ')} |`; }

function protectiveSection() {
  const lines = [];
  lines.push('## "What to do" — 7 novel protective candidates');
  lines.push('');
  lines.push('Top non-polyphenol hits from the 127-molecule sweep with no published α-syn direct-binding data. Ranked by `delta_activity_gated` (more negative = stronger predicted destabilisation of the toxic shape). All are dietary / endogenous / supplemental, not drug candidates.');
  lines.push('');
  lines.push(row(['#', 'Candidate', 'Category', 'Δact (gated)', 'Affinity (kcal/mol)', 'Reaches brain?', 'Supply / reg.']));
  lines.push(row(['---', '---', '---', '---:', '---:', '---', '---']));
  PROTECTIVE.forEach((id, i) => {
    const m = metaById.get(id), r = byId.get(id), n = NOTES[id], ba = BRAIN[id];
    const cat = ROUTE_CATEGORY[m?.delivery?.route] || m?.delivery?.route || '—';
    const reach = ba?.verdict ?? 'unknown';
    const star = n.starter ? ' ⭐' : '';
    lines.push(row([
      String(i + 1) + star,
      `**${m.name}**`,
      cat,
      Number(r.delta_activity_gated).toFixed(3),
      r.vina_top_affinity_kcal_per_mol,
      reach,
      n.starter ? 'first-pass trio' : '—',
    ]));
  });
  lines.push('');
  lines.push('⭐ = **recommended first-pass trio** (best supply + lowest regulatory friction): DHEA, piperine, urolithin A. All OTC and cheap; per issue #11 step "pick 3 most actionable based on supply / regulatory ease".');
  lines.push('');
  PROTECTIVE.forEach((id, i) => {
    const m = metaById.get(id), n = NOTES[id], ba = BRAIN[id];
    lines.push(`### ${i + 1}. ${m.name}`);
    lines.push('');
    lines.push(`- **Provenance:** ${provenance(m)}`);
    lines.push(`- **Predicted mechanism:** ${n.mech}`);
    lines.push(`- **Suggested assay:** ${n.assay}`);
    lines.push(`- **Delivery / actionable lever:** ${m.delivery?.notes ?? '—'}${ba?.lever ? ` — ${ba.lever}` : ''}`);
    lines.push(`- **Brain access:** ${ba?.route ?? '—'} (typical brain level ${brainLevel(m)})`);
    lines.push('');
  });
  return lines.join('\n');
}

function harmfulSection() {
  const lines = [];
  lines.push('## "What to avoid" — 4 reactive metabolites');
  lines.push('');
  lines.push('Top hits on the covalent / adduct channel (`aspr_score`). The −6 kcal/mol affinity gate correctly collapsed their non-covalent `delta_activity_gated` to ~zero — Vina cannot see covalent damage, so the covalent channel is the right axis. Ranked by `aspr_score` (higher = more reactive).');
  lines.push('');
  lines.push(row(['#', 'Metabolite', 'Category', 'aspr_score', 'Exposure lever']));
  lines.push(row(['---', '---', '---', '---:', '---']));
  HARMFUL.forEach((id, i) => {
    const m = metaById.get(id), r = byId.get(id), n = NOTES[id];
    const cat = ROUTE_CATEGORY[m?.delivery?.route] || m?.delivery?.route || '—';
    lines.push(row([
      String(i + 1),
      `**${m.name}**`,
      cat,
      `+${Number(r.aspr_score).toFixed(2)}`,
      n.supply.replace(/^Lever:\s*/, ''),
    ]));
  });
  lines.push('');
  HARMFUL.forEach((id, i) => {
    const m = metaById.get(id), n = NOTES[id], ba = BRAIN[id];
    lines.push(`### ${i + 1}. ${m.name}`);
    lines.push('');
    lines.push(`- **Provenance:** ${provenance(m)}`);
    lines.push(`- **Predicted mechanism:** ${n.mech}`);
    lines.push(`- **Suggested assay:** ${n.assay}`);
    lines.push(`- **Source / exposure:** ${ba?.route ?? m.delivery?.notes ?? '—'}`);
    lines.push('');
  });
  return lines.join('\n');
}

const today = new Date().toISOString().slice(0, 10);

const md = `# Wet-lab handoff package

> Generated by \`node docs/build_handoff.mjs\` from the live sweep CSV — do not edit by hand.
> Last generated: ${today}. Tracking issue: [#11](${REPO}/issues/11).

**This is computational research output. Nothing here is medical, dietary, or clinical advice. Every ranking is an unvalidated hypothesis pending experimental testing.**

## Purpose

This is the first external-validation moment for the screen. The pre-registered hold-out ([#6](${REPO}/issues/6)) shows the framework recovers *known* α-syn modulators; it does **not** show that novel predictions are real destabilisers in cells. This package hands 11 molecules — 7 novel protective candidates and 4 reactive metabolites — to a wet-lab partner for orthogonal-mechanism assays, turning sweep rankings into the project's two-sided deliverable: confirmed **what to do** (protective dietary / endogenous items) and confirmed **what to avoid** (reactive metabolites linked to diet, cooking, lifestyle).

**Receptor:** \`fusco_parallel_3mer_core70-88_relaxed.pdb\` (apo activity +17.78). **Channels:** \`delta_activity_gated\` (non-covalent rearrangement, −6 kcal/mol affinity-gated) for the protective set; \`aspr_score\` (covalent adduct propensity) for the reactive set.

${protectiveSection()}
${harmfulSection()}

## Assay guidance

- **Orthogonal mechanism is required (caveat [#29](${REPO}/issues/29)).** A direct-binding-only assay cannot confirm a destabiliser whose mechanism is conformational or autophagic. Use readouts that detect aggregation/oligomer-state change (ThT kinetics, DLS), not just binding.
- **Protective set:** ThT aggregation kinetics and/or DLS oligomer-size distribution on recombinant α-syn, dose–response across the reported brain-concentration band.
- **Reactive set:** ThT **plus** LC-MS adduct mapping to localise the covalent modification (predicted Lys-CEL / Lys-Schiff / Michael adducts on the disordered N-tail).
- **Start with the first-pass trio** (DHEA, piperine, urolithin A) for supply / regulatory ease, then expand.

## Caveats carried into the handoff

- **Polyphenol bias is partly framework-coupled** (caveat [#24](${REPO}/issues/24)) — the novel non-polyphenol hits carry more signal per rank than the polyphenol confirmatory band.
- **Scores are ordinal, not magnitude-interpretable** (caveat [#27](${REPO}/issues/27)) — read the rank order, not the absolute Δact.
- **Δact is sign-bound by static-pose features** (caveat [#23](${REPO}/issues/23)) — the screen currently sees destabilisers only.

## Pending: anti-target axis

The "what to avoid" half here is the **covalent** axis only. The symmetric **shape-stabiliser** anti-target axis (ligands that stabilise toxic shapes) is blocked on the multi-replica short-MD dwell-time channel ([#14](${REPO}/issues/14)) and its post-processing ([#30](${REPO}/issues/30)). When that lands, this package gains an anti-target flag column. Until then, the covalent reactive metabolites are the only confirmed-mechanism "avoid" entries.

## Reproduce

\`\`\`
python screen/sweep_oligomer.py --skip-existing   # refresh the sweep CSV
node docs/build_handoff.mjs                        # regenerate this file
\`\`\`

Full method, candidate list, and every caveat: ${REPO}
`;

writeFileSync(OUT_PATH, md);
console.log(`wrote ${OUT_PATH}`);
console.log(`  protective: ${PROTECTIVE.length}  reactive: ${HARMFUL.length}`);
const missing = [...PROTECTIVE, ...HARMFUL].filter((id) => !byId.get(id) || !metaById.get(id) || !NOTES[id]);
if (missing.length) console.log(`  ⚠ missing data for: ${missing.join(', ')}`);
else console.log('  ✓ every handoff molecule has sweep + metadata + notes');
