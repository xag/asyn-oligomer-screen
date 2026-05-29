// Build a self-contained interactive results page that joins the Stage 3 sweep
// ranking (results/sweep/*.csv) with the candidate-list delivery / CNS-reach
// metadata (data/vicinity_molecules.js).
//
// The sweep CSV carries the ranking signal (delta_activity_gated, aspr_score);
// it does NOT carry how a molecule reaches the brain. That lives only in the
// candidate list. This script joins the two on `id` and emits a single
// self-contained HTML file (data embedded inline — opens over file://, no
// server, no fetch).
//
//   node docs/build_display.mjs            # writes docs/index.html
//
// Run from the repo root. Output lives in docs/ so GitHub Pages can serve it.

import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '..');

const CSV_PATH = join(
  root,
  'results',
  'sweep',
  'fusco_parallel_3mer_core70-88_relaxed_sweep.csv',
);
const MOL_PATH = join(root, 'data', 'vicinity_molecules.js');
const OUT_PATH = join(here, 'index.html');

// --- load candidate metadata (the JS module is the source of delivery info) ---
function loadMolecules() {
  const text = readFileSync(MOL_PATH, 'utf8');
  // The file is an ES module with several `export const` declarations and no
  // imports. Strip the export keyword and evaluate in a function scope, then
  // return the array we care about.
  const body = text.replace(/^export\s+/gm, '');
  // eslint-disable-next-line no-new-func
  const fn = new Function(`${body}\nreturn VICINITY_MOLECULES;`);
  return fn();
}

// --- minimal CSV parser (handles quoted SMILES fields with commas) ---
function parseCsv(text) {
  const rows = [];
  let field = '';
  let row = [];
  let inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (inQuotes) {
      if (c === '"') {
        if (text[i + 1] === '"') {
          field += '"';
          i++;
        } else inQuotes = false;
      } else field += c;
    } else if (c === '"') inQuotes = true;
    else if (c === ',') {
      row.push(field);
      field = '';
    } else if (c === '\n' || c === '\r') {
      if (c === '\r' && text[i + 1] === '\n') i++;
      if (field !== '' || row.length) {
        row.push(field);
        rows.push(row);
        row = [];
        field = '';
      }
    } else field += c;
  }
  if (field !== '' || row.length) {
    row.push(field);
    rows.push(row);
  }
  const header = rows.shift();
  return rows.map((r) => Object.fromEntries(header.map((h, i) => [h, r[i]])));
}

const mols = loadMolecules();
const metaById = new Map(mols.map((m) => [m.id, m]));

const sweep = parseCsv(readFileSync(CSV_PATH, 'utf8'));

// reachability tier from delivery.feasibility (+ microbiome route = conditional)
function reachTier(meta) {
  const f = meta?.delivery?.feasibility ?? 'unknown';
  const route = meta?.delivery?.route ?? 'none-known';
  if (route === 'microbiome') return 'conditional';
  if (f === 'native' || f === 'achievable') return 'reaches';
  if (f === 'low-bioavailability' || f === 'invasive-only') return 'barrier';
  return 'unknown';
}

// provenance flag: how much prior evidence backs this entry
function provenance(meta) {
  if (meta?.validation_holdout) return 'holdout';
  if (meta?.role === 'anchor' || meta?.role === 'both') return 'known';
  return 'novel';
}

const records = [];
for (const r of sweep) {
  if (r.status !== 'ok') continue;
  const meta = metaById.get(r.mol_id);
  if (!meta) continue;
  const dActGated = Number(r.delta_activity_gated);
  const aspr = Number(r.aspr_score);
  const aff = Number(r.vina_top_affinity_kcal_per_mol);
  records.push({
    id: r.mol_id,
    name: meta.name ?? r.mol_name,
    group: meta.group,
    dActGated: Number.isFinite(dActGated) ? dActGated : 0,
    aspr: Number.isFinite(aspr) ? aspr : 0,
    aff: Number.isFinite(aff) ? aff : null,
    route: meta?.delivery?.route ?? 'none-known',
    feasibility: meta?.delivery?.feasibility ?? 'unknown',
    deliveryNotes: meta?.delivery?.notes ?? '',
    cnsLow: meta?.cns_conc?.low ?? null,
    cnsHigh: meta?.cns_conc?.high ?? null,
    cnsNote: meta?.cns_conc?.note ?? '',
    evidence: meta?.evidence ?? '',
    reach: reachTier(meta),
    prov: provenance(meta),
  });
}

// harmful = positive adduct reactivity; protective = destabilising (Δact < 0)
const harmful = records
  .filter((r) => r.aspr > 0)
  .sort((a, b) => b.aspr - a.aspr);

const harmfulIds = new Set(harmful.map((r) => r.id));
const protective = records
  .filter((r) => r.dActGated < 0 && !harmfulIds.has(r.id))
  .sort((a, b) => a.dActGated - b.dActGated);

// rank within each panel
protective.forEach((r, i) => (r.rank = i + 1));
harmful.forEach((r, i) => (r.rank = i + 1));

const payload = JSON.stringify({ protective, harmful });

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>α-synuclein oligomer screen — results</title>
<style>
  :root {
    --bg:#0e1117; --card:#161b22; --line:#2b3240; --txt:#e3e8ef; --dim:#9aa5b4;
    --reaches:#3fb950; --conditional:#d29922; --barrier:#f0883e; --unknown:#6e7681;
    --novel:#58a6ff; --known:#8b95a5; --holdout:#bc8cff;
    --prot:#3fb950; --harm:#f85149;
  }
  * { box-sizing:border-box; -webkit-text-size-adjust:100%; }
  body { margin:0; background:var(--bg); color:var(--txt);
    font:16px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; }
  .wrap { max-width:680px; margin:0 auto; padding:18px 14px 70px; }
  h1 { font-size:18px; margin:0 0 6px; line-height:1.3; }
  .disc { color:var(--dim); font-size:13px; margin:0 0 4px; }
  .disc b { color:var(--conditional); }

  /* segmented toggle */
  .seg { display:flex; gap:6px; margin:18px 0 6px; position:sticky; top:0;
    background:var(--bg); padding:8px 0; z-index:5; }
  .seg button { flex:1; min-height:46px; border:1px solid var(--line); background:var(--card);
    color:var(--dim); border-radius:10px; font-size:14px; font-weight:600; cursor:pointer; }
  .seg button.on.prot { color:#0e1117; background:var(--prot); border-color:var(--prot); }
  .seg button.on.harm { color:#0e1117; background:var(--harm); border-color:var(--harm); }
  .intro { color:var(--dim); font-size:13px; margin:2px 0 10px; }

  .card { background:var(--card); border:1px solid var(--line); border-radius:12px;
    padding:13px 14px; margin:0 0 11px; }
  .top { display:flex; align-items:center; gap:8px; flex-wrap:wrap; }
  .rk { color:var(--dim); font-variant-numeric:tabular-nums; font-size:13px; }
  .nm { font-weight:700; font-size:16px; flex:1; min-width:130px; }
  .prov { font-size:11px; padding:2px 8px; border-radius:999px; white-space:nowrap; }
  .prov.novel { color:var(--novel); border:1px solid var(--novel); }
  .prov.known { color:var(--known); border:1px solid var(--known); }
  .prov.holdout { color:var(--holdout); border:1px solid var(--holdout); }

  .eff { display:flex; align-items:center; gap:9px; margin:9px 0 11px; }
  .bar { height:8px; border-radius:4px; flex:0 0 auto; }
  .bar.prot { background:var(--prot); } .bar.harm { background:var(--harm); }
  .efftxt { color:var(--dim); font-size:12.5px; }

  .how { font-size:13.5px; line-height:1.5; border-top:1px solid var(--line); padding-top:10px; }
  .how .lbl { display:flex; align-items:center; gap:7px; margin-bottom:3px; }
  .how .lbl b { font-size:11px; text-transform:uppercase; letter-spacing:.04em; color:var(--dim); }
  .dot { width:10px; height:10px; border-radius:50%; flex:0 0 auto; }
  .dot.reaches{background:var(--reaches);} .dot.conditional{background:var(--conditional);}
  .dot.barrier{background:var(--barrier);} .dot.unknown{background:var(--unknown);}
  .route { display:inline-block; padding:1px 9px; border-radius:999px; font-size:12px;
    background:var(--bg); border:1px solid var(--line); color:var(--dim); }
  .reachline { margin:2px 0; }
  .reachline.barrier { color:var(--barrier); }
  .reachline.conditional { color:var(--conditional); }
  .reachline { margin-top:4px; }
  .mini { font-size:10.5px; text-transform:uppercase; letter-spacing:.04em; color:var(--dim); }
  .notes { color:var(--dim); }
  .lever { color:var(--reaches); }

  .more { display:block; width:100%; min-height:46px; margin:6px 0 0; border:1px solid var(--line);
    background:var(--card); color:var(--txt); border-radius:10px; font-size:14px; cursor:pointer; }
  .legend { color:var(--dim); font-size:12px; margin:18px 0 0; line-height:1.9; }
  .legend .dot { display:inline-block; vertical-align:middle; margin:0 3px 0 6px; }
  a { color:var(--novel); }
</style>
</head>
<body>
<div class="wrap">
  <h1>α-synuclein oligomer screen — results</h1>
  <p class="disc">Molecules ranked by a computational model of the toxic α-synuclein
    oligomer in Parkinson's, shown with how each typically reaches the brain.
    <b>Unvalidated hypotheses — not medical or dietary advice.</b> Ranks order
    molecules; they do not measure effect size.</p>

  <div class="seg">
    <button id="bProt" class="prot on" data-panel="protective">Worth testing</button>
    <button id="bHarm" class="harm" data-panel="harmful">Reduce exposure</button>
  </div>
  <p class="intro" id="intro"></p>

  <div id="list"></div>
  <button class="more" id="more"></button>

  <div class="legend" id="legend"></div>
</div>

<script>
const DATA = ${payload};
const TOPN = 20;
let panel = 'protective', showAll = false;

const ROUTE_LABEL = {
  endogenous:'made by the body', diet:'diet', supplement:'supplement',
  precursor:'dietary precursor', microbiome:'gut microbiome', environmental:'environmental',
  injection:'injection', 'none-known':'route unknown',
};

// concrete achievable brain concentration, straight from the curated CNS data
function brainLevel(r) {
  const range = (r.cnsLow && r.cnsHigh) ? r.cnsLow + '–' + r.cnsHigh
    : r.cnsLow ? r.cnsLow : null;
  if (!range && !r.cnsNote) return '';
  return (range || '') + (r.cnsNote ? (range ? ' — ' : '') + r.cnsNote : '');
}

// actionable lever for an anti-target, keyword-derived from its source notes
// only (not the mechanism text, which can mention "oxidative" etc. spuriously)
function harmLever(r) {
  const t = r.deliveryNotes.toLowerCase();
  const out = [];
  if (/glyc|\\bage\\b|ages|browning|maillard|hyperglyc/.test(t)) out.push('lower glycaemic load; limit browned & ultra-processed foods');
  if (/lipid peroxidation|polyunsaturat|oxidative|peroxid/.test(t)) out.push('antioxidant status; avoid oxidised / rancid fats');
  if (/smok|tobacco/.test(t)) out.push('avoid tobacco smoke');
  if (/heated|overheat|cooking oil|frying/.test(t)) out.push('avoid overheated cooking oils');
  return out.length ? out.join('; ') : null;
}

function shortName(n){ return n.replace(/\\s*\\(.*\\)/, ''); }

function cardProtective(r, max) {
  const w = Math.max(4, Math.round((-r.dActGated / max) * 140));
  const reachCls = r.reach === 'barrier' ? 'barrier' : r.reach === 'conditional' ? 'conditional' : '';
  return \`<div class="card">
    <div class="top">
      <span class="rk">#\${r.rank}</span>
      <span class="nm">\${shortName(r.name)}</span>
      <span class="prov \${r.prov}">\${r.prov === 'novel' ? 'new lead' : r.prov === 'holdout' ? 'sanity-check' : 'known'}</span>
    </div>
    <div class="eff"><span class="bar prot" style="width:\${w}px"></span>
      <span class="efftxt">predicted to destabilise — rank \${r.rank} of \${DATA.protective.length}</span></div>
    <div class="how">
      <div class="lbl"><span class="dot \${r.reach}"></span><b>Getting it to the brain</b>
        <span class="route">\${ROUTE_LABEL[r.route]}</span></div>
      \${r.deliveryNotes ? '<div class="notes">' + r.deliveryNotes + '.</div>' : ''}
      \${brainLevel(r) ? '<div class="reachline ' + reachCls + '"><span class="mini">Brain level reached</span> ' + brainLevel(r) + '</div>' : ''}
    </div>
  </div>\`;
}

function cardHarmful(r, max) {
  const w = Math.max(4, Math.round((r.aspr / max) * 140));
  const lever = harmLever(r);
  return \`<div class="card">
    <div class="top">
      <span class="rk">#\${r.rank}</span>
      <span class="nm">\${shortName(r.name)}</span>
      <span class="prov \${r.prov}">\${r.prov === 'known' ? 'documented' : 'suspected'}</span>
    </div>
    <div class="eff"><span class="bar harm" style="width:\${w}px"></span>
      <span class="efftxt">reactivity toward α-synuclein — +\${r.aspr.toFixed(2)}</span></div>
    <div class="how">
      <div class="lbl"><b>Where it comes from</b></div>
      \${r.deliveryNotes ? '<div class="notes">' + r.deliveryNotes + '.</div>' : ''}
      \${lever ? '<div class="reachline"><b style="color:var(--dim);font-size:11px;text-transform:uppercase;letter-spacing:.04em">Lower it by</b><br><span class="lever">' + lever + '</span></div>' : ''}
    </div>
  </div>\`;
}

function render() {
  const prot = panel === 'protective';
  const all = DATA[panel];
  const max = prot ? Math.max(...all.map(r => -r.dActGated)) : Math.max(...all.map(r => r.aspr));

  document.getElementById('intro').textContent = prot
    ? 'Candidates the model predicts could break up the toxic oligomer. Top of the list = strongest predicted effect.'
    : 'Reactive molecules the model flags as likely to damage α-synuclein and accelerate aggregation. These are exposures to reduce, not things to take.';

  const shown = showAll ? all : all.slice(0, TOPN);
  document.getElementById('list').innerHTML =
    shown.map(r => prot ? cardProtective(r, max) : cardHarmful(r, max)).join('');

  const more = document.getElementById('more');
  if (all.length > TOPN) {
    more.style.display = 'block';
    more.textContent = showAll ? 'Show fewer' : ('Show all ' + all.length);
  } else more.style.display = 'none';
}

document.querySelectorAll('.seg button').forEach(b => b.onclick = () => {
  panel = b.dataset.panel; showAll = false;
  document.getElementById('bProt').classList.toggle('on', panel === 'protective');
  document.getElementById('bHarm').classList.toggle('on', panel === 'harmful');
  window.scrollTo(0, 0);
  render();
});
document.getElementById('more').onclick = () => { showAll = !showAll; render(); };

document.getElementById('legend').innerHTML =
  'Reaches brain readily <span class="dot reaches"></span> · conditional <span class="dot conditional"></span> · '
  + 'delivery barrier <span class="dot barrier"></span>. '
  + 'Tags: <span class="prov novel" style="font-size:11px">new lead</span> no prior α-syn evidence · '
  + '<span class="prov known" style="font-size:11px">known</span> published modulator · '
  + '<span class="prov holdout" style="font-size:11px">sanity-check</span> blind validation entry. '
  + 'Method &amp; caveats: <a href="https://github.com/xag/asyn-oligomer-screen">repository README</a>.';

render();
</script>
</body>
</html>
`;

writeFileSync(OUT_PATH, html);
console.log(
  `wrote ${OUT_PATH}\n  protective: ${protective.length} molecules\n  harmful:    ${harmful.length} molecules`,
);
