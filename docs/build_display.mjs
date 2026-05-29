// Build a self-contained interactive results page that joins the Stage 3 sweep
// ranking (results/sweep/*.csv) with the candidate-list delivery metadata
// (data/vicinity_molecules.js) and the plain-English brain-access justification
// layer (data/brain_access.js).
//
//   node docs/build_display.mjs            # writes docs/index.html
//
// Run from the repo root. Output lives in docs/ so GitHub Pages can serve it.

import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '..');

const CSV_PATH = join(root, 'results', 'sweep', 'fusco_parallel_3mer_core70-88_relaxed_sweep.csv');
const MOL_PATH = join(root, 'data', 'vicinity_molecules.js');
const BRAIN_PATH = join(root, 'data', 'brain_access.js');
const OUT_PATH = join(here, 'index.html');
const REPO = 'https://github.com/xag/asyn-oligomer-screen';

// --- load an ES-module data file by stripping `export` and returning a const ---
function loadConst(path, name) {
  const body = readFileSync(path, 'utf8').replace(/^export\s+/gm, '');
  // eslint-disable-next-line no-new-func
  return new Function(`${body}\nreturn ${name};`)();
}

// --- minimal CSV parser (handles quoted SMILES fields with commas) ---
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
const sweep = parseCsv(readFileSync(CSV_PATH, 'utf8'));

function provenance(meta) {
  if (meta?.validation_holdout) return 'holdout';
  if (meta?.role === 'anchor' || meta?.role === 'both') return 'known';
  return 'novel';
}

const records = [];
const missingBrain = [];
for (const r of sweep) {
  if (r.status !== 'ok') continue;
  const meta = metaById.get(r.mol_id);
  if (!meta) continue;
  const ba = BRAIN[r.mol_id];
  if (!ba) missingBrain.push(r.mol_id);
  const dActGated = Number(r.delta_activity_gated);
  const aspr = Number(r.aspr_score);
  records.push({
    id: r.mol_id,
    name: meta.name ?? r.mol_name,
    dActGated: Number.isFinite(dActGated) ? dActGated : 0,
    aspr: Number.isFinite(aspr) ? aspr : 0,
    deliveryNotes: meta?.delivery?.notes ?? '',
    cnsLow: meta?.cns_conc?.low ?? null,
    cnsHigh: meta?.cns_conc?.high ?? null,
    verdict: ba?.verdict ?? 'unknown',
    brainRoute: ba?.route ?? '',
    lever: ba?.lever ?? '',
    prov: provenance(meta),
  });
}

const harmful = records.filter((r) => r.aspr > 0).sort((a, b) => b.aspr - a.aspr);
const harmfulIds = new Set(harmful.map((r) => r.id));
const protective = records
  .filter((r) => r.dActGated < 0 && !harmfulIds.has(r.id))
  .sort((a, b) => a.dActGated - b.dActGated);

const payload = JSON.stringify({ protective, harmful });

const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>α-synuclein oligomer screen — results</title>
<style>
  :root {
    --bg:#0e1117; --card:#161b22; --line:#2b3240; --txt:#e3e8ef; --dim:#9aa5b4; --track:#222a36;
    --crosses:#3fb950; --boost:#2dd4bf; --limited:#d29922; --marker:#8b95a5;
    --sub:#f0883e; --none:#f85149; --prot:#3fb950; --harm:#f85149;
    --known:#8b95a5; --novel:#58a6ff; --holdout:#bc8cff;
  }
  * { box-sizing:border-box; -webkit-text-size-adjust:100%; }
  body { margin:0; background:var(--bg); color:var(--txt);
    font:16px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; }
  .wrap { max-width:680px; margin:0 auto; padding:18px 14px 70px; }
  h1 { font-size:18px; margin:0 0 6px; line-height:1.3; }
  .disc { color:var(--dim); font-size:13px; margin:0 0 10px; }
  .disc b { color:var(--limited); }
  .contribute { display:inline-block; font-size:13px; margin:0 0 4px; color:var(--novel); text-decoration:none; }
  .contribute:hover { text-decoration:underline; }

  .seg { display:flex; gap:6px; margin:14px 0 6px; position:sticky; top:0;
    background:var(--bg); padding:8px 0; z-index:5; }
  .seg button { flex:1; min-height:46px; border:1px solid var(--line); background:var(--card);
    color:var(--dim); border-radius:10px; font-size:14px; font-weight:600; cursor:pointer; }
  .seg button.on.prot { color:#0e1117; background:var(--prot); border-color:var(--prot); }
  .seg button.on.harm { color:#0e1117; background:var(--harm); border-color:var(--harm); }
  .intro { color:var(--dim); font-size:13px; margin:2px 0 12px; }

  details.mol { background:var(--card); border:1px solid var(--line); border-radius:10px;
    margin:0 0 7px; overflow:hidden; }
  details.mol > summary { list-style:none; cursor:pointer; display:flex; flex-direction:column;
    align-items:stretch; gap:8px; padding:11px 13px; }
  details.mol > summary::-webkit-details-marker { display:none; }
  .nm { font-weight:600; font-size:15.5px; line-height:1.3; }
  .meta { display:flex; align-items:center; gap:10px; }
  .track { flex:1 1 auto; height:9px; border-radius:5px; background:var(--track); overflow:hidden; }
  .track > i { display:block; height:100%; border-radius:5px; }
  .chipwrap { flex:0 0 132px; display:flex; justify-content:flex-end; }
  .track > i.prot { background:var(--prot); } .track > i.harm { background:var(--harm); }
  .vchip { font-size:11px; padding:2px 9px; border-radius:999px; white-space:nowrap; flex:0 0 auto; }
  .v-crosses { color:var(--crosses); border:1px solid var(--crosses); }
  .v-boost { color:var(--boost); border:1px solid var(--boost); }
  .v-limited { color:var(--limited); border:1px solid var(--limited); }
  .v-marker { color:var(--marker); border:1px solid var(--marker); }
  .v-sub { color:var(--sub); border:1px solid var(--sub); }
  .v-none { color:var(--none); border:1px solid var(--none); }
  .v-unknown { color:var(--dim); border:1px solid var(--dim); }

  .body { padding:2px 13px 13px; font-size:13.5px; }
  .pchip { display:inline-block; font-size:11px; padding:2px 9px; border-radius:999px; margin:4px 0 2px; }
  .p-known { color:var(--known); border:1px solid var(--known); }
  .p-novel { color:var(--novel); border:1px solid var(--novel); }
  .p-holdout { color:var(--holdout); border:1px solid var(--holdout); }
  .field { margin:9px 0; }
  .field .lbl { font-size:11px; text-transform:uppercase; letter-spacing:.04em; color:var(--dim);
    display:block; margin-bottom:1px; }
  .field .val { line-height:1.5; }
  .raise { color:var(--boost); }
  .reduce { color:var(--harm); }

  .grouphdr { color:var(--dim); font-size:12.5px; margin:20px 0 8px; display:flex; align-items:center; gap:8px; }
  .grouphdr::before, .grouphdr::after { content:""; flex:1; height:1px; background:var(--line); }
  details.group { margin:8px 0 0; }
  details.group > summary { list-style:none; cursor:pointer; color:var(--dim); font-size:13px;
    padding:12px 13px; min-height:46px; background:var(--card); border:1px dashed var(--line);
    border-radius:10px; }
  details.group > summary::-webkit-details-marker { display:none; }
  details.group[open] > summary { margin-bottom:8px; }
  .gcount { color:var(--txt); }

  .foot { color:var(--dim); font-size:12px; margin:22px 0 0; line-height:2; }
  .foot a { color:var(--novel); }
</style>
</head>
<body>
<div class="wrap">
  <h1>α-synuclein oligomer screen — results</h1>
  <a class="contribute" href="${REPO}">Contribute, critique, or test a candidate on GitHub →</a>
  <p class="disc">Molecules ranked by a computational model of the toxic α-synuclein
    oligomer in Parkinson's, each with a plain-English account of how it actually
    reaches the brain. <b>Unvalidated hypotheses — not medical or dietary advice.</b>
    Tap any molecule for detail.</p>

  <div class="seg">
    <button id="bProt" class="prot on" data-panel="protective">Worth testing</button>
    <button id="bHarm" class="harm" data-panel="harmful">Reduce exposure</button>
  </div>
  <p class="intro" id="intro"></p>

  <div id="list"></div>
  <div class="foot" id="foot"></div>
</div>

<script>
const DATA = ${payload};
let panel = 'protective';

const VERDICT = {
  crosses:        { label:'reaches the brain',      cls:'v-crosses' },
  boost:          { label:'raise it naturally',     cls:'v-boost' },
  limited:        { label:'a little gets in',       cls:'v-limited' },
  marker:         { label:'not a target',          cls:'v-marker' },
  subtherapeutic: { label:'too little gets in',     cls:'v-sub' },
  'does-not-reach':{ label:'can’t get in',          cls:'v-none' },
  unknown:        { label:'unclear',                cls:'v-unknown' },
};
const REACH = new Set(['crosses', 'boost', 'limited', 'unknown']);
const MARKER = new Set(['marker']);

function provLabel(prov, prot) {
  if (prot) {
    if (prov === 'known') return ['Known α-synuclein modulator', 'p-known'];
    if (prov === 'holdout') return ['Known modulator (used as a blind test)', 'p-holdout'];
    return ['New lead — no prior α-syn evidence', 'p-novel'];
  }
  if (prov === 'known') return ['Documented α-synuclein-damaging agent', 'p-known'];
  return ['Suspected — weaker evidence', 'p-novel'];
}
function brainLevel(r) {
  return (r.cnsLow && r.cnsHigh) ? r.cnsLow + '–' + r.cnsHigh : (r.cnsLow || '');
}
function harmLever(r) {
  const t = r.deliveryNotes.toLowerCase(); const o = [];
  if (/glyc|\\bage\\b|ages|browning|maillard|hyperglyc/.test(t)) o.push('keep blood sugar down; cut back on browned & ultra-processed foods');
  if (/lipid peroxidation|polyunsaturat|oxidative|peroxid/.test(t)) o.push('eat more antioxidants; avoid oxidised / rancid fats');
  if (/smok|tobacco/.test(t)) o.push('avoid tobacco smoke');
  if (/heated|overheat|cooking oil|frying/.test(t)) o.push('don’t overheat cooking oils');
  return o.join('; ');
}
function shortName(n){ return n.replace(/\\s*\\(.*\\)/, ''); }
function field(lbl, val, cls){ return val ? '<div class="field"><span class="lbl">' + lbl + '</span><span class="val ' + (cls||'') + '">' + val + '</span></div>' : ''; }
function pct(v, max){ return Math.max(3, Math.round((v / max) * 100)); }

function molProtective(r, max) {
  const v = VERDICT[r.verdict];
  const [plabel, pcls] = provLabel(r.prov, true);
  return \`<details class="mol">
    <summary>
      <span class="nm">\${shortName(r.name)}</span>
      <div class="meta">
        <span class="track"><i class="prot" style="width:\${pct(-r.dActGated, max)}%"></i></span>
        <span class="chipwrap"><span class="vchip \${v.cls}">\${v.label}</span></span>
      </div>
    </summary>
    <div class="body">
      <span class="pchip \${pcls}">\${plabel}</span>
      \${field('How it gets into the brain', r.brainRoute)}
      \${field('How to make more', r.lever, 'raise')}
      \${field('Typical brain level', brainLevel(r))}
    </div>
  </details>\`;
}

function molHarmful(r, max) {
  const [plabel, pcls] = provLabel(r.prov, false);
  const lever = harmLever(r);
  return \`<details class="mol">
    <summary>
      <span class="nm">\${shortName(r.name)}</span>
      <div class="meta">
        <span class="track"><i class="harm" style="width:\${pct(r.aspr, max)}%"></i></span>
      </div>
    </summary>
    <div class="body">
      <span class="pchip \${pcls}">\${plabel}</span>
      \${field('Where it comes from', r.brainRoute)}
      \${field('How to reduce your exposure', lever, 'reduce')}
    </div>
  </details>\`;
}

function render() {
  const prot = panel === 'protective';
  const all = DATA[panel];
  const max = prot ? Math.max(...all.map(r => -r.dActGated)) : Math.max(...all.map(r => r.aspr));
  document.getElementById('intro').textContent = prot
    ? 'Candidates the model predicts could break up the toxic clump, grouped by whether you can actually get more of them into the brain. Longer bar = stronger predicted effect.'
    : 'Reactive molecules the model flags as damaging to α-synuclein. These are exposures to reduce — longer bar = more reactive.';

  const list = document.getElementById('list');
  if (!prot) { list.innerHTML = all.map(r => molHarmful(r, max)).join(''); return; }

  const reach = all.filter(r => REACH.has(r.verdict));
  const marker = all.filter(r => MARKER.has(r.verdict));
  const barrier = all.filter(r => !REACH.has(r.verdict) && !MARKER.has(r.verdict));
  const group = (title, sub, arr) => arr.length
    ? '<details class="group"><summary><span class="gcount">' + arr.length + '</span> ' + title
      + ' — ' + sub + '</summary>' + arr.map(r => molProtective(r, max)).join('') + '</details>'
    : '';
  list.innerHTML =
    reach.map(r => molProtective(r, max)).join('')
    + group('made in the brain', 'by-products or markers — nothing useful to raise, or you’d want less', marker)
    + group('can’t get into the brain', 'bind in the model but don’t reach the brain from outside — a delivery problem, not benefits', barrier);
}

document.querySelectorAll('.seg button').forEach(b => b.onclick = () => {
  panel = b.dataset.panel;
  document.getElementById('bProt').classList.toggle('on', panel === 'protective');
  document.getElementById('bHarm').classList.toggle('on', panel === 'harmful');
  window.scrollTo(0, 0);
  render();
});

document.getElementById('foot').innerHTML =
  '<b>Tags:</b> '
  + '<span class="vchip v-crosses">reaches the brain</span> eat/supplement it · '
  + '<span class="vchip v-boost">raise it naturally</span> the brain makes it — diet/activity/sleep raises it · '
  + '<span class="vchip v-limited">a little gets in</span> modest/uncertain · '
  + '<span class="vchip v-marker">by-product</span> not a target · '
  + '<span class="vchip v-sub">too little gets in</span> / <span class="vchip v-none">can’t get in</span>.<br>'
  + 'Every result is an unvalidated computational hypothesis. Method, candidate list, and caveats: '
  + '<a href="${REPO}">github.com/xag/asyn-oligomer-screen</a> — contributions, critiques, and wet-lab tests welcome.';

render();
</script>
</body>
</html>
`;

writeFileSync(OUT_PATH, html);
console.log(`wrote ${OUT_PATH}`);
console.log(`  protective: ${protective.length}  harmful: ${harmful.length}`);
if (missingBrain.length) console.log(`  ⚠ missing brain_access for: ${missingBrain.join(', ')}`);
else console.log('  ✓ every candidate has a brain-access justification');
