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
//   node site/build_display.mjs            # writes site/index.html
//
// Run from the repo root.

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
<title>asyn-oligomer-screen — results</title>
<style>
  :root {
    --bg: #0e1117; --panel: #161b22; --line: #2b3240; --txt: #d7dde5;
    --dim: #8b95a5; --reaches: #2ea043; --conditional: #d29922; --barrier: #f85149;
    --unknown: #6e7681; --novel: #58a6ff; --known: #8b95a5; --holdout: #bc8cff;
    --bar-prot: #3fb950; --bar-harm: #f85149;
  }
  * { box-sizing: border-box; }
  body { margin: 0; background: var(--bg); color: var(--txt);
    font: 14px/1.5 -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; }
  header { padding: 22px 28px 8px; }
  h1 { font-size: 19px; margin: 0 0 4px; }
  .sub { color: var(--dim); font-size: 13px; max-width: 860px; }
  .warn { color: var(--conditional); }
  .wrap { padding: 0 28px 60px; }
  .tabs { display: flex; gap: 8px; margin: 18px 0 6px; }
  .tab { padding: 7px 16px; border: 1px solid var(--line); background: var(--panel);
    color: var(--dim); border-radius: 7px 7px 0 0; cursor: pointer; font-weight: 600; }
  .tab.active { color: var(--txt); border-bottom-color: var(--bg); background: var(--bg); }
  .tab.harm.active { color: var(--barrier); }
  .tab.prot.active { color: var(--bar-prot); }
  .controls { display: flex; flex-wrap: wrap; gap: 14px; align-items: center;
    padding: 12px 14px; background: var(--panel); border: 1px solid var(--line);
    border-radius: 0 7px 7px 7px; margin-bottom: 14px; }
  .controls label { color: var(--dim); font-size: 12px; cursor: pointer; user-select: none; }
  .controls .grp { display: flex; gap: 9px; align-items: center; }
  .controls .grp b { color: var(--txt); font-size: 11px; text-transform: uppercase;
    letter-spacing: .04em; }
  select { background: var(--bg); color: var(--txt); border: 1px solid var(--line);
    border-radius: 5px; padding: 4px 7px; }
  table { width: 100%; border-collapse: collapse; }
  th, td { text-align: left; padding: 9px 10px; border-bottom: 1px solid var(--line);
    vertical-align: top; }
  th { color: var(--dim); font-size: 11px; text-transform: uppercase; letter-spacing: .04em;
    cursor: pointer; white-space: nowrap; user-select: none; }
  th.sorted::after { content: " \\25BC"; font-size: 9px; }
  th.sorted.asc::after { content: " \\25B2"; }
  td.rank { color: var(--dim); width: 34px; }
  td.name { font-weight: 600; }
  .badge { display: inline-block; padding: 1px 7px; border-radius: 999px; font-size: 11px;
    border: 1px solid var(--line); color: var(--dim); white-space: nowrap; }
  .prov { font-size: 10px; padding: 1px 6px; border-radius: 999px; margin-left: 6px;
    vertical-align: middle; }
  .prov.novel { color: var(--novel); border: 1px solid var(--novel); }
  .prov.known { color: var(--known); border: 1px solid var(--known); }
  .prov.holdout { color: var(--holdout); border: 1px solid var(--holdout); }
  .dot { display: inline-block; width: 9px; height: 9px; border-radius: 50%; margin-right: 6px;
    vertical-align: middle; }
  .dot.reaches { background: var(--reaches); }
  .dot.conditional { background: var(--conditional); }
  .dot.barrier { background: var(--barrier); }
  .dot.unknown { background: var(--unknown); }
  .reach-txt.barrier { color: var(--barrier); }
  .reach-txt.conditional { color: var(--conditional); }
  .barcell { width: 190px; }
  .bar { height: 9px; border-radius: 3px; }
  .bar.prot { background: var(--bar-prot); }
  .bar.harm { background: var(--bar-harm); }
  .barnum { color: var(--dim); font-size: 11px; font-variant-numeric: tabular-nums; }
  tr.mol { cursor: pointer; }
  tr.mol:hover td { background: #1b2230; }
  tr.detail td { background: #10141c; color: var(--dim); font-size: 12.5px; padding: 4px 14px 14px 44px; }
  tr.detail .k { color: var(--txt); }
  .hidden { display: none; }
  .matrix { margin: 26px 0 6px; }
  .matrix h2 { font-size: 14px; margin: 0 0 4px; }
  .matrix .msub { color: var(--dim); font-size: 12px; margin-bottom: 10px; }
  .grid { display: grid; grid-template-columns: 130px 1fr 1fr; gap: 8px; }
  .grid .h { color: var(--dim); font-size: 11px; text-transform: uppercase; letter-spacing: .04em;
    align-self: end; padding-bottom: 4px; }
  .cell { background: var(--panel); border: 1px solid var(--line); border-radius: 8px;
    padding: 11px 13px; min-height: 64px; }
  .cell.first { background: #11261a; border-color: #1f6f3a; }
  .rowlabel { color: var(--dim); font-size: 12px; display: flex; align-items: center; }
  .chip { display: inline-block; margin: 3px 5px 0 0; padding: 2px 9px; border-radius: 999px;
    font-size: 12px; background: var(--bg); border: 1px solid var(--line); }
  .legend { color: var(--dim); font-size: 11.5px; margin: 14px 0 0; line-height: 1.9; }
  .legend .dot { margin-left: 2px; }
  code { background: #1b2230; padding: 1px 5px; border-radius: 4px; font-size: 12px; }
  a { color: var(--novel); }
</style>
</head>
<body>
<header>
  <h1>α-synuclein oligomer screen — candidate display</h1>
  <div class="sub">
    Ranked outputs of the Stage 3 perturbation screen against the toxic-oligomer model,
    joined with each molecule's route to the brain. <span class="warn">These are
    unvalidated computational hypotheses, not medical or dietary advice.</span>
    Rankings order molecules; they do not measure effect size. Reachability is
    <b>not</b> in the score — it is shown alongside so a tightly-binding but
    poorly-bioavailable hit is not mistaken for an actionable one.
  </div>
</header>
<div class="wrap">
  <div class="tabs">
    <div class="tab prot active" data-panel="protective">Protective — candidates to test</div>
    <div class="tab harm" data-panel="harmful">Suspected harmful — anti-targets</div>
  </div>
  <div class="controls">
    <div class="grp"><b>Route</b><select id="routeFilter"></select></div>
    <div class="grp"><b>Reach</b><select id="reachFilter"></select></div>
    <div class="grp"><b>Evidence</b><select id="provFilter"></select></div>
    <label><input type="checkbox" id="showAll"> show all rows</label>
  </div>
  <table id="tbl"><thead></thead><tbody></tbody></table>

  <div class="matrix" id="matrix"></div>

  <div class="legend" id="legend"></div>
</div>

<script>
const DATA = ${payload};
const ROUTE_LABEL = {
  endogenous: 'endogenous', diet: 'diet', supplement: 'supplement',
  precursor: 'precursor', microbiome: 'microbiome', environmental: 'environmental',
  injection: 'injection', 'none-known': 'none known',
};
const REACH_LABEL = {
  reaches: 'reaches brain', conditional: 'conditional', barrier: 'delivery barrier',
  unknown: 'unknown',
};
const PROV_LABEL = { novel: 'novel hypothesis', known: 'known modulator', holdout: 'hold-out (sanity check)' };

let state = { panel: 'protective', route: '', reach: '', prov: '', showAll: false,
  sortKey: 'rank', sortAsc: true };

const tbl = document.getElementById('tbl');
const TOPN = 25;

function rows() {
  let rs = DATA[state.panel].slice();
  if (state.route) rs = rs.filter(r => r.route === state.route);
  if (state.reach) rs = rs.filter(r => r.reach === state.reach);
  if (state.prov) rs = rs.filter(r => r.prov === state.prov);
  const k = state.sortKey, asc = state.sortAsc ? 1 : -1;
  rs.sort((a, b) => {
    let av = a[k], bv = b[k];
    if (typeof av === 'string') return av.localeCompare(bv) * asc;
    return (av - bv) * asc;
  });
  if (!state.showAll && !state.route && !state.reach && !state.prov) rs = rs.slice(0, TOPN);
  return rs;
}

function fmtConc(r) {
  if (r.cnsLow && r.cnsHigh) return r.cnsLow + '–' + r.cnsHigh;
  if (r.cnsLow) return r.cnsLow;
  return '—';
}

function render() {
  const prot = state.panel === 'protective';
  const all = DATA[state.panel];
  const maxMag = prot
    ? Math.max(...all.map(r => -r.dActGated))
    : Math.max(...all.map(r => r.aspr));

  const cols = prot
    ? [['rank','#'],['name','molecule'],['dActGated','predicted destabilisation'],
       ['route','route to brain'],['reach','reachability'],['cns','CNS conc.']]
    : [['rank','#'],['name','molecule'],['aspr','adduct reactivity'],
       ['route','origin'],['cns','CNS conc.']];

  tbl.querySelector('thead').innerHTML = '<tr>' + cols.map(([k,l]) => {
    const sorted = state.sortKey === k;
    const cls = sorted ? 'sorted ' + (state.sortAsc ? 'asc' : '') : '';
    return '<th class="'+cls+'" data-k="'+k+'">'+l+'</th>';
  }).join('') + '</tr>';

  const body = rows().map(r => {
    const mag = prot ? -r.dActGated : r.aspr;
    const w = Math.max(3, Math.round((mag / maxMag) * 150));
    const num = prot ? r.dActGated.toFixed(3) : '+' + r.aspr.toFixed(2);
    const barCls = prot ? 'prot' : 'harm';
    const provBadge = '<span class="prov '+r.prov+'" title="'+PROV_LABEL[r.prov]+'">'
      + (r.prov==='novel'?'★ novel':r.prov==='holdout'?'⚑ hold-out':'known') + '</span>';
    const routeCell = '<span class="badge">'+ROUTE_LABEL[r.route]+'</span>';
    const reachCell = '<span class="dot '+r.reach+'"></span>'
      + '<span class="reach-txt '+r.reach+'">'+REACH_LABEL[r.reach]+'</span>';
    const barCell = '<div style="display:flex;gap:8px;align-items:center">'
      + '<div class="bar '+barCls+'" style="width:'+w+'px"></div>'
      + '<span class="barnum">'+num+'</span></div>';

    const cells = prot
      ? ['<td class="rank">'+r.rank+'</td>',
         '<td class="name">'+r.name+provBadge+'</td>',
         '<td class="barcell">'+barCell+'</td>',
         '<td>'+routeCell+'</td>',
         '<td>'+reachCell+'</td>',
         '<td class="barnum">'+fmtConc(r)+'</td>']
      : ['<td class="rank">'+r.rank+'</td>',
         '<td class="name">'+r.name+provBadge+'</td>',
         '<td class="barcell">'+barCell+'</td>',
         '<td>'+routeCell+'</td>',
         '<td class="barnum">'+fmtConc(r)+'</td>'];

    const span = cols.length;
    const detail = '<tr class="detail hidden"><td colspan="'+span+'">'
      + (r.evidence ? '<div><span class="k">why in scope:</span> '+r.evidence+'</div>' : '')
      + (r.deliveryNotes ? '<div><span class="k">'+(prot?'route detail':'source / reduce via')+':</span> '+r.deliveryNotes+'</div>' : '')
      + (r.cnsNote ? '<div><span class="k">compartment:</span> '+r.cnsNote+'</div>' : '')
      + (r.aff!=null && prot ? '<div><span class="k">top dock affinity:</span> '+r.aff.toFixed(2)+' kcal/mol</div>' : '')
      + '</td></tr>';
    return '<tr class="mol">'+cells.join('')+'</tr>'+detail;
  }).join('');
  tbl.querySelector('tbody').innerHTML = body;

  tbl.querySelectorAll('th').forEach(th => th.onclick = () => {
    const k = th.dataset.k;
    if (state.sortKey === k) state.sortAsc = !state.sortAsc;
    else { state.sortKey = k; state.sortAsc = (k === 'rank' || k === 'name'); }
    render();
  });
  tbl.querySelectorAll('tr.mol').forEach(tr => tr.onclick = () => {
    const d = tr.nextElementSibling;
    if (d && d.classList.contains('detail')) d.classList.toggle('hidden');
  });

  renderMatrix();
  renderLegend();
}

function renderMatrix() {
  const m = document.getElementById('matrix');
  if (state.panel !== 'protective') { m.innerHTML = ''; return; }
  const rs = DATA.protective;
  const mags = rs.map(r => -r.dActGated).sort((a,b)=>a-b);
  const med = mags[Math.floor(mags.length/2)];
  const reaches = r => r.reach === 'reaches';
  const strong = r => (-r.dActGated) >= med;
  const bucket = (s, re) => rs.filter(r => strong(r)===s && reaches(r)===re)
    .sort((a,b)=>a.dActGated-b.dActGated)
    .map(r => '<span class="chip">'+r.name.replace(/\\s*\\(.*\\)/,'')+'</span>').join('');
  m.innerHTML = '<h2>Actionability matrix — predicted effect × reachability</h2>'
    + '<div class="msub">Protective candidates split at the median predicted destabilisation. '
    + 'Top-left = strong predicted effect <i>and</i> reaches the brain natively or through diet — '
    + 'the first tier worth testing.</div>'
    + '<div class="grid">'
    + '<div class="h"></div><div class="h">reaches brain readily</div><div class="h">delivery barrier / conditional</div>'
    + '<div class="rowlabel">stronger effect</div>'
    + '<div class="cell first">'+bucket(true,true)+'</div>'
    + '<div class="cell">'+bucket(true,false)+'</div>'
    + '<div class="rowlabel">weaker effect</div>'
    + '<div class="cell">'+bucket(false,true)+'</div>'
    + '<div class="cell">'+bucket(false,false)+'</div>'
    + '</div>';
}

function renderLegend() {
  document.getElementById('legend').innerHTML =
    'Reachability: <span class="dot reaches"></span>reaches brain (native / dietary) '
    + '<span class="dot conditional"></span>conditional (e.g. microbiome producer status) '
    + '<span class="dot barrier"></span>delivery barrier (low bioavailability / invasive only)<br>'
    + 'Evidence: <span class="prov novel">★ novel</span> no prior α-syn binding evidence · '
    + '<span class="prov known">known</span> published modulator · '
    + '<span class="prov holdout">⚑ hold-out</span> blind sanity-check entry (not a discovery). '
    + 'Click any row for the delivery / source detail.';
}

// build filter dropdowns
function fillSelect(id, label, values, labels) {
  const sel = document.getElementById(id);
  sel.innerHTML = '<option value="">all</option>' +
    values.map(v => '<option value="'+v+'">'+(labels[v]||v)+'</option>').join('');
}
function refreshFilters() {
  const rs = DATA[state.panel];
  fillSelect('routeFilter','route',[...new Set(rs.map(r=>r.route))],ROUTE_LABEL);
  fillSelect('reachFilter','reach',[...new Set(rs.map(r=>r.reach))],REACH_LABEL);
  fillSelect('provFilter','prov',[...new Set(rs.map(r=>r.prov))],PROV_LABEL);
}

document.querySelectorAll('.tab').forEach(t => t.onclick = () => {
  document.querySelectorAll('.tab').forEach(x=>x.classList.remove('active'));
  t.classList.add('active');
  state = { ...state, panel: t.dataset.panel, route:'', reach:'', prov:'',
    sortKey:'rank', sortAsc:true };
  refreshFilters();
  document.getElementById('routeFilter').value='';
  document.getElementById('reachFilter').value='';
  document.getElementById('provFilter').value='';
  render();
});
document.getElementById('routeFilter').onchange = e => { state.route=e.target.value; render(); };
document.getElementById('reachFilter').onchange = e => { state.reach=e.target.value; render(); };
document.getElementById('provFilter').onchange = e => { state.prov=e.target.value; render(); };
document.getElementById('showAll').onchange = e => { state.showAll=e.target.checked; render(); };

refreshFilters();
render();
</script>
</body>
</html>
`;

writeFileSync(OUT_PATH, html);
console.log(
  `wrote ${OUT_PATH}\n  protective: ${protective.length} molecules\n  harmful:    ${harmful.length} molecules`,
);
