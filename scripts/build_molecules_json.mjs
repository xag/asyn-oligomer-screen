// Emit molecules.json — the candidate-molecule registry health's GET /molecules
// reads — from the canonical seed list data/vicinity_molecules.js. Primary
// entries only (source:"primary"); the broker appends source:"contributed"
// proposals to the published file at runtime, and `hf_store publish-molecules`
// merges so a reseed never clobbers them.
//
//   node scripts/build_molecules_json.mjs [out.json]   # default: molecules.json
//   .venv/bin/python -m hf_store publish-molecules --repo <repo> --file molecules.json
//
// Carries only the fields a client needs for discovery + the docking prep that
// turns an entry into runnable work (smiles); the heavier provenance stays in
// the registry source file.

import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

// The asyn repo has no package.json `type:module`, so a plain ESM import of the
// registry .js fails under node's CJS interop. The file is self-contained data
// (no imports), so we read it, strip the `export` keywords, and eval it as a
// script to capture the array — independent of module system.
const here = dirname(fileURLToPath(import.meta.url));
const src = readFileSync(join(here, '..', 'data', 'vicinity_molecules.js'), 'utf8')
  .replace(/\bexport\s+/g, '');
// eslint-disable-next-line no-new-func
const { VICINITY_MOLECULES } = new Function(`${src}\nreturn { VICINITY_MOLECULES };`)();

const out = VICINITY_MOLECULES.map((m) => ({
  id: m.id,
  name: m.name,
  source: 'primary',
  group: m.group ?? null,
  role: m.role ?? null,
  smiles: m.smiles ?? null,
  pdb_ligand: m.pdb_ligand ?? null,
  mw_da: m.mw_da ?? null,
  cns_conc: m.cns_conc ?? null,
  delivery: m.delivery ?? null,
  evidence: m.evidence ?? null,
  refs: m.refs ?? [],
}));

const path = process.argv[2] || 'molecules.json';
writeFileSync(path, JSON.stringify(out, null, 2));
console.error(`wrote ${out.length} primary molecules -> ${path}`);
