# Glossary

Plain-English definitions of the technical terms that appear in this repository.
Skim this once and the rest of the documentation should read more easily.
Linked from [README.md](README.md).

If something here is unclear or missing, open an issue — or paste the
relevant passage from a README into an AI assistant and ask it to explain.

---

## Biology and disease

**α-synuclein (alpha-synuclein, "α-syn")** — a small protein found in
brain cells, especially at the connection points between neurons
(synapses). Each molecule is a chain of 140 amino acids. In its healthy
form it helps regulate the release of neurotransmitters. In Parkinson's
disease and related diseases (Lewy body dementia, multiple system
atrophy), it misbehaves: it clumps together, and the clumps damage
brain cells.

**Aggregation** — the process of individual α-syn molecules sticking
together into larger assemblies.

**Monomer** — a single, free, unclumped α-syn molecule.

**Oligomer** — a small clump of a few (2–50ish) α-syn molecules stuck
together. The toxic oligomer is the small clump that actually damages
brain cells. The atomic structure of the toxic oligomer has never been
deposited in any public database.

**Fibril** — a large, mature clump of many α-syn molecules wound
together into a rope-like fibre. Fibrils are what you see in autopsy
images of Parkinson's brains (in structures called "Lewy bodies"). Mature
fibrils are *not* the main cause of cell damage; the small oligomers are.

**Type B\* oligomer** — a specific toxic oligomer characterised by
Fusco and colleagues in a 2017 paper. They described its general shape
(three α-syn molecules, β-sheet around residues 70–88, the rest
floppy) but did not release atomic coordinates. This work builds a model
of it.

**Lewy body** — a cellular inclusion (visible under a microscope) full
of α-syn fibrils. The pathological hallmark of Parkinson's disease.

**Parkinson's disease (PD)** — a neurodegenerative disease whose main
visible cause is the loss of dopamine-producing neurons in a brain region
called the *substantia nigra*. α-syn aggregation is causally upstream of
the cell loss.

**Familial Parkinson's** — rare inherited forms of PD caused by point
mutations in the α-syn gene (A53T, A30P, E46K, H50Q, etc.). The mutant
proteins aggregate more readily than the normal version. We use deposited
structures of mutant fibrils as the "more toxic" calibration class
because the field has measured them to be more dangerous.

**Substantia nigra** — the brain region where dopamine-producing
neurons die in Parkinson's disease.

**Dopaminergic neurons** — neurons that release dopamine. Their loss
in the substantia nigra causes Parkinson's motor symptoms.

**Blood-brain barrier (BBB)** — the filter that controls which
molecules in the bloodstream can reach the brain. A candidate molecule
that can't cross the BBB at meaningful concentrations is unlikely to be
useful against α-syn in vivo, even if it binds well in a test tube.

---

## Protein structure and chemistry

**Amino acid** — the building block of proteins. There are 20 standard
ones; each has a one-letter code (K = lysine, R = arginine, C = cysteine,
H = histidine, Y = tyrosine, etc.).

**Residue** — a single amino acid within a protein chain. "Residues
70–88" means amino acids 70 through 88 along the chain.

**β-sheet (beta-sheet)** — a flat, extended folding pattern where
multiple stretches of protein chain ("β-strands") line up side by
side, held together by hydrogen bonds. β-sheets are the structural
core of both α-syn fibrils and the toxic oligomer.

**β-strand** — a single extended stretch of protein chain that
participates in a β-sheet.

**β-core** — the part of α-syn (around residues 70–88 in the Fusco
oligomer) that folds into a β-sheet and holds the oligomer together.

**Parallel vs antiparallel β-sheet** — two β-strands can run in the
same direction (parallel) or opposite directions (antiparallel) within
a sheet. The Fusco oligomer is built as parallel by default; we test
antiparallel as an alternative.

**NAC region** — the "non-amyloid β component" region of α-syn,
residues 60–95. The most hydrophobic ("greasy") part of α-syn and the
core driver of aggregation.

**Hydrophobic / hydrophilic** — hydrophobic = water-repelling
("greasy"); hydrophilic = water-loving. Hydrophobic surfaces stick to
fats — including cell membranes, which is how the toxic oligomer damages
neurons.

**SASA (solvent-accessible surface area)** — the surface area of a
protein that water can touch (the not-buried part). Measured in
square Ångströms (Å²).

**Cα (alpha carbon)** — the central carbon atom of each amino acid;
one Cα per residue. Used as a reference point in structural calculations.

**PDB (Protein Data Bank)** — the public database of experimentally
determined atomic-resolution protein structures. PDB IDs (e.g., 6PEO)
are short codes for specific entries.

**REMARK 350** — the section of a PDB file that records the
biologically meaningful assembly of protein chains (e.g., "this fibril
has 6 chains arranged in two protofilaments"). Different from the chains
the experimentalist happened to put in the file.

**SMILES** — a compact text representation of a small molecule's
chemical structure. Used as the input format for the candidate molecules.

**Polyphenol** — a class of plant-derived molecules with multiple
phenol (benzene-ring-with-OH) groups. Includes EGCG, curcumin,
quercetin, silibinin, fisetin, etc. Most published α-syn-disrupting
molecules are polyphenols.

**Steroid hormone** — a class of molecules with a four-ring carbon
skeleton, made by the body from cholesterol. Includes DHEA,
allopregnanolone, testosterone, cortisol.

**Neurosteroid** — a steroid produced or active in the brain. DHEA and
allopregnanolone are neurosteroids.

---

## Methods and statistics

**Docking** — the computational process of predicting where a small
molecule (the *ligand*) will sit on a protein (the *receptor*), and how
tightly it will bind. AutoDock Vina is the tool we use.

**AutoDock Vina** — the open-source docking program we use. Given a
receptor and a ligand, it tries many possible binding poses and reports
the most favourable ones plus their predicted binding affinities.

**Binding affinity** — how tightly a small molecule sticks to a
protein. Reported in kcal/mol; *more negative* = *tighter binding*. A
value around −6 kcal/mol is the rough "drug-like" threshold; anything
weaker is unlikely to bind well enough in the body to matter.

**Pose** — one specific predicted geometry of a ligand sitting on a
receptor. Vina returns several poses per ligand; we average over them.

**Boltzmann weighting** — a standard physics formula that gives more
weight to favourable poses (lower energy) and less weight to unfavourable
ones, in proportion to how often they would occur at body temperature.

**kcal/mol** — kilocalories per mole. The standard energy unit in
biochemistry. Binding affinities are reported in kcal/mol.

**Δactivity (delta-activity)** — the change in our toxicity score
when a ligand binds. Δactivity = activity(complex) − activity(empty
oligomer). A *negative* Δactivity means the ligand made the oligomer
look *less* toxic — what we want.

**Δact_gated** — Δactivity with the binding-affinity penalty applied.
The primary ranking metric in this work. Negative = candidate. The more
negative, the better the candidate.

**Covalent vs reversible binding** — reversible binders sit on the
protein and can leave again (the normal case for drugs). Covalent binders
form permanent chemical bonds. The four reactive aldehydes
(malondialdehyde, acrolein, 4-HNE, methylglyoxal) bind covalently and
*damage* α-syn — they are "anti-targets" to be minimised.

**Covalent adduct** — the chemical product formed when a reactive
molecule attaches itself permanently to a protein.

**Schiff base** — a specific type of covalent adduct formed between an
aldehyde and the side chain of lysine.

**Michael adduct** — a specific type of covalent adduct formed between
an electron-poor double bond (e.g., in acrolein or 4-HNE) and the side
chain of cysteine, histidine, or lysine.

**Molecular dynamics (MD)** — a computer simulation that models how
the atoms of a protein move over time under the laws of physics. Used
here to relax the oligomer model into a physically plausible shape.

**Implicit solvent / OBC2** — a way of representing the water around a
protein in MD simulations as a smooth background rather than as
individual water molecules. Faster than explicit water; less accurate
for fine details. OBC2 is one specific implicit-solvent model.

**Topology prior** — a "starting belief" about the general shape of
the structure being modelled, used to constrain the build process. We
use the published topology of the Fusco oligomer (chain count, β-sheet
location, sheet direction) as our topology prior.

**Anchor structures** — the 14 deposited α-syn structures used to
calibrate the toxicity score (nine mature fibrils as "less toxic", five
familial-mutant fibrils as "more toxic").

**Pre-registration / blind hold-out** — writing down predictions
*before* running the analysis, so you can't unconsciously tune the
framework to pass its own test. We pre-registered which of five
literature-validated α-syn modulators should rank where, then ran the
docking without modifying anything.

**AUC (Area Under the Curve)** — a statistic that measures how well a
score separates two classes. 0.5 = chance (no separation), 1.0 = perfect
separation. Our calibration AUC of 0.84 means the score correctly
ranks a more-toxic structure above a less-toxic one in 84% of pairs.

**z-score** — a value re-scaled to be "how many standard deviations
above or below the average". Used here to put the five features on a
common scale before combining them.

**Spearman ρ (rho)** — a statistic that measures rank correlation.
+1 = perfect agreement on ordering, 0 = no correlation, −1 = perfect
disagreement.

---

## Bench-side assays (lab techniques the screen would feed into)

**ThT (thioflavin T) assay** — a standard lab test for α-syn
aggregation. ThT is a fluorescent dye that glows brightly when it binds
to amyloid (β-sheet-rich) structures. If a candidate molecule slows
the rise of ThT fluorescence over time, it slows aggregation.

**DLS (dynamic light scattering)** — a technique that measures the
size distribution of particles (here, α-syn aggregates) in solution.
Used to detect oligomers vs fibrils.

**LC-MS (liquid chromatography–mass spectrometry)** — an analytical
technique that identifies chemical modifications on a protein. Used to
verify covalent adducts (e.g., does methylglyoxal actually attach to
α-syn lysines in a cell?).

**In vitro / in vivo** — *in vitro* = in a test tube or dish (purified
protein, cultured cells). *In vivo* = in a living organism (mouse,
human).

---

## Project-specific terms

**Anti-target** — a molecule we want to *reduce* exposure to, not
increase. The reactive aldehydes (MDA, acrolein, 4-HNE, MGO) are the
anti-targets in this work. The project has equal weight on what to seek
(destabilisers) and what to avoid (anti-targets).

**Destabiliser** — a molecule that, by binding to the toxic
oligomer, makes it less stable / more likely to fall apart. The desired
property for a candidate.

**Vicinity molecules** — the candidate list (`data/vicinity_molecules.js`).
Named because the inclusion rule is "plausibly reaches the substantia
nigra at non-trivial concentrations" — they live in the "vicinity" of
the target by way of diet, supplementation, endogenous synthesis, or
lifestyle exposure.

**Validation hold-out** — five molecules (silibinin, EGCG, fisetin,
rosmarinic acid, CAPE) the framework had not seen during development.
Used as the pre-registered blind test.

---

## When in doubt

Paste the unclear passage into an AI assistant (Claude, ChatGPT, Copilot,
Gemini) along with this glossary, and ask for an explanation tailored to
your background. Then send an issue if the glossary should be expanded
to cover the term that confused you.
