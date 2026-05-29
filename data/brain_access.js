// Per-candidate brain-access justification layer for the candidate list.
//
// The original inclusion gate for data/vicinity_molecules.js ("plausibly reaches
// the substantia nigra") was never enforced per molecule. This file makes that
// gate explicit and auditable: every docked candidate must carry a concrete,
// real reason it reaches the brain, and a verdict. Molecules that do not reach
// the brain at useful levels are demoted out of the "worth testing" view (they
// stay visible, flagged as a delivery problem).
//
// `route` is written in PLAIN ENGLISH — it is shown to readers, so it avoids
// jargon (no "BBB", "LAT1", "glucuronidation"). The pharmacology underneath is
// standard absorption / blood-brain-barrier knowledge, not a claim about
// α-synuclein activity. "Too little reaches" and "doesn't reach" are honest
// findings, not failures.
//
// verdict (drives display grouping):
//   crosses        taking it in actually raises the level in the brain
//   endogenous     the body already makes it inside the brain; the lever is
//                  metabolic, not swallowing it (often can't be raised from outside)
//   limited        some gets in, but only a little / uncertain
//   subtherapeutic gets in only far below the amount the effect would need
//   does-not-reach doesn't get into the brain intact by any practical route

export const BRAIN_ACCESS = {
  // ---- fat-soluble small molecules that genuinely get in ----
  thc: { verdict: 'crosses', route: 'Fat-soluble, so it slips straight into the brain — the fact that it gets you high is direct proof it arrives there.' },
  cannabidiol: { verdict: 'crosses', route: 'Fat-soluble; passes easily into the brain after oil or capsule.' },
  caffeine: { verdict: 'crosses', route: 'Small and fat-soluble; reaches the brain almost completely within minutes of a coffee or tea.' },
  nicotine: { verdict: 'crosses', route: 'Small and fat-soluble; in the brain within seconds when inhaled, minutes from a patch or gum.' },
  theobromine: { verdict: 'crosses', route: 'The caffeine-like compound in cocoa; small enough to pass into the brain after eating chocolate.' },
  theophylline: { verdict: 'crosses', route: 'A caffeine relative found in tea; passes into the brain the same way caffeine does.' },
  honokiol: { verdict: 'crosses', route: 'A fat-soluble magnolia compound — one of the few plant polyphenols that genuinely gets into the brain.' },
  pterostilbene: { verdict: 'limited', route: 'A more stable, fattier cousin of resveratrol, so a modest amount reaches the brain where resveratrol itself essentially does not.' },
  piperine: { verdict: 'crosses', route: 'The pungent compound in black pepper; gets into the brain itself, and also slows the gut/liver from clearing other compounds taken with it.' },
  sulforaphane: { verdict: 'crosses', route: 'A small fat-soluble compound from broccoli sprouts; passes into the brain and switches on the cell’s antioxidant defences.' },
  melatonin: { verdict: 'crosses', route: 'Small and fat-soluble; passes freely into the brain, so a bedtime dose reaches it directly.' },
  'retinoic-acid': { verdict: 'crosses', route: 'The active form of vitamin A; the body makes it from dietary vitamin A and it is fat-soluble enough to enter the brain.' },
  thiamine: { verdict: 'crosses', route: 'Vitamin B1; carried into the brain by a dedicated uptake system, and high-dose or the fat-soluble benfotiamine form raises brain levels.' },
  nicotinamide: { verdict: 'crosses', route: 'A B3 vitamin form; gets into the brain and tops up cellular energy molecules when supplemented.' },
  ascorbate: { verdict: 'limited', route: 'Vitamin C; pumped into the brain by a dedicated carrier, but brain levels are held steady and max out — hard to push much higher.' },
  bmaa: { verdict: 'crosses', route: 'An algae-derived amino-acid mimic; sneaks into the brain on the carrier meant for normal amino acids — which is exactly why it is neurotoxic.' },
  hydroxytyrosol: { verdict: 'limited', route: 'A small olive-oil compound; absorbed better than most polyphenols and partly enters the brain, but levels stay low.' },

  // ---- polyphenols that rank well on binding but do NOT reach useful brain levels ----
  silibinin: { verdict: 'subtherapeutic', route: 'From milk thistle; very little survives the gut and liver, and only trace amounts reach the brain — far below what works in a dish. Special absorption-boosted forms help only a little.' },
  egcg: { verdict: 'subtherapeutic', route: 'The main green-tea catechin; barely absorbed and quickly broken down. Brain levels end up roughly 100× below the amount that works in lab tests.' },
  curcumin: { verdict: 'subtherapeutic', route: 'From turmeric; cleared by the body almost as fast as it is absorbed. Even with black-pepper, liposomal, or nanoparticle tricks, brain levels stay too low to matter.' },
  demethoxycurcumin: { verdict: 'subtherapeutic', route: 'A turmeric compound with the same poor absorption as curcumin; very little reaches the brain.' },
  baicalein: { verdict: 'subtherapeutic', route: 'A Scutellaria (skullcap) flavone; poorly absorbed and quickly cleared, leaving only low brain levels.' },
  fisetin: { verdict: 'subtherapeutic', route: 'A flavonoid from strawberries; gets into the brain a bit better than most of its class, but absorption is still low and brain levels likely fall short.' },
  naringenin: { verdict: 'subtherapeutic', route: 'A citrus flavonoid; poorly absorbed and quickly cleared.' },
  luteolin: { verdict: 'subtherapeutic', route: 'A flavonoid in herbs and vegetables; poorly absorbed, only low brain levels.' },
  myricetin: { verdict: 'subtherapeutic', route: 'A widespread flavonoid; very poorly absorbed and quickly cleared.' },
  hesperetin: { verdict: 'subtherapeutic', route: 'A citrus flavonoid; low absorption keeps brain levels low.' },
  genistein: { verdict: 'subtherapeutic', route: 'A soy isoflavone; absorbed but largely tagged for removal, so brain levels are low.' },
  daidzein: { verdict: 'subtherapeutic', route: 'A soy isoflavone; mostly cleared before reaching the brain (its gut-made product equol matters more — and only some people make it).' },
  kaempferol: { verdict: 'subtherapeutic', route: 'A common flavonoid; poorly absorbed, low brain levels.' },
  epicatechin: { verdict: 'subtherapeutic', route: 'A cocoa/tea flavanol; circulates mostly in a deactivated form, so little active compound reaches the brain.' },
  quercetin: { verdict: 'subtherapeutic', route: 'A very common flavonoid; almost all of it is tagged and cleared, so free quercetin in the brain is very low.' },
  apigenin: { verdict: 'subtherapeutic', route: 'A flavonoid from parsley and chamomile; poorly soluble and poorly absorbed.' },
  'rosmarinic-acid': { verdict: 'subtherapeutic', route: 'A rosemary/herb compound; poorly absorbed and quickly cleared, with very little reaching the brain.' },
  resveratrol: { verdict: 'subtherapeutic', route: 'The red-wine compound; almost entirely deactivated within minutes of absorption, so virtually none reaches the brain in active form.' },
  'caffeic-acid': { verdict: 'subtherapeutic', route: 'A coffee/plant acid; poorly absorbed and quickly cleared.' },
  'ferulic-acid': { verdict: 'subtherapeutic', route: 'A whole-grain plant acid; some is absorbed but little reaches the brain.' },
  'gallic-acid': { verdict: 'subtherapeutic', route: 'A water-loving tannin acid; struggles to cross into the brain.' },
  cape: { verdict: 'subtherapeutic', route: 'A propolis (bee glue) compound; broken down to caffeic acid and poorly absorbed by mouth, so little reaches the brain.' },

  // ---- does not reach intact ----
  trehalose: { verdict: 'does-not-reach', route: 'A sugar that the gut splits into glucose before absorption — almost none gets into the blood intact, and it cannot enter the brain. Its anti-clumping effect in animals works through the gut and body, not by reaching the brain.' },
  mannitol: { verdict: 'does-not-reach', route: 'A sugar-alcohol the brain barrier actively keeps out — it is used in hospitals precisely because it stays out of the brain.' },

  // ---- body-made molecules that DO get in from outside ----
  dhea: { verdict: 'crosses', route: 'A fat-soluble hormone the body makes; passes into the brain, and a supplement raises brain levels. It naturally declines with age.' },
  allopregnanolone: { verdict: 'crosses', route: 'A fat-soluble brain steroid; gets in, and its precursor pregnenolone (taken by mouth) raises how much the brain makes.' },
  'l-dopa': { verdict: 'crosses', route: 'The textbook example: dopamine itself cannot get into the brain, but L-DOPA rides an amino-acid shuttle across and is turned into dopamine inside — the basis of Parkinson’s medication.' },
  inosine: { verdict: 'crosses', route: 'Taken by mouth, it raises the antioxidant urate in the brain and spinal fluid — actually tested this way in Parkinson’s trials.' },
  tryptophan: { verdict: 'crosses', route: 'A dietary amino acid that rides an amino-acid shuttle into the brain (competing with others in a protein meal); it is the raw material for serotonin.' },
  kynurenine: { verdict: 'crosses', route: 'Rides the amino-acid shuttle into the brain; diet and inflammation set how much arrives.' },
  glutamine: { verdict: 'crosses', route: 'A dietary amino acid taken up into the brain as a major fuel and building block.' },
  creatine: { verdict: 'limited', route: 'Carried into the brain by its own shuttle; a supplement raises brain creatine, but slowly and only modestly.' },
  'acetyl-l-carnitine': { verdict: 'crosses', route: 'A fattier form of carnitine that gets into the brain (plain carnitine barely does); taken as a supplement and studied for the brain.' },
  carnitine: { verdict: 'limited', route: 'Water-loving, so it barely enters the brain — the acetyl form (ALCAR) is the practical way in.' },
  taurine: { verdict: 'limited', route: 'Carried into the brain by its own shuttle; a supplement raises brain levels modestly.' },
  'choline-endo': { verdict: 'crosses', route: 'Carried into the brain on a dedicated shuttle; dietary choline (eggs, soy) feeds the brain’s acetylcholine supply.' },
  lactate: { verdict: 'crosses', route: 'Carried into the brain as a fuel; rises with exercise.' },
  'beta-hydroxybutyrate': { verdict: 'crosses', route: 'A ketone carried into the brain as fuel; a ketogenic diet or ketone drinks raise brain levels substantially — a real dietary lever.' },
  betaine: { verdict: 'crosses', route: 'A nutrient (from beets) carried into the brain; reaches it after a supplement or food.' },
  pyruvate: { verdict: 'limited', route: 'A fuel molecule carried into the brain in small amounts.' },
  succinate: { verdict: 'limited', route: 'An energy-cycle molecule that barely crosses into the brain; the brain mostly makes its own.' },
  citrate: { verdict: 'limited', route: 'An energy-cycle molecule; only some crosses in, and the brain makes its own.' },
  'alpha-ketoglutarate': { verdict: 'limited', route: 'An energy-cycle molecule with limited entry into the brain.' },
  guanosine: { verdict: 'limited', route: 'A building-block molecule taken into the brain in modest amounts.' },
  urate: { verdict: 'limited', route: 'A natural antioxidant; the brain level is a fraction of the blood level and can be nudged up via oral inosine, but not freely controlled.' },
  spermidine: { verdict: 'limited', route: 'A compound in wheat germ and aged cheese, studied for cellular "self-cleaning"; only a little reaches the brain, but some effect is reported.' },
  spermine: { verdict: 'limited', route: 'A natural cell compound present in the brain; little crosses in from outside.' },
  putrescine: { verdict: 'limited', route: 'A natural cell compound; only a little enters the brain from outside.' },
  homocysteine: { verdict: 'endogenous', route: 'Not something you take — its brain level rises when B-vitamins (B12, folate, B6) run low, so the lever is B-vitamin status.' },
  glutathione: { verdict: 'subtherapeutic', route: 'The body’s main antioxidant, but swallowed glutathione is mostly destroyed in the gut and barely enters the brain. The realistic lever is its building block NAC, which raises brain levels only a little.' },

  // ---- body-made, but you CANNOT raise them from outside (don't cross) ----
  dopamine: { verdict: 'does-not-reach', route: 'Cannot get into the brain at all — the whole reason Parkinson’s patients take L-DOPA instead. Only the brain’s own supply counts.' },
  norepinephrine: { verdict: 'does-not-reach', route: 'Cannot cross into the brain; the brain makes its own from dopamine.' },
  epinephrine: { verdict: 'does-not-reach', route: 'Cannot cross into the brain; it is a body hormone, with only a small brain supply made locally.' },
  'glutamate-endo': { verdict: 'does-not-reach', route: 'Kept out of the brain; the brain makes and recycles its own — dietary glutamate (MSG) does not raise it.' },
  'gaba-endo': { verdict: 'does-not-reach', route: 'Barely crosses into the brain; the brain makes its own. Swallowed GABA acts mostly in the body, not the brain.' },
  hva: { verdict: 'endogenous', route: 'A dopamine breakdown product made in the brain — a meter of dopamine turnover, not something you take.' },
  dopac: { verdict: 'endogenous', route: 'A dopamine breakdown product formed inside neurons — a turnover marker, not something taken.' },
  mhpg: { verdict: 'endogenous', route: 'A noradrenaline breakdown product made in the brain — a turnover marker.' },
  '3-methoxytyramine': { verdict: 'endogenous', route: 'A dopamine breakdown product formed in the brain — not something taken.' },
  '3-hydroxykynurenine': { verdict: 'endogenous', route: 'Made inside the brain (more so during inflammation); the lever is the pathway upstream, not intake.' },
  'kynurenic-acid': { verdict: 'endogenous', route: 'Barely crosses into the brain — the brain makes its own from kynurenine, so swallowing it is not the route.' },
  'quinolinic-acid': { verdict: 'endogenous', route: 'Made inside the brain by activated immune cells; does not get in well from outside.' },
  aminochrome: { verdict: 'endogenous', route: 'Forms inside dopamine neurons as dopamine oxidises — it cannot be taken; it arises on the spot.' },
  xanthine: { verdict: 'endogenous', route: 'A normal breakdown molecule present in the brain — a metabolic step, not a delivery target.' },
  hypoxanthine: { verdict: 'endogenous', route: 'A normal breakdown molecule recycled within the brain.' },

  // ---- signalling molecules ----
  adenosine: { verdict: 'limited', route: 'Already abundant in the brain; you influence it indirectly (caffeine blocks its receptors) rather than by swallowing it.' },
  serotonin: { verdict: 'does-not-reach', route: 'Cannot cross into the brain; the brain makes its own from tryptophan or 5-HTP — those precursors are the lever, not serotonin itself.' },
  histamine: { verdict: 'does-not-reach', route: 'Cannot cross into the brain; the brain makes its own.' },
  acetylcholine: { verdict: 'does-not-reach', route: 'Cannot cross into the brain and is destroyed within moments; the lever is its building block choline, not acetylcholine itself.' },
  gaba: { verdict: 'does-not-reach', route: 'Barely crosses into the brain; swallowed GABA acts mostly on the gut, not the brain.' },
  glycine: { verdict: 'limited', route: 'An amino acid that enters the brain, but levels are tightly controlled, so high doses raise it only a little.' },
  anandamide: { verdict: 'endogenous', route: 'A natural "bliss" molecule the brain makes on demand and destroys within seconds; raised indirectly by blocking its breakdown, not by swallowing it.' },
  '2-ag': { verdict: 'endogenous', route: 'A natural cannabis-like molecule the brain makes on demand — not something you take.' },
  '5-hiaa': { verdict: 'endogenous', route: 'A serotonin breakdown product made in the brain — a turnover marker.' },
  pea: { verdict: 'crosses', route: 'A natural fat-like calming molecule; the micronised supplement form reaches the brain and is used that way.' },
  'hydrogen-sulfide': { verdict: 'endogenous', route: 'A gas the brain makes itself; the lever is foods that feed its production (garlic, cysteine, NAC), not the gas directly.' },

  // ---- gut-bacteria molecules ----
  'urolithin-a': { verdict: 'limited', route: 'Made from pomegranate/walnut compounds by gut bacteria in only about 40% of people; a supplement raises blood levels, but only a little reaches the brain — and only if your gut makes it.' },
  equol: { verdict: 'limited', route: 'Made from soy by gut bacteria in only some people (~30–50%); a little reaches the brain in those who produce it.' },
  'indole-3-propionate': { verdict: 'crosses', route: 'A fat-soluble molecule gut bacteria make from tryptophan; it gets into the brain (studied as protective) — if you carry the right gut bacteria.' },
  'indole-3-acetate': { verdict: 'limited', route: 'A gut-bacteria tryptophan product; some reaches the brain, depending on your microbiome.' },
  tryptamine: { verdict: 'crosses', route: 'A gut/trace molecule that is fat-soluble and gets into the brain, but is broken down within minutes — so exposure is brief.' },
  skatole: { verdict: 'crosses', route: 'A fat-soluble molecule from gut protein breakdown; gets into the brain.' },
  indole: { verdict: 'crosses', route: 'A small fat-soluble gut-bacteria molecule; gets into the brain.' },
  butyrate: { verdict: 'limited', route: 'A fibre-fermentation fat made by gut bacteria; a small fraction reaches the brain, as most is used up by the gut and liver first.' },
  propionate: { verdict: 'limited', route: 'A fibre-fermentation fat from gut bacteria; partly reaches the brain.' },
  acetate: { verdict: 'crosses', route: 'A fibre-fermentation fat from gut bacteria; carried into the brain and used as fuel.' },
  valerate: { verdict: 'limited', route: 'A gut-bacteria fat; a modest amount reaches the brain.' },
  hippurate: { verdict: 'subtherapeutic', route: 'A water-loving by-product of gut polyphenol breakdown; barely crosses into the brain — mostly seen in urine.' },
  'p-cresol-sulfate': { verdict: 'limited', route: 'A gut by-product that reaches the brain mainly when blood levels run high (e.g. with kidney problems).' },
  phenylacetate: { verdict: 'limited', route: 'A gut/metabolic acid; a modest amount reaches the brain.' },
  cadaverine: { verdict: 'limited', route: 'A gut-bacteria amine; little reaches the brain.' },
  tmao: { verdict: 'limited', route: 'A gut by-product (from red meat/egg compounds); a small amount reaches the brain and has been found in spinal fluid.' },

  // ---- dietary fats (fats do get into the brain) ----
  'arachidonic-acid': { verdict: 'crosses', route: 'A membrane fat taken into the brain from the blood; a major brain building and signalling fat.' },
  'oleic-acid': { verdict: 'crosses', route: 'The main olive-oil fat; taken into the brain and also made there.' },
  'palmitic-acid': { verdict: 'crosses', route: 'A common dietary fat taken into the brain and also made there.' },

  // ---- minerals (tightly controlled) ----
  zinc: { verdict: 'limited', route: 'Carried into the brain by dedicated shuttles, but brain zinc is kept within a narrow range, so diet shifts it only a little.' },
  magnesium: { verdict: 'limited', route: 'Ordinary magnesium pills barely raise brain magnesium; the special "magnesium L-threonate" form was made specifically to get into the brain.' },
  manganese: { verdict: 'crosses', route: 'Enters the brain on iron shuttles; too much (diet or industrial exposure) builds up and is toxic — a real but double-edged route.' },
  calcium: { verdict: 'limited', route: 'Brain calcium is tightly controlled; eating more does not meaningfully raise it.' },

  // ---- environmental ----
  paraquat: { verdict: 'crosses', route: 'A weedkiller that sneaks into the brain on normal nutrient shuttles and damages dopamine neurons — an exposure to avoid, not a benefit.' },

  // ====================================================================
  // anti-targets (covalent channel) — made on the spot inside the brain
  // ====================================================================
  malondialdehyde: { verdict: 'endogenous', route: 'Not eaten as such — it forms inside brain tissue when fats are damaged by oxidation; the lever is lowering that oxidative damage.' },
  '4-hne': { verdict: 'endogenous', route: 'Forms inside the brain when omega-6 fats are oxidised; rises with oxidative stress, so the lever is antioxidant status and fat quality.' },
  acrolein: { verdict: 'endogenous', route: 'Forms inside the brain from fat and amine breakdown, and is also delivered directly by tobacco smoke and overheated cooking oils.' },
  methylglyoxal: { verdict: 'endogenous', route: 'A sugar-handling by-product made throughout the body, including the brain; it rises with high blood sugar, so the lever is glucose control.' },
  'nitric-oxide': { verdict: 'endogenous', route: 'A gas the brain makes itself; dietary nitrate (beets, leafy greens) feeds the body’s supply, but the brain controls its own.' },
};
