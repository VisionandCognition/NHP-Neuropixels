# NHP-Neuropixels

Macaque Neuropixels ephys analysis. Main active pipeline: `SaccadeDirection/` (subject m032 / "Kwibus", 8-target saccade task). See `SaccadeDirection/PIPELINE_REPORT.md` for full pipeline detail — don't duplicate that here.

Per-channel anatomy labels come from `SaccadeDirection/anatomy/anatomy_template_clusters_hierarchical.csv` (hand-assigned against CHARM/SARM atlas by eye, not verified histology), looked up via `SaccadeDirection/code/getChannelArea.m`. Hierarchy columns `L6` (finest) → `L3` (coarsest) are CHARM (cortex) / SARM (subcortex) parcel names.

Don't run git commands in this project — user handles git themselves.

## VisualizeOnAreas — anatomical visualization of saccade metrics (in progress, 2026-07-22/23)

Goal: plot saccade-related metrics (tuning %, decoding accuracy, planning/execution recruitment) across recorded areas in an anatomically-organized 2D layout, similar in spirit to `VisualizeOnAreas/modha-singh-network-layout/` — eventually also as a volumetric CHARM/SARM parcel map in the animal's own aligned atlas space.

### The Modha-Singh box layout (what already exists there)
`modha-singh-network-layout/make_modha_mds_box_layout.m` builds a grid-box layout of 383 macaque brain regions from **tracer-based** connectivity (Modha & Singh 2010, PNAS — CoCoMac-collated axonal tract-tracing data, not diffusion MRI; `sd01.txt`=nodes, `sd02.txt`=6602 directed edges, `sd03.txt`=hierarchy). Node neighbors in the box grid ≈ real anatomical connectivity neighbors, via MDS on connectivity profiles + soft AP/DV anatomical anchors (`region_constraints.csv`) + a refinement pass. Includes both cortical and subcortical nodes (thalamic/basal_ganglia/diencephalon/"lower" major groups), not just cortex.

### Crosswalk status: our CHARM/SARM labels → Modha node acronyms
Compared our recorded-area finest labels (`L6` column) against Modha's 383 node acronyms. Result: **6 exact matches, ~5 mechanically splittable/aggregatable (Modha is finer- or coarser-grained for the same region), ~5 fuzzy (need manual disambiguation), 8 true gaps** (region has no tracer-network node at all — most importantly **superior colliculus (SC, 13 recorded rows)** and **habenula (Hb, 15 rows)**; confirmed absent by grep, not a naming issue). True gaps need an "off-network" placement policy (anatomical-position anchor, not connectivity-derived) since there's no tracer data to place them by.

**Gotcha found:** `Pu`/`Pu_r`/`Pu_c` in Modha's node list are **Putamen** (basal_ganglia major_group), not Pulvinar — a false-friend acronym collision. The real pulvinar cluster is `Pul#1`, `Pul.o`, `PIl`, `PIl-s`, `PLvl`, `PLvm`, `PLd`, `PLa#1`, `PI#3`, `PIc`, `PIp`, `PIm` (all spatially clustered next to V1/MT/LGN, as expected). Don't trust acronym string-matching alone for subcortex — cross-check `major_group` and spatial position too.

**Resolved (2026-07-27):** legend found — PNAS SI Appendix (`sapp1.pdf`, in user's ~/Downloads; not the 6-page main-text PDF already in this folder, which has no Table S1). Table S1 = "Hierarchical brain map in depth-first order" with full names + merged-region lists; extracted to `table_s1_acronym_index.txt` in this folder. Confirms:
- Pulvinar cluster (Pul#1, PI#3, PIl/PIl-s/PIp/PIm/PIc, PM#3/PMm/PMl, Pul.o, PL#3/PLa#1/PLd/PLvl/PLvm) — all genuine pulvinar subdivisions, full names all "Nucleus pulvinaris ...".
- `Pu`/`Pu_r`/`Pu_c` confirmed = Putamen (rostral/caudal), not pulvinar — matches the false-friend gotcha above.
- `CM#2` = "nucleus centrum medianum (thalamus)" — clean, unambiguous match for our recorded CM.
- **SC (superior colliculus) is not its own vertex** — it's one of many regions merged into the coarse `MB#2` "Mid Brain" node (along with APT, CG, PPN, VTA, etc.). So SC has an approximate coarse-level anchor (MB#2's position) if needed, rather than being a total gap.
- **Habenula confirmed a total gap** — zero mentions anywhere in the 57-page SI Appendix (checked "habenula" and "epithalam*"), not even merged into a coarser node. No tracer-derived placement possible at any granularity; needs a hand-anchored anatomical position.

### MDS re-run gotcha
Re-running `make_modha_mds_box_layout.m` from scratch (no cached MDS) can converge to a **degenerate near-1D embedding** (mds_x/mds_y collapse to ~constant, stress deceptively near-zero) when the first `mdscale` attempt fails ("points co-located") and falls back to a random-init retry that finds a bad local optimum. Don't trust freshly-computed `mds_x`/`mds_y`/`refined_x`/`refined_y` without checking their spread first (real runs should show real 2D variance, e.g. compare against V1 vs FEF vs Amyg positions).

Workaround used this session: the tracked reference `figures/modha_boxed_layout_refined.svg` has good (non-degenerate) box positions from a prior run. Its label positions are recoverable by regex despite Batik's split-tag SVG formatting:
```python
# <g transform="translate(X,Y)" ...> ... <text ...>ACRONYM</text  (note: closing '>' of <text> tag is on the NEXT line, don't require it in the regex)
re.finditer(r'<g transform="translate\(([\-\d.]+),([\-\d.]+)\)"', svg)
# then search next ~500 chars for: r'<text[^>]*>([^<]*)</text'
```

### Atlas files (for the later volumetric/3D step)
User will point to the correct/current atlas files directly — don't re-search the filesystem for CHARM/SARM atlas volumes. (For reference only, not confirmed current/correct: CHARM/SARM NIfTI + key files were seen under `/home/chris/MATLAB Add-Ons/Apps/AtlasQuery/` and a NIN Sharepoint "VandC-NHP Electrode Planning" share with per-monkey warped atlas folders — none of the monkey folders there were obviously "Kwibus"/m032 by name.)

### Full crosswalk re-check against Table S1 (2026-07-27)
Went back through all 43 distinct L6 labels against Table S1 full names now that the legend is available. Several previously-"fuzzy"/"gap" calls tightened up:
- **Newly clean/mechanical matches**: APul→`Pul.o`, LPul→`PL#3`, CM-pf→`CM#2`, MD→`MD`, MG→`MG`, SG→`SG`, VPI→`VPI`, VPM-VPL→`VPM`+`VPL`, BA29→`29`, PCC/BA23C→`23`, PE[BA5]→`PEm`, PEa[BA5]→`MIP`, S1[BA1-2]→`1#1`+`2#1`, S1[BA3]→`3a`+`3b`, VIP→`VIP`, SN→`SN`, Cd(tail)→`Cd`/`Cd t`, **ZI→merged into `Sub.Th`** (not its own node).
- **Another false-friend acronym collision caught**: `CL-pc Thalamus` (our recorded label) ≠ `CL#4` (Modha's `CL#4` = caudal-lateral auditory belt cortex). Real match is `Cl#2` = "Nucleus centralis lateralis thalami."
- **SC is not a true gap** — merged into coarse `MB#2` "Mid Brain" node (see Resolved note above), can anchor there.
- **Actual true gaps remaining** (zero mentions anywhere in the 57-page SI, not even merged into a coarser node): **Habenula (Hb)**, Suprafascicular nucleus (SPFC), mesencephalic reticular formation (mRt — plausibly foldable into `MB#2` like SC/reticular structures, but not explicitly listed there, so treat as unresolved), precommissural nucleus / posterior commissure (PCR), Fornix (f — white matter tract, expected absent, not a placement target). So down from "8 true gaps" to really just Hb as the one that matters for anatomical anchoring (SPFC/mRt/PCR are each ≤1 recorded row; Fornix isn't a placement target).

### Gap-area occurrence check (2026-07-27) — decided how much placement effort each gap deserves
Counted rows/channels in `anatomy_template_clusters_hierarchical.csv` (122 clusters, 4608 channels total) matching each gap label. Many of these rows are themselves uncertain/disjunctive labels ("SC or mRt", "Hb or PCR or CM"), so true counts may be lower once resolved:
- **Hb: 15 rows, 600 ch (13.0%)** — big enough to need a real placement decision.
- Fornix: 7 rows, 184 ch (4.0%) — not a placement target anyway, it's white matter.
- mRt: 2 rows, 140 ch (3.0%); PCR: 2 rows, 40 ch (0.9%); SPFC: 1 row, 20 ch (0.4%) — small, and low-confidence labels to begin with. **Decision: leave these as gaps, not worth building placement logic for ~8% combined of mostly-uncertain rows.**

### Habenula placement decision (2026-07-27)
Habenula has zero tracer-network representation at any granularity (confirmed absent from the entire 57-page SI, not merged into any node's bracket list — unlike SC, which is explicitly listed inside `MB#2`'s merged-region list). So there's no data-driven anchor for it, coarse or fine.

**Decision: anchor Hb to the `Tha` (thalamus) box position** from `results/modha_boxed_layout_coordinates.csv`. Rationale: `Tha` is a real plotted box in the layout (hierarchy_depth=2, `is_anatomically_anchored=1`, position fixed via `region_constraints.csv`, not MDS-derived — see its row: `anchor_notes="thalamus root weak"`), one level up from all the individual thalamic nuclei (MD, CM#2, Pul.o, etc.) that Hb sits adjacent to anatomically. This keeps Hb visually next to its real neighbors instead of dropped in a generic diencephalon blob.

**Caveat user flagged (2026-07-27): this is a coarse approximation, acknowledged as imperfect** — habenula is epithalamus, not thalamus proper, and `Tha` itself has weak/placeholder connectivity data (`in_degree=0`). Accepted only because there is no better option (no tracer data exists for Hb at any level). **Any figure/output using this placement must visually flag Hb as coarse-level-anchored, not a real leaf-node match** (e.g. different box style/hatching from areas with genuine Table S1 matches), so it isn't mistaken for a resolved mapping.

### Metrics joined to crosswalk (2026-07-27)
Pulled the 3 pooled-by-area outputs (labelField=L6) from `saccadeDataDir()` (v7.3 .mat, needed `h5py` not `scipy.io`): `tuning_pooled_by_area_L6.mat`, `decode_by_area_L6.mat`, `planexec_pooled_by_area_L6.mat`. Joined all three against the crosswalk into `VisualizeOnAreas/area_metrics_by_modha_node_L6.csv` (26 pooled-area rows, one per distinct sanitized L6 label — e.g. `AREA_23`, `LPUL`, `HB`).

New finding while joining: the pooled tuning data has a **`CM` bucket distinct from `CMN-PF`** — "central medial thalamus", not centromedian/parafascicular. Added to the crosswalk: `central medial thalamus (CM)` → `C#4` ("Nucleus centralis thalami", leaf_mechanical_aggregate — `CeM#2` is in its merged-region bracket list). Don't confuse with `CM-pf Thalamus`→`CM#2`.

Also excluded pooled label **`CD`** (distinct from `TAIL`) from the crosswalk — its pooled `FullName` is "Ventricle or Caudate", i.e. itself an unresolved disjunctive label, not a clean single-region match. `TAIL` (Caudate tail proper) maps fine to `Cd t`.

Result: 20 of 26 pooled areas have real Modha nodes (leaf_exact / leaf_mechanical_aggregate / leaf_split / coarse_merged), 1 is the flagged coarse_anchor (Hb→`Tha`), 4 are gap_unplaced/gap_excluded (F, MRT, PCR, SPFC) and correctly carry no Modha node, 1 (CD) excluded for label ambiguity.

`VisualizeOnAreas/area_to_modha_crosswalk.csv` updated with the `central medial thalamus (CM)`→`C#4` row so the two files stay in sync.

### Renderer (2026-07-27, first draft)
`VisualizeOnAreas/render_area_metrics.py` → `VisualizeOnAreas/figures/area_metrics_layout.png`. Reads `results/modha_boxed_layout_coordinates.csv` (box positions) + `area_metrics_by_modha_node_L6.csv` (metrics), draws 3 side-by-side panels (% directionally tuned, decoding accuracy execution-epoch, % execution-recruited), each its own single-hue sequential colorbar (light→dark blue, not a preset multi-hue map), cropped to the bounding box of matched areas (not the full 383-node grid) for legibility.

Visual encoding (per user's requirement that placement confidence stay visible, not just the metric):
- solid black border = real leaf-node tracer match
- light diagonal hatch = coarse-merged (SC→`MB#2`) — real tracer data, coarser granularity
- heavy cross-hatch + dashed border = coarse-anchored (Hb→`Tha`) — **no tracer data at all**, per the 2026-07-27 placement decision above
- faint unfilled outline = other Modha boxes nearby with no recorded metric (context only)
- excluded gaps (CD, F, MRT, PCR, SPFC) listed in a caption line, not drawn anywhere

User reviewed the first draft (2026-07-27): labels/fonts are small since matched areas are sparse across the full grid — acceptable for now, "will probably want to make more edits later." Not polished, just functional.

### Next steps (not yet done)
1–6. Legend lookup, full re-check, gap-placement decisions, crosswalk table, metrics join, and first-draft renderer are all done — see sections above.
7. Iterate on the renderer's visual polish (label size/placement, layout) per user request, next time this is picked up.
8. Later: volumetric CHARM/SARM parcel painting in the animal's native atlas space, once atlas files are confirmed.
