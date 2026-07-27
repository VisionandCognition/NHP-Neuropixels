"""Render saccade-direction metrics onto the Modha-Singh box layout.

Reads:
  - modha-singh-network-layout/results/modha_boxed_layout_coordinates.csv  (383 boxes, positions)
  - area_metrics_by_modha_node_L6.csv                                      (26 pooled areas -> metrics)
Writes:
  - figures/area_metrics_layout.png

Color: single-hue sequential ramp (light->dark), one per metric panel.
Coarse-anchored areas (Hb->Tha): colored by metric but diagonal-hatched + dashed border.
Gap areas (no Modha node): excluded, never drawn colored (they have no box to draw on).
Unmatched boxes (no recorded data): thin light-gray outline only, no fill.
"""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

VIS_DIR = '/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels/VisualizeOnAreas'
COORDS_CSV = f'{VIS_DIR}/modha-singh-network-layout/results/modha_boxed_layout_coordinates.csv'
METRICS_CSV = f'{VIS_DIR}/area_metrics_by_modha_node_L6.csv'
OUT_PNG = f'{VIS_DIR}/figures/area_metrics_layout.png'

BOX_W, BOX_H = 1.0, 0.45

# single-hue sequential ramp: light -> dark blue (not a preset multi-hue colormap)
SEQ_CMAP = LinearSegmentedColormap.from_list('single_hue_blue', ['#dfeaf7', '#0b3d6e'])
GAP_FACE = '#e5e5e5'      # unmatched / no-data box fill
GAP_EDGE = '#b0b0b0'
BG_EDGE = '#cfcfcf'       # boxes with no recorded data at all (context only)

METRICS = [
    ('tuning_PctDirectional', '% Directionally Tuned'),
    ('decode_execution', 'Decoding Accuracy (execution epoch)'),
    ('planexec_PctExecOnly', '% Execution-Recruited (of plan+exec cells)'),
]

def load_coords():
    boxes = {}
    with open(COORDS_CSV) as f:
        for row in csv.DictReader(f):
            if row['is_plotted'] != '1':
                continue
            boxes[row['acronym']] = {
                'x': float(row['box_x']), 'y': float(row['box_y']),
                'hierarchy_depth': int(row['hierarchy_depth']),
            }
    return boxes

def load_metrics():
    rows = []
    with open(METRICS_CSV) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows

def node_list(modha_node):
    return [n.strip() for n in modha_node.split(';') if n.strip()]

def draw_panel(ax, boxes, metrics, metric_key, metric_label):
    values = []
    for row in metrics:
        if not row['modha_node']:
            continue
        v = row[metric_key]
        if v:
            values.append(float(v))
    vmin, vmax = (min(values), max(values)) if values else (0, 1)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    # background: every plotted box, faint outline only (anatomical context)
    for acr, b in boxes.items():
        ax.add_patch(Rectangle((b['x'], b['y']), BOX_W, BOX_H,
                                facecolor='none', edgecolor=BG_EDGE, linewidth=0.3, zorder=1))

    # colored boxes: our mapped areas
    for row in metrics:
        nodes = node_list(row['modha_node'])
        if not nodes:
            continue
        v = row[metric_key]
        if not v:
            continue
        val = float(v)
        color = SEQ_CMAP(norm(val))
        for node in nodes:
            b = boxes.get(node)
            if b is None:
                continue
            rect = Rectangle((b['x'], b['y']), BOX_W, BOX_H,
                              facecolor=color, edgecolor='black', linewidth=0.6, zorder=3)
            if row['match_type'] == 'coarse_anchor':
                # no tracer data at any level (Hb): heavy hatch, dashed edge
                rect.set_hatch('xxxx')
                rect.set_linestyle('--')
                rect.set_edgecolor('#333333')
            elif row['match_type'] == 'coarse_merged':
                # real tracer data, but at coarser granularity (SC->MB#2): light hatch
                rect.set_hatch('///')
                rect.set_edgecolor('#333333')
            ax.add_patch(rect)
            label = node if row['match_type'] != 'coarse_anchor' else f"{node}\n({row['pooled_label']})"
            ax.text(b['x'] + BOX_W / 2, b['y'] + BOX_H / 2, label,
                    ha='center', va='center', fontsize=3.2, zorder=4,
                    color='white' if norm(val) > 0.55 else 'black')

    ax.set_title(metric_label, fontsize=9)
    matched_xs, matched_ys = [], []
    for row in metrics:
        for node in node_list(row['modha_node']):
            b = boxes.get(node)
            if b is not None:
                matched_xs.append(b['x']); matched_ys.append(b['y'])
    pad = 1.5
    ax.set_xlim(min(matched_xs) - pad, max(matched_xs) + BOX_W + pad)
    ax.set_ylim(min(matched_ys) - pad, max(matched_ys) + BOX_H + pad)
    ax.set_aspect('equal')
    ax.axis('off')

    sm = ScalarMappable(norm=norm, cmap=SEQ_CMAP)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.ax.tick_params(labelsize=6)

def main():
    boxes = load_coords()
    metrics = load_metrics()

    fig, axes = plt.subplots(1, len(METRICS), figsize=(27, 6.5))
    for ax, (key, label) in zip(axes, METRICS):
        draw_panel(ax, boxes, metrics, key, label)

    # legend for gap/coarse encoding
    from matplotlib.patches import Patch
    legend_elems = [
        Patch(facecolor='#dfeaf7', edgecolor='black', label='Leaf-node match (real tracer node)'),
        Patch(facecolor='#dfeaf7', edgecolor='#333333', hatch='///',
              label='Coarse-merged (e.g. SC→MB#2) — real tracer data, coarser granularity'),
        Patch(facecolor='#dfeaf7', edgecolor='#333333', hatch='xxxx', linestyle='--',
              label='Coarse-anchored (Hb→Tha) — NO tracer data, anatomical placement only'),
        Patch(facecolor='none', edgecolor=BG_EDGE, label='Unmatched Modha node (no recorded data)'),
    ]
    fig.legend(handles=legend_elems, loc='lower center', bbox_to_anchor=(0.5, 0.02),
               ncol=4, fontsize=8, frameon=False)

    excluded = [r['pooled_label'] for r in metrics if not r['modha_node']]
    fig.suptitle('Saccade-direction metrics on Modha-Singh anatomical box layout (m032/Kwibus, labelField=L6)',
                 fontsize=12, y=0.98)
    fig.text(0.5, -0.02,
              f"Gap areas excluded from layout (no Modha tracer node, see area_to_modha_crosswalk.csv): {', '.join(excluded)}",
              ha='center', fontsize=7, color='#555555')

    fig.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
    print('wrote', OUT_PNG)

if __name__ == '__main__':
    main()
