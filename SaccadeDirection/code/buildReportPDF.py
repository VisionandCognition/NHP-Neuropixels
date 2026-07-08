#!/usr/bin/env python3
"""
Builds the full collaborator-facing PDF from PIPELINE_REPORT.md + all
per-session figures in <saccadeDataDir()>/figures/, for all 12 sessions.

Pipeline: markdown (report + generated appendix) -> self-contained HTML
(pandoc, images embedded as base64) -> PDF (headless Chrome print-to-pdf).
No LaTeX engine needed. Requires: pandoc, google-chrome (or
google-chrome-stable) on PATH.

Usage: python3 buildReportPDF.py [output.pdf]
Default output: <saccadeDataDir()>/../PIPELINE_REPORT_<date>.pdf
"""
import subprocess
import sys
import datetime
import shutil
import tempfile
from pathlib import Path

REPO_ROOT = Path('/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels/SaccadeDirection')
REPORT_MD = REPO_ROOT / 'PIPELINE_REPORT.md'
FIGDIR = Path('/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_analysis/SaccadeDirection/data/figures')
DATADIR = FIGDIR.parent
OUTDIR = DATADIR.parent  # .../SaccadeDirection/ on the server, sibling to data/ and info/

DAYS = ['20260224','20260226','20260303','20260310','20260312','20260317',
        '20260318','20260320','20260323','20260324','20260325','20260326']
RUN_OVERRIDE = {'20260224': '001'}  # all others are run-003

FIG_ORDER = [
    ('overview', 'Overview: full-probe MUA heatmap, real-depth ruler, per-epoch significance map, anatomy annotation'),
    ('eye_summary', 'Eye summary: trajectories, fixation, main sequence, latency, pupil, accuracy by direction'),
    ('tuning_examples', 'All responsive channels: tuning curves, category-colored (red=directional, orange=complex, gray=untuned)'),
    ('tuning_matrix', 'Tuning matrix: all responsive channels sorted by preferred direction, with strength and category strips'),
    ('population_summary', 'Population summary (this session only): preferred direction, tuning strength, vs. real depth'),
    ('artifact_check', 'Shared-artifact check: trial-residual correlation matrix and vs. distance'),
    ('shared_kinematics', 'Shared component vs. kinematics: PC1 loadings and correlation with peak velocity/amplitude'),
    ('planexec_gocue', 'Planning vs execution (go-cue-anchored): effect map, tuning-strength and preferred-direction scatter'),
    ('planexec_gocue_curves', 'Planning vs execution (go-cue-anchored): per-channel tuning curves, planning (top) vs execution (bottom)'),
    ('planexec_saclocked', 'Planning vs execution (saccade-onset-locked): effect map, tuning-strength and preferred-direction scatter'),
    ('planexec_saclocked_curves', 'Planning vs execution (saccade-onset-locked): per-channel tuning curves, planning (top) vs execution (bottom)'),
]

CSS = """
body { font-family: -apple-system, Helvetica, Arial, sans-serif; max-width: 900px; margin: 2em auto; padding: 0 1em; line-height: 1.45; color: #111; }
h1 { font-size: 1.8em; border-bottom: 2px solid #333; padding-bottom: 0.2em; }
h2 { font-size: 1.4em; margin-top: 2em; border-bottom: 1px solid #999; padding-bottom: 0.15em; page-break-before: always; }
h3 { font-size: 1.15em; margin-top: 1.5em; }
table { border-collapse: collapse; margin: 1em 0; font-size: 0.9em; }
th, td { border: 1px solid #999; padding: 4px 8px; text-align: left; }
th { background: #eee; }
code { background: #f2f2f2; padding: 1px 4px; border-radius: 3px; font-size: 0.9em; }
pre { background: #f2f2f2; padding: 0.7em; border-radius: 4px; overflow-x: auto; }
img { max-width: 100%; display: block; margin: 0.5em 0 1em 0; border: 1px solid #ccc; }
strong { color: #000; }
hr { border: none; border-top: 1px solid #ccc; margin: 2em 0; }
"""


def build_appendix() -> str:
    lines = ['## 9. Appendix: Per-Session Figures\n',
             'All figures for all 12 sessions, in full. See sec 5.2 for a description of each figure type.\n']
    for day in DAYS:
        run = RUN_OVERRIDE.get(day, '003')
        lines.append(f'### {day} (run-{run})\n')
        for suffix, caption in FIG_ORDER:
            fn = FIGDIR / f'{day}_run-{run}_{suffix}.png'
            if not fn.exists():
                lines.append(f'*(missing: {fn.name})*\n')
                continue
            lines.append(f'**{caption}**\n')
            lines.append(f'![{day} run-{run} {suffix}]({fn})\n')
        lines.append('\\newpage\n')
    return '\n'.join(lines)


def main():
    for tool in ('pandoc',):
        if shutil.which(tool) is None:
            sys.exit(f'{tool} not found on PATH')
    chrome = shutil.which('google-chrome') or shutil.which('google-chrome-stable')
    if chrome is None:
        sys.exit('google-chrome not found on PATH')

    date_str = datetime.date.today().isoformat()
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else OUTDIR / f'PIPELINE_REPORT_{date_str}.pdf'

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        combined_md = tmp / 'combined.md'
        with open(combined_md, 'w') as f:
            f.write(REPORT_MD.read_text())
            f.write('\n\n')
            f.write(build_appendix())

        css_path = tmp / 'style.css'
        css_path.write_text(CSS)

        html_path = tmp / 'combined.html'
        subprocess.run(['pandoc', str(combined_md), '-o', str(html_path),
                         '--standalone', '--embed-resources',
                         '--metadata', 'title=Saccade Direction Pipeline Report — m032',
                         '-c', str(css_path)], check=True)

        subprocess.run([chrome, '--headless', '--disable-gpu', '--no-sandbox',
                         f'--print-to-pdf={out_path}', '--print-to-pdf-no-header',
                         '--no-pdf-header-footer', f'file://{html_path}'], check=True)

    print(f'Wrote {out_path} ({out_path.stat().st_size/1e6:.1f} MB)')


if __name__ == '__main__':
    main()
