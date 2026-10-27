"""Phase 0: inventory every OpenEphys recording under Data_collection.

Standard library only. Port of the recursive search idea in
SaccadeDirection/code/matchRunRecording.m (recordings are not always at
<day>/experimentN/recordingM; e.g. 20260224 is nested in a 'Record Node' folder).

A sort unit is (monkey, day, node, experiment). `node` is the folder that directly
contains experimentN/ and settings.xml, relative to the day folder ('.' for the flat
layout). This matters: on 20260224 both <day>/experiment1 and
<day>/<...>/Record Node 110/experiment1 exist and are different recordings.

Usage:  python find_recordings.py [--root PATH] [--out inventory.csv]
"""
import argparse
import csv
import json
import re
from pathlib import Path

import config

FIELDS = ["monkey", "day", "node", "experiment", "recording", "path", "settings_xml",
          "ap_folder", "ap_num_channels", "ap_fs", "ap_bytes", "duration_s",
          "has_ap_dat", "sortable", "note"]


def find_settings(exp_dir):
    n = exp_dir.name.replace("experiment", "")
    for name in ("settings_%s.xml" % n, "settings.xml"):
        p = exp_dir.parent / name
        if p.is_file():
            return p
    return None


def inspect_recording(rec, monkey, day_dir):
    exp = rec.parent
    node = exp.parent.relative_to(day_dir)
    row = dict.fromkeys(FIELDS, "")
    row.update(monkey=monkey, day=day_dir.name, node=str(node), experiment=exp.name,
               recording=rec.name, path=str(rec))
    settings = find_settings(exp)
    row["settings_xml"] = str(settings) if settings else ""
    row["sortable"] = 0
    oebin = rec / "structure.oebin"
    if not oebin.is_file():
        row["note"] = "no structure.oebin"
        return row
    try:
        cont = json.loads(oebin.read_text())["continuous"]
    except Exception as e:  # corrupt/partial recording
        row["note"] = "oebin unreadable: %s" % e
        return row
    ap = next((c for c in cont if config.AP_STREAM_KEY in c["folder_name"]), None)
    if ap is None:
        row["note"] = "no AP stream"
        return row
    dat = rec / "continuous" / ap["folder_name"].strip("/\\") / "continuous.dat"
    row.update(ap_folder=ap["folder_name"].strip("/\\"), ap_num_channels=ap["num_channels"],
               ap_fs=ap["sample_rate"], has_ap_dat=int(dat.is_file()))
    if dat.is_file():
        nbytes = dat.stat().st_size
        row["ap_bytes"] = nbytes
        row["duration_s"] = round(nbytes / 2 / ap["num_channels"] / ap["sample_rate"], 1)
    row["sortable"] = int(dat.is_file() and settings is not None)
    if settings is None:
        row["note"] = "no settings.xml (geometry unknown)"
    return row


def scan_day(day_dir, monkey):
    """All recordings under one day folder (recursive)."""
    rows = []
    for rec in sorted(day_dir.rglob("recording*")):
        if rec.is_dir() and re.fullmatch(r"recording\d+", rec.name) \
                and re.fullmatch(r"experiment\d+", rec.parent.name):
            rows.append(inspect_recording(rec, monkey, day_dir))
    return rows


def scan_all(root):
    rows = []
    for monkey in sorted(p for p in root.iterdir() if p.is_dir() and re.fullmatch(r"m\d+", p.name)):
        for day in sorted(p for p in monkey.iterdir() if p.is_dir() and re.fullmatch(r"\d{8}", p.name)):
            rows.extend(scan_day(day, monkey.name))
    return rows


def rec_number(row):
    return int(re.sub(r"\D", "", row["recording"]))


def exp_number(row):
    return int(re.sub(r"\D", "", row["experiment"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(config.DATA_COLLECTION))
    ap.add_argument("--out", default=str(config.INVENTORY))
    args = ap.parse_args()

    rows = scan_all(Path(args.root))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(str(out), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    n_ok = sum(int(r["sortable"] or 0) for r in rows)
    units = {(r["monkey"], r["day"], r["node"], r["experiment"]) for r in rows if r["sortable"]}
    hours = sum(float(r["duration_s"] or 0) for r in rows if r["sortable"]) / 3600
    print("%d recordings (%d sortable, %.1f h) in %d sort units -> %s"
          % (len(rows), n_ok, hours, len(units), out))
    bad = [r for r in rows if not r["sortable"]]
    for r in bad:
        print("  NOT sortable: %s/%s/%s/%s: %s" % (r["monkey"], r["day"], r["experiment"], r["recording"], r["note"]))


if __name__ == "__main__":
    main()
