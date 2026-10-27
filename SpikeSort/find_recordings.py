"""Phase 0: inventory every OpenEphys recording under Data_collection.

Standard library only. Port of the recursive search idea in
SaccadeDirection/code/matchRunRecording.m (recordings are not always at
<day>/experimentN/recordingM; e.g. 20260224 is nested in a 'Record Node' folder).

Usage:  python find_recordings.py [--root PATH] [--out inventory.csv]
Writes one row per recording with probe/channel info and a 'sortable' flag.
"""
import argparse
import csv
import json
import re
from pathlib import Path

import config

FIELDS = ["monkey", "day", "experiment", "recording", "path", "settings_xml",
          "ap_folder", "ap_num_channels", "ap_fs", "ap_bytes", "duration_s",
          "has_ap_dat", "sortable", "note"]


def find_settings(exp_dir: Path):
    for name in (f"settings_{exp_dir.name.replace('experiment', '')}.xml", "settings.xml"):
        p = exp_dir.parent / name
        if p.is_file():
            return p
    return None


def inspect_recording(rec: Path, monkey: str, day: str):
    exp = rec.parent
    row = dict.fromkeys(FIELDS, "")
    row.update(monkey=monkey, day=day, experiment=exp.name, recording=rec.name, path=str(rec))
    settings = find_settings(exp)
    row["settings_xml"] = str(settings) if settings else ""
    oebin = rec / "structure.oebin"
    if not oebin.is_file():
        row["note"] = "no structure.oebin"
        row["sortable"] = 0
        return row
    try:
        cont = json.loads(oebin.read_text())["continuous"]
    except Exception as e:  # corrupt/partial recording
        row["note"] = f"oebin unreadable: {e}"
        row["sortable"] = 0
        return row
    ap = next((c for c in cont if config.AP_STREAM_KEY in c["folder_name"]), None)
    if ap is None:
        row["note"] = "no AP stream"
        row["sortable"] = 0
        return row
    dat = rec / "continuous" / ap["folder_name"].strip("/\\") / "continuous.dat"
    row.update(ap_folder=ap["folder_name"], ap_num_channels=ap["num_channels"],
               ap_fs=ap["sample_rate"], has_ap_dat=int(dat.is_file()))
    if dat.is_file():
        nbytes = dat.stat().st_size
        row["ap_bytes"] = nbytes
        row["duration_s"] = round(nbytes / 2 / ap["num_channels"] / ap["sample_rate"], 1)
    row["sortable"] = int(dat.is_file() and settings is not None)
    if settings is None:
        row["note"] = "no settings.xml (geometry unknown)"
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(config.DATA_COLLECTION))
    ap.add_argument("--out", default=str(config.INVENTORY))
    args = ap.parse_args()
    root = Path(args.root)

    rows = []
    for monkey in sorted(p for p in root.iterdir() if p.is_dir() and re.fullmatch(r"m\d+", p.name)):
        for day in sorted(p for p in monkey.iterdir() if p.is_dir() and re.fullmatch(r"\d{8}", p.name)):
            for rec in sorted(day.rglob("recording*")):
                if rec.is_dir() and re.fullmatch(r"recording\d+", rec.name) \
                        and re.fullmatch(r"experiment\d+", rec.parent.name):
                    rows.append(inspect_recording(rec, monkey.name, day.name))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    n_ok = sum(int(r["sortable"] or 0) for r in rows)
    print(f"{len(rows)} recordings, {n_ok} sortable -> {out}")


if __name__ == "__main__":
    main()
