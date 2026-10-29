"""Stage the raw AP data of one experiment from the (slow) server to local disk.

The server link is the bottleneck (see PLAN.md), so each session is copied once, sequentially,
to config.STAGE_ROOT and sorted from there. Only what sorting needs is copied: the AP
continuous.dat plus all small files (oebin, sync messages, events, timestamps, settings*.xml).
The reader opens every stream's continuous.dat, so the LFP / NI-DAQ ones are sparse
placeholders of the right size (no disk used, never read).

    python stage_raw.py m032 20251010 1 [--node NODE] [--remove]

Resumable: rsync --partial --inplace; a recording counts as staged once `.staged_ok` exists.
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import config
from find_recordings import exp_number

log = logging.getLogger("stage_raw")


def _node_slug(node):
    return "" if node in (".", "") else node.replace("\\", "/").replace("/", "__").replace(" ", "_")


def staged_root(row):
    """Local mirror of the Record-Node folder (holds experimentN/ and settings*.xml)."""
    # One folder per experiment: the reader scans every experimentN/ under its root, so a sibling
    # experiment (or Record-Node folder) that is still being copied ahead (prefetch) would break it.
    slug = _node_slug(row["node"])
    name = "exp%d" % exp_number(row) + ("__" + slug if slug else "")
    return config.STAGE_ROOT / row["monkey"] / row["day"] / name


def _rsync(src, dst, *extra):
    cmd = ["rsync", "-a", "--partial", "--inplace", *extra, str(src), str(dst)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("rsync failed (%d): %s\n%s" % (r.returncode, " ".join(cmd), r.stderr[-500:]))


def _stage_one(row):
    src = Path(row["path"])
    src_root = Path(row["settings_xml"]).parent
    dst_root = staged_root(row)
    dst = dst_root / src.relative_to(src_root)
    ok = dst / ".staged_ok"
    ap_name = row["ap_folder"]
    dst_ap = dst / "continuous" / ap_name / "continuous.dat"
    expected = int(row["ap_bytes"])
    if ok.exists() and dst_ap.is_file() and dst_ap.stat().st_size == expected:
        return dst
    ok.unlink(missing_ok=True)
    dst.mkdir(parents=True, exist_ok=True)
    # everything small first (also creates the continuous/<stream>/ folders)
    _rsync(str(src) + "/", str(dst) + "/", "--exclude=continuous.dat")
    # the AP data itself
    _rsync(src / "continuous" / ap_name / "continuous.dat", dst_ap)
    if dst_ap.stat().st_size != expected:
        raise RuntimeError("%s: staged %d bytes, expected %d" % (dst_ap, dst_ap.stat().st_size, expected))
    # sparse placeholders for the other streams
    for s in (src / "continuous").iterdir():
        if s.name == ap_name or not (s / "continuous.dat").is_file():
            continue
        ph = dst / "continuous" / s.name / "continuous.dat"
        size = (s / "continuous.dat").stat().st_size
        if not ph.is_file() or ph.stat().st_size != size:
            with open(str(ph), "wb") as f:
                f.truncate(size)
    ok.write_text("ok")
    return dst


def check_link():
    """Warn loudly if the NIC negotiated below gigabit (staging would take ~8x longer)."""
    try:
        speed = int(Path("/sys/class/net/%s/speed" % config.NET_IFACE).read_text())
    except Exception:
        return None
    if speed < 1000:
        log.warning("network link %s is only %d Mb/s (expected 1000): copies will be ~8x slower. "
                    "Try `sudo ethtool -r %s`", config.NET_IFACE, speed, config.NET_IFACE)
    return speed


def stage_rows(rows, streams=None):
    """Copy every recording of the experiment; return rows pointing at the local copy."""
    streams = streams or config.STAGE_STREAMS
    check_link()
    first = rows[0]
    src_root = Path(first["settings_xml"]).parent
    dst_root = staged_root(first)
    dst_root.mkdir(parents=True, exist_ok=True)
    for x in src_root.glob("settings*.xml"):
        shutil.copy2(str(x), str(dst_root / x.name))
    total = sum(int(r["ap_bytes"]) for r in rows) / 1e9
    log.info("staging %d recordings (%.1f GB) -> %s with %d streams", len(rows), total, dst_root, streams)
    with ThreadPoolExecutor(max_workers=streams) as ex:
        dests = list(ex.map(_stage_one, rows))
    out = []
    for r, d in zip(rows, dests):
        r = dict(r)
        r["path"] = str(d)
        r["settings_xml"] = str(dst_root / Path(r["settings_xml"]).name)
        out.append(r)
    return out


def remove_staged(staged_rows):
    """Delete the staged recordings of this experiment (never touches anything outside STAGE_ROOT)."""
    root = config.STAGE_ROOT.resolve()
    for r in staged_rows:
        p = Path(r["path"]).resolve()
        if root in p.parents and p.name.startswith("recording"):
            shutil.rmtree(str(p), ignore_errors=True)
            exp = p.parent
            if exp.name.startswith("experiment") and not any(exp.iterdir()):
                exp.rmdir()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("monkey")
    ap.add_argument("day")
    ap.add_argument("experiment", type=int)
    ap.add_argument("--node", default=".")
    ap.add_argument("--remove", action="store_true", help="delete the staged copy instead")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    import sort_session as ss
    rows = ss.select_recordings(args.monkey, args.day, args.experiment, args.node)
    if args.remove:
        remove_staged([dict(r, path=str(staged_root(r) / Path(r["path"]).relative_to(Path(r["settings_xml"]).parent)))
                       for r in rows])
    else:
        stage_rows(rows)


if __name__ == "__main__":
    main()
