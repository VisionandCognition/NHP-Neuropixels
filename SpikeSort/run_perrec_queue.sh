#!/usr/bin/env bash
# Per-recording sorting (sort_perrec.py) for sessions whose concatenated sort has low presence ratios.
# One session after the other; `touch SpikeSort/STOP` ends the queue after the current session.
cd "$(dirname "$0")"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate kilosort
for s in "20260915 1" "20260916 1" "20260918 1" "20260320 1" "20260903 1"; do
  [ -e STOP ] && { echo "STOP file found, queue ended"; break; }
  set -- $s
  echo "$(date '+%F %T') START $1 exp$2"
  python sort_perrec.py m032 "$1" "$2" >> "logs/perrec_$1.log" 2>&1
  echo "$(date '+%F %T') END $1 exp$2 (exit $?)"
done
echo "$(date '+%F %T') QUEUE DONE"
