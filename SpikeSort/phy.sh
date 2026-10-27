#!/bin/bash
# Open a Phy folder written by sort_session.py --phy.
#   ./phy.sh                      # m032 20251010 exp2 (default)
#   ./phy.sh <path/to/phy>        # any phy folder
# Needs a desktop display (DWService session); over plain ssh set DISPLAY=:0 first.
# Flags: QtWebEngine's GPU thread segfaults on this NVIDIA driver (kernel log: Chrome_InProcGp in
# libQt5WebEngineCore), so run it software-rendered.
D="${1:-/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_analysis/SpikeSorting/m032/20251010/exp2/phy}"
export DISPLAY="${DISPLAY:-:0}"
export QTWEBENGINE_CHROMIUM_FLAGS="--disable-gpu --no-sandbox"
cd "$D" && exec /home/chris/anaconda3/envs/phy2/bin/phy template-gui params.py
