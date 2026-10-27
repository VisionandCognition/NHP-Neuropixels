"""Central paths and parameters for the SpikeSort pipeline.

Runs on the Linux desktop that mounts the server at /media/NETDISKS.
"""
from pathlib import Path

SERVER_ROOT = Path("/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP")
DATA_COLLECTION = SERVER_ROOT / "Data_collection"           # raw: <monkey>/<YYYYMMDD>/...
OUTPUT_ROOT = SERVER_ROOT / "Data_analysis" / "SpikeSorting"  # spike-sorting results only
MANIFEST = OUTPUT_ROOT / "_manifest.csv"
INVENTORY = OUTPUT_ROOT / "_inventory.csv"

# Local fast scratch for preprocessed binary (Kilosort needs a binary file); deleted after sort.
SCRATCH = Path.home() / "spikesort_scratch"

AP_STREAM_KEY = "ProbeA-AP"
NIDAQ_STREAM_KEY = "NI-DAQmx"
