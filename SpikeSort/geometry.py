"""Classify a session's electrode configuration and sanity-check the sorting result.

Electrode selection varies a lot across sessions (coverage was prioritised over density on the
long NHP probe): Bank A (dense), Single Column, 4-bank-even, and hand-picked custom mixes.
`describe` turns the channel positions + settings.xml into a profile name and a few numbers;
config.KS4_PROFILES maps profile -> Kilosort parameter overrides (empty until piloted).
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

import config


def _preset_and_banks(settings_xml):
    try:
        np_ = ET.parse(str(settings_xml)).getroot().find(".//NP_PROBE")
        ch = np_.find("CHANNELS")
        banks = Counter(ch.attrib.values()) if ch is not None else Counter()
        return np_.get("electrodeConfigurationPreset") or "", dict(banks)
    except Exception:
        return "", {}


def describe(locations_um, settings_xml):
    """locations_um: (n_channels, 2) array of (x, y) in um (channels that actually get sorted)."""
    loc = np.asarray(locations_um, dtype=float)
    y = np.sort(np.unique(loc[:, 1]))
    gaps = np.diff(y)
    contig = float(np.mean(gaps <= 25.0)) if len(gaps) else 1.0
    span = float(np.ptp(loc[:, 1]))
    preset, banks = _preset_and_banks(settings_xml)
    if preset == "Bank A":
        profile = "dense_bankA"
    elif preset == "Single Column":
        profile = "single_column"
    elif sorted(banks.values()) == [96, 96, 96, 96] and sorted(banks) == ["0", "1", "2", "3"]:
        profile = "four_bank_even"
    else:
        profile = "custom_contiguous" if contig >= 0.9 else "custom_sparse"
    return dict(profile=profile, preset=preset, banks=banks, n_channels=int(len(loc)),
                n_x=int(len(np.unique(loc[:, 0]))), y_span_um=span,
                median_y_gap_um=float(np.median(gaps)) if len(gaps) else 0.0,
                contiguous_frac=contig)


def sanity_flags(n_units, n_good, n_channels, n_spikes, duration_s):
    """Cheap flags for sessions worth a look. Thresholds are provisional (config.SANITY)."""
    s = config.SANITY
    flags = []
    if n_units < s["min_units"]:
        flags.append("few_units")
    if n_channels and n_units / n_channels > s["max_units_per_channel"]:
        flags.append("many_units")
    if n_units and n_good / n_units < s["min_good_frac"]:
        flags.append("low_good_frac")
    if n_channels and duration_s and n_spikes / duration_s / n_channels < s["min_spikes_per_s_per_ch"]:
        flags.append("low_spike_rate")
    return flags


def patch_labels(locations_um, gap_um=100.0):
    """Label channels by contiguous patch along the shank: a new patch starts where the sorted
    y positions jump by more than gap_um. Returns an int label per channel (input order)."""
    y = np.asarray(locations_um, dtype=float)[:, 1]
    order = np.argsort(y, kind="stable")
    new = np.r_[0, (np.diff(y[order]) > gap_um).astype(int)]
    labels = np.empty(len(y), dtype=int)
    labels[order] = np.cumsum(new)
    return labels
