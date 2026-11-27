# -*- coding: utf-8 -*-
"""
Figure 2B — Divergence-time tree with gene-family expansion/contraction
========================================================================
Reproduces panel 2B: calibrated ultrametric tree (green), geological-era
background bands, node labels "point (interval)" in Ma, and gene-family
expansion (+, red) / contraction (-, blue) counts at internal nodes and
terminal branches.

Inputs  (data/):
  tree_calibrated.nwk    calibrated ultrametric tree, branch lengths in Ma
  node_annotations.tsv   node labels, expansion/contraction, label offsets
  leaf_annotations.tsv   species display names, expansion/contraction
  geological_eras.tsv    background bands (name, span, color, label row)

Outputs (output/):
  Figure2B.pdf / .svg (vector) and Figure2B.png (300 dpi)
"""
import csv
import os

import matplotlib as mpl
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "output")
os.makedirs(OUT, exist_ok=True)

# ---------------- style ----------------
mpl.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",          # keep text editable in SVG
    "pdf.fonttype": 42,              # embed TrueType (editable text in PDF)
    "axes.linewidth": 0.8,
})

COL_BRANCH = "#5b8c5a"
COL_NODE = "#3c7a3c"
COL_TIME = "#8b3a3a"
COL_EXP = "#c0504d"
COL_CON = "#4a7a9b"
COL_LEAF = "#3a3a3a"
COL_ERA = "#6e6e6e"

# ---------------- IO helpers ----------------
def read_tsv(path):
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

# ---------------- Newick parser ----------------
class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self):
        self.name, self.length, self.children, self.parent = None, 0.0, [], None

    @property
    def leaves(self):
        if not self.children:
            return [self.name]
        out = []
        for c in self.children:
            out.extend(c.leaves)
        return out


def parse_newick(s):
    s = s.split("=", 1)[1].strip().rstrip(";") if "=" in s else s.strip().rstrip(";")
    pos = 0

    def parse():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                c = parse()
                c.parent = n
                n.children.append(c)
                if s[pos] == ",":
                    pos += 1
                elif s[pos] == ")":
                    pos += 1
                    break
        name = ""
        while pos < len(s) and s[pos] not in ",():;":
            name += s[pos]
            pos += 1
        n.name = name or None
        if pos < len(s) and s[pos] == ":":
            pos += 1
            num = ""
            while pos < len(s) and (s[pos].isdigit() or s[pos] in ".eE+-"):
                num += s[pos]
                pos += 1
            n.length = float(num)
        return n

    return parse()


def mrca(root, leafset):
    """ Deepest node whose leaf descendant set contains `leafset` (exact match
    of the smallest such clade). """
    target = set(leafset)
    best = [None]

    def walk(n):
        ls = set(n.leaves)
        if target <= ls:
            best[0] = n
            for c in n.children:
                walk(c)

    walk(root)
    return best[0]

# ---------------- load data ----------------
with open(os.path.join(DATA, "tree_calibrated.nwk"), encoding="utf-8") as f:
    root = parse_newick(f.read().strip())

nodes_ann = read_tsv(os.path.join(DATA, "node_annotations.tsv"))
leaves_ann = {r["taxon"]: r for r in read_tsv(os.path.join(DATA, "leaf_annotations.tsv"))}
eras = read_tsv(os.path.join(DATA, "geological_eras.tsv"))

# ---------------- layout ----------------
leaf_y, n_leaves = {}, 0

def assign_y(n):
    global n_leaves
    if not n.children:
        leaf_y[n.name] = n_leaves
        n_leaves += 1
    else:
        for c in n.children:
            assign_y(c)

def depth_of(n):
    """ Cumulative branch length from root (root depth = 0). """
    d, p = 0.0, n
    while p.parent is not None:
        d += p.length
        p = p.parent
    return d

def collect_leaves(n, acc):
    if not n.children:
        acc.append(n)
    for c in n.children:
        collect_leaves(c, acc)

_leaves = []
collect_leaves(root, _leaves)
ROOT_AGE = max(depth_of(n) for n in _leaves)

def age_of(n):
    """ Node age before present = root height - depth from root. """
    return ROOT_AGE - depth_of(n)

assign_y(root)

def y_of(n):
    if not n.children:
        return leaf_y[n.name]
    return sum(y_of(c) for c in n.children) / len(n.children)

# ---------------- figure ----------------
fig, ax = plt.subplots(figsize=(7.6, 4.6))
fig.subplots_adjust(left=0.02, right=0.99, top=0.99, bottom=0.11)

ax.set_xlim(203.0, -30.0)          # age decreases to the right
ax.set_ylim(9.7, -2.7)             # row 0 on top
BAND_TOP, BAND_BOT = -0.12, 9.45

# --- era bands ---
for e in eras:
    lo, hi = float(e["end_ma"]), float(e["start_ma"])   # end (young) left, start (old) right
    ax.fill_between([lo, hi], BAND_TOP, BAND_BOT, color=e["color"], zorder=0, lw=0)
    if e["label_row"]:
        y = -1.45 if e["label_row"] == "1" else -0.75
        ha = e.get("label_ha") or "center"
        x = (lo + hi) / 2 if ha == "center" else lo + 0.4
        ax.text(x, y, e["name"], ha=ha, va="center",
                fontsize=7.5, fontweight="bold", color=COL_ERA, zorder=1)

# --- gridlines ---
for t in range(0, 201, 20):
    ax.plot([t, t], [BAND_TOP, BAND_BOT], color="#b8b8b8", lw=0.7,
            ls=(0, (3, 3)), zorder=0.5)

# --- tree ---
def draw_tree(n):
    if not n.children:
        return
    a = age_of(n)
    ys = [y_of(c) for c in n.children]
    ax.plot([a, a], [min(ys), max(ys)], color=COL_BRANCH, lw=1.6,
            solid_capstyle="round", zorder=2)
    for c in n.children:
        ca = age_of(c)
        cy = y_of(c)
        ax.plot([a, ca], [cy, cy], color=COL_BRANCH, lw=1.6,
                solid_capstyle="round", zorder=2)
        draw_tree(c)

draw_tree(root)

# --- leaves: names + expansion/contraction ---
for taxon, row in leaves_ann.items():
    y = leaf_y[taxon]
    ax.text(-2.2, y, row["display_name"], ha="left", va="center",
            fontsize=8.5, color=COL_LEAF, zorder=3)
    if row["expansion"] not in ("", "NA"):
        lx = float(row["label_age_ma"])
        ax.text(lx, y - 0.38, f"+{row['expansion']}", ha="right", va="center",
                fontsize=7.5, color=COL_EXP, zorder=3)
        ax.text(lx, y - 0.38, f"/-{row['contraction']}", ha="left", va="center",
                fontsize=7.5, color=COL_CON, zorder=3)

# --- internal nodes: dots, time labels, expansion/contraction ---
for ann in nodes_ann:
    leafset = ann["leaf_set"].split(";")
    n = mrca(root, leafset)
    a = age_of(n)
    y = y_of(n)
    ax.plot([a], [y], "o", ms=4.2, color=COL_NODE, zorder=4)
    # time label: "point (interval)"
    ax.text(a + float(ann["dx_time"] or 0), y + float(ann["dy_time"] or 0),
            ann["time_label"], ha="center", va="center",
            fontsize=7.5, fontweight="bold", color=COL_TIME, zorder=4)
    # expansion / contraction
    if ann["expansion"] not in ("", "NA"):
        ex = a + float(ann["dx_ec"] or 0)
        ey = y + float(ann["dy_ec"] or 0)
        ax.text(ex, ey, f"+{ann['expansion']}", ha="right", va="center",
                fontsize=7.5, color=COL_EXP, zorder=4)
        ax.text(ex, ey, f"/-{ann['contraction']}", ha="left", va="center",
                fontsize=7.5, color=COL_CON, zorder=4)

# --- axes ---
ax.set_xticks(range(0, 201, 20))
ax.tick_params(axis="x", labelsize=8, length=3, width=0.8)
ax.set_yticks([])
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.set_xlabel("Time (Ma)", fontsize=9.5)

# --- panel label ---
ax.text(201, -2.35, "(B)", fontsize=12, fontweight="bold", color="#555555",
        ha="left", va="center")

# ---------------- save ----------------
for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(OUT, f"Figure2B.{ext}"), dpi=300)
print("Saved:", ", ".join(f"Figure2B.{e}" for e in ("pdf", "svg", "png")),
      "->", OUT)
