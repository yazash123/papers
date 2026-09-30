"""Digitise the figure data used in session 2 (REPORT2.md, Part 0 and Part C).

Every value produced here is READ OFF A FIGURE and is flagged as such in
REPORT2.md.  Two methods are used:

* vector figures (Carter & Cross 2005, Fig. 2): marker centres are taken
  directly from the PDF drawing commands, and the axes are calibrated from
  tick marks / axis lines in the same drawing.  The residual error is the
  calibration (~0.05 pN, ~1 nm/s, ~1% in log axes).
* raster figures (Visscher 1999 Fig. 3a and 4; Block 2003 Fig. 4; Sozanski
  2015 Fig. 2b): the page is rendered at 600 dpi, the axes are located from
  tick marks, and markers are found as connected dark (or coloured) blobs.
  Error bars are read along the marker's column.  The residual error is a few
  pixels (~0.3% of the axis range) plus marker overlap, which is flagged.

Output: data/digitized/*.csv (one file per panel) and a printed summary.
Run:  cd analysis && PYTHONPATH=.. python digitize_figures.py
Needs PyMuPDF (pip install pymupdf).
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import ndimage

try:
    import pymupdf
except ImportError:  # older PyMuPDF
    import fitz as pymupdf

ROOT = Path(__file__).resolve().parents[1]
PAPERS = ROOT / "papers"
OUT = ROOT / "data" / "digitized"
OUT.mkdir(parents=True, exist_ok=True)


def write(name, header, rows, note):
    path = OUT / name
    with path.open("w", newline="") as fh:
        fh.write(f"# {note}\n")
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow([f"{v:.4g}" if isinstance(v, float) else v for v in r])
    return path


def render(fname, page, clip, dpi=600, rgb=False):
    doc = pymupdf.open(str(PAPERS / fname))
    pg = doc[page - 1]
    pix = pg.get_pixmap(dpi=dpi, clip=pymupdf.Rect(*clip),
                        colorspace=pymupdf.csRGB if rgb else pymupdf.csGRAY)
    a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    a = a.astype(float) / 255.0
    return a if rgb else a[:, :, 0]


def runs(idx):
    idx = list(idx)
    out = []
    if not idx:
        return out
    s = p = idx[0]
    for i in idx[1:]:
        if i != p + 1:
            out.append((s, p))
            s = i
        p = i
    out.append((s, p))
    return out


# ---------------------------------------------------------------------------
# Carter & Cross 2005, Nature 435:308, Fig. 2 (p.310) -- vector graphics
# ---------------------------------------------------------------------------
CC = "Carter 2005 - Mechanics of the kinesin step.pdf"
NAVY, RED = (0.2, 0.21, 0.58), (0.91, 0.21, 0.18)


def _cc_markers(y_lo, y_hi, x_hi=300.0):
    doc = pymupdf.open(str(PAPERS / CC))
    out = []
    for d in doc[2].get_drawings():
        r = d["rect"]
        if r.x0 > x_hi or not (y_lo < r.y0 < y_hi):
            continue
        if d["type"] in ("f", "fs") and 1.5 < r.width < 4.0 and 1.5 < r.height < 4.0 and d.get("fill"):
            out.append((tuple(round(c, 2) for c in d["fill"]), (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2, r.width, d["items"][0][0]))
    return out


def carter_cross():
    # panel c: velocity vs load.  Zero-load line x = 175.095 pt, axis -15..15 pN
    # between 91.795 and 258.715 pt; v = 0 at y = 712.785 pt, -200 at 740.745 pt.
    sx = (258.715 - 91.795) / 30.0
    Xc = lambda x: (x - 175.095) / sx
    Yc = lambda y: (712.785 - y) / (27.96 / 200.0)
    rows = []
    for fill, x, y, w, kind in _cc_markers(560, 740):
        atp = {NAVY: "1mM", RED: "10uM"}.get(fill)
        if atp and 2.0 < w < 3.2:
            rows.append((atp, round(Xc(x), 3), round(Yc(y), 2), "trap-off square" if kind == "re" else "bin"))
    rows = sorted(set(rows))
    write("carter2005_fig2c_velocity.csv", ["atp", "load_pN", "velocity_nm_s", "marker"], rows,
          "Carter & Cross 2005 Fig.2c (p.310), vector-extracted. Velocity = mean amplitude / mean dwell per 1-pN bin; "
          "squares = trap-off bead velocity. Load > 0 hindering.")
    # panel b: dwell vs load (log axis 1 ms..10 s between y = 564.08 and 397.16 pt)
    Xb = lambda x: (x - (175.095 - 0.36)) / sx
    Yb = lambda y: 10 ** ((564.08 - y) / ((564.08 - 397.16) / 4.0) - 3)
    marks = _cc_markers(395, 566)
    whites = {(round(x, 1), round(y, 1)) for f, x, y, w, k in marks if f == (1.0, 1.0, 1.0)}
    rowsb = []
    for fill, x, y, w, kind in marks:
        atp = {NAVY: "1mM", RED: "10uM", (0.92, 0.19, 0.18): "10uM"}.get(fill)
        if not atp:
            continue
        step = "backward" if (round(x, 1), round(y, 1)) in whites else "forward"
        rowsb.append((atp, step, round(Xb(x), 3), round(1e3 * Yb(y), 3), "trap-off square" if kind == "re" else "bin"))
    rowsb = sorted(set(rowsb))
    write("carter2005_fig2b_dwell.csv", ["atp", "step", "load_pN", "mean_dwell_ms", "marker"], rowsb,
          "Carter & Cross 2005 Fig.2b (p.310), vector-extracted. Filled = forward-step dwells, open = backward-step dwells "
          "(1-pN bins, s.e.m. bars not extracted).")
    # inset of panel a: forward/backward ratio (log axis, 16.79 pt/decade, ratio 1 at y = 334.98;
    # load 0 at x = 117.47 pt, 10 pN at 173.97 pt)
    Xi = lambda x: (x - 117.47) / 5.65
    Yi = lambda y: 10 ** ((334.98 - y) / 16.79)
    rowsi = []
    for fill, x, y, w, kind in _cc_markers(270, 385, x_hi=300):
        atp = {(0.1, 0.34, 0.65): "1mM", (0.92, 0.19, 0.18): "10uM"}.get(fill)
        if atp and 2.0 < w < 2.6 and x < 180:
            rowsi.append((atp, round(Xi(x), 3), round(Yi(y), 4)))
    rowsi = sorted(set(rowsi))
    write("carter2005_fig2a_ratio.csv", ["atp", "load_pN", "forward_to_backward_ratio"], rowsi,
          "Carter & Cross 2005 Fig.2a inset (p.310), vector-extracted.")
    return rows, rowsb, rowsi


# ---------------------------------------------------------------------------
# Visscher, Schnitzer & Block 1999, Nature 400:184 -- raster figures
# ---------------------------------------------------------------------------
VIS = "Visscher 1999 - Single kinesin molecules studied with a molecular force clamp.pdf"


def _errbar(dark, c0, r0, rmax):
    col = dark[:rmax, int(round(c0))]
    top = int(r0)
    while top > 0 and col[top - 1]:
        top -= 1
    bot = int(r0)
    while bot < rmax - 1 and col[bot + 1]:
        bot += 1
    return bot, top


def visscher_fig4b():
    a = render(VIS, 4, (40, 250, 275, 460))
    dark = a < 0.5
    # calibration from tick marks (see REPORT2.md): x 1 pN at col 584.0 ... 6 pN at 1780.0;
    # y 1.0 at row 203, 0.5 at row 589
    X = lambda c: (c - 344.8) / 239.2
    Y = lambda r: (975.0 - r) / 772.0
    lab, _ = ndimage.label(ndimage.binary_opening(dark, structure=np.ones((15, 15)))[:960, 360:])
    rows = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        m = lab[sl] == i
        cy, cx = ndimage.center_of_mass(m)
        r0, c0 = sl[0].start + cy, sl[1].start + cx + 360
        bot, top = _errbar(dark, c0, r0, 960)
        rows.append((round(X(c0), 3), round(Y(r0), 4), round((Y(top) - Y(bot)) / 2, 4)))
    rows.sort()
    write("visscher1999_fig4b_randomness_vs_load.csv", ["load_pN", "r", "sem_approx"], rows,
          "Visscher 1999 Fig.4b (p.187): randomness vs load at 2 mM ATP (squid kinesin), raster-digitised; "
          "sem = half the error-bar extent (includes cap thickness ~0.006).")
    return rows


def visscher_fig4a():
    a = render(VIS, 4, (40, 45, 290, 205))
    dark = a < 0.5
    Xc = lambda c: 10 ** ((c - 499.0) / 386.8)
    Y = lambda r: (1251.0 - r) / 662.0
    in_legend = lambda c, r: (1500 < c < 1660 and r < 440) or (c > 1660 and r < 460)
    rows = []
    filled = ndimage.binary_opening(dark, structure=np.ones((13, 13)))
    lab, _ = ndimage.label(filled[:1240, 360:])
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        m = lab[sl] == i
        cy, cx = ndimage.center_of_mass(m)
        r0, c0 = sl[0].start + cy, sl[1].start + cx + 360
        if in_legend(c0, r0):
            continue
        bot, top = _errbar(dark, c0, r0, 1245)
        rows.append(("filled", round(Xc(c0), 2), round(Y(r0), 4), round((Y(top) - Y(bot)) / 2, 4)))
    holes = ndimage.binary_fill_holes(dark) & ~dark
    lab, _ = ndimage.label(ndimage.binary_opening(holes, structure=np.ones((7, 7)))[:1240, 360:])
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        m = lab[sl] == i
        h, w = m.shape
        if not (15 < h < 40 and 15 < w < 40):
            continue
        cy, cx = ndimage.center_of_mass(m)
        r0, c0 = sl[0].start + cy, sl[1].start + cx + 360
        if in_legend(c0, r0):
            continue
        bot, top = _errbar(dark, c0, r0 + h / 2 + 4, 1245)
        rows.append(("open_5.69pN", round(Xc(c0), 2), round(Y(r0), 4), round(Y(bot) - Y(r0), 4)))
    rows.sort(key=lambda t: (t[0], t[1], t[2]))
    write("visscher1999_fig4a_randomness_vs_atp_raw.csv", ["marker", "atp_uM", "r", "err_approx"], rows,
          "Visscher 1999 Fig.4a (p.187), raster-digitised. 'filled' merges triangles (1.05 pN) and filled circles "
          "(3.59 pN): assignment by eye in REPORT2.md; overlapping markers flagged there. err for open circles = "
          "distance from centre to the lower cap.")
    return rows


def _colour_blobs(a, target, tol=0.18, min_area=30, max_area=4000, open_px=5):
    d = np.sqrt(((a - np.array(target)) ** 2).sum(axis=2))
    mask = ndimage.binary_opening(d < tol, structure=np.ones((open_px, open_px)))
    lab, _ = ndimage.label(mask)
    out = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        m = lab[sl] == i
        if min_area <= m.sum() <= max_area:
            cy, cx = ndimage.center_of_mass(m)
            out.append((sl[1].start + cx, sl[0].start + cy, int(m.sum()), m.shape))
    return out, d < tol


def _bar(mask, x, y):
    col = mask[:, int(round(x))]
    top = bot = int(round(y))
    while top > 0 and col[top - 1]:
        top -= 1
    while bot < mask.shape[0] - 1 and col[bot + 1]:
        bot += 1
    return top, bot


def visscher_fig3a():
    """Force-velocity at 2 mM (blue circles, right axis) and 5 uM (red triangles, left axis).
    Calibration from ticks: 0.5-pN ticks from col 487.0 (0.5 pN) to 1713.0 (7 pN) -> 188.4 px/pN,
    0 at 392.8; left axis 15.08 px per nm/s (0 at row 1000.4); right axis 103.0 px per 100 nm/s
    (0 at row 1000.0)."""
    a = render(VIS, 3, (30, 55, 280, 200), rgb=True)
    F = lambda c: (c - 392.8) / 188.4
    conv = {"2mM": lambda r: (1000.0 - r) * 100 / 103.0, "5uM": lambda r: (1000.4 - r) / 15.08}
    rows = []
    for atp, tgt in (("2mM", (0.20, 0.185, 0.53)), ("5uM", (0.886, 0.10, 0.12))):
        blobs, mask = _colour_blobs(a, tgt, tol=0.2, min_area=400)
        for x, y, area, shape in blobs:
            if not (393 < x < 1807 and 22 < y < 998) or (1250 < x < 1400 and y < 280):
                continue  # outside the plot, or the legend
            top, bot = _bar(mask, x, y)
            c = conv[atp]
            rows.append((atp, round(F(x), 3), round(c(y), 1), round((c(top) - c(bot)) / 2, 1)))
    rows.sort()
    write("visscher1999_fig3a_force_velocity.csv", ["atp", "load_pN", "velocity_nm_s", "sem_approx"], rows,
          "Visscher 1999 Fig.3a (p.186), raster-digitised (colour). Force clamp, squid kinesin. "
          "sem = half the error-bar extent incl. caps.")
    return rows


# ---------------------------------------------------------------------------
# Block et al. 2003, PNAS 100:2351, Fig. 4 (p.2354) -- raster (JPEG)
# Block's sign convention: negative longitudinal force = hindering (backward) load.
# Stored here with load_pN = hindering load = -(Block's force).
# ---------------------------------------------------------------------------
BLK = "Block 2003 - Probing the kinesin reaction cycle with a 2D optical force clamp.pdf"
BLUE_B, RED_B = (0.18, 0.33, 0.62), (0.83, 0.18, 0.20)


def block_fig4():
    rows = []
    # panel C (randomness): x ticks every 2 pN, -8 at col 324.5 ... 8 at 1595.5 (79.4 px/pN,
    # 0 at 959.5); y: 0.25 per 180 px, 0 at row 1121.
    a = render(BLK, 4, (315, 390, 543, 545), rgb=True)
    Fx = lambda c: -(c - 959.5) / 79.44
    Yr = lambda r: (1121.0 - r) / 720.0
    for atp, tgt in (("1.6mM", BLUE_B), ("4.2uM", RED_B)):
        blobs, mask = _colour_blobs(a, tgt, tol=0.25, min_area=250, open_px=9)
        for x, y, area, shape in blobs:
            if not (280 < x < 1640 and 40 < y < 1115) or (y < 260 and x > 1150) or (y > 880 and x > 1150):
                continue  # plot area only; skip the coloured labels
            top, bot = _bar(mask, x, y)
            rows.append(("randomness", atp, round(Fx(x), 3), round(Yr(y), 4), round((Yr(top) - Yr(bot)) / 2, 4)))
    # panel A (velocity).  Ticks: x every 2 pN, -8 at col 316.5 ... 8 at 1589.5 (79.55 px/pN,
    # 0 at 953.5); left axis (1.6 mM, blue) 133.9 px per 100 nm/s with 0 at row 1109.4;
    # right axis (4.2 uM, red) 102.2 px per 5 nm/s with 0 at row 1109.7.
    a = render(BLK, 4, (315, 52, 543, 215), rgb=True)
    Fa = lambda c: -(c - 953.5) / 79.55
    conv = {"1.6mM": lambda r: (1109.4 - r) / 1.339, "4.2uM": lambda r: (1109.7 - r) / 20.44}
    pl, pr = 1.339, 20.44
    for atp, tgt in (("1.6mM", BLUE_B), ("4.2uM", RED_B)):
        blobs, mask = _colour_blobs(a, tgt, tol=0.25, min_area=250, open_px=9)
        for x, y, area, shape in blobs:
            if not (275 < x < 1630 and 30 < y < 1104):
                continue
            if (atp == "1.6mM" and 350 < x < 780 and 230 < y < 330) or (atp == "4.2uM" and 950 < x < 1330 and 600 < y < 690):
                continue  # coloured labels
            top, bot = _bar(mask, x, y)
            c = conv[atp]
            rows.append(("velocity_nm_s", atp, round(Fa(x), 3), round(c(y), 1), round((c(top) - c(bot)) / 2, 1)))
    rows.sort()
    write("block2003_fig4_velocity_randomness.csv", ["quantity", "atp", "load_pN_hindering", "value", "sem_approx"], rows,
          "Block et al. 2003 Fig.4A,C (p.2354), raster-digitised (colour). load_pN_hindering = -(Block's longitudinal force). "
          f"Velocity axes calibrated from ticks: {pl:.4f} px/(nm/s) left, {pr:.4f} px/(nm/s) right.")
    return rows


# ---------------------------------------------------------------------------
# Sozanski et al. 2015, PRL 115:218102, Fig. 2b (p.218102-2) -- raster
# ---------------------------------------------------------------------------
SOZ = "Sozanski 2015 - Small crowders slow down kinesin-1 stepping by hindering motor domain diffusion.pdf"
SOZ_SERIES = {
    "PEG 6 kg/mol": (0.0, 0.0, 0.0), "PEG 18 kg/mol": (1.0, 0.0, 0.0), "PEG 1000 kg/mol": (0.0, 0.0, 1.0),
    "Dextran 10 kg/mol": (0.5, 0.5, 0.0), "Dextran 500 kg/mol": (0.0, 0.0, 0.5), "TetraEG": (0.0, 0.5, 0.5),
    "Sucrose": (1.0, 0.0, 1.0), "BSA": (0.5, 0.0, 0.0),
}


def sozanski_fig2b():
    """x = eta_eff/eta0: long ticks 1..8 at cols 294.5 + 190.43*(x-1); y = v (um/s): 0 at row
    1215.5, 1004 px per um/s."""
    a = render(SOZ, 2, (66, 215, 283, 387), rgb=True)
    X = lambda c: 1.0 + (c - 294.5) / 190.43
    Y = lambda r: (1215.5 - r) / 1004.0
    rows = []
    for name, tgt in SOZ_SERIES.items():
        blobs, mask = _colour_blobs(a, tgt, tol=0.22, min_area=150, open_px=11)
        for x, y, area, shape in blobs:
            if not (262 < x < 1790 and 0 < y < 1240):
                continue
            if 1060 < x < 1790 and 20 < y < 700:
                continue  # legend box
            top, bot = _bar(mask, x, y)
            rows.append((name, round(X(x), 3), round(Y(y), 4), area))
    rows.sort()
    write("sozanski2015_fig2b_velocity_vs_viscosity.csv", ["crowder", "eta_eff_over_eta0", "velocity_um_s", "blob_area_px"], rows,
          "Sozanski et al. 2015 Fig.2b (p.218102-2), raster-digitised by colour; 297 +/- 2 K, 1 mM ATP. "
          "Markers overlapping in colour/position may be merged or missed; check against the figure.")
    return rows


def main():
    rows, rowsb, rowsi = carter_cross()
    print("C&C Fig.2c:", len(rows), "points; 2b:", len(rowsb), "; inset:", len(rowsi))
    r4b = visscher_fig4b()
    print("Visscher Fig.4b:", r4b)
    r4a = visscher_fig4a()
    print("Visscher Fig.4a:", len(r4a))
    r3a = visscher_fig3a()
    print("Visscher Fig.3a:", r3a)
    b4 = block_fig4()
    print("Block Fig.4:", b4)
    s2 = sozanski_fig2b()
    print("Sozanski Fig.2b:", s2)


if __name__ == "__main__":
    main()
