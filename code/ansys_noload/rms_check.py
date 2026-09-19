# -*- coding: utf-8 -*-
"""Fundamental rms against true rms of the reference currents (no load and rated load), window [1, 2) s."""
import numpy as np, re, os
W = 2 * np.pi * 50.0
def read_tab(f):
    with open(f, encoding='utf-8', errors='replace') as fh:
        names = re.findall(r'"([^"]+)"', fh.readline())
    return names, np.loadtxt(f, skiprows=1)
for lab, d in [('no load', r"<home>\Desktop\ANSYS résultat 18.5KW\transitoire\a vide"), ('rated load', r"<home>\Desktop\ANSYS résultat 18.5KW\transitoire\en charge")]:
    names, data = read_tab(os.path.join(d, "Winding Plot 4.tab")); t = data[:, 0]; sel = (t >= 1.0 - 1e-9) & (t < 2.0 - 1e-9)
    out = []
    for ph in 'ABC':
        j = [k for k, n in enumerate(names) if n.startswith("Current(Phase_%s)" % ph)][0]; y = data[:, j]
        c = (2.0 / sel.sum()) * np.sum(y[sel] * np.exp(-1j * W * t[sel])); f1 = abs(c) / np.sqrt(2); rms = np.sqrt(np.mean(y[sel] ** 2))
        out.append((f1, rms))
    f1 = np.mean([o[0] for o in out]); rms = np.mean([o[1] for o in out])
    print("%-10s : fundamental rms %.4f A | true rms %.4f A | ratio %.4f | THD %.1f %%" % (lab, f1, rms, f1 / rms, 100 * np.sqrt(max(rms ** 2 - f1 ** 2, 0)) / f1))
