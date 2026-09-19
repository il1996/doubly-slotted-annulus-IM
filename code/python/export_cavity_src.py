# -*- coding: utf-8 -*-
"""Export the rotor-cavity bar-current source vector for the MATLAB chain (T5).
s = Phi_J - Q_r [ramp; 0] per ampere of bar current, with Q_r taken from the
SAME cavity_nO16.mat the network uses (so that the homogeneous and source parts
are consistent), Phi_J from cavity_graded.cavity_source_rotor."""
import sys
sys.dont_write_bytecode = True
import numpy as np, scipy.io as sio
from dsop import Machine
from cavity_graded import cavity_source_rotor
M = Machine(); L = M.L; b0 = 2.0
src = r"<home>\AppData\Local\Temp\claude\C--Users-hp-Desktop-claude\fb8f7acf-d703-40c5-b60c-b7b9dd7398fa\scratchpad\repo\doubly-slotted-annulus-IM\code\MEC_IM"
for nO in [16]:
    cav = sio.loadmat(f"{src}\\cavity_nO{nO}.mat")
    Qr = cav['Qr']; Qs = cav['Qs']
    xc = -b0 / 2 + (np.arange(nO) + 0.5) * b0 / nO
    PhiJ, info = cavity_source_rotor(xc, L, lc_mouth=0.01)
    ramp = np.zeros(nO + 1); ramp[:nO] = -(xc + b0 / 2) / b0
    s = PhiJ - Qr[:nO + 1, :nO + 1] @ ramp
    print(f"nO={nO}: Abar {info['Abar_mm2']:.3f} mm2 | sum f_J {info['tot_src']:.6f} | sum f_b {info['tot_bnd']:.6f}")
    print("  Phi_J (Wb/A):", np.array2string(PhiJ, precision=3))
    print("  s     (Wb/A):", np.array2string(s, precision=3), " sum %.2e" % s.sum())
    print("  |Phi_J|max / |Q ramp|max = %.2e" % (np.abs(PhiJ).max() / np.abs(Qr[:nO+1,:nO+1] @ ramp).max()))
    sio.savemat(f"{src}\\cavity_src_nO{nO}.mat", {'sr': s.reshape(-1, 1), 'PhiJ': PhiJ.reshape(-1, 1), 'ramp': ramp.reshape(-1, 1),
                                                   'xc_mm': xc.reshape(-1, 1), 'nO': nO, 'Abar_mm2': info['Abar_mm2'],
                                                   'note': 'flux (Wb) per ampere of bar current (+z), dofs = nO mouth columns then the wall; U_right - U_left = -I_z'})
    print("  saved cavity_src_nO%d.mat" % nO)
