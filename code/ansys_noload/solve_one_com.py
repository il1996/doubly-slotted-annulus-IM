# -*- coding: utf-8 -*-
"""Solve ONE design of the no-load sweep by COM, in its own process, and
return.  The results are written to the .aedtresults folder by AEDT as the
solve proceeds; the export is done afterwards by noload_export_com.py in a
fresh session, so a crash of the desktop after the solve loses nothing.
  python solve_one_com.py NL_V550"""
import os, sys, time
for v in ("PYTHONPATH", "PYTHONHOME"):
    os.environ.pop(v, None)
import win32com.client
HERE = r"C:\Users\hp\Desktop\claude\T2_ansys"
PROJ = os.path.join(HERE, "IM_18kW_690V_noload_sweep.aedt")
LOG = os.path.join(HERE, "sweep_log.txt")
name = sys.argv[1]
pname = os.path.splitext(os.path.basename(PROJ))[0]


def log(s):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + s + "\n")
    print(s, flush=True)


log("=== solve_one_com %s ===" % name)
app = win32com.client.DispatchEx("Ansoft.ElectronicsDesktop.2023.1")
oDesktop = app.GetAppDesktop()
openp = list(oDesktop.GetProjectList())
log("desktop %s ; open projects %s" % (oDesktop.GetVersion(), openp))
if pname in openp:
    oDesktop.SetActiveProject(pname); oProject = oDesktop.GetActiveProject()
else:
    lock = PROJ + ".lock"
    if os.path.exists(lock):
        os.remove(lock); log("removed stale lock")
    oProject = oDesktop.OpenProject(PROJ)
    if oProject is None:
        oProject = oDesktop.GetActiveProject()
log("project %s" % oProject.GetName())
oDesign = oProject.SetActiveDesign(name)
t0 = time.time()
log("  Analyze Setup1 of %s ..." % name)
try:
    oDesign.Analyze("Setup1")
    log("  solved in %.0f s" % (time.time() - t0))
except Exception as e:
    log("  Analyze raised after %.0f s: %s" % (time.time() - t0, str(e)))
try:
    oProject.Save(); log("  saved")
except Exception as e:
    log("  save failed: %s" % str(e))
try:
    oDesktop.CloseProject(pname); log("  closed")
except Exception as e:
    log("  close failed: %s" % str(e))
log("=== solve_one_com %s end ===" % name)
