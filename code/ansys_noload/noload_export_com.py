# -*- coding: utf-8 -*-
"""Export the no-load waveforms of solved designs by COM.  Attaches to a
running AEDT that holds the project if there is one, else opens the project.
  python noload_export_com.py NL_V690 NL_V550 ...   [--close]
--close closes the project at the end (needed before a -batchsolve)."""
import os, sys, time
for v in ("PYTHONPATH", "PYTHONHOME"):
    os.environ.pop(v, None)
import win32com.client
HERE = r"<home>\Desktop\claude\T2_ansys"
PROJ = os.path.join(HERE, "IM_18kW_690V_noload_sweep.aedt")
OUT = os.path.join(HERE, "exports"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(HERE, "sweep_log.txt")
names = [a for a in sys.argv[1:] if not a.startswith('--')]
CLOSE = '--close' in sys.argv
pname = os.path.splitext(os.path.basename(PROJ))[0]


def log(s):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + s + "\n")
    print(s, flush=True)


log("=== export_com start %s ===" % names)
app = win32com.client.Dispatch("Ansoft.ElectronicsDesktop.2023.1")
oDesktop = app.GetAppDesktop()
openp = list(oDesktop.GetProjectList())
log("desktop %s ; open projects %s" % (oDesktop.GetVersion(), openp))
if pname in openp:
    oDesktop.SetActiveProject(pname)
    oProject = oDesktop.GetActiveProject()
else:
    lock = PROJ + ".lock"
    if os.path.exists(lock):
        os.remove(lock); log("removed stale lock")
    oProject = oDesktop.OpenProject(PROJ)
    if oProject is None:
        oProject = oDesktop.GetActiveProject()
log("project %s ; designs %s" % (oProject.GetName(), list(oProject.GetTopDesignList())))
YQ = ["%s(%s)" % (q, ph) for ph in ["Phase_A", "Phase_B", "Phase_C"] for q in ["Current", "InputVoltage", "InducedVoltage", "FluxLinkage"]]
for name in names:
    try:
        oDesign = oProject.SetActiveDesign(name)
        oR = oDesign.GetModule("ReportSetup")
        for rep, ys in [("NL_wave", YQ), ("NL_misc", ["CoreLoss", "Moving1.Torque", "Moving1.Speed", "SolidLoss"])]:
            try:
                oR.DeleteReports([rep])
            except Exception:
                pass
            oR.CreateReport(rep, "Transient", "Data Table", "Setup1 : Transient",
                            ["Domain:=", "Sweep"], ["Time:=", ["All"]],
                            ["X Component:=", "Time", "Y Component:=", ys], [])
            fn = os.path.join(OUT, "%s_%s.tab" % (name, rep))
            oR.ExportToFile(rep, fn)
            log("  %s: exported %s (%d bytes)" % (name, os.path.basename(fn), os.path.getsize(fn) if os.path.exists(fn) else -1))
    except Exception as e:
        log("  %s: EXPORT ERROR %s" % (name, str(e)))
try:
    oProject.Save(); log("project saved")
except Exception as e:
    log("save failed: %s" % str(e))
if CLOSE:
    try:
        oDesktop.CloseProject(pname); log("project closed")
    except Exception as e:
        log("close failed: %s" % str(e))
log("=== export_com end ===")
