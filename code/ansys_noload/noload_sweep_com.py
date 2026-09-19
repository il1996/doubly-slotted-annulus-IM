# -*- coding: utf-8 -*-
"""
Numerical no-load test on a COPY of IM_18kW_690V.aedt, driven by COM from
CPython (pywin32).  The design 'a vide' is duplicated once per line voltage,
the winding voltage amplitude is changed, the rotor is driven at constant
synchronous speed (1500 rpm, mechanical transient off); everything else
(Rs = $Rs, external end-winding inductance $Ls, dt = 1 ms, 2 s, nonlinear
residual 1e-4, one slice, post-processed core loss, mesh) is inherited
unchanged from the reference design.  After each solve the phase currents,
input voltages, induced voltages and flux linkages are exported to .tab.

  python noload_sweep_com.py            -> dry run (designs created and saved, no solve)
  python noload_sweep_com.py solve      -> solve and export
"""
import os, sys, time
for v in ("PYTHONPATH", "PYTHONHOME"):
    os.environ.pop(v, None)
import win32com.client

HERE = r"<home>\Desktop\claude\T2_ansys"
PROJ = os.path.join(HERE, "IM_18kW_690V_noload_sweep.aedt")
OUT = os.path.join(HERE, "exports")
LOG = os.path.join(HERE, "sweep_log.txt")
SOLVE = len(sys.argv) > 1 and sys.argv[1] == "solve"
VOLTS = [int(a) for a in sys.argv[2:]] if len(sys.argv) > 2 else [690, 550, 410, 275, 140]
SRC = "a vide"
os.makedirs(OUT, exist_ok=True)


def log(s):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + s + "\n")
    print(s, flush=True)


log("=== noload_sweep_com start (solve=%s, volts=%s) ===" % (SOLVE, VOLTS))
app = win32com.client.DispatchEx("Ansoft.ElectronicsDesktop.2023.1")     # always a fresh server process
oDesktop = app.GetAppDesktop()
log("desktop %s" % oDesktop.GetVersion())
pname = os.path.splitext(os.path.basename(PROJ))[0]
openp = list(oDesktop.GetProjectList())
log("projects already open in this desktop: %s" % openp)
if pname in openp:
    oProject = oDesktop.SetActiveProject(pname)
else:
    oProject = oDesktop.OpenProject(PROJ)
    if oProject is None:
        oProject = oDesktop.SetActiveProject(pname)
log("project %s ; designs %s" % (oProject.GetName(), list(oProject.GetTopDesignList())))

YQ = []
for ph in ["Phase_A", "Phase_B", "Phase_C"]:
    for q in ["Current", "InputVoltage", "InducedVoltage", "FluxLinkage"]:
        YQ.append("%s(%s)" % (q, ph))

for V in VOLTS:
    name = "NL_V%d" % V
    try:
        if name in list(oProject.GetTopDesignList()):
            log("design %s exists, reusing" % name)
            oDesign = oProject.SetActiveDesign(name)
        else:
            before = set(oProject.GetTopDesignList())
            oProject.CopyDesign(SRC)
            oProject.Paste()
            after = [n for n in oProject.GetTopDesignList() if n not in before]
            newname = after[0]
            oDesign = oProject.SetActiveDesign(newname)
            oDesign.RenameDesignInstance(newname, name)
            oDesign = oProject.SetActiveDesign(name)
            log("design %s created from '%s' (pasted as '%s')" % (name, SRC, newname))
            oB = oDesign.GetModule("BoundarySetup")
            for ph, shift in [("Phase_A", ""), ("Phase_C", "+2*pi/3"), ("Phase_B", "-2*pi/3")]:
                expr = "%d*sqrt(2/3)*sin(time*2*$f*pi%s)" % (V, shift)
                oB.EditWindingGroup(ph, ["NAME:" + ph, "Type:=", "Voltage", "IsSolid:=", False,
                                         "Current:=", "0mA", "Resistance:=", "$Rs", "Inductance:=", "$Ls",
                                         "Voltage:=", expr, "ParallelBranchesNum:=", "2"])
                log("  %s voltage := %s" % (ph, expr))
            oM = oDesign.GetModule("ModelSetup")
            oM.EditMotionSetup("MotionSetup1", ["NAME:Data", "Move Type:=", "Rotate", "Coordinate System:=", "Global",
                                                "Axis:=", "Z", "Is Positive:=", True, "InitPos:=", "0deg",
                                                "HasRotateLimit:=", False, "NonCylindrical:=", False,
                                                "Consider Mechanical Transient:=", False,
                                                "Angular Velocity:=", "1500rpm"])
            log("  motion := 1500 rpm constant, mechanical transient off")
            oProject.Save()
            log("  saved")
        if not SOLVE:
            continue
        t0 = time.time()
        log("  Analyze Setup1 of %s ..." % name)
        oDesign.Analyze("Setup1")
        log("  solved in %.0f s" % (time.time() - t0))
        oProject.Save()
        oR = oDesign.GetModule("ReportSetup")
        for rep, ys in [("NL_wave", YQ), ("NL_misc", ["CoreLoss", "Moving1.Torque", "Moving1.Speed", "SolidLoss"])]:
            try:
                oR.DeleteReports([rep])
            except Exception:
                pass
            try:
                oR.CreateReport(rep, "Transient", "Data Table", "Setup1 : Transient",
                                ["Domain:=", "Sweep"], ["Time:=", ["All"]],
                                ["X Component:=", "Time", "Y Component:=", ys], [])
                fn = os.path.join(OUT, "%s_%s.tab" % (name, rep))
                oR.ExportToFile(rep, fn)
                log("  exported " + fn)
            except Exception as e:
                log("  report %s failed: %s" % (rep, str(e)))
        oProject.Save()
    except Exception as e:
        log("ERROR on %s: %s" % (name, str(e)))

oProject.Save()
log("=== noload_sweep_com end ===")
if SOLVE:
    try:
        oProject.Close()
        oDesktop.QuitApplication()
    except Exception as e:
        log("quit: %s" % str(e))
