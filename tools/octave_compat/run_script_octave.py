#!/usr/bin/env python3
"""Run a MATLAB RUN_*.m script of the MEC_IM chain under Octave.

Octave needs local functions to be defined before use (or in their own
files), so the trailing local functions of the script are split into
separate files in a temporary directory that is put on the path.  The main
body is written as <name>_main.m and executed with the MEC_IM folder, the
compatibility shims and the function directory on the path.
"""
import re, sys, os, subprocess, pathlib, shutil

HERE = pathlib.Path(__file__).resolve().parent
MEC = HERE / "MEC_IM"
COMPAT = HERE / "compat"


def octave_fixes(src):
    """Syntax MATLAB accepts but Octave does not: indexing the result of a
    function call (f(x).field) -> getfield(f(x),'field')."""
    pat = re.compile(r'(mec\.\w+\((?:[^()]|\([^()]*\))*\))\.(\w+)')
    return pat.sub(lambda m: f"getfield({m.group(1)},'{m.group(2)}')", src)


def split_script(src):
    src = octave_fixes(src)
    lines = src.splitlines()
    # first column-0 'function' line after which everything is local functions
    idx = None
    for i, l in enumerate(lines):
        if re.match(r'^function\b', l):
            idx = i; break
    if idx is None:
        return src, {}
    main = "\n".join(lines[:idx])
    funcs = {}
    cur = None; buf = []
    for l in lines[idx:]:
        if re.match(r'^function\b', l):
            if cur is not None:
                funcs[cur] = "\n".join(buf)
            m = re.match(r'^function\s+(?:\[?[^=\]]*\]?\s*=\s*)?([A-Za-z_]\w*)', l)
            cur = m.group(1); buf = [l]
        else:
            buf.append(l)
    if cur is not None:
        funcs[cur] = "\n".join(buf)
    return main, funcs


def run(script, outdir, extra_prelude="", timeout=3600):
    name = pathlib.Path(script).stem
    src = open(MEC / script, encoding="utf-8", errors="replace").read()
    main, funcs = split_script(src)
    out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
    fdir = out / f"{name}_fn"; fdir.mkdir(exist_ok=True)
    for fn, body in funcs.items():
        (fdir / f"{fn}.m").write_text(body + "\n", encoding="utf-8")
    # the scripts start with clear; clc which would wipe the path? No: clear does not clear the path.
    (out / f"{name}_main.m").write_text(extra_prelude + "\n" + main + "\n", encoding="utf-8")
    cmd = ["octave", "--no-gui", "--eval",
           f"addpath('{MEC}'); addpath('{COMPAT}'); addpath('{fdir}'); cd('{out}'); "
           f"try, {name}_main; catch e, disp(['ERROR: ' e.message]); for k=1:min(4,numel(e.stack)), disp(e.stack(k)); end; end"]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    txt = p.stdout + p.stderr
    txt = "\n".join(l for l in txt.splitlines() if "X11" not in l and "GUI" not in l)
    (out / f"{name}_octave.log").write_text(txt, encoding="utf-8")
    return txt


if __name__ == "__main__":
    script = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else str(MEC / "octave_out")
    print(run(script, outdir))
