# -*- coding: utf-8 -*-
"""
verify_invariants.py — re-check every physical and numerical invariant this
study relies on, against the artifacts actually on disk.

check_claims.py guards the NUMBERS QUOTED IN PROSE. This guards the PROPERTIES
THE RESULTS MUST HAVE, which is the other half: a value can be quoted correctly
and still be wrong.

Every check here corresponds to a defect that was actually found and fixed in
this project, so each one is a regression test rather than a hypothetical:

  * setup thermodynamics    - flow_conditions.csv was once internally
                              over-determined, with rho and U disagreeing.
  * geometry                - the section is checked against the analytic NACA
                              0012 polynomial and its exact enclosed area.
  * mesh validity           - 22 folded cells were once hidden by an abs() in
                              the area metric; inverted_cells must be zero.
  * mesh spacing precision  - the published spacing table was once rounded by
                              decimal places, so the geometric growth law it
                              documents could not be recovered from it.
  * reconstruction closure  - integrating the surface Cp must return the C_L the
                              reconstruction was given (Blasius). This read
                              -12.4 % before three evaluation errors were fixed.
  * Kutta reference         - the trailing-edge jump must vanish when the
                              inviscid attached circulation is imposed.
  * solver edge cases       - zero, negative and extreme inputs must stay finite.
  * published metrics       - re-derived from the raw time histories.
  * drag                    - instantaneous C_D goes negative (real unsteady
                              thrust), but the CYCLE MEAN must stay positive.
  * response surface        - no published point may sit outside the range the
                              separation law was calibrated on.
  * field bounds            - Cp <= 1 everywhere.
  * experimental provenance - the frame conditions must match the .mat files.
  * dead code               - no unused imports.
  * stale prose             - no comment may describe removed behaviour.

Run from the repository root, after the pipeline:  python3 verify_invariants.py
Exits non-zero on any failure, so run_all.py stops.
"""
import sys, os, glob, json, re, subprocess, ast
import numpy as np, pandas as pd
sys.path.insert(0,'04_solver'); import unistall_solver as us
FAIL=[]
def ck(name, cond, detail=""):
    (print(f"  OK   {name}") if cond else (FAIL.append(name), print(f"  FAIL {name}  {detail}")))

fl=pd.read_csv('03_model_setup/flow_conditions.csv').set_index('parameter')
mA=pd.read_csv('05_solution/metrics_A_validation.csv').set_index('metric')['value']
mB=pd.read_csv('05_solution/metrics_B_application.csv').set_index('metric')['value']
G='01_geometry/naca0012_coordinates.csv'
c,U,M=[float(fl['case_A_validation'][k]) for k in ('chord_c','freestream_velocity_U','freestream_mach_M')]

# --- physics of the setup
R,g=287.05,1.4
for case in ('case_A_validation','case_B_application'):
    col=fl[case]; T=float(col['static_temperature_T_inf']); p=float(col['static_pressure_p_inf'])
    ck(f"{case[5:6]} rho=p/RT", abs(p/(R*T)-float(col['air_density_rho']))/float(col['air_density_rho'])<2e-4)
    ck(f"{case[5:6]} a=sqrt(gRT)", abs(np.sqrt(g*R*T)-float(col['speed_of_sound_a']))/float(col['speed_of_sound_a'])<2e-4)
    ck(f"{case[5:6]} M=U/a", abs(float(col['freestream_velocity_U'])/float(col['speed_of_sound_a'])-float(col['freestream_mach_M']))<1e-4)

# --- geometry
gg=pd.read_csv(G); xc,yc=gg['x_over_c'].values,gg['y_over_c'].values
yt=lambda x:0.6*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
ck("geometry matches NACA0012 polynomial", np.max(np.abs(np.abs(yc[yc>0])-yt(np.clip(xc[yc>0],0,1))))<1e-5)
X,Y=np.append(xc,xc[0]),np.append(yc,yc[0])
ck("enclosed area vs analytic", abs(0.5*abs(np.sum(X[:-1]*Y[1:]-X[1:]*Y[:-1]))-0.082210)<1e-4)

# --- mesh
q=pd.read_csv('02_mesh/mesh_quality_metrics.csv').set_index('metric')['value']
ck("mesh inverted_cells == 0", int(float(q['inverted_cells']))==0)
r=pd.read_csv('02_mesh/mesh_radial_spacing.csv'); yv=r['normal_coord_chords'].values
gr=np.diff(yv)[1:]/np.diff(yv)[:-1]
ck("growth law recoverable from published file", gr.max()-gr.min()<1e-4, f"spread {gr.max()-gr.min():.2e}")
ck("growth ratio matches metric", abs(gr.mean()-float(q['wall_normal_growth_ratio']))<1e-3)

# --- the mesh's own quality metrics must be recomputable from its published nodes.
#     They were not: x_m/y_m were rounded to 6 dp, a 1e-06 m quantum against a
#     3.94e-06 m first cell, so recomputing gave max AR 1183 against the true 952.
_nd=pd.read_csv('02_mesh/mesh_nodes.csv')
_I,_J=int(float(q['i_nodes_wrap'])),int(float(q['j_nodes_normal']))
_X=np.full((_J,_I),np.nan); _Y=np.full((_J,_I),np.nan)
_X[_nd.j.values,_nd.i.values]=_nd.x_m.values; _Y[_nd.j.values,_nd.i.values]=_nd.y_m.values
_x1,_y1=_X[:-1,:-1],_Y[:-1,:-1]; _x2,_y2=_X[:-1,1:],_Y[:-1,1:]
_x3,_y3=_X[1:,1:],_Y[1:,1:];     _x4,_y4=_X[1:,:-1],_Y[1:,:-1]
_ar_area=0.5*((_x1*_y2-_x2*_y1)+(_x2*_y3-_x3*_y2)+(_x3*_y4-_x4*_y3)+(_x4*_y1-_x1*_y4))
_e=[np.hypot(_x2-_x1,_y2-_y1),np.hypot(_x3-_x2,_y3-_y2),np.hypot(_x4-_x3,_y4-_y3),np.hypot(_x1-_x4,_y1-_y4)]
_arr=np.maximum.reduce(_e)/np.maximum(np.minimum.reduce(_e),1e-30)
ck("mesh max aspect ratio recomputes from published nodes",
   abs(_arr.max()-float(q['max_aspect_ratio']))<0.5, f"{_arr.max():.1f} vs {float(q['max_aspect_ratio']):.1f}")
ck("mesh inverted cells recompute from published nodes",
   int(np.sum(_ar_area*np.sign(np.median(_ar_area))<=0))==int(float(q['inverted_cells'])))
ck("near-wall growth ratio recovers from published nodes",
   abs(np.hypot(np.diff(_X[:6,128]),np.diff(_Y[:6,128]))[1]/
       np.hypot(np.diff(_X[:6,128]),np.diff(_Y[:6,128]))[0]
       -float(q['wall_normal_growth_ratio']))<1e-3)

# --- reconstruction invariants
for a,CL in ((2.,0.22),(10.,1.10),(17.5,1.91)):
    _,pct,_=us.surface_load_closure(G,c,U,M,a,CL,0.,0.)
    ck(f"closure |err|<1.5% at alpha={a}", abs(pct)<1.5, f"{pct:+.2f}%")
ckl=us.kutta_reference_CL(G,c,U,M,10.0)
_,_,tj=us.surface_load_closure(G,c,U,M,10.0,ckl,0.,0.)
ck("TE jump vanishes at CL_kutta", tj<5e-3, f"{tj:.4f}")
for a,CL,CNv,tau in ((0.,0.,0.,0.),(-10.,-1.1,0.,0.),(45.,2.5,0.,0.),(17.5,1.91,0.24,5.0)):
    x,cp,_=us.surface_cp(G,c,U,M,a,CL,CNv,tau)
    ck(f"surface_cp finite at edge case a={a},CL={CL}", np.all(np.isfinite(cp)))

# --- metrics re-derived
for case,m in (('A_validation',mA),('B_application',mB)):
    th=pd.read_csv(f'05_solution/time_history_{case}.csv').iloc[-720:]
    a_,CL_,CM_,CD_=(th[k].values for k in ('alpha_deg','CL','CM_c4','CD'))
    ck(f"{case} CL_max", abs(CL_.max()-float(m['CL_max_dynamic']))<1e-3)
    ck(f"{case} CM_min", abs(CM_.min()-float(m['CM_min(c/4)']))<1e-3)
    ck(f"{case} CD_cycle_mean", abs(CD_.mean()-float(m['CD_cycle_mean']))<1e-3)
    ar=np.radians(np.append(a_,a_[0])); CMc=np.append(CM_,CM_[0])
    ck(f"{case} damping Xi", abs(-np.trapz(CMc,ar)-float(m['aero_damping_Xi']))<2e-4)
    ck(f"{case} CD cycle mean positive (no net propulsion)", CD_.mean()>0)

# --- response surface within calibration
rs=pd.read_csv('05_solution/response_surface.csv')
ck("response surface fully calibrated", bool(rs['within_calibration'].all()))
ck("peak alpha <= polar range", rs['peak_alpha_deg'].max()<=float(pd.read_csv('03_model_setup/static_polar_reference.csv')['alpha_deg'].max())+1e-9)

# --- fields
for f in sorted(glob.glob('05_solution/field_*.csv')):
    d=pd.read_csv(f); ck(f"{f.split('/')[-1]} Cp<=1", bool((d['Cp'].dropna()<=1+1e-9).all()))

# --- provenance
import scipy.io
res=pd.read_csv('06_postprocessing/validation/validation_nasa_real.csv')
bad=0
for _,row in res.iterrows():
    dd=scipy.io.loadmat(f"06_postprocessing/validation/experimental/nasa_frames/{row['frame']}.mat")
    gv=lambda k: float(np.asarray(dd[k]).ravel()[0])
    if abs(np.degrees(gv('alpha_0'))-float(row['alpha0_deg']))>0.06 or abs(gv('M')-float(row['M']))>5e-4: bad+=1
ck("all 6 frame provenances match the .mat", bad==0)

# --- dead code / imports
tot=0
for f in sorted(glob.glob('0*/**/*.py',recursive=True)+glob.glob('*.py')):
    src=open(f,encoding='utf-8').read(); tree=ast.parse(src)
    imp={}
    for n in ast.walk(tree):
        if isinstance(n,ast.Import):
            for al in n.names: imp[(al.asname or al.name).split('.')[0]]=1
        elif isinstance(n,ast.ImportFrom):
            for al in n.names:
                if al.name!='*': imp[al.asname or al.name]=1
    used={n.id for n in ast.walk(tree) if isinstance(n,ast.Name)}|{n.attr for n in ast.walk(tree) if isinstance(n,ast.Attribute)}
    tot+=len([k for k in imp if k not in used and len(re.findall(r'\b'+re.escape(k)+r'\b',src))<2])
ck("no unused imports", tot==0, f"{tot} found")

# --- stale prose
# This file is EXCLUDED from its own scan. It necessarily quotes the phrases it
# forbids, so including it makes every stale-prose check match itself and fail --
# which is exactly what happened the first time this file was committed, since
# it had passed while still untracked and git ls-files could not see it.
_SELF = os.path.basename(__file__)
allsrc="".join(open(f,encoding='utf-8',errors='replace').read()
               for f in subprocess.run(['git','ls-files'],capture_output=True,text=True).stdout.split()
               if f.endswith(('.py','.md','.json')) and os.path.basename(f) != _SELF)
for pat in (r'surface_cp probes at', r'demands it be zero', r'the 0\.015c\s*used', r'clip is kept'):
    ck(f"no stale prose: {pat}", len(re.findall(pat,allsrc))==0)

print()
print(f"  ==> {len(FAIL)} FAILURES" if FAIL else "  ==> ALL CHECKS PASS")
sys.exit(1 if FAIL else 0)
