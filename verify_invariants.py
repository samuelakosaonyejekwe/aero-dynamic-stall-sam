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
                              over-determined, with rho and U disagreeing; rho,
                              a, M and Re are all re-derived from the
                              independent quantities, and the kinematics are
                              re-derived from the conditions AS PUBLISHED.
  * geometry                - the section is checked against the analytic NACA
                              0012 polynomial and its exact enclosed area.
  * mesh validity           - 22 folded cells were once hidden by an abs() in
                              the area metric; inverted_cells must be zero, and
                              the worst cells must be where the header says
                              they are (it named the wrong place for a while).
  * mesh spacing precision  - the published spacing table was once rounded by
                              decimal places, so the geometric growth law it
                              documents could not be recovered from it.
  * panel solver vs exact   - a circle's Cp is 1 - 4 sin^2(theta) in closed
                              form, so the panel routines are checked against
                              something outside themselves, and the second-order
                              convergence is pinned too. Symmetry and the
                              Kirchhoff round-trip likewise.
  * reconstruction closure  - integrating the surface Cp must return the C_L the
                              reconstruction was given (Blasius). This read
                              -12.4 % before three evaluation errors were
                              fixed, and the three figures the solver docstring
                              publishes are parsed out of it and re-measured,
                              because two further copies of them had gone stale.
  * Kutta reference         - the trailing-edge jump must vanish when the
                              inviscid attached circulation is imposed.
  * solver edge cases       - zero, negative and extreme inputs must stay finite.
  * physical limits         - zero amplitude must give exactly the static
                              polar, and k -> 0 must collapse the hysteresis
                              loop onto it. Both check the march against
                              physics rather than against its own output.
  * published metrics       - re-derived from the raw time histories.
  * drag                    - instantaneous C_D goes negative (real unsteady
                              thrust), but the CYCLE MEAN must stay positive.
  * response surface        - no published point may sit outside the range the
                              separation law was calibrated on.
  * reconstruction scaling  - its non-dimensional outputs must not depend on
                              chord, speed or Mach, as the march's do not.
  * dynamic-stall vortex    - its circulation and core are DERIVED from the
                              shear-layer vorticity flux and the model's own
                              vortex clock; every published number recomputes
                              from that, the core lands in the measured -3 to
                              -6 band, and the field shows the flow reversal
                              beneath it.
  * field bounds            - Cp <= 1 everywhere, the flow stays subsonic and
                              the static temperature stays positive, on the
                              grid the config declares the fields are written
                              on.
  * closure bound           - the cycle-wide closure metric must bound the
                              peak-lift one; the peak-lift figure was once
                              quoted as though it bounded the whole cycle.
  * report equations        - the report must state the Cp and impulsive-lag
                              expressions the solver actually evaluates.
  * dossier completeness    - a matplotlib page that overruns is CLIPPED, in
                              silence: the solver config printed 82 of its 97
                              lines and every long table 33 of the 38 rows its
                              caption claimed.
  * report completeness     - every paragraph, table cell and image of
                              case.docx must reach the rendered PDF, and no
                              page may overlap its own content or run off it.
  * experimental provenance - the frame conditions must match the .mat files.
  * config vs code          - solver_config.json must still describe the
                              solver: f_min, n_panels, the default grid, the
                              domain and the masked ring were duplicated
                              literals in both, and are now imported.
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
    # Re is derived and published to four significant figures, like rho and a
    # above; it was the one derived free-stream quantity nothing re-checked.
    _re=float(col['air_density_rho'])*float(col['freestream_velocity_U'])*float(col['chord_c'])/float(col['dynamic_viscosity_mu'])
    ck(f"{case[5:6]} Re=rho*U*c/mu", abs(_re-float(col['reynolds_number_Re_c']))/_re<5e-4,
       f"{_re:.4g} vs {col['reynolds_number_Re_c']}")

# --- the kinematics must recompute from the flow conditions AS PUBLISHED.
#     They did not: omega was derived from the full-precision U while
#     flow_conditions.csv published U rounded to 2 dp, so kinematics.csv said
#     68.058 rad/s and the solver -- which reads the published U and recomputes
#     omega = 2kU/c itself -- marched at 68.060. Two published files, two
#     different motions, neither recomputable from the other.
_kin=pd.read_csv('03_model_setup/kinematics.csv').set_index('case_id')
for _cid,_col in (('A_validation_rig','case_A_validation'),
                  ('B_application_rotor','case_B_application')):
    _k=_kin.loc[_cid]; _U=float(fl[_col]['freestream_velocity_U']); _c=float(fl[_col]['chord_c'])
    _w=2.0*float(_k['reduced_freq_k'])*_U/_c
    ck(f"{_cid} omega recomputes from the published U and c",
       abs(_w-float(_k['omega_rad_s']))<5e-4, f"{_w:.5f} vs {_k['omega_rad_s']}")
    ck(f"{_cid} frequency and period agree with that omega",
       abs(_w/(2*np.pi)-float(_k['freq_Hz']))<5e-4
       and abs(2*np.pi/_w-float(_k['period_s']))<5e-6)

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
# The header names WHERE the worst cells are, after having named the wrong place
# for a long time ("near the trailing edge"). Re-measure it, so the explanation
# cannot go stale again: every one of the twenty worst-aspect-ratio cells must
# be in the first wall-normal layer, and none of them at the trailing edge.
_top=np.dstack(np.unravel_index(np.argsort(_arr.ravel())[-20:], _arr.shape))[0]
ck("the worst cells are all in the first wall-normal layer",
   bool(all(int(t[0])==0 for t in _top)), f"layers {sorted({int(t[0]) for t in _top})}")
ck("the worst cells are at mid-chord, not at the trailing edge",
   bool(all(0.15*_I < int(t[1]) < 0.85*_I for t in _top)),
   f"wrap indices {sorted({int(t[1]) for t in _top})[:4]}")
# and the far-field azimuthal spacing the metrics now publish must recompute
_thf=np.unwrap(np.arctan2(_Y[-1,:]-0.0, _X[-1,:]-0.5*c))
_dth=np.abs(np.diff(_thf)); _dth=_dth[_dth>0]
ck("far-field angular spacing ratio recomputes from published nodes",
   abs(_dth.max()/_dth.min()-float(q['farfield_angular_spacing_ratio']))<0.5,
   f"{_dth.max()/_dth.min():.1f} vs {q['farfield_angular_spacing_ratio']}")
ck("near-wall growth ratio recovers from published nodes",
   abs(np.hypot(np.diff(_X[:6,128]),np.diff(_Y[:6,128]))[1]/
       np.hypot(np.diff(_X[:6,128]),np.diff(_Y[:6,128]))[0]
       -float(q['wall_normal_growth_ratio']))<1e-3)

# --- THE PANEL SOLVER AGAINST AN EXACT SOLUTION. Every other check on the
#     reconstruction is internal -- closure against the circulation it was
#     handed, the Kutta reference, convergence under refinement -- so all of
#     them would pass a solver that was self-consistently wrong. A circle in
#     uniform flow has Cp = 1 - 4 sin^2(theta) in closed form, so feeding the
#     panel routines a circle tests them against something outside themselves.
#     It also pins the ORDER: refining 160 -> 640 panels must cut the error
#     roughly fourfold, which a first-order bug would not do.
import tempfile as _tf
_th=np.linspace(0,2*np.pi,401)
_cd=pd.DataFrame({"x_over_c":0.5+0.5*np.cos(_th),"y_over_c":0.5*np.sin(_th)})
_cf=_tf.NamedTemporaryFile(suffix='.csv',delete=False,mode='w')
_cd.to_csv(_cf.name,index=False); _cf.close()
_orig_af=us._airfoil_surface
_cerr={}
for _np_ in (160,640):
    us._airfoil_surface=lambda csv,cc,n_panel=_np_,_o=_orig_af,_n=_np_: _o(csv,cc,_n)
    _,_cp,_=us.surface_cp(_cf.name,1.0,100.0,0.0,0.0,0.0,0.0,0.0)
    us._airfoil_surface=_orig_af
    _xp,_yp=_orig_af(_cf.name,1.0,_np_)
    _xcc=0.5*(_xp[:-1]+_xp[1:]); _ycc=0.5*(_yp[:-1]+_yp[1:])
    _exact=1.0-4.0*np.sin(np.arctan2(_ycc,_xcc-0.5))**2
    _cerr[_np_]=float(np.abs(_cp-_exact).max())
os.unlink(_cf.name)
ck("panel Cp matches the exact circle solution at 160 panels",
   _cerr[160]<0.06, f"max|dCp| {_cerr[160]:.4f}")
ck("panel Cp matches the exact circle solution at 640 panels",
   _cerr[640]<0.02, f"max|dCp| {_cerr[640]:.4f}")
ck("the panel discretisation is second order (4x refinement -> ~4x less error)",
   3.0 < _cerr[160]/_cerr[640] < 5.5, f"ratio {_cerr[160]/_cerr[640]:.2f}")

# --- the symmetric section must behave symmetrically. A sign or index error in
#     the panel normals, the sheet, or the surface self-terms would break this
#     while leaving every self-referential check intact.
for _asym,_clsym in ((5.,0.55),(12.,1.32)):
    _,_c1,_=us.surface_cp(G,c,U,M, _asym, _clsym,0,0)
    _,_c2,_=us.surface_cp(G,c,U,M,-_asym,-_clsym,0,0)
    ck(f"surface Cp mirrors between alpha=+-{_asym}", float(np.abs(_c1-_c2[::-1]).max())<1e-9,
       f"{np.abs(_c1-_c2[::-1]).max():.2e}")
    _l1,_,_=us.surface_load_closure(G,c,U,M, _asym, _clsym,0,0)
    _l2,_,_=us.surface_load_closure(G,c,U,M,-_asym,-_clsym,0,0)
    ck(f"integrated CL is antisymmetric at alpha=+-{_asym}", abs(_l1+_l2)<1e-9,
       f"{_l1:+.6f} vs {_l2:+.6f}")

# --- the Kirchhoff inversion must round-trip: the f it recovers from the
#     measured polar must reproduce that polar through the forward relation.
_sref=pd.read_csv('03_model_setup/static_polar_reference.csv')
_CNa=json.load(open('03_model_setup/solver_config.json'))['lift_curve_slope_CNalpha_per_rad']
_fst=us.calibrate_separation(_sref.alpha_deg,_sref.Cl,_sref.Cd,_CNa)
_ar=np.radians(_sref.alpha_deg.values)
_CNref=_sref.Cl.values*np.cos(_ar)+_sref.Cd.values*np.sin(_ar)
_fv=_fst(_sref.alpha_deg.values)
_CNfwd=_CNa*((1+np.sqrt(_fv))/2)**2*_ar
_mk=(_fv>us.F_MIN+1e-9)&(_sref.alpha_deg.values>0.5)
ck("Kirchhoff inversion round-trips to the measured polar",
   float(np.abs(_CNfwd-_CNref)[_mk].max())<1e-12,
   f"max |dCN| {np.abs(_CNfwd-_CNref)[_mk].max():.2e} over {int(_mk.sum())} points")

# --- reconstruction invariants
# The three closure figures are PARSED OUT of the solver's own docstring and
# checked against a fresh measurement, rather than restated here. Two other
# copies of these numbers, in surface_cp and _surface_velocity, had gone stale
# (-0.18 % against the header's -0.37 % for the same condition) because nothing
# re-measured them. The docstring is the source; this is what keeps it true.
_hdr=re.search(r'([-\d.]+) % at alpha 2 deg \(C_L ([\d.]+)\).*?'
               r'([-\d.]+) % at 10 deg \(C_L ([\d.]+)\).*?'
               r'([-\d.]+) % at 17\.5 deg \(C_L ([\d.]+)\)',
               us.__doc__, re.S)
ck("solver docstring states its three closure figures with their C_L", _hdr is not None)
_trip=([(2.,float(_hdr.group(2)),float(_hdr.group(1))),
        (10.,float(_hdr.group(4)),float(_hdr.group(3))),
        (17.5,float(_hdr.group(6)),float(_hdr.group(5)))] if _hdr else
       [(2.,0.22,None),(10.,1.10,None),(17.5,1.91,None)])
for a,CL,want in _trip:
    _,pct,_=us.surface_load_closure(G,c,U,M,a,CL,0.,0.)
    ck(f"closure |err|<1.5% at alpha={a}", abs(pct)<1.5, f"{pct:+.2f}%")
    if want is not None:
        ck(f"closure at alpha={a} matches the figure the docstring publishes",
           abs(pct-want)<0.01, f"measured {pct:+.3f}% vs documented {want:+.2f}%")
ckl=us.kutta_reference_CL(G,c,U,M,10.0)
_,_,tj=us.surface_load_closure(G,c,U,M,10.0,ckl,0.,0.)
ck("TE jump vanishes at CL_kutta", tj<5e-3, f"{tj:.4f}")
for a,CL,CNv,tau in ((0.,0.,0.,0.),(-10.,-1.1,0.,0.),(45.,2.5,0.,0.),(17.5,1.91,0.24,5.0)):
    x,cp,_=us.surface_cp(G,c,U,M,a,CL,CNv,tau)
    ck(f"surface_cp finite at edge case a={a},CL={CL}", np.all(np.isfinite(cp)))

# --- PHYSICAL LIMITS OF THE MARCH. Everything else here checks the march
#     against its own published output; these check it against physics. A
#     reduced-order unsteady model has two limits it must hit exactly, and
#     nothing tested either.
_cfgm=json.load(open('03_model_setup/solver_config.json'))
_CNam=_cfgm['lift_curve_slope_CNalpha_per_rad']
_cst=dict(**_cfgm['indicial_circulatory'],**_cfgm['time_constants_semichords'])
_cst.update({k:v for k,v in _cfgm['calibrated_constants'].items() if k!='comment'})
_srf=pd.read_csv('03_model_setup/static_polar_reference.csv')
_fsm=us.calibrate_separation(_srf.alpha_deg,_srf.Cl,_srf.Cd,_CNam)
# the indicial deficiency must relax to the full circulatory load: A1+A2 = 1,
# or the step response never reaches its steady value
_ai=_cfgm['indicial_circulatory']
ck("indicial A1+A2 = 1 (the step response relaxes to unity)",
   abs(_ai['A1']+_ai['A2']-1.0) < 1e-12, f"A1+A2 = {_ai['A1']+_ai['A2']}")
# ZERO AMPLITUDE: no motion, so no unsteady content at all -- the loads must be
# constant and equal the static polar the separation law was fitted to.
_o0=us.solve_dynamic_stall(8.0,0.0,0.10,0.30,0.30,102.09,_fsm,CNalpha=_CNam,
                           consts=_cst,n_per_cycle=360,n_cycles=3)
_f8=float(_fsm(8.0)); _a8=np.radians(8.0)
_cl8=(_CNam*((1+np.sqrt(_f8))/2)**2*_a8)*np.cos(_a8) \
     + _cst['eta']*_CNam*_a8**2*np.sqrt(_f8)*np.sin(_a8)
ck("zero-amplitude motion produces no unsteady content",
   float(_o0['CL'].max()-_o0['CL'].min())<1e-12, f"spread {_o0['CL'].max()-_o0['CL'].min():.2e}")
ck("zero-amplitude load equals the static polar",
   abs(float(_o0['CL'][0])-_cl8)<1e-9, f"{_o0['CL'][0]:.6f} vs {_cl8:.6f}")
# QUASI-STEADY LIMIT: as the reduced frequency goes to zero the hysteresis loop
# must collapse onto that same static polar.
_oq=us.solve_dynamic_stall(8.0,2.0,1e-3,0.30,0.30,102.09,_fsm,CNalpha=_CNam,
                           consts=_cst,n_per_cycle=720,n_cycles=6)
_aq=_oq['alpha_deg']; _fq=_fsm(_aq); _arq=np.radians(_aq)
_clq=(_CNam*((1+np.sqrt(_fq))/2)**2*_arq)*np.cos(_arq) \
     + _cst['eta']*_CNam*_arq**2*np.sqrt(_fq)*np.sin(_arq)
ck("k -> 0 collapses the loop onto the static polar",
   float(np.abs(_oq['CL']-_clq).max())<5e-3, f"max |dCL| {np.abs(_oq['CL']-_clq).max():.5f}")
ck("k -> 0 closes the hysteresis loop",
   abs(float(us._trapz(_oq['CL'],_arq)))<2e-3, f"loop area {abs(us._trapz(_oq['CL'],_arq)):.2e}")

# --- metrics re-derived
for case,m in (('A_validation',mA),('B_application',mB)):
    th=pd.read_csv(f'05_solution/time_history_{case}.csv').iloc[-720:]
    a_,CL_,CM_,CD_=(th[k].values for k in ('alpha_deg','CL','CM_c4','CD'))
    ck(f"{case} CL_max", abs(CL_.max()-float(m['CL_max_dynamic']))<1e-3)
    ck(f"{case} CM_min", abs(CM_.min()-float(m['CM_min(c/4)']))<1e-3)
    ck(f"{case} CD_cycle_mean", abs(CD_.mean()-float(m['CD_cycle_mean']))<1e-3)
    ar=np.radians(np.append(a_,a_[0])); CMc=np.append(CM_,CM_[0])
    # us._trapz, not np.trapz: np.trapz was REMOVED in NumPy 2.0 and
    # requirements.txt admits numpy<3, so this line was the one place in the
    # pipeline that would have raised AttributeError on a permitted NumPy. The
    # solver already binds whichever name exists.
    ck(f"{case} damping Xi", abs(-us._trapz(CMc,ar)-float(m['aero_damping_Xi']))<2e-4)
    ck(f"{case} CD cycle mean positive (no net propulsion)", CD_.mean()>0)

# --- the PUBLISHED cp_distribution must integrate back to the solver's own C_L at
#     the same phase. This guards the written artifact rather than the code that
#     wrote it: the two mesh files that failed this class of check were both cases
#     of data published in a form that could not reproduce its own metric.
_geo=pd.read_csv(G); _gu=_geo[_geo.y_over_c>=0].sort_values('x_over_c'); _gl=_geo[_geo.y_over_c<=0].sort_values('x_over_c')
for _case in ('A_validation','B_application'):
    _cc=float(fl['case_'+_case[0]+('_validation' if _case[0]=='A' else '_application')]['chord_c'])
    _d=pd.read_csv(f'05_solution/cp_distribution_{_case}.csv')
    _th=pd.read_csv(f'05_solution/time_history_{_case}.csv').iloc[-720:]
    _ad=np.gradient(_th.alpha_deg.values)
    _worst=0.0
    for _ph,_g in _d.groupby('phase_tag'):
        _a=float(_g['alpha_deg'].iloc[0]); _dn='down' in _ph
        _sub=_th[(_ad<0) if _dn else (_ad>0)]
        _CLs=float(_th.loc[(_sub.alpha_deg-_a).abs().idxmin(),'CL'])
        _u=_g[_g.surface=='upper'].sort_values('x_c'); _l=_g[_g.surface=='lower'].sort_values('x_c',ascending=False)
        _X=np.concatenate([_u.x_c.values,_l.x_c.values])*_cc
        _Y=np.concatenate([np.interp(_u.x_c.values,_gu.x_over_c,_gu.y_over_c),
                           np.interp(_l.x_c.values,_gl.x_over_c,_gl.y_over_c)])*_cc
        _CP=np.concatenate([_u.Cp.values,_l.Cp.values])
        _Xc,_Yc,_CPc=np.append(_X,_X[0]),np.append(_Y,_Y[0]),np.append(_CP,_CP[0])
        _cpm=0.5*(_CPc[1:]+_CPc[:-1]); _ar=np.radians(_a)
        _sg=1.0 if 0.5*np.sum(_Xc[:-1]*_Yc[1:]-_Xc[1:]*_Yc[:-1])>0 else -1.0
        _CN=_sg*np.sum(_cpm*np.diff(_Xc))/_cc; _CA=-_sg*np.sum(_cpm*np.diff(_Yc))/_cc
        _cl=_CN*np.cos(_ar)-_CA*np.sin(_ar)
        _worst=max(_worst, abs(100*(_cl-_CLs)/_CLs))
    ck(f"{_case} published Cp integrates to published CL (<2%)", _worst<2.0, f"worst {_worst:.2f}%")

# --- the model static polar must say which of its rows are calibrated. It runs
#     past the reference's last tabulated incidence, and the calibration figures
#     drew the extrapolated tail exactly like the fitted part.
_mp=pd.read_csv('05_solution/model_static_polar.csv')
_acal=float(pd.read_csv('03_model_setup/static_polar_reference.csv')['alpha_deg'].max())
ck("model polar flags its calibrated rows", 'within_calibration' in _mp.columns)
if 'within_calibration' in _mp.columns:
    ck("model polar within_calibration agrees with the reference range",
       bool((_mp['within_calibration'] == (_mp['alpha_deg'] <= _acal + 1e-9)).all()))
    ck("the model polar does extend past calibration (so the flag is not vacuous)",
       bool((~_mp['within_calibration']).any()), "no extrapolated rows")

# --- response surface within calibration
rs=pd.read_csv('05_solution/response_surface.csv')
ck("response surface fully calibrated", bool(rs['within_calibration'].all()))
ck("peak alpha <= polar range", rs['peak_alpha_deg'].max()<=float(pd.read_csv('03_model_setup/static_polar_reference.csv')['alpha_deg'].max())+1e-9)

# --- fields
_fcfg=json.load(open('03_model_setup/solver_config.json'))['field_reconstruction']
_nxy=_fcfg['grid_nx_solution']*_fcfg['grid_ny_solution']
for f in sorted(glob.glob('05_solution/field_*.csv')):
    d=pd.read_csv(f); ck(f"{f.split('/')[-1]} Cp<=1", bool((d['Cp'].dropna()<=1+1e-9).all()))
    # The reconstruction is subsonic-only: the compressibility it carries is the
    # Prandtl-Glauert beta of the indicial march, and the thermal module is the
    # isentropic stagnation relation. A published field that went sonic would
    # invalidate both, and the report prints local-Mach contours as if it could
    # not. Measured worst over the eight fields: M = 0.78, T = 261 K.
    ck(f"{f.split('/')[-1]} stays subsonic", bool((d['Mach_local'].dropna()<1.0).all()),
       f"max M {d['Mach_local'].max():.3f}")
    ck(f"{f.split('/')[-1]} static temperature stays positive",
       bool((d['T_static_K'].dropna()>0).all()), f"min T {d['T_static_K'].min():.1f} K")
    # the config DECLARES grid_*_solution as "the sizes written to
    # 05_solution/field_*.csv"; the writer used to carry its own literals, and
    # the DSV-core metric was measured on the coarser DEFAULT grid instead, so
    # it described no field this study ships.
    ck(f"{f.split('/')[-1]} is on the declared solution grid", len(d)==_nxy,
       f"{len(d)} rows vs {_nxy}")

# --- the dynamic-stall vortex is DERIVED, and every number it publishes must
#     recompute from that derivation. It used to carry two chosen constants,
#     which gave a core an order of magnitude too shallow, no surface footprint
#     and no flow reversal -- while three places in the study asserted one.
for _case,_m in (('A_validation',mA),('B_application',mB)):
    _col='case_'+('A_validation' if _case[0]=='A' else 'B_application')
    _cc=float(fl[_col]['chord_c']); _uu=float(fl[_col]['freestream_velocity_U'])
    _mm=float(fl[_col]['freestream_mach_M'])
    _thc=pd.read_csv(f'05_solution/time_history_{_case}.csv')
    _rc_=_thc.iloc[int(_thc.CN_vortex.idxmax())]
    _tvl=json.load(open('03_model_setup/solver_config.json'))['calibrated_constants']['Tvl']
    _st=us.dsv_vortex_state(G,_cc,_uu,_mm,_rc_.alpha_deg,_rc_.CL,_rc_.CN_vortex,
                            _rc_.tau_v_semichords/_tvl,_tvl)
    for _k,_row in (('Gamma_over_Uc','DSV_circulation_over_Uc'),
                    ('rc_chords','DSV_core_radius_chords'),
                    ('peak_swirl_over_U','DSV_peak_swirl_over_U'),
                    ('induced_at_wall_over_U','DSV_induced_at_wall_over_U'),
                    ('edge_speed_over_U','DSV_edge_speed_over_U')):
        ck(f"{_case} {_row} recomputes from the derivation",
           abs(_st[_k]-float(_m[_row]))<1e-3, f"{_st[_k]:.4f} vs {_m[_row]}")
    # the core must be self-consistent: peak swirl = the edge speed it was
    # derived from, and rc = LAMB_OSEEN_PEAK*Gamma/(2 pi Ve)
    ck(f"{_case} peak swirl equals the shear-layer edge speed",
       abs(_st['peak_swirl_over_U']-_st['edge_speed_over_U'])<2e-3)
    ck(f"{_case} core radius follows from the circulation and that speed",
       abs(_st['rc_chords']-us.LAMB_OSEEN_PEAK*_st['Gamma_over_Uc']
           /(2*np.pi*_st['edge_speed_over_U']))<1e-4)
    # the circulation budget must close: the vortex carries exactly the share
    # of the total circulation the model attributes to it, CNv/CN, because
    # Gamma_v = 0.5*CNv*U*c is the same Kutta-Joukowski relation the bound sheet
    # uses. A vortex sized from the shear-layer vorticity flux instead --
    # 1.79*U*c, tried during this audit -- is nearly twice the whole
    # circulation and inflates the body's integrated lift by 60 %.
    _share = float(_rc_.CN_vortex/_rc_.CN)
    ck(f"{_case} the vortex carries the CNv/CN share of the circulation",
       abs(_st['Gamma_over_Uc']/(0.5*float(_rc_.CL)) - _share) < 2e-3,
       f"{_st['Gamma_over_Uc']/(0.5*float(_rc_.CL)):.4f} vs {_share:.4f}")
    ck(f"{_case} the vortex barely moves the body's integrated lift",
       abs(float(_m['DSV_induced_lift_dCL'])) < 0.05, _m['DSV_induced_lift_dCL'])
    # it does NOT turn the flow over at the wall, and that is a measured
    # comparison now -- what it induces there against what it must overcome --
    # rather than the assertion three places in this study used to make
    ck(f"{_case} the vortex cannot reverse the flow at the wall beneath it",
       _st['induced_at_wall_over_U'] < _st['edge_speed_over_U'],
       f"induced {_st['induced_at_wall_over_U']:.3f} U vs local {_st['edge_speed_over_U']:.2f} U")
    # the published core depth must land where measurement puts a real one
    ck(f"{_case} core suction is in the measured -3 to -6 band",
       -6.5 <= float(_m['Cp_DSV_core_min']) <= -2.5, f"{_m['Cp_DSV_core_min']}")
    # the field does not resolve that core, and publishes by how much
    ck(f"{_case} the under-resolution of the core is published",
       0.0 < float(_m['DSV_core_radius_cells']) < 3.0, _m['DSV_core_radius_cells'])
    # THE COUPLING MUST BE CONVERGED. The core radius depends on the edge speed
    # and the panel solution depends on the vortex, so a single sweep publishes
    # an edge speed belonging to the solution BEFORE the vortex entered the
    # boundary condition -- 1.098*U was shipped that way against a converged
    # 1.057*U, with the core radius 3.9 % out. Rebuild the state from the
    # solution the reconstruction actually uses and require it to be unchanged.
    _a_ = np.radians(float(_rc_.alpha_deg))
    _sol = us._solve_with_dsv(G, _cc, _uu, _a_, float(_rc_.CL),
                              float(_rc_.CN_vortex), float(_rc_.tau_v_semichords)/_tvl)
    _re = us._dsv_state(_sol[2], _sol[3], _sol[4], _sol[5], _sol[6], _uu, _a_, _cc,
                        float(_rc_.CN_vortex), float(_rc_.tau_v_semichords)/_tvl)
    ck(f"{_case} the vortex/panel coupling is converged, not a single sweep",
       abs(_re[4]/_uu - _st['edge_speed_over_U']) < 1e-9
       and abs(_re[3]/_cc - _st['rc_chords']) < 1e-9,
       f"resolved Ve/U {_re[4]/_uu:.6f} vs published {_st['edge_speed_over_U']:.6f}")
# the Lamb-Oseen peak-swirl coefficient must be the swirl one (0.6382), not the
# enclosed-circulation one (0.7153) -- they differ by 12 % and were confused
ck("LAMB_OSEEN_PEAK is the peak-swirl coefficient",
   abs(us.LAMB_OSEEN_PEAK-0.638173)<1e-5, f"{us.LAMB_OSEEN_PEAK:.6f}")
ck("config Tvl literature default comes from the solver",
   json.load(open('03_model_setup/solver_config.json'))['time_constants_semichords']['Tvl']
   == us.TVL_DEFAULT)

# --- the RECONSTRUCTION must be chord- and speed-independent in its
#     non-dimensional outputs, the way the march already is (the validation
#     stage asserts that for the march). Nothing checked it for the
#     reconstruction, and the vortex derivation added two more places where a
#     stray dimensional term could hide.
_ref=None
for _cc2,_uu2,_mm2 in ((0.30,102.0,0.30),(1e-4,102.0,0.30),(50.0,102.0,0.30),
                       (0.30,1e-3,0.30),(0.30,102.0,0.90)):
    _s2=us.dsv_vortex_state(G,_cc2,_uu2,_mm2,17.5,1.9,0.24,1.0,6.0)
    _,_,_c2=us.dsv_core_cp(G,_cc2,_uu2,_mm2,17.5,1.9,0.24,1.0,6.0)
    _v2=(round(_s2['rc_chords'],9), round(_s2['Gamma_over_Uc'],9),
         round(_s2['peak_swirl_over_U'],9), round(_c2,7))
    if _ref is None: _ref=_v2
    ck(f"reconstruction is chord/speed independent (c={_cc2}, U={_uu2}, M={_mm2})",
       _v2==_ref, f"{_v2} vs {_ref}")

# --- the published vortex-core depth must be grid-independent. It was read off
#     the field at the nearest grid node, which drifted -0.425 / -0.381 / -0.358
#     / -0.368 over 110x85 .. 880x680 for one unchanged instant. dsv_core_cp
#     evaluates it at the centre, so it must now equal the published number and
#     must not move when a grid is refined around it.
_thA=pd.read_csv('05_solution/time_history_A_validation.csv')
_rA=_thA.iloc[int(_thA.CN_vortex.idxmax())]
_tvl=json.load(open('03_model_setup/solver_config.json'))['calibrated_constants']['Tvl']
_,_,_core=us.dsv_core_cp(G,c,U,M,_rA.alpha_deg,_rA.CL,_rA.CN_vortex,
                         _rA.tau_v_semichords/_tvl)
ck("published Cp_DSV_core_min is the value at the vortex centre",
   abs(_core-float(mA['Cp_DSV_core_min']))<5e-4, f"{_core:.4f} vs {mA['Cp_DSV_core_min']}")

# --- solver_config.json must still be DESCRIBING the solver. 03_model_setup
#     imports these five from unistall_solver rather than restating them; if that
#     import is ever unwired back into literals, this fails. They agreed before
#     the import existed, but nothing made them agree.
for _k,_v in (('separation_model.f_min', us.F_MIN),
              ('field_reconstruction.n_panels', us.N_PANELS_DEFAULT),
              ('field_reconstruction.grid_nx_default', us.GRID_NX_DEFAULT),
              ('field_reconstruction.grid_ny_default', us.GRID_NY_DEFAULT),
              ('field_reconstruction.near_wall_cells_masked', us.NEAR_WALL_CELLS_MASKED)):
    _a,_b=_k.split('.'); _got=json.load(open('03_model_setup/solver_config.json'))[_a][_b]
    ck(f"config {_k} matches the solver", _got==_v, f"{_got} vs {_v}")
_dom=json.load(open('03_model_setup/solver_config.json'))['field_reconstruction']['domain_chords']
ck("config domain_chords matches the solver", tuple(_dom)==tuple(us.DOMAIN_CHORDS),
   f"{_dom} vs {list(us.DOMAIN_CHORDS)}")

# --- the cycle-wide closure metric must actually bound the instant it contains.
#     The peak-lift figure was once quoted as though it were a cycle-wide bound;
#     these two rows exist so that reading can never be made again silently.
for _case,_m in (('A_validation',mA),('B_application',mB)):
    _clm=float(_m['CL_max_dynamic'])
    _pk=abs(float(_m['Cp_closure_error_pct'])/100.0*_clm)
    ck(f"{_case} worst-cycle closure bounds the peak-lift closure",
       float(_m['Cp_closure_worst_dCL_cycle'])>=_pk-5e-4,
       f"worst {_m['Cp_closure_worst_dCL_cycle']} vs peak-lift dCL {_pk:.4f}")
    ck(f"{_case} worst-cycle closure percentage is consistent",
       abs(100.0*float(_m['Cp_closure_worst_dCL_cycle'])/_clm
           -float(_m['Cp_closure_worst_dCL_pct_of_CLmax']))<0.02)

# --- the dynamic validation must say which frames extrapolate. Two of the four
#     held-out frames peak 5 deg past the incidence the separation law is fitted
#     to, and nothing said so, in a study that flags the same boundary in the
#     response surface, the model polar, the field sampling and two figures.
_nrv=pd.read_csv('06_postprocessing/validation/validation_nasa_real.csv')
_acal2=float(pd.read_csv('03_model_setup/static_polar_reference.csv')['alpha_deg'].max())
ck("validation frames carry a calibration-range flag",
   'within_static_calibration' in _nrv.columns and 'peak_alpha_deg' in _nrv.columns)
if 'within_static_calibration' in _nrv.columns:
    ck("peak_alpha_deg equals alpha0 + amplitude",
       bool((( _nrv.alpha0_deg+_nrv.amp_deg-_nrv.peak_alpha_deg).abs()<0.051).all()))
    ck("within_static_calibration agrees with the polar range",
       bool((_nrv.within_static_calibration ==
             (_nrv.peak_alpha_deg <= _acal2+1e-9)).all()))
    ck("the flag is not vacuous (some frame does extrapolate)",
       bool((~_nrv.within_static_calibration).any()))

# --- provenance
import scipy.io
res=pd.read_csv('06_postprocessing/validation/validation_nasa_real.csv')
bad=0
for _,row in res.iterrows():
    dd=scipy.io.loadmat(f"06_postprocessing/validation/experimental/nasa_frames/{row['frame']}.mat")
    gv=lambda k: float(np.asarray(dd[k]).ravel()[0])
    if abs(np.degrees(gv('alpha_0'))-float(row['alpha0_deg']))>0.06 or abs(gv('M')-float(row['M']))>5e-4: bad+=1
ck("all 6 frame provenances match the .mat", bad==0)

# --- every generated artifact must be DECLARED in some stage's docstring manifest.
#     run_case.py's manifest listed a field-file pattern that no file has matched
#     since the phase tags were introduced, and omitted summary_all_cases.csv
#     entirely; make_3d_plots.py declared a wildcard for a file that has no
#     suffix. Nothing checked the manifests against the artifacts, so they drifted.
_docs=""
for _f in sorted(glob.glob('0*/**/*.py',recursive=True)+glob.glob('*.py')):
    _docs += (ast.get_docstring(ast.parse(open(_f,encoding='utf-8').read())) or "") + "\n"
_decl=set(re.findall(r'[A-Za-z0-9_<>\-]+\.(?:csv|png|json|pdf|docx)', _docs))
# Placeholders are bound to their real vocabularies. They used to expand to a
# greedy [A-Za-z0-9_.]+, which let <case>_<phase> split "A_validation" as
# case="A", phase="validation" -- so temperature_profile_<case>_<phase>.png
# "matched" files that have no phase at all, and the wrong pattern went unnoticed.
_CASE = r'(?:A_validation|B_application)'
_PHASE = r'(?:rise|peak|fall|dsv)'
_DEG  = r'\d+'
def _declared(name):
    for d in _decl:
        if d == name:
            return True
        if '<' in d:
            pat = (re.escape(d).replace('<case>', _CASE)
                              .replace('<phase>', _PHASE)
                              .replace('<deg>', _DEG))
            # unbound placeholders may contain underscores (<frame> is
            # "frame_10118"); the bound ones above stop this being greedy
            # enough to hide a wrong pattern.
            pat = re.sub(r'<[^>]+>', r'[A-Za-z0-9_]+', pat)
            if re.fullmatch(pat, name):
                return True
    return False
_gen=[f for f in subprocess.run(['git','ls-files'],capture_output=True,text=True).stdout.split()
      if re.match(r'0[1-8]_', f) and f.endswith(('.csv','.png'))
      and 'experimental' not in f]
_undeclared=[f for f in _gen if not _declared(os.path.basename(f))]
ck("every generated artifact is declared in a stage manifest",
   not _undeclared, f"{len(_undeclared)} undeclared, e.g. {_undeclared[:3]}")

# --- every repository path the prose points a reader at must exist. The licence
#     carve-out added in this audit names four of them (the .mat directory, its
#     PROVENANCE.txt, the exp_frame_* extracts and LICENSE); nothing checked that
#     a named path was real, and a licence notice that points at a file which is
#     not there is worse than none.
_refsrc=open('NOTICE',encoding='utf-8').read()+"\n"+open('README.md',encoding='utf-8').read()
_refs=set(re.findall(r'(?<![\w/])(0[0-8]_[A-Za-z0-9_./*{},\-]*[A-Za-z0-9_*}])', _refsrc))
_badref=[]
for _r in sorted(_refs):
    _cands=[_r]
    _m=re.search(r'\{([^}]*)\}', _r)               # {CL,CM} alternation
    if _m: _cands=[_r[:_m.start()]+_alt+_r[_m.end():] for _alt in _m.group(1).split(',')]
    if not any(glob.glob(_c) or os.path.exists(_c) for _c in _cands): _badref.append(_r)
ck("every repository path named in NOTICE/README exists", not _badref, f"missing {_badref}")

# --- every shipped top-level entry must be described in the README's structure
#     table. Making 00_overview/ and 07_report/ ship in an earlier pass left both
#     undocumented there, so the table described a repository that no longer
#     existed. Site furniture and the README itself are excluded.
_rm=open('README.md',encoding='utf-8').read()
_sec=_rm[_rm.index('## Repository structure'):]
_sec=_sec[:_sec.index('\n## ',5)] if '\n## ' in _sec[5:] else _sec
_tops=sorted({t.split('/')[0] if '/' in t else t
              for t in subprocess.run(['git','ls-files'],capture_output=True,text=True).stdout.split()})
_skip={'_config.yml','_includes','_layouts','assets','favicon.ico','README.md'}
_undoc=[t for t in _tops if not t.startswith('.') and t not in _skip
        and f'`{t}/`' not in _sec and f'`{t}`' not in _sec]
ck("README structure table covers everything that ships", not _undoc, f"missing {_undoc}")

# --- the consolidated report must be navigable, and its outline must point at the
#     right pages. It shipped with no contents page and no PDF outline at all.
try:
    import fitz as _fz
    _rp=_fz.open('aero_dynamic_stall_report.pdf'); _toc=_rp.get_toc()
    ck("report has a PDF outline", len(_toc)>=18, f"{len(_toc)} entries")
    _nums={int(t[1].split('.')[0]) for t in _toc if t[0]==1 and t[1][:1].isdigit()}
    ck("outline covers every top-level section", not (set(range(1,17))-_nums),
       f"missing {sorted(set(range(1,17))-_nums)}")
    _bad=0
    for _lvl,_ti,_pg in _toc:
        _t=_rp[_pg-1].get_text()
        _key=_ti.replace('Appendix \u2014 ','').split('\u2014')[0].strip()[:24]
        if _key not in _t and _ti.split('.')[0]+'.' not in _t: _bad+=1
    ck("every outline entry lands on its own heading", _bad==0, f"{_bad} wrong targets")
    # No table cell may break a filename mid-word. Widening the data inventory's
    # folder column in this audit squeezed the file column until three filenames
    # broke as "...peak_a19.pn" with an orphaned "g" on the next line -- no data
    # lost, but it reads as a typo in a table whose purpose is exact filenames.
    # Exact test, no heuristics: take the column names of every published CSV.
    # If a name appears in the report with the newlines stripped out but NOT on a
    # single line, the renderer broke it. Adjacent table cells cannot produce a
    # false positive this way, which a line-adjacency heuristic could not avoid.
    _pages=[_p.get_text() for _p in _rp]
    _ids=set()
    for _cf in subprocess.run(['git','ls-files'],capture_output=True,text=True).stdout.split():
        if _cf.endswith('.csv') and re.match(r'0[1-8]_',_cf):
            try: _ids.update(str(_c) for _c in pd.read_csv(_cf, nrows=0).columns)
            except Exception: pass
            # metrics_*.csv is a key/value table: its IDENTIFIERS are the row
            # keys, not the column headings ("metric", "value"), so scanning
            # headers alone never looked at the names the prose sends readers to
            # -- Cp_closure_worst_dCL_pct_of_CLmax among them.
            if os.path.basename(_cf).startswith('metrics_'):
                try: _ids.update(str(_v) for _v in pd.read_csv(_cf)['metric'])
                except Exception: pass
    _ids={i for i in _ids if len(i)>6 and re.fullmatch(r'[A-Za-z0-9_]+', i)}
    # PER PAGE. Checking the whole document at once was useless: a name broken on
    # one page but printed intact on another looked fine, and every one of these
    # names appears in the section 15 inventory as well as in its own table.
    # A split is ACCEPTABLE when it falls on a separator: "CLmax_" + "exp" is the
    # deliberate behaviour. It is a defect only when a piece ends mid-word, e.g.
    # "alpha_mean_de" + "g". So find the consecutive lines that reconstruct the
    # name and check where the joins land.
    def _split_ok(lines, ident):
        for _a in range(len(lines)):
            acc, parts = "", []
            for _b in range(_a, min(_a + 6, len(lines))):
                acc += lines[_b]; parts.append(lines[_b])
                if acc == ident:
                    return all(x.endswith(("_", "/")) for x in parts[:-1])
                if not ident.startswith(acc):
                    break
        return None                     # not reconstructible: not this page's problem
    _broken=set()
    for _pt in _pages:
        _lines=[l.strip() for l in _pt.splitlines() if l.strip()]
        _pf=re.sub(r'\s+','',_pt)
        for _i in _ids:
            if _i in _pf and _i not in _pt:
                if _split_ok(_lines, _i) is False:
                    _broken.add(_i)
    _broken=sorted(_broken)
    ck("no CSV column name is broken across lines in the report", not _broken,
       f"{len(_broken)} broken, e.g. {_broken[:4]}")
    _rp.close()
except ImportError:
    pass

# --- NO PAGE of the shipped report may overlap its own content. The house rule
#     the figures follow ("text never overlaps a figure") was only ever enforced
#     inside individual matplotlib figures; nothing checked the assembled
#     document, where the renderer lays out text, tables and images together.
#     Text on text, text on an image, anything outside the page rectangle, or a
#     body page left near-empty by a bad break -- all of it, on every page.
try:
    import fitz as _fz6
    if os.path.exists('aero_dynamic_stall_report.pdf'):
        def _ov(a, b, tol=1.5):
            return (a[0] < b[2]-tol and b[0] < a[2]-tol
                    and a[1] < b[3]-tol and b[1] < a[3]-tol)
        _rp6 = _fz6.open('aero_dynamic_stall_report.pdf')
        _nb6 = _rp6.page_count
        for _cv in ('UNISTALL_data_dossier.pdf', 'UNISTALL_plots_album.pdf'):
            if os.path.exists(_cv):
                with _fz6.open(_cv) as _c6: _nb6 -= _c6.page_count
        _bad6 = []
        for _i6, _p6 in enumerate(_rp6):
            _bl = [b for b in _p6.get_text("blocks") if b[4].strip()]
            _im = [im['bbox'] for im in _p6.get_image_info()]
            for _j in range(len(_bl)):
                for _k in range(_j+1, len(_bl)):
                    if _ov(_bl[_j][:4], _bl[_k][:4]): _bad6.append((_i6+1, 'text/text'))
            for _b in _bl:
                if any(_ov(_b[:4], _q) for _q in _im): _bad6.append((_i6+1, 'text/image'))
                if (_b[0] < -1 or _b[1] < -1
                        or _b[2] > _p6.rect.x1+1 or _b[3] > _p6.rect.y1+1):
                    _bad6.append((_i6+1, 'outside page'))
            if _i6 < _nb6 and not _im and len("".join(b[4] for b in _bl).strip()) < 40:
                _bad6.append((_i6+1, 'near-empty body page'))
        ck("no page of the report overlaps its own content or runs off it",
           not _bad6, f"{len(_bad6)} on pages {sorted({b[0] for b in _bad6})[:6]}")
        _rp6.close()
except ImportError:
    pass

# --- the rendered report must contain EVERYTHING case.docx does. The docx is
#     assembled from the solver outputs and then rendered by build_report_pdf;
#     nothing between the two may silently drop content, and the dossier showed
#     this study can lose content without any trace (a clipped page leaves
#     nothing behind to notice). Paragraphs, table cells and images, all three.
#
#     Page numbers are stripped before the comparison, and a paragraph that
#     spans a page break is joined across it: PyMuPDF puts the footer number at
#     the START of some pages' text, so a naive concatenation drops a page
#     number into the middle of such a paragraph and reports three false
#     absences.
try:
    import fitz as _fz5
    from docx import Document as _Doc
    if os.path.exists('07_report/case.docx') and os.path.exists('aero_dynamic_stall_report.pdf'):
        _rp5=_fz5.open('aero_dynamic_stall_report.pdf')
        _nb=_rp5.page_count
        for _cv in ('UNISTALL_data_dossier.pdf','UNISTALL_plots_album.pdf'):
            if os.path.exists(_cv):
                with _fz5.open(_cv) as _c5: _nb-=_c5.page_count
        _txt=[]
        for _i in range(_nb):
            _t=_rp5[_i].get_text()
            _t=re.sub(r'^\s*%d\s*\n' % (_i+1), '', _t)      # leading page number
            _t=re.sub(r'\n\s*%d\s*\n?$' % (_i+1), '\n', _t)  # trailing page number
            _txt.append(_t)
        _flatb=re.sub(r'\s+','', "".join(_txt))
        _dx=_Doc('07_report/case.docx')
        _ps=[x.text.strip() for x in _dx.paragraphs if x.text.strip()]
        _miss=[x for x in _ps if re.sub(r'\s+','',x) not in _flatb]
        ck("every case.docx paragraph reaches the rendered report", not _miss,
           f"{len(_miss)} missing, e.g. {[m[:60] for m in _miss[:2]]}")
        _cl=[" ".join(y.text for y in cc.paragraphs).strip()
             for t in _dx.tables for r in t.rows for cc in r.cells]
        _cl=[x for x in _cl if x]
        _mc=[x for x in _cl if re.sub(r'\s+','',x) not in _flatb]
        ck("every case.docx table cell reaches the rendered report", not _mc,
           f"{len(_mc)} of {len(_cl)} missing")
        _ni=sum(len(_rp5[_i].get_images()) for _i in range(_nb))
        ck("every case.docx image reaches the rendered report",
           _ni == len(_dx.inline_shapes), f"{_ni} placed vs {len(_dx.inline_shapes)} in the docx")
        _rp5.close()
except ImportError:
    pass

# --- section 15 promises "all generated files". It listed 135 of 137: the folder
#     list was hardcoded and omitted validation/experimental/, so the two files
#     the digitizer stage writes there were absent from an inventory whose title
#     says it is complete.
try:
    import fitz as _fz2
    _rp2=_fz2.open('aero_dynamic_stall_report.pdf'); _t2=_rp2.get_toc()
    _st=[pg for lv,ti,pg in _t2 if ti.startswith('15.')]
    _en=[pg for lv,ti,pg in _t2 if ti.startswith('16.')]
    if _st and _en:
        _inv="".join(_rp2[_i].get_text() for _i in range(_st[0]-1,_en[0]))
        # normalise wrapping: a filename split across two lines by the table
        # renderer is present, not missing, so newlines are removed first.
        _flat=re.sub(r'\s+','',_inv)
        _named=set(re.findall(r'[A-Za-z0-9_\-]+\.(?:csv|png|json)', _flat))
        _gen=[f for f in subprocess.run(['git','ls-files'],capture_output=True,text=True).stdout.split()
              if re.match(r'0[1-8]_',f) and f.endswith(('.csv','.png','.json'))]
        _uncov=[f for f in _gen if os.path.basename(f) not in _named
                and os.path.basename(f) not in _flat]
        ck("data inventory lists every generated file", not _uncov,
           f"{len(_uncov)} missing, e.g. {[os.path.basename(x) for x in _uncov[:3]]}")
    _rp2.close()
except ImportError:
    pass

# --- the DATA DOSSIER must print what it says it prints. Two defects of the
#     same kind were found here, and neither left any trace to notice: a
#     matplotlib figure that overruns its page is CLIPPED at the media box, so
#     the content that falls off does not appear anywhere -- not off the page,
#     not in a warning, nowhere.
#       * solver_config.json wrapped to 97 lines and 82 were rendered, so the
#         dossier printed JSON that ended mid-string and never reached
#         calibration_state or the closing brace.
#       * every table with more than 33 rows printed 33 under a caption that
#         said 38, response_surface.csv (36 rows) included -- the one table
#         max_rows was raised to 38 in order to show whole.
try:
    import fitz as _fz4
    if os.path.exists('UNISTALL_data_dossier.pdf'):
        _dos=_fz4.open('UNISTALL_data_dossier.pdf')
        _dostxt=[_p.get_text() for _p in _dos]
        _flatdos=re.sub(r'\s+','', "".join(_dostxt))
        _cfgtxt=json.load(open('03_model_setup/solver_config.json'))
        _tail=re.sub(r'\s+','', json.dumps({"calibration_state": _cfgtxt["calibration_state"]},
                                           indent=2).strip("{} \n"))
        ck("dossier prints solver_config.json to its last key", _tail in _flatdos,
           "calibration_state missing -- the config page is being clipped")
        # Counted EXACTLY, against the file the page names: take the first N
        # values of the table's key column and see how many appear on the page
        # as whole words. A baseline count would do, but a wrapped title or a
        # wrapped cell shifts it either way, and a row-count check that can pass
        # by miscounting is the thing this check exists to replace.
        _csvs={os.path.basename(_f):_f for _f in
               subprocess.run(['git','ls-files'],capture_output=True,text=True).stdout.split()
               if _f.endswith('.csv')}
        _short=[]
        for _pi,_pt in enumerate(_dostxt):
            _one=re.sub(r'\s+',' ',_pt)
            _m=re.search(r'first (\d+) of (\d+) rows', _one)
            _hit=[_n for _n in re.findall(r'[A-Za-z0-9_\-]+\.csv', _one) if _n in _csvs]
            if not _m or not _hit: continue
            _fn=_hit[0]
            _want=int(_m.group(1))
            _df=pd.read_csv(_csvs[_fn]).head(_want)
            _words={_w[4] for _w in _dos[_pi].get_text("words")}
            _key=_df.columns[0]
            _got=sum(1 for _v in _df[_key] if str(_v) in _words)
            if _got < _want: _short.append((_fn,_got,_want))
        ck("every sampled dossier table prints the rows its caption claims",
           not _short, f"{_short[:3]}")
        _dos.close()
except ImportError:
    pass

# --- a drawing sheet states a SCALE, so its paper units must be millimetres on
#     the page. They were not: the figure was 8.5 in tall while the sheet was 210
#     paper units, and savefig.bbox="tight" then grew it further, so a sheet whose
#     title block read "1:270" printed at 1:263.8.
try:
    from PIL import Image as _Img
    _sh='08_engineering_drawings/sheet1_general_arrangement_3view.png'
    if os.path.exists(_sh):
        _im=_Img.open(_sh); _dpi=_im.info.get('dpi',(150,150))[1]
        _mm=_im.size[1]/_dpi*25.4
        ck("drawing sheet is 210 mm tall, so its stated scale is true",
           abs(_mm-210.0)<0.5, f"{_mm:.2f} mm")
except ImportError:
    pass

# --- the report must not restate model constants that live in solver_config.json.
#     The indicial pair (A1,A2,b1,b2) was written into an equation as literals, a
#     second copy of numbers the solver reads from the config.
_cfg=json.load(open('03_model_setup/solver_config.json'))
_icc=_cfg['indicial_circulatory']
_bd=open('07_report/build_docx.py',encoding='utf-8').read() if os.path.exists('07_report/build_docx.py') else ''
if _bd:
    _lit=re.search(r'\(A_\{1\},A_\{2\},b_\{1\},b_\{2\}\)=\(([-0-9.]+)', _bd)
    ck("report does not hardcode the indicial constants", _lit is None,
       f"found literal {_lit.group(1) if _lit else ''}")
    # --- the report's Cp equation must NOT carry a Prandtl-Glauert factor. The
    #     solver deliberately does not apply one (the circulation it is handed
    #     already carries compressibility through beta), and section 4.9 of the
    #     same report calls applying it one of three errors -- while section 4.7
    #     printed C_p = 1/sqrt(1-M^2)[...] for two revisions.
    ck("report's Cp equation carries no Prandtl-Glauert factor",
       re.search(r'C_\{p\}\s*=\s*\\frac\{1\}\{\\sqrt\{1-M', _bd) is None)
    # --- and the impulsive time constant must be the one the solver computes.
    #     It read T_I = K_alpha c/a, which is 1/M times what solve_dynamic_stall
    #     has ever used (T_I = K_alpha c/U).
    ck("report states the impulsive time constant the solver uses",
       ('T_{I}=\\frac{K_{\\alpha}c}{U}' in _bd
        and 'T_{I}=\\frac{K_{\\alpha}c}{a}' not in _bd))
try:
    import fitz as _fz3
    _rp3=_fz3.open('aero_dynamic_stall_report.pdf'); _rp3.close()
except ImportError:
    pass

# --- dead code / imports
# The name an import binds is a NAME. It was being looked for among attribute
# names as well, and then again as a bare word anywhere in the source -- both of
# which let a real dead import hide behind an unrelated method of the same name:
# "import glob" survived in build_docx.py and build_pdfs.py because those files
# call Path.glob(), so the attribute "glob" and dozens of textual matches made
# the import look used. An attribute is never how an imported name is READ
# (a.b.c parses as Attribute(Attribute(Name('a')))), so Name is the whole test.
_dead=[]
for f in sorted(glob.glob('0*/**/*.py',recursive=True)+glob.glob('*.py')):
    src=open(f,encoding='utf-8').read(); tree=ast.parse(src)
    imp={}
    for n in ast.walk(tree):
        if isinstance(n,ast.Import):
            for al in n.names: imp[(al.asname or al.name).split('.')[0]]=1
        elif isinstance(n,ast.ImportFrom):
            for al in n.names:
                if al.name!='*': imp[al.asname or al.name]=1
    used={n.id for n in ast.walk(tree) if isinstance(n,ast.Name)}
    _dead+=[f"{f}:{k}" for k in imp if k not in used]
ck("no unused imports", not _dead, f"{len(_dead)} found: {_dead[:4]}")

# --- stale prose
# This file is EXCLUDED from its own scan. It necessarily quotes the phrases it
# forbids, so including it makes every stale-prose check match itself and fail --
# which is exactly what happened the first time this file was committed, since
# it had passed while still untracked and git ls-files could not see it.
_SELF = os.path.basename(__file__)
allsrc="".join(open(f,encoding='utf-8',errors='replace').read()
               for f in subprocess.run(['git','ls-files'],capture_output=True,text=True).stdout.split()
               if f.endswith(('.py','.md','.json')) and os.path.basename(f) != _SELF)
for pat in (r'surface_cp probes at', r'demands it be zero', r'the 0\.015c\s*used', r'clip is kept',
            # the DSV's two constants are derived now; an assignment to either
            # of the old chosen names, or the wrong Lamb-Oseen coefficient,
            # would mean the derivation had been unwired
            r'DSV_GAMMA_FACTOR\s*=\s*[0-9]', r'DSV_CORE_RADIUS_CHORDS\s*=\s*[0-9]',
            r'0\.7152\*Gv'):
    ck(f"no stale prose: {pat}", len(re.findall(pat,allsrc))==0)

print()
print(f"  ==> {len(FAIL)} FAILURES" if FAIL else "  ==> ALL CHECKS PASS")
sys.exit(1 if FAIL else 0)
