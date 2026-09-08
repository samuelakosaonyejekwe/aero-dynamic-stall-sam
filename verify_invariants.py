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
