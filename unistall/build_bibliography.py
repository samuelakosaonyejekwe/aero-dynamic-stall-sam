# -*- coding: utf-8 -*-
"""
unistall / build_bibliography.py
----------------------------
Build the manuscript's bibliography from identifiers RESOLVED ONLINE, not from
memory. Each candidate below is only a search request: (first author, a few
title words, year). The entry written to references.bib takes its title,
authors, journal, volume, pages and year from the registry's own record
(Crossref for DOIs, the NASA Technical Reports Server for NASA and NACA
reports). A candidate that cannot be matched is listed as unresolved in
bibliography_resolution.csv and is NOT written to the bibliography. A work
that neither registry holds (a laboratory report, a journal without registry
records, a thesis, a repository) is entered in MANUAL with the page its
details were read from; the resolution file marks it "other".

Author: Akosa Samuel Onyejekwe (independent)

A Crossref match is accepted only if the first author's family name matches,
the year is within one of the requested year, and at least 60 % of the
requested title words appear in the registered title.
"""
import json, os, re, time, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
from unistall.paths import DOCS
# The registries ask for a contact address in the request header. It is read
# from the environment (UNISTALL_CONTACT) and left out if that is not set.
CONTACT = os.environ.get("UNISTALL_CONTACT", "")
UA = {"User-Agent": "unistall bibliography builder" + (f" (mailto:{CONTACT})" if CONTACT else "")}

# key, kind, first author family name, year, title words, topic
C = [
 ("leishman1989", "doi", "Leishman", 1989, "A Semi-Empirical Model for Dynamic Stall", "model"),
 ("leishman1990", "doi", "Leishman", 1990, "State-space representation of unsteady airfoil behavior", "model"),
 ("leishman1988", "doi", "Leishman", 1988, "Validation of approximate indicial aerodynamic functions for two-dimensional subsonic flow", "indicial"),
 ("leishman1993", "doi", "Leishman", 1993, "Indicial lift approximations for two-dimensional subsonic flow as obtained from oscillatory measurements", "indicial"),
 ("sheng2006", "doi", "Sheng", 2006, "A new stall-onset criterion for low speed dynamic-stall", "onset"),
 ("sheng2008", "doi", "Sheng", 2008, "A modified dynamic stall model for low Mach numbers", "lowmach"),
 ("goman1994", "doi", "Goman", 1994, "State-space representation of aerodynamic characteristics of an aircraft at high angles of attack", "model"),
 ("mccroskey1982arfm", "doi", "McCroskey", 1982, "Unsteady airfoils", "review"),
 ("carr1988", "doi", "Carr", 1988, "Progress in analysis and prediction of dynamic stall", "review"),
 ("mulleners2012", "doi", "Mulleners", 2012, "The onset of dynamic stall revisited", "onset"),
 ("mulleners2013", "doi", "Mulleners", 2013, "Dynamic stall development", "experiment"),
 ("corke2015", "doi", "Corke", 2015, "Dynamic stall in pitching airfoils: aerodynamic damping and compressibility effects", "review"),
 ("ekaterinaris1998", "doi", "Ekaterinaris", 1998, "Computational prediction of airfoil dynamic stall", "cfd"),
 ("richter2011", "doi", "Richter", 2011, "Improved two-dimensional dynamic stall prediction with structured and hybrid numerical methods", "cfd"),
 ("dossantos2021", "doi", "Santos", 2021, "Improvements on the Beddoes-Leishman dynamic stall model for low speed applications", "lowmach"),
 ("larsen2007", "doi", "Larsen", 2007, "Dynamic stall model for wind turbine airfoils", "tabulated"),
 ("gupta2006", "doi", "Gupta", 2006, "Dynamic stall modelling of the S809 aerofoil and comparison with experiments", "model"),
 ("pierce1995", "doi", "Pierce", 1995, "Prediction of wind turbine rotor loads using the Beddoes-Leishman model for dynamic stall", "model"),
 ("ericsson1988", "doi", "Ericsson", 1988, "Fluid mechanics of dynamic stall part I. Unsteady flow concepts", "review"),
 ("johnson1972", "doi", "Johnson", 1972, "On the mechanism of dynamic stall", "rotor"),
 ("ham1968", "doi", "Ham", 1968, "Dynamic stall considerations in helicopter rotors", "rotor"),
 ("bousman1998", "doi", "Bousman", 1998, "A qualitative examination of dynamic stall from flight test data", "rotor"),
 ("gardner2023", "doi", "Gardner", 2023, "Review of rotating wing dynamic stall: experiments and flow control", "review"),
 ("visbal2018", "doi", "Visbal", 2018, "Analysis of dynamic stall on a pitching airfoil using high-fidelity large-eddy simulations", "cfd"),
 ("bangga2020", "doi", "Bangga", 2020, "An improved second-order dynamic stall model for wind turbine airfoils", "model"),
 ("elgammi2016", "doi", "Elgammi", 2016, "A modified Beddoes-Leishman model for unsteady aerodynamic blade load computations on wind turbine blades", "model"),
 ("boutet2020", "doi", "Boutet", 2020, "A modified Leishman-Beddoes model for airfoil sections undergoing dynamic stall at low Reynolds numbers", "lowmach"),
 ("wagner1925", "doi", "Wagner", 1925, "Über die Entstehung des dynamischen Auftriebes von Tragflügeln", "theory"),
 ("mccroskey1976", "doi", "McCroskey", 1976, "Dynamic stall experiments on oscillating airfoils", "experiment"),
 ("mccroskey1980forum", "doi", "McCroskey", 1980, "Dynamic stall on advanced airfoil sections", "experiment"),
 ("leishman2002", "doi", "Leishman", 2002, "Challenges in modelling the unsteady aerodynamics of wind turbines", "review"),
 ("holierhoek2013", "doi", "Holierhoek", 2013, "Comparing different dynamic stall models", "comparison"),
 ("mcalister1978tp", "ntrs", "McAlister", 1978, "Dynamic stall experiments on the NACA 0012 airfoil", "experiment"),
 ("mccroskey1981tm", "ntrs", "McCroskey", 1981, "The phenomenon of dynamic stall", "review"),
 ("carr1977", "ntrs", "Carr", 1977, "Analysis of the development of dynamic stall based on oscillating airfoil experiments", "experiment"),
 ("theodorsen1935", "ntrs", "Theodorsen", 1935, "General theory of aerodynamic instability and the mechanism of flutter", "theory"),
 ("jones1940", "ntrs", "Jones", 1940, "The unsteady lift of a wing of finite aspect ratio", "theory"),
 ("lomax1952", "ntrs", "Lomax", 1952, "Two- and three-dimensional unsteady lift problems in high-speed flight", "theory"),
 ("garrick1936", "ntrs", "Garrick", 1936, "Propulsion of a flapping and oscillating airfoil", "theory"),
 ("allen1944", "ntrs", "Allen", 1944, "Wall interference in a two-dimensional-flow wind tunnel, with consideration of the effect of compressibility", "experiment"),
 ("mccroskey1982v1", "ntrsid", "19820024438", 1982, "", "data"),
 ("mcalister1982v2", "ntrsid", "19830003778", 1982, "", "data"),
 ("carr1982v3", "ntrsid", "19830009234", 1982, "", "data"),
 ("mccroskey1987", "ntrsid", "19880002254", 1987, "", "data"),
 ("bousman2003", "ntrsid", "20040081236", 2003, "", "rotor"),
]


# Entries whose identifier is not in Crossref or NTRS but was resolved on the
# registry or host named in `where`. The metadata below is what that page gives.
MANUAL = [
 dict(key="damiani2019", topic="model", identifier="doi:10.2172/1576488",
      where="https://www.osti.gov/biblio/1576488 (U.S. DOE Office of Scientific and Technical Information), read 2026-10-06",
      bib="""@techreport{damiani2019,
  author = {Damiani, Rick R. and Hayman, Gregory},
  title = {The Unsteady Aerodynamics Module for {FAST} 8},
  institution = {National Renewable Energy Laboratory},
  number = {NREL/TP-5000-66347},
  year = {2019},
  doi = {10.2172/1576488},
}
"""),
 dict(key="pancini2021data", topic="data", identifier="https://github.com/luizpancini/BL-DSM-JFS-2021",
      where="GitHub repository record (api.github.com/repos/luizpancini/BL-DSM-JFS-2021), read 2026-10-06",
      bib="""@misc{pancini2021data,
  author = {Pancini, Luiz},
  title = {{BL-DSM-JFS-2021}: source code for a modified {B}eddoes--{L}eishman dynamic stall model, with digitised experimental data},
  howpublished = {GitHub repository},
  year = {2022},
  url = {https://github.com/luizpancini/BL-DSM-JFS-2021},
  note = {Last updated 17 May 2022; no licence stated},
}
"""),
 dict(key="green2017", topic="data", identifier="doi:10.5525/gla.researchdata.464",
      where="https://researchdata.gla.ac.uk/464/ (University of Glasgow research data record), read 2026-10-07",
      bib="""@misc{green2017,
  author = {Green, R. B. and Giuni, M.},
  title = {Dynamic Stall Database {R} and {D} 1570-{AM}-01: Final Report},
  howpublished = {University of Glasgow, research data},
  year = {2017},
  doi = {10.5525/gla.researchdata.464},
}
"""),
 dict(key="leishman1989crouse", topic="model", identifier="doi:10.2514/6.1989-1319",
      where="the registry record of this DOI names the first author only; second author and paper number as cited in The Aeronautical Journal and Journal of the Brazilian Society of Mechanical Sciences and Engineering reference lists, read 2026-10-07",
      bib="""@inproceedings{leishman1989crouse,
  author = {Leishman, J. G. and Crouse, G. L.},
  title = {State-Space Model for Unsteady Airfoil Behavior and Dynamic Stall},
  booktitle = {Proceedings of the 30th AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics and Materials Conference, Mobile, AL, Paper AIAA 89-1319},
  year = {1989},
  doi = {10.2514/6.1989-1319},
}
"""),
 dict(key="hansen2004", topic="tabulated", identifier="Risø-R-1354(EN); ISBN 87-550-3090-4",
      where="https://orbit.dtu.dk/en/publications/a-beddoes-leishman-type-dynamic-stall-model-in-state-space-and-in/ (DTU Research Database record), read 2026-10-07",
      bib="""@techreport{hansen2004,
  author = {Hansen, Morten H. and Gaunaa, Mac and Madsen, Helge Aagaard},
  title = {A {B}eddoes--{L}eishman Type Dynamic Stall Model in State-Space and Indicial Formulations},
  institution = {Ris{\\o} National Laboratory},
  address = {Roskilde, Denmark},
  number = {Ris{\\o}-R-1354(EN)},
  year = {2004},
  url = {https://orbit.dtu.dk/en/publications/a-beddoes-leishman-type-dynamic-stall-model-in-state-space-and-in/},
}
"""),
 dict(key="beddoes1983", topic="model", identifier="Vertica 7(2):183-197",
      where="the journal Vertica has no registry record; volume, number and pages as cited in The Aeronautical Journal (Cambridge University Press) reference lists, read 2026-10-07",
      bib="""@article{beddoes1983,
  author = {Beddoes, T. S.},
  title = {Representation of Airfoil Behaviour},
  journal = {Vertica},
  volume = {7},
  number = {2},
  pages = {183--197},
  year = {1983},
}
"""),
 dict(key="beddoes1984", topic="indicial", identifier="Vertica 8(1):55-71",
      where="the journal Vertica has no registry record; volume, number and pages as cited in The Aeronautical Journal (Cambridge University Press) reference lists, read 2026-10-07",
      bib="""@article{beddoes1984,
  author = {Beddoes, T. S.},
  title = {Practical Computation of Unsteady Lift},
  journal = {Vertica},
  volume = {8},
  number = {1},
  pages = {55--71},
  year = {1984},
}
"""),
 dict(key="tran1981", topic="onera", identifier="Vertica 5(1):35-53",
      where="the journal Vertica has no registry record; volume and pages as cited in Aerospace 4(2):21 (doi:10.3390/aerospace4020021) and other reference lists, read 2026-10-07",
      bib="""@article{tran1981,
  author = {Tran, C. T. and Petot, D.},
  title = {Semi-Empirical Model for the Dynamic Stall of Airfoils in View of the Application to the Calculation of Responses of a Helicopter Blade in Forward Flight},
  journal = {Vertica},
  volume = {5},
  number = {1},
  pages = {35--53},
  year = {1981},
}
"""),
 dict(key="khan2018", topic="review", identifier="https://repository.tudelft.nl/record/uuid:f1ee9368-ca44-47ca-abe2-b816f64a564f",
      where="TU Delft repository record, read 2026-10-07",
      bib="""@mastersthesis{khan2018,
  author = {Khan, Muhammad Arsalan},
  title = {Dynamic Stall Modeling for Wind Turbines},
  school = {Delft University of Technology},
  year = {2018},
  url = {https://repository.tudelft.nl/record/uuid:f1ee9368-ca44-47ca-abe2-b816f64a564f},
}
"""),
]


def _get(url, data=None):
    req = urllib.request.Request(url, data=data, headers={**UA, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


def _words(s):
    return {w for w in re.findall(r"[a-zäöüß0-9]+", s.lower()) if len(w) > 2}


def crossref(author, year, title):
    q = urllib.parse.quote(f"{author} {title} {year}")
    d = _get(f"https://api.crossref.org/works?query.bibliographic={q}&rows=6")
    want = _words(title)
    for it in d["message"]["items"]:
        fam = [a.get("family", "") for a in it.get("author", [])]
        # the year of the printed issue where the registry gives one; "issued"
        # is often the earlier date of online posting
        yr = ((it.get("published-print") or it.get("issued") or {}).get("date-parts") or [[None]])[0][0]
        t = (it.get("title") or [""])[0]
        if not fam or yr is None:
            continue
        ok_author = any(author.lower() in f.lower() for f in fam[:1]) or author.lower() in fam[0].lower()
        overlap = len(want & _words(t))/max(len(want), 1)
        if ok_author and abs(yr - year) <= 1 and overlap >= 0.6:
            return dict(id=it["DOI"], title=t.title() if t.isupper() else t, year=yr,
                        authors=" and ".join(f"{a.get('family','')}, {a.get('given','')}".strip(", ")
                                             for a in it.get("author", [])),
                        journal=(it.get("container-title") or [""])[0], volume=it.get("volume", ""),
                        number=it.get("issue", ""), pages=it.get("page") or it.get("article-number", ""),
                        publisher=it.get("publisher", ""), type=it.get("type", ""))
    return None


def _surname_case(name):
    """The registry stores 'Mccroskey, W. J.'; restore the capital after Mc."""
    return re.sub(r"\bMc([a-z])", lambda m: "Mc" + m.group(1).upper(), name)


def ntrs_record(r):
    yr = int((r.get("publications") or [{}])[0].get("publicationDate", "0000")[:4] or 0)
    nums = [n for n in r.get("otherReportNumbers", []) if not n.startswith("Report Number")]
    return dict(id=str(r["id"]), title=r["title"].strip(), year=yr,
                authors=" and ".join(_surname_case(a["meta"]["author"]["name"]) for a in r.get("authorAffiliations", [])),
                journal="", volume="", number="; ".join(nums[:3]), pages="", publisher="NASA", type="report")


def ntrs_search(author, year, title):
    body = json.dumps({"q": title, "page": {"size": 8}}).encode()
    d = _get("https://ntrs.nasa.gov/api/citations/search", body)
    want = _words(title)
    for r in d.get("results", []):
        rec = ntrs_record(r)
        if author.lower() in rec["authors"].lower() and len(want & _words(rec["title"]))/max(len(want), 1) >= 0.7:
            return rec                      # NTRS dates are often re-issue dates, so the year is not tested
    return None


def bibtex(key, kind, rec, year_requested):
    def f(k, v):
        return f"  {k} = {{{v}}},\n" if v else ""
    if kind == "doi":
        typ = {"journal-article": "article", "proceedings-article": "inproceedings"}.get(
            rec["type"], "techreport" if "report" in rec["type"] else "misc")
        venue = {"article": "journal", "inproceedings": "booktitle"}.get(typ, "howpublished")
        return (f"@{typ}{{{key},\n" + f("author", rec["authors"]) + f("title", rec["title"])
                + f(venue, rec["journal"] or rec["publisher"])
                + f("volume", rec["volume"]) + f("number", rec["number"]) + f("pages", rec["pages"].replace("-", "--"))
                + f("year", rec["year"]) + f("doi", rec["id"]) + "}\n")
    return (f"@techreport{{{key},\n" + f("author", rec["authors"]) + f("title", rec["title"])
            + f("institution", "National Aeronautics and Space Administration" if year_requested >= 1958
                else "National Advisory Committee for Aeronautics")
            + f("number", rec["number"]) + f("year", year_requested)
            + f("note", f"NASA Technical Reports Server document {rec['id']}")
            + f("url", f"https://ntrs.nasa.gov/citations/{rec['id']}") + "}\n")


def _previous():
    """Entries and resolution rows written by the last successful run, by key,
    so that a registry that does not answer today does not remove a reference
    it resolved before."""
    bibfile, resfile = DOCS/"references.bib", DOCS/"bibliography_resolution.csv"
    if not (bibfile.exists() and resfile.exists()):
        return {}, {}
    entries = {m.group(1): m.group(0).strip() + "\n" for m in
               re.finditer(r"@\w+\{([^,]+),\n.*?\n\}", bibfile.read_text(encoding="utf-8"), re.S)}
    res = pd.read_csv(resfile).fillna("")
    return entries, {r.key: r._asdict() for r in res.itertuples(index=False) if r.resolved}


def _lookup(kind: str, a: str, year: int, title: str) -> tuple:
    """(record or None, note) for one candidate: the registry's record, or
    None with the reason it could not be had."""
    try:
        if kind == "doi":
            rec = crossref(a, year, title)
        elif kind == "ntrs":
            rec = ntrs_search(a, year, title)
        else:
            rec = ntrs_record(_get(f"https://ntrs.nasa.gov/api/citations/{a}"))
    except Exception as e:                                    # network or registry failure
        return None, f"lookup failed: {e}"
    return rec, "" if rec is not None else "no acceptable match"


def main():
    rows, bib = [], []
    old_entries, old_rows = _previous()
    for key, kind, a, year, title, topic in C:
        rec, note = _lookup(kind, a, year, title)
        time.sleep(0.25)
        if rec is None and key in old_entries and key in old_rows:
            rows.append(dict(old_rows[key]))                  # resolved by an earlier run: kept
            bib.append(old_entries[key])
        elif rec is None:
            rows.append(dict(key=key, topic=topic, resolved=False, identifier="", registry=kind,
                             title_requested=title, title_registered="", note=note))
        else:
            ident = ("doi:" + rec["id"]) if kind == "doi" else ("NTRS " + rec["id"])
            rows.append(dict(key=key, topic=topic, resolved=True, identifier=ident, registry=kind,
                             title_requested=title, title_registered=rec["title"], note=""))
            bib.append(bibtex(key, "doi" if kind == "doi" else "ntrs", rec, year))
    for m in MANUAL:
        rows.append(dict(key=m["key"], topic=m["topic"], resolved=True, identifier=m["identifier"],
                         registry="other", title_requested="", title_registered="", note=m["where"]))
        bib.append(m["bib"])
    df = pd.DataFrame(rows)
    df.to_csv(DOCS/"bibliography_resolution.csv", index=False)
    (DOCS/"references.bib").write_text("% Generated by unistall/build_bibliography.py from registry records. Do not edit by hand.\n\n"
                                      + "\n".join(bib), encoding="utf-8")
    print(f"[bibliography] {int(df.resolved.sum())} of {len(df)} candidates resolved; "
          f"unresolved: {list(df[~df.resolved].key)}")


if __name__ == "__main__":
    main()
