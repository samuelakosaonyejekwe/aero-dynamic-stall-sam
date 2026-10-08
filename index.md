---
layout: default
title: "Prediction of Dynamic Stall on a Pitching NACA 0012 Section"
---

# Prediction of Dynamic Stall on a Pitching NACA 0012 Section

*A case study at a tunnel condition and at a rotor-blade-station condition*

Akosa Samuel Onyejekwe (independent)

This repository holds UNISTALL (Python package `unistall`), a Leishman-Beddoes dynamic-stall load model
for a pitching NACA 0012 section, with a tabulated separation point. The model
is calibrated on part of the oscillating-aerofoil measurements of NASA
TM-84245 and assessed against the rest. A case study applies it at two
conditions:

- **Case A, tunnel condition:** a condition of the NASA TM-84245 experiment,
  so the result can be set beside a measured loop. That loop is a calibration
  loop, so the comparison is an illustration and not an independent test.
- **Case B, rotor-blade-station condition:** an aerofoil at the condition of a
  rotor-blade station, as an illustration. It is not a rotor calculation: the
  stream is steady and the pitch is prescribed.

The conditions of the two cases are given in
[00_overview/case_definition.md](00_overview/case_definition.md).

The agreement with measurement, and where the model falls short, are reported
in the tables and figures linked below.

## Start here

- [README](README.md): what the model is, how to run it, and what it does and does not show
- [Formulation](docs/formulation.md): the equations and constants
- [Tables of results](results/tables.md)
- [Figures]({{ site.repo_url }}/tree/main/results/figures)

## Case study folders

| Folder | Contents |
| --- | --- |
| [00_overview]({{ site.repo_url }}/tree/main/00_overview) | Definition of the two cases |
| [01_geometry]({{ site.repo_url }}/tree/main/01_geometry) | Section geometry |
| [02_mesh]({{ site.repo_url }}/tree/main/02_mesh) | O-grid round the section, an illustration at whose nodes the reconstruction is evaluated; the load model does not use a mesh and no flow equation is solved on it |
| [03_model_setup]({{ site.repo_url }}/tree/main/03_model_setup) | Flow conditions, kinematics and model settings |
| [04_solver]({{ site.repo_url }}/tree/main/04_solver) | Running the load model for the two cases |
| [05_solution]({{ site.repo_url }}/tree/main/05_solution) | Loads predicted by the model; flow fields that are a reconstruction drawn round the predicted lift, not a flow solution |
| [06_postprocessing]({{ site.repo_url }}/tree/main/06_postprocessing) | Load loops (Case A with its calibration loop) and maps of the reconstruction |
| [07_report]({{ site.repo_url }}/tree/main/07_report) | Report |
| [08_engineering_drawings]({{ site.repo_url }}/tree/main/08_engineering_drawings) | Illustrative drawings; every printed dimension is listed with its source |

## Sources of the measurements and parameters

- McCroskey, McAlister, Carr & Pucci (1982), NASA TM-84245 Vol. 1; McAlister, Pucci, McCroskey & Carr (1982), NASA TM-84245 Vol. 2: the
  oscillating-aerofoil experiment.
- Bousman, NASA/TP-2003-212265: the main-rotor parameters used to set Case B.
