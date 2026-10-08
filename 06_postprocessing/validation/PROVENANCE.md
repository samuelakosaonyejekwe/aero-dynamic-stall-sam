# Provenance of the measured data

The measured loops are third-party data. They are fetched by `unistall/fetch_frames.py` from the repository
named in that script at the fixed commit `84e7945c5b16c6da4075c7fb467c01125dace36b`, each file is checked against the
SHA-256 recorded in `data/frames_inventory.csv`, and they are kept in a local cache
(`unistall/cache/frames`, 149 frame files at the time of writing). They are not redistributed in
this repository, and none is copied into this folder: the `loop_<frame>_*.csv` tables here hold the measured
points of 4 loops only, beside the model's values.

The measurements themselves are those of the oscillating-aerofoil experiment reported in
McCroskey, McAlister, Carr & Pucci (1982), NASA TM-84245 Vol. 1; McAlister, Pucci, McCroskey & Carr (1982), NASA TM-84245 Vol. 2 (entries `mccroskey1982v1` and `mcalister1982v2` of `docs/references.bib`).

The static points are digitised or transcribed from that report; each row of `data/static_naca0012_M030_*.csv`
names the figure or table and the page it comes from.
