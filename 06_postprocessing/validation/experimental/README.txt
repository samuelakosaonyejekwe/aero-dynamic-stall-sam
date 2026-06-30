Digitised-experiment validation harness
=======================================
1. Open the target figure in McAlister NASA TP-1100 or McCroskey NASA
   TM-84245 (e.g. a C_L vs alpha or C_M vs alpha dynamic-stall loop).
2. Digitise it (e.g. WebPlotDigitizer https://automeris.io) into a CSV
   with columns: alpha_deg, CL [, CM] [, CD] [, stroke=up/down].
3. Add a matching row to conditions.csv (file, mean, amp, k, M, c, U,
   source) describing the exact test point.
4. Re-run:  python3 validate_digitized.py
   -> writes validation_digitized_<file>.csv + overlay figure with
      true point-by-point RMS / peak / loop-area errors.
No data are fabricated; results appear only for files you provide.
