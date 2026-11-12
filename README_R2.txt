Physics manuscript physics-4642818: data and code for revision 2.
All members from the previous data.zip are preserved byte for byte.
Added: simulation and analysis code in src/, completed results/S3_past/,
config/S3_past.json, validation logs, residual_drift.csv, macro provenance,
and the script that reads these saved results to generate the revised S3.
No simulations were run during this editorial revision.

Python dependencies: NumPy, Numba, Matplotlib, Pillow. Scripts also use the
standard library. Keep src/, config/, results/, and qa/ under one root.
Do not launch run_*.py to inspect existing results: those scripts launch
simulations. model.py implements the main model, model_s3.py proportional
selection, and model_s3_past.py permits a=0 for the auxiliary comparison.
The first manuscript's numerical macros come from src/build_paper_assets.py.
It targets artigo/ and should only be run in a separate extracted workspace.

To regenerate only the added macros, S3, and drift table from saved data:
  python "submissão/03 Physics R2/revision_work/build_revision_assets.py"
This writes inside submissão/03 Physics R2/ and does not run simulations.
The accompanying manuscript_sources.zip contains the complete clean article.
Extract that ZIP into submissão/03 Physics R2/sources/ if compiling there.

In residual_drift.csv, drifting_realizations counts runs whose absolute
half_difference exceeds 0.02. The table includes points above 10 percent.
source identifies the selected raw result file. PP(0,m) data use w=1,
b=1.15, epsilon=0, eta=1e-6, L=100, and 30 independent runs per length.
They do not provide an a=0 comparison at the main value w=0.5.
