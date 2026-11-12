"""Regenerate supplementary outputs from the saved simulation data."""
from pathlib import Path
import runpy

if __name__ == '__main__':
    root=Path(__file__).resolve().parent
    script=root/'submissão/03 Physics R2/revision_work/build_revision_assets.py'
    runpy.run_path(str(script),run_name='__main__')
