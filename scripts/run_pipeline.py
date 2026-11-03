"""Run the complete documented workflow, stopping at the first failed step."""
from pathlib import Path
import subprocess,sys,os
ROOT=Path(__file__).resolve().parent.parent
STEPS=['prepare_data.py','test_analysis.py','analyse_graz.py','make_figures.py','additional_analysis.py','native_knot_analysis.py','uncertainty_analysis.py','make_supplemental_figures.py','qc_sensitivity.py']
def main():
    env=os.environ.copy();env['PYTHONUTF8']='1';env['MPLBACKEND']='Agg'
    for name in STEPS:
        print('RUN '+name,flush=True)
        subprocess.run([sys.executable,str(ROOT/'scripts'/name)],cwd=ROOT,env=env,check=True)
    print('PIPELINE COMPLETE',flush=True)
if __name__=='__main__':main()
