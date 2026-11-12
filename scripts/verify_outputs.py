"""Compare regenerated outputs with the distributed result archives.

python scripts/verify_outputs.py --rerun  # complete pipeline, then numeric checks
python scripts/verify_outputs.py          # check already regenerated output folders
"""
from pathlib import Path
from zipfile import ZipFile
import argparse,io,json,hashlib,subprocess,sys
import pandas as pd,numpy as np
R=Path(__file__).resolve().parent.parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rerun',action='store_true');args=ap.parse_args()
    if args.rerun:subprocess.run([sys.executable,str(R/'scripts/run_pipeline.py')],cwd=R,check=True)
    manifest=json.loads((R/'verification/SHA256_manifest.json').read_text());hashes=[]
    for rel,expected in manifest.items():
        actual=hashlib.sha256((R/rel).read_bytes()).hexdigest()
        if actual!=expected:raise AssertionError('Release file hash differs: '+rel)
        hashes.append(rel)
    checked=[]
    for archive in ['analysis_results.zip','supplemental_results.zip','revision_results.zip']:
        with ZipFile(R/'results/archives'/archive) as z:
            for rel in z.namelist():
                if not rel.endswith(('.csv','.csv.gz')):continue
                compression='gzip' if rel.endswith('.gz') else None
                expected=pd.read_csv(io.BytesIO(z.read(rel)),compression=compression)
                actual=pd.read_csv(R/rel,compression=compression)
                pd.testing.assert_frame_equal(expected,actual,check_exact=False,atol=1e-10,rtol=1e-10)
                checked.append({'path':rel,'rows':len(actual)})
    report={'pipeline_rerun':args.rerun,'tables_matched':checked,'release_hashes_matched':hashes,'rtol':1e-10,'atol':1e-10,'scope':'Numerical reproducibility, not independent validation of field measurements.'}
    (R/'verification/verification_report.json').write_text(json.dumps(report,indent=2))
    print('VERIFIED',len(checked),'tables;',len(hashes),'release file hashes',flush=True)
if __name__=='__main__':main()
