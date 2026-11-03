"""Small numerical and design checks; no field-data truth is assumed."""
import unittest
import numpy as np,pandas as pd
import analyse_graz as b
from qc_sensitivity import nearest_native

class AnalysisChecks(unittest.TestCase):
    def test_existing_integrals(self):b.verify_calculations()
    def test_subminute_interpolation(self):
        idx=pd.to_datetime(['2010-01-01T00:00:30Z','2010-01-01T00:01:30Z'])
        vals,exact=b.at_minutes(pd.Series([0.,60.],index=idx),pd.date_range('2010-01-01',periods=3,freq='min',tz='UTC'),3)
        self.assertEqual(vals[1],30.);self.assertFalse(exact.any())
    def test_schedules_penalties_and_counts(self):
        n=182;grid=pd.date_range('2010-01-01',periods=n+121,freq='min',tz='UTC')
        e={'event_id':'test','duration_min':n,'start':grid[60]}
        d={'Q':np.linspace(1.,2.,len(grid)),'COD':np.linspace(2.,5.,len(grid)),'TSS':np.linspace(1.,3.,len(grid))}
        rows=b.run_event(e,grid,d);self.assertEqual(len(rows),3456)
        self.assertTrue(all(r[10]==0 for r in rows if r[3]=='none'))
        for peak in [0,1,30,n-1,n]:
            left=int(np.clip(peak-30,0,n-60))
            self.assertEqual(len(np.arange(n+1)[(np.arange(n+1)>=left)&(np.arange(n+1)<left+60)]),60)
    def test_flow_regression(self):
        np.testing.assert_allclose(b.flow_predict(np.ones(4),np.array([1.,2.,3.,4.]),np.arange(5)),2.5)
        pred=b.flow_predict(np.array([1.,2.,3.]),np.array([3.,2.,1.]),np.array([0.,10.]))
        np.testing.assert_allclose(pred,[4.,0.],atol=1e-12)
    def test_native_nearest(self):
        np.testing.assert_array_equal(nearest_native(np.array([1,2,4,9]),np.array([0,3,6])),[0,3])
    def test_native_quadrature(self):
        t=np.array([0.,1.,4.,6.]);q=2+t;c=4+3*t
        dense=np.arange(7.)
        self.assertAlmostEqual(b.mass_kg(q,c,dt=np.diff(t)),b.mass_kg(2+dense,4+3*dense,dt=1),places=12)

if __name__=='__main__':unittest.main(verbosity=2)
