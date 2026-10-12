"""Validate the numerical solver against an independent analytic energy decay."""
from pathlib import Path
import sys, unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'models'/'python'))
from bve import run_simulation

class AnalyticVortexTest(unittest.TestCase):
    def test_viscous_energy_decay(self):
        result=run_simulation(1)
        expected=0.25*np.exp(-4*result['energy_times'])
        np.testing.assert_allclose(result['energy'],expected,rtol=1e-8,atol=0)
        x,y=np.meshgrid(result['x'],result['y'],indexing='ij')
        exact=-2*np.cos(x)*np.cos(y)*np.exp(-2*result['time'])
        np.testing.assert_allclose(result['vorticity'],exact,rtol=1e-7,atol=1e-13)
if __name__=='__main__':unittest.main()
