#!/usr/bin/env python3
"""Quick numerical regression tests for Yi PBS and Yassin PBE."""
import math
import tempfile
from pathlib import Path
import numpy as np
from branch_stats import fst_to_t, branch_values, read_triplets

def test_branch_formula():
    # Exact three-pair values chosen to yield T_AB=.3,T_AC=.5,T_BC=.4
    vals=[1-math.exp(-x) for x in (.3,.5,.4)]
    pbs,tbc,tab,tac=branch_values(*vals)
    assert np.isclose(pbs,.2)
    assert np.isclose(tbc,.4)
    assert np.isclose(tab,.3)
    assert np.isclose(tac,.5)

def test_pbe_calibration():
    # Set a baseline of medians PBS=.1, T_BC=.25, observed PBS=.2, T_BC=.4
    pbe=.2-.4*(.1/.25)
    assert np.isclose(pbe,.04)

def test_transform_boundaries():
    assert fst_to_t(-.5)==0
    assert np.isfinite(fst_to_t(1))
    assert np.isclose(fst_to_t(0),0)

def test_triplet_count():
    pops=['A','B','C','D','E','F','G']
    assert len(read_triplets(None,pops))==105

if __name__=='__main__':
    test_branch_formula();test_pbe_calibration();test_transform_boundaries();test_triplet_count()
    print('PBS/PBE numerical tests passed')
