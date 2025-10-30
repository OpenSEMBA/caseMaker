from utils import *
from src.caseMaker import *

def test_casemaker_setGridFromVTK():
    cm = CaseMaker()
    cm.setGridFromVTK(GEOMETRIES_FOLDER + 'sphere.grid.vtp')

    grid = cm['mesh']['grid']
    
    steps = grid['steps']
    assert np.allclose(np.array(steps['x']), 2.0)
    assert np.allclose(np.array(steps['y']), 2.0)
    assert np.allclose(np.array(steps['z']), 2.0)
    assert np.array(steps['x']).size == 100
    assert np.array(steps['y']).size == 100
    assert np.array(steps['z']).size == 100

    assert grid['numberOfCells'] == [100, 100, 100]
    assert grid['origin'] == [-100.0, -100.0, -100.0]
