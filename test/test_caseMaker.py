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

def test_casemaker_sphere_rcs(tmp_path):
    os.chdir(tmp_path)
    cm = CaseMaker()

    cm.setNumberOfTimeSteps(1)
    cm.setTimeStep(1e-12)

    cm.setAllBoundaries("mur")

    cm.setGridFromVTK(GEOMETRIES_FOLDER + "sphere.grid.vtp")
    sphereId = cm.addCellElementsFromVTK(
        GEOMETRIES_FOLDER + "buggy_sphere.str.vtp")

    pecId = cm.addPECMaterial()
    cm.addMaterialAssociation(pecId, [sphereId])

    planewaveBoxId = cm.addCellElementBox(
        [[-75.0, -75.0, -75.0], [75.0, 75.0, 75.0]])
    direction = {"theta": np.pi/2, "phi": 0.0}
    polarization = {"theta": np.pi/2, "phi": np.pi/2}
    dt = 1e-12
    w0 = 0.1e-9    # ~ 2 GHz bandwidth
    t0 = 10*w0
    t = np.arange(0, t0+20*w0, dt)
    data = np.empty((len(t), 2))
    data[:,0] = t
    data[:,1] = np.exp( -np.power(t-t0,2)/ w0**2 )
    np.savetxt('gauss.exc', data)
    cm.addPlanewaveSource(planewaveBoxId, 'gauss.exc', direction, polarization)

    pointProbeNodeId = cm.addNodeElement([-65.0, 0.0, 0.0])
    cm.addPointProbe(pointProbeNodeId, name="front")

    n2ffBoxId = cm.addCellElementBox(
        [[-85.0, -85.0, -85.0], [85.0, 85.0, 85.0]])
    theta = {"initial": np.pi/2, "final": np.pi/2, "step": 0.0}
    phi = {"initial": np.pi, "final": np.pi, "step": 0.0}
    domain = {
        "type": "frequency",
        "initialFrequency": 10e6,
        "finalFrequency": 1e9,
        "numberOfFrequencies": 10,
        "frequencySpacing": "logarithmic"
    }
    cm.addFarFieldProbe(n2ffBoxId, "n2ff", theta, phi, domain)

    case_name = 'sphere_rcs'
    cm.exportCase(case_name)
    
    with open(case_name+'.fdtd.json', 'r') as f:
        a = json.load(f)
    with open(OUTPUTS_FOLDER+'sphere_rcs.fdtd.json') as f:
        b = json.load(f)
        
    assert  (sorted(a.items()) == sorted(b.items())) == True


def test_casemaker_material_fragment(tmp_path):
    grid = {
        "numberOfCells": [4, 4, 2],
        "origin": [0.0, 0.0, 0.0],
        "steps": {"x": [0.1], "y": [0.1], "z": [0.1]},
    }
    fragment = {
        "material": {
            "id": 1,
            "name": "plasma",
            "type": "lorentz",
            "epsInf": 1.0,
            "omega1": 0.0,
            "levels": [
                {"omegaP": 1.0e6, "gamma": 1.0e5},
                {"omegaP": 2.0e6, "gamma": 2.0e5},
            ],
        },
        "mesh": {
            "grid": grid,
            "elements": [
                {"id": 1, "name": "fp_0_nu_0", "type": "cell",
                 "intervals": [[[0, 0, 0], [4, 1, 1]]]},
                {"id": 2, "name": "fp_1_nu_1", "type": "cell",
                 "intervals": [[[0, 1, 0], [4, 2, 1]]]},
            ],
        },
        "materialAssociations": [
            {"materialId": 1, "level": 0, "elementIds": [1]},
            {"materialId": 1, "level": 1, "elementIds": [2]},
        ],
    }
    fragment_path = tmp_path / "fragment.json"
    fragment_path.write_text(json.dumps(fragment))

    cm = CaseMaker()
    cm.input["mesh"] = {"grid": grid, "elements": []}
    pec_id = cm.addPECMaterial()
    box_id = cm.addCellElementIntervals([[[0, 0, 0], [1, 1, 1]]], name="box")
    cm.addMaterialAssociation(pec_id, [box_id])

    fragment_id = cm.addMaterialFragment(str(fragment_path))

    assert fragment_id == 2
    assert cm.input["materials"][1]["id"] == 2
    assert len(cm.input["materials"][1]["levels"]) == 2
    assert [element["id"] for element in cm.input["mesh"]["elements"]] == [1, 2, 3]
    associations = cm.input["materialAssociations"]
    assert associations[1]["materialId"] == 2
    assert associations[1]["level"] == 0
    assert associations[1]["elementIds"] == [2]
    assert associations[2]["level"] == 1
    assert associations[2]["elementIds"] == [3]

    bad = dict(fragment)
    bad_mesh = dict(fragment["mesh"])
    bad_grid = dict(grid)
    bad_grid["numberOfCells"] = [4, 4, 3]
    bad_mesh["grid"] = bad_grid
    bad["mesh"] = bad_mesh
    bad_path = tmp_path / "bad_fragment.json"
    bad_path.write_text(json.dumps(bad))
    try:
        cm.addMaterialFragment(str(bad_path))
    except ValueError:
        pass
    else:
        raise AssertionError("a mismatched fragment grid must be rejected")
