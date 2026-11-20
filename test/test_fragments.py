import json
import os
import subprocess

import numpy as np
import pytest

from utils import *
from src.caseMaker import *


def build_grid(nx=4, ny=4, nz=2, step=0.1, origin=(0.0, 0.0, 0.0),
               single_step=False):
    def steps(count):
        return [step] if single_step else [step] * count

    return {
        "numberOfCells": [nx, ny, nz],
        "origin": list(origin),
        "steps": {
            "x": steps(nx),
            "y": steps(ny),
            "z": steps(nz),
        },
    }


def build_fragment(grid=None, levels=None, elements=None, associations=None):
    if grid is None:
        grid = build_grid()
    if levels is None:
        levels = [
            {"omegaP": 1.0e6, "gamma": 1.0e5,
             "f_p_level": 0, "nu_c_level": 0},
            {"omegaP": 2.0e6, "gamma": 2.0e5,
             "f_p_level": 1, "nu_c_level": 1},
        ]
    if elements is None:
        elements = [
            {"id": 1, "name": "fp_0_nu_0", "type": "cell",
             "intervals": [[[0, 0, 0], [4, 1, 1]]]},
            {"id": 2, "name": "fp_1_nu_1", "type": "cell",
             "intervals": [[[0, 1, 0], [4, 2, 1]]]},
        ]
    if associations is None:
        associations = [
            {"materialId": 1, "level": 0, "elementIds": [1]},
            {"materialId": 1, "level": 1, "elementIds": [2]},
        ]
    return {
        "material": {
            "id": 1,
            "name": "plasma",
            "type": "lorentz",
            "epsInf": 1.0,
            "omega1": 0.0,
            "levels": levels,
        },
        "mesh": {
            "grid": grid,
            "elements": elements,
        },
        "materialAssociations": associations,
    }


def write_fragment(path, fragment):
    path.write_text(json.dumps(fragment))
    return path


def build_case_with_grid(grid):
    cm = CaseMaker()
    cm.input["mesh"] = {"grid": grid, "elements": []}
    return cm


def assert_fragment_merged(cm, material_id=1):
    assert cm['materials'][material_id - 1]['type'] == 'lorentz'
    assert cm['materials'][material_id - 1]['name'] == 'plasma'
    assert len(cm['materials'][material_id - 1]['levels']) == 2
    assert len(cm['mesh']['elements']) == 2
    assert len(cm['materialAssociations']) == 2
    assert [a['level'] for a in cm['materialAssociations']] == [0, 1]


def test_casemaker_fragment_grid_roundtrip(tmp_path):
    grid = build_grid()
    grid_path = tmp_path / "grid.json"
    grid_path.write_text(json.dumps(grid))

    cm = CaseMaker()
    cm.setGridFromJSON(grid_path)
    assert cm['mesh']['grid'] == grid

    output = tmp_path / "written_grid.json"
    cm.writeGridJSON(output)
    assert json.loads(output.read_text()) == grid


def test_casemaker_fragment_grid_accepts_single_and_expanded_steps(tmp_path):
    full_grid = build_grid()
    single_grid = build_grid(single_step=True)

    fragment = build_fragment(grid=single_grid)
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)
    cm = build_case_with_grid(full_grid)
    cm.addMaterialFragment(fragment_path)

    fragment = build_fragment(grid=full_grid)
    fragment_path = write_fragment(tmp_path / "fragment2.json", fragment)
    cm = build_case_with_grid(single_grid)
    cm.addMaterialFragment(fragment_path)


def test_casemaker_fragment_grid_rejects_number_of_cells(tmp_path):
    fragment = build_fragment(grid=build_grid(nx=5))
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)
    cm = build_case_with_grid(build_grid())

    with pytest.raises(ValueError):
        cm.addMaterialFragment(fragment_path)


def test_casemaker_fragment_grid_rejects_origin(tmp_path):
    fragment = build_fragment(grid=build_grid(origin=(0.01, 0.0, 0.0)))
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)
    cm = build_case_with_grid(build_grid())

    with pytest.raises(ValueError):
        cm.addMaterialFragment(fragment_path)


def test_casemaker_fragment_grid_rejects_steps(tmp_path):
    fragment = build_fragment(grid=build_grid(step=0.2))
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)
    cm = build_case_with_grid(build_grid())

    with pytest.raises(ValueError):
        cm.addMaterialFragment(fragment_path)


def test_casemaker_fragment_grid_tolerance(tmp_path):
    case_grid = build_grid()

    tiny = build_grid(origin=(5.0e-13, 0.0, 0.0), step=0.1 + 5.0e-13)
    fragment = build_fragment(grid=tiny)
    fragment_path = write_fragment(tmp_path / "tiny.json", fragment)
    build_case_with_grid(case_grid).addMaterialFragment(fragment_path)

    large = build_grid(origin=(1.0e-9, 0.0, 0.0))
    fragment = build_fragment(grid=large)
    fragment_path = write_fragment(tmp_path / "large.json", fragment)
    with pytest.raises(ValueError):
        build_case_with_grid(case_grid).addMaterialFragment(fragment_path)


def test_casemaker_fragment_adopts_grid_when_case_has_none(tmp_path):
    fragment = build_fragment()
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    cm = CaseMaker()
    cm.addMaterialFragment(fragment_path)

    assert cm['mesh']['grid'] == fragment['mesh']['grid']
    assert_fragment_merged(cm)


def test_casemaker_fragment_remaps_ids(tmp_path):
    fragment = build_fragment()
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    cm = build_case_with_grid(build_grid())
    first = cm.addPECMaterial()
    second = cm.addPECMaterial()
    cm.addCellElementIntervals([[[0, 0, 0], [1, 1, 1]]], name="first")
    cm.addCellElementIntervals([[[1, 1, 1], [2, 2, 2]]], name="second")
    cm.addMaterialAssociation(first, [1])
    cm.addMaterialAssociation(second, [2])

    fragment_id = cm.addMaterialFragment(fragment_path)

    assert fragment_id == 3
    assert [m['id'] for m in cm['materials']] == [1, 2, 3]
    assert [e['id'] for e in cm['mesh']['elements']] == [1, 2, 3, 4]
    fragment_associations = cm['materialAssociations'][2:]
    assert [a['materialId'] for a in fragment_associations] == [3, 3]
    assert [a['level'] for a in fragment_associations] == [0, 1]
    assert [a['elementIds'] for a in fragment_associations] == [[3], [4]]


def test_casemaker_fragment_preserves_levels(tmp_path):
    fragment = build_fragment()
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)
    before = json.loads(fragment_path.read_text())

    cm = build_case_with_grid(build_grid())
    cm.addMaterialFragment(fragment_path)

    merged = cm['materials'][0]
    assert merged['id'] == 1
    assert merged['name'] == 'plasma'
    assert merged['type'] == 'lorentz'
    assert merged['epsInf'] == 1.0
    assert merged['omega1'] == 0.0
    assert merged['levels'] == before['material']['levels']
    assert cm['mesh']['elements'][0]['name'] == 'fp_0_nu_0'
    assert cm['materialAssociations'][0]['level'] == 0
    assert cm['materialAssociations'][1]['level'] == 1


def test_casemaker_fragment_supports_multiple_elements(tmp_path):
    fragment = build_fragment(
        associations=[
            {"materialId": 1, "level": 0, "elementIds": [1, 2]},
        ],
    )
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    cm = build_case_with_grid(build_grid())
    cm.addMaterialFragment(fragment_path)
    assert cm['materialAssociations'][0]['elementIds'] == [1, 2]


def test_casemaker_fragment_does_not_mutate_input(tmp_path):
    fragment = build_fragment()
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)
    before = fragment_path.read_text()

    build_case_with_grid(build_grid()).addMaterialFragment(fragment_path)
    assert fragment_path.read_text() == before


def test_casemaker_fragment_rejects_missing_material(tmp_path):
    fragment = build_fragment()
    del fragment['material']
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    with pytest.raises(ValueError):
        build_case_with_grid(build_grid()).addMaterialFragment(fragment_path)


def test_casemaker_fragment_rejects_missing_mesh(tmp_path):
    fragment = build_fragment()
    del fragment['mesh']
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    with pytest.raises(ValueError):
        build_case_with_grid(build_grid()).addMaterialFragment(fragment_path)


def test_casemaker_fragment_rejects_missing_associations(tmp_path):
    fragment = build_fragment()
    del fragment['materialAssociations']
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    with pytest.raises(ValueError):
        build_case_with_grid(build_grid()).addMaterialFragment(fragment_path)


def test_casemaker_fragment_rejects_missing_grid(tmp_path):
    fragment = build_fragment()
    del fragment['mesh']['grid']
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    with pytest.raises(ValueError):
        build_case_with_grid(build_grid()).addMaterialFragment(fragment_path)


def test_casemaker_fragment_rejects_missing_elements(tmp_path):
    fragment = build_fragment()
    del fragment['mesh']['elements']
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    with pytest.raises(ValueError):
        build_case_with_grid(build_grid()).addMaterialFragment(fragment_path)


def test_casemaker_fragment_rejects_unknown_material(tmp_path):
    fragment = build_fragment()
    fragment['materialAssociations'][0]['materialId'] = 2
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    with pytest.raises(ValueError):
        build_case_with_grid(build_grid()).addMaterialFragment(fragment_path)


def test_casemaker_fragment_rejects_unknown_element(tmp_path):
    fragment = build_fragment()
    fragment['materialAssociations'][0]['elementIds'] = [99]
    fragment_path = write_fragment(tmp_path / "fragment.json", fragment)

    with pytest.raises(ValueError):
        build_case_with_grid(build_grid()).addMaterialFragment(fragment_path)


def test_casemaker_add_cell_element_intervals():
    cm = CaseMaker()
    first = cm.addCellElementIntervals(
        [[[0, 0, 0], [2, 1, 1]]], name="first")
    second = cm.addCellElementIntervals(
        [[[1, 1, 1], [3, 2, 2]]])

    assert first == 1
    assert second == 2
    assert cm['mesh']['elements'][0] == {
        "id": 1,
        "name": "first",
        "type": "cell",
        "intervals": [[[0, 0, 0], [2, 1, 1]]],
    }
    assert 'name' not in cm['mesh']['elements'][1]


def test_casemaker_add_movie_probe():
    cm = CaseMaker()
    cm.addMovieProbe(3, "movie_default")
    cm.addMovieProbe(
        4, "movie_magnetic",
        field="magnetic",
        component="x",
        domain={"type": "time", "samplingPeriod": 1e-9},
    )

    assert cm['probes'][0] == {
        "name": "movie_default",
        "type": "movie",
        "field": "electric",
        "component": "magnitude",
        "elementIds": [3],
    }
    assert cm['probes'][1]['field'] == 'magnetic'
    assert cm['probes'][1]['component'] == 'x'
    assert cm['probes'][1]['domain']['samplingPeriod'] == 1e-9


def build_fragment_case(tmp_path, monkeypatch):
    grid = build_grid(nx=8, ny=8, nz=6, step=0.01)
    grid_path = tmp_path / "grid.json"
    grid_path.write_text(json.dumps(grid))

    fragment = build_fragment(grid=grid)
    fragment_path = tmp_path / "fragment.json"
    fragment_path.write_text(json.dumps(fragment))

    monkeypatch.chdir(tmp_path)

    cm = CaseMaker()
    cm.setGridFromJSON(grid_path)
    cm.setNumberOfTimeSteps(1)
    cm.setTimeStep(1.0e-12)
    cm.setAllBoundaries("pec")

    tfsf = cm.addCellElementIntervals(
        [[[2, 2, 2], [6, 6, 5]]], name="tfsf_box")
    ntff = cm.addCellElementIntervals(
        [[[1, 1, 1], [7, 7, 6]]], name="ntff_box")
    movie = cm.addCellElementIntervals(
        [[[2, 3, 3], [6, 4, 4]]], name="movie_z0")
    body = cm.addCellElementIntervals(
        [[[3, 3, 2], [5, 5, 4]]], name="body_pec")
    pec = cm.addPECMaterial()
    cm.addMaterialAssociation(pec, [body])

    cm.addMaterialFragment(fragment_path)

    t = np.arange(0.0, 2.0e-10, 1.0e-12)
    data = np.column_stack(
        (t, np.exp(-((t - 5.0e-11) / 1.0e-11) ** 2)))
    np.savetxt(tmp_path / "gauss.exc", data)
    cm.addPlanewaveSource(
        tfsf, "gauss.exc",
        {"theta": np.pi / 2.0, "phi": 0.0},
        {"theta": 0.0, "phi": 0.0})
    cm.addFarFieldProbe(
        ntff, "backscatter",
        {"initial": 90.0, "final": 90.0, "step": 1.0},
        {"initial": 180.0, "final": 180.0, "step": 1.0},
        {"type": "frequency", "initialFrequency": 0.2e9,
         "finalFrequency": 1.5e9, "numberOfFrequencies": 5,
         "frequencySpacing": "linear"})
    cm.addMovieProbe(
        movie, "movie_electric", domain={
            "type": "time",
            "initialTime": 0.0,
            "finalTime": 1.0e-11,
            "samplingPeriod": 1.0e-12,
        })
    cm.exportCase("fragment_case")
    return cm, tmp_path / "fragment_case.fdtd.json"


def test_casemaker_fragment_case_exports_expected_json(tmp_path, monkeypatch):
    cm, case_path = build_fragment_case(tmp_path, monkeypatch)
    exported = json.loads(case_path.read_text())
    assert exported == cm.input

    assert exported['general']['numberOfSteps'] == 1
    assert exported['boundary']['all']['type'] == 'pec'
    assert len(exported['materials']) == 2
    assert exported['materials'][0]['type'] == 'pec'
    assert exported['materials'][1]['type'] == 'lorentz'
    assert len(exported['materials'][1]['levels']) == 2
    assert len(exported['mesh']['elements']) == 6
    assert exported['mesh']['elements'][3]['name'] == 'body_pec'
    assert exported['mesh']['elements'][4]['name'] == 'fp_0_nu_0'
    assert exported['mesh']['elements'][5]['name'] == 'fp_1_nu_1'
    assert len(exported['materialAssociations']) == 3
    assert exported['materialAssociations'][1]['level'] == 0
    assert exported['materialAssociations'][2]['level'] == 1
    assert [p['type'] for p in exported['probes']] == [
        'farField', 'movie']
    assert exported['sources'][0]['type'] == 'planewave'


def test_casemaker_fragment_case_launches_fdtd_preprocessing(
        tmp_path, monkeypatch):
    semba = find_fdtd_executable()
    if semba is None:
        pytest.skip("semba-fdtd not found; build fdtd or set SEMBA_EXE.")

    _, case_path = build_fragment_case(tmp_path, monkeypatch)
    result = subprocess.run(
        [semba, "-i", case_path.name, "-n", "-1"],
        cwd=tmp_path,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
        capture_output=True,
        text=True,
        timeout=300)
    output = result.stdout + result.stderr

    assert result.returncode == 0, output
    assert "ERROR" not in output, output
    assert "STOP" not in output, output
    assert "ENDED NFDE" in output or "BEGUN" in output, output
