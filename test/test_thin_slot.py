import json
import os
import subprocess

import numpy as np
import pytest

from utils import *
from src.caseMaker import *


def build_thin_slot_line_case():
    cm = CaseMaker()
    cm.setGridFromVTK(
        THIN_SLOT_FOLDER + 'squareInSquare.tessellator.grid.vtk')

    plate_id = cm.addCellElementsFromVTK(
        THIN_SLOT_FOLDER + 'SquareInSquare.tessellator.str.vtk',
        name='plate')
    slot_id = cm.addThinSlotElementsFromVTK(
        THIN_SLOT_FOLDER + 'Discretize_Wire.tessellator.str.vtk',
        name='slot')

    pec_id = cm.addPECMaterial()
    cm.addMaterialAssociation(pec_id, [plate_id])

    slot_material_id = cm.addThinSlotMaterial(1e-3, name='slot')
    cm.addMaterialAssociation(slot_material_id, [slot_id])

    return cm, plate_id, slot_id, slot_material_id


def build_thin_slot_surface_case():
    cm = CaseMaker()
    cm.setGridFromVTK(
        THIN_SLOT_SURFACE_FOLDER
        + 'squareInSquare.embeddedSurface.tessellator.grid.vtk')

    plate_id = cm.addCellElementsFromVTK(
        THIN_SLOT_SURFACE_FOLDER + 'SquareInSquare.tessellator.str.vtk',
        name='plate')
    slot_id = cm.addThinSlotElementsFromVTK(
        THIN_SLOT_SURFACE_FOLDER + 'Discretize_Wire.tessellator.str.vtk',
        name='slot')

    pec_id = cm.addPECMaterial()
    cm.addMaterialAssociation(pec_id, [plate_id])

    slot_material_id = cm.addThinSlotMaterial(1e-3, name='slot')
    cm.addMaterialAssociation(slot_material_id, [slot_id])

    return cm, plate_id, slot_id, slot_material_id


def assert_thin_slot_materials(cm, plate_id, slot_id, slot_material_id):
    assert cm['materials'] == [
        {'id': 1, 'type': 'pec'},
        {'id': 2, 'type': 'thinSlot', 'width': 1e-3, 'name': 'slot'}
    ]
    assert cm['materialAssociations'] == [
        {'materialId': 1, 'elementIds': [plate_id]},
        {'materialId': slot_material_id, 'elementIds': [slot_id],
         'type': 'thinSlot'}
    ]


def test_casemaker_thin_slot_line_from_vtk(tmp_path, monkeypatch):
    cm, plate_id, slot_id, slot_material_id = build_thin_slot_line_case()

    assert cm['mesh']['grid']['numberOfCells'] == [40, 40, 23]
    assert cm['mesh']['grid']['origin'] == [-20.0, -20.0, -10.0]
    assert_thin_slot_materials(cm, plate_id, slot_id, slot_material_id)

    plate = cm['mesh']['elements'][plate_id - 1]
    assert plate['name'] == 'plate'
    assert len(plate['intervals']) == 124
    for start, end in plate['intervals']:
        assert np.count_nonzero(np.array(start) != np.array(end)) == 2

    slot = cm['mesh']['elements'][slot_id - 1]
    assert slot['name'] == 'slot'
    assert len(slot['intervals']) == 44
    for start, end in slot['intervals']:
        assert np.count_nonzero(np.array(start) != np.array(end)) == 1

    monkeypatch.chdir(tmp_path)
    cm.exportCase('thin_slot_line')
    with open('thin_slot_line.fdtd.json') as f:
        exported = json.load(f)
    assert exported == cm.input


def test_casemaker_thin_slot_surface_from_vtk(tmp_path, monkeypatch):
    cm, plate_id, slot_id, slot_material_id = build_thin_slot_surface_case()

    assert cm['mesh']['grid']['numberOfCells'] == [40, 40, 23]
    assert cm['mesh']['grid']['origin'] == [-20.0, -20.0, -10.0]
    assert_thin_slot_materials(cm, plate_id, slot_id, slot_material_id)

    plate = cm['mesh']['elements'][plate_id - 1]
    assert plate['name'] == 'plate'
    assert len(plate['intervals']) == 124

    slot = cm['mesh']['elements'][slot_id - 1]
    assert slot['name'] == 'slot'
    assert len(slot['intervals']) == 43
    for start, end in slot['intervals']:
        assert np.count_nonzero(np.array(start) != np.array(end)) == 2

    monkeypatch.chdir(tmp_path)
    cm.exportCase('thin_slot_surface')
    with open('thin_slot_surface.fdtd.json') as f:
        exported = json.load(f)
    assert exported == cm.input


def test_casemaker_thin_slot_surface_launches_fdtd(tmp_path, monkeypatch):
    semba = find_fdtd_executable()
    if semba is None:
        pytest.skip("semba-fdtd not found; build fdtd or set SEMBA_EXE.")

    cm, plate_id, slot_id, slot_material_id = build_thin_slot_surface_case()
    cm.setNumberOfTimeSteps(1)
    cm.setTimeStep(1e-12)
    cm.setAllBoundaries('mur')

    box_id = cm.addCellElementBox([[-19.0, -19.0, -9.0], [19.0, 19.0, 12.0]])
    t = np.arange(0.0, 3e-9, 1e-12)
    excitation = np.column_stack(
        (t, np.exp(-((t - 1e-9) / 0.2e-9) ** 2)))

    monkeypatch.chdir(tmp_path)
    np.savetxt('gauss.exc', excitation)
    cm.addPlanewaveSource(
        box_id, 'gauss.exc',
        {"theta": np.pi / 2, "phi": 0.0},
        {"theta": np.pi / 2, "phi": np.pi / 2})
    cm.addPointProbe(cm.addNodeElement([0.0, 0.0, 0.0]), 'center')
    cm.exportCase('thin_slot_surface')

    result = subprocess.run(
        [semba, '-i', 'thin_slot_surface.fdtd.json'],
        cwd=tmp_path,
        env={**os.environ, 'OMP_NUM_THREADS': '1'},
        capture_output=True,
        text=True,
        timeout=300)
    output = result.stdout + result.stderr

    assert result.returncode == 0, output
    assert 'ERROR' not in output, output
    assert 'STOP' not in output, output
    assert 'DONE :' in output or 'ENDED' in output, output
    assert list(tmp_path.glob('thin_slot_surface.fdtd_center_Ex_*.dat'))


def test_casemaker_grid_alignment_with_non_integer_steps():
    cm = CaseMaker()
    cm.input['mesh'] = {
        'grid': {
            'origin': [-20.0, -20.0, -10.0],
            'numberOfCells': [40, 40, 20],
            'steps': {
                'x': [1.0] * 40,
                'y': [1.0] * 40,
                'z': [1.15] * 20,
            },
        }
    }

    box_id = cm.addCellElementBox(
        [[-20.0, -20.0, 0.3499999940395355], [20.0, 20.0, 13.0]])

    assert cm['mesh']['elements'][box_id - 1]['intervals'] == [
        [[0, 0, 9], [40, 40, 20]]
    ]


def test_casemaker_thin_slot_rejects_unsupported_cells():
    cm = CaseMaker()
    cm.setGridFromVTK(
        THIN_SLOT_FOLDER + 'squareInSquare.tessellator.grid.vtk')

    with pytest.raises(ValueError):
        cm.addThinSlotElementsFromVTK(
            THIN_SLOT_FOLDER + 'SquareInSquare.vtu')


def test_casemaker_thin_slot_rejects_non_positive_width():
    cm = CaseMaker()

    with pytest.raises(ValueError):
        cm.addThinSlotMaterial(0.0)


def test_casemaker_conformal_volume_from_legacy_vtu():
    cm = CaseMaker()
    cm.setGridFromVTK(
        THIN_SLOT_FOLDER + 'squareInSquare.tessellator.grid.vtk')

    element_id = cm.addConformalVolumeFromVTK(
        THIN_SLOT_FOLDER + 'SquareInSquare.vtu', name='surface')

    element = cm['mesh']['elements'][element_id - 1]
    assert element['type'] == 'cell'
    assert element['subtype'] == 'volume'
    assert element['name'] == 'surface'
    assert len(element['triangles']) == 54
