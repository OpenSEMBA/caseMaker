import json

import pytest

from utils import *
from src.caseMaker import *


def build_thin_slot_case():
    cm = CaseMaker()
    cm.setGridFromVTK(
        THIN_SLOT_FOLDER + 'squareInSquare.tessellator.grid.vtk')

    plateId = cm.addCellElementsFromVTK(
        THIN_SLOT_FOLDER + 'SquareInSquare.tessellator.str.vtk',
        name='plate')
    slotId = cm.addThinSlotElementsFromVTK(
        THIN_SLOT_FOLDER + 'Discretize_Wire.tessellator.str.vtk',
        name='slot')

    pecId = cm.addPECMaterial()
    cm.addMaterialAssociation(pecId, [plateId])

    slotMaterialId = cm.addThinSlotMaterial(1e-3, name='slot')
    cm.addMaterialAssociation(slotMaterialId, [slotId])

    return cm, plateId, slotId, slotMaterialId


def test_casemaker_thin_slot_from_vtk(tmp_path, monkeypatch):
    cm, plateId, slotId, slotMaterialId = build_thin_slot_case()

    assert cm['mesh']['grid']['numberOfCells'] == [40, 40, 23]
    assert cm['mesh']['grid']['origin'] == [-20.0, -20.0, -10.0]

    assert cm['materials'] == [
        {'id': 1, 'type': 'pec'},
        {'id': 2, 'type': 'thinSlot', 'width': 1e-3, 'name': 'slot'}
    ]
    assert cm['materialAssociations'] == [
        {'materialId': 1, 'elementIds': [plateId]},
        {'materialId': slotMaterialId, 'elementIds': [slotId],
         'type': 'thinSlot'}
    ]

    plate = cm['mesh']['elements'][plateId - 1]
    assert plate['name'] == 'plate'
    assert len(plate['intervals']) == 124
    for start, end in plate['intervals']:
        assert np.count_nonzero(np.array(start) != np.array(end)) == 2

    slot = cm['mesh']['elements'][slotId - 1]
    assert slot['name'] == 'slot'
    assert len(slot['intervals']) == 44
    for start, end in slot['intervals']:
        assert np.count_nonzero(np.array(start) != np.array(end)) == 1

    monkeypatch.chdir(tmp_path)
    cm.exportCase('thin_slot')
    with open('thin_slot.fdtd.json') as f:
        exported = json.load(f)
    assert exported == cm.input


def test_casemaker_thin_slot_rejects_non_line_cells():
    cm = CaseMaker()
    cm.setGridFromVTK(
        THIN_SLOT_FOLDER + 'squareInSquare.tessellator.grid.vtk')

    with pytest.raises(ValueError):
        cm.addThinSlotElementsFromVTK(
            THIN_SLOT_FOLDER + 'SquareInSquare.tessellator.str.vtk')


def test_casemaker_thin_slot_rejects_non_positive_width():
    cm = CaseMaker()

    with pytest.raises(ValueError):
        cm.addThinSlotMaterial(0.0)


def test_casemaker_conformal_volume_from_legacy_vtu():
    cm = CaseMaker()
    cm.setGridFromVTK(
        THIN_SLOT_FOLDER + 'squareInSquare.tessellator.grid.vtk')

    elementId = cm.addConformalVolumeFromVTK(
        THIN_SLOT_FOLDER + 'SquareInSquare.vtu', name='surface')

    element = cm['mesh']['elements'][elementId - 1]
    assert element['type'] == 'cell'
    assert element['subtype'] == 'volume'
    assert element['name'] == 'surface'
    assert len(element['triangles']) == 54
