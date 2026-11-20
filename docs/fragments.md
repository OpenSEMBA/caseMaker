# External material fragments

CaseMaker can merge a material/mesh fragment produced by an external
preprocessing tool.  CaseMaker only consumes the JSON file; it does not run,
import, or inspect the tool that produced it.

## Fragment format

The fragment has three top-level fields:

```json
{
  "material": {
    "id": 1,
    "name": "plasma",
    "type": "lorentz",
    "epsInf": 1.0,
    "omega1": 0.0,
    "levels": [
      {"omegaP": 27497.99, "gamma": 95598295.23,
       "f_p_level": 0, "nu_c_level": 0}
    ]
  },
  "mesh": {
    "grid": {
      "numberOfCells": [128, 128, 128],
      "origin": [-0.51516904, -1.09476239, -1.09476239],
      "steps": {
        "x": [0.0164875],
        "y": [0.0150898],
        "z": [0.0150898]
      }
    },
    "elements": [
      {
        "id": 1,
        "name": "fp_0_nu_0",
        "type": "cell",
        "intervals": [
          [[24, 30, 64], [31, 31, 65]]
        ]
      }
    ]
  },
  "materialAssociations": [
    {"materialId": 1, "level": 0, "elementIds": [1]}
  ]
}
```

- `material` is the material definition.  Level-based Lorentz materials use a
  `levels` array; other existing material types use their normal fields.
- `mesh.elements` contains one `cell` element per material level.  Each
  element can contain many closed-open intervals.
- `mesh.grid.steps` may contain one value per axis or one value per cell.
- `materialAssociations[].level`, when present, is the zero-based index into
  `material.levels`.

The fragment does not contain sources, probes, boundaries or the PEC body.
Those remain caseMaker inputs.

## Merging

```python
cm.setGridFromJSON("grid.json")
cm.addMaterialFragment("fragment.json")
```

`addMaterialFragment()`:

- validates that `material`, `mesh`, `materialAssociations`, `mesh.grid` and
  `mesh.elements` are present;
- validates `numberOfCells` exactly and `origin`/`steps` within `1e-12`;
- adopts the fragment grid if the case has no grid yet;
- appends the material to `materials`;
- appends its elements to `mesh.elements`;
- rewrites `materialId` and `elementIds` to the next free ids;
- preserves each association’s `level` index.

The fragment is not modified.

## Grid export

`setGridFromJSON(path)` loads a grid file, `writeGridJSON(path)` writes the
current case grid, and the grid can also come from tessellator through
`setGridFromVTK()`.

## Supporting caseMaker helpers

Two helpers are useful when assembling a case around a fragment:

- `addCellElementIntervals(intervals, name=None)` adds a `cell` element
  directly from closed-open relative node intervals.  This is used for TFSF
  and NTFF boxes, movie slices, or any geometry already expressed as grid
  intervals.
- `addMovieProbe(elementId, name, field="electric",
  component="magnitude", domain=None)` adds a ParaView movie probe.  The
  referenced element must be one cell region with one interval.

## PEC body

The fragment is geometry-agnostic.  Supply the PEC body separately, for
example from tessellator:

```python
body = cm.addCellElementsFromVTK(
    "body_pec_surface.vtk", name="body_pec")
pec = cm.addPECMaterial()
cm.addMaterialAssociation(pec, [body])
```

Overlaps between fragment elements and PEC are resolved by fdtd material
priority; PEC wins.

## Example

```python
from caseMaker import CaseMaker

cm = CaseMaker()
cm.setGridFromJSON("grid.json")
cm.setNumberOfTimeSteps(3300)
cm.setTimeStep(2.390423293775277e-11)
cm.setAllBoundaries("pml")

tfsf = cm.addCellElementIntervals(
    [[[16, 16, 16], [112, 112, 112]]], name="tfsf_box")
ntff = cm.addCellElementIntervals(
    [[[8, 8, 8], [120, 120, 120]]], name="ntff_box")

body = cm.addCellElementsFromVTK(
    "body_pec_surface.vtk", name="body_pec")
pec = cm.addPECMaterial()
cm.addMaterialAssociation(pec, [body])

cm.addMaterialFragment("fragment.json")

cm.addPlanewaveSource(
    tfsf, "gauss.exc",
    {"theta": 1.5707963267948966, "phi": 0.0},
    {"theta": 0.0, "phi": 0.0})
cm.addFarFieldProbe(
    ntff, "backscatter",
    {"initial": 90.0, "final": 90.0, "step": 1.0},
    {"initial": 180.0, "final": 180.0, "step": 1.0},
    {"type": "frequency", "initialFrequency": 0.2e9,
     "finalFrequency": 1.5e9, "numberOfFrequencies": 141,
     "frequencySpacing": "linear"})

cm.exportCase("case")
```

## Tests

`test/test_fragments.py` covers:

- grid round-trip and validation;
- fragment ID remapping and level preservation;
- malformed and unknown associations;
- `addCellElementIntervals` and `addMovieProbe`;
- a small end-to-end case that is passed to `semba-fdtd -n -1` when the
  solver is available.

Run with:

```shell
python3 -m pytest test/test_fragments.py
```
