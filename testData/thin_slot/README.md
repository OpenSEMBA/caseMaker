# Thin-slot test data

- `SquareInSquare.vtu`: staircased conductive plate.
- `Discretize_Wire.vtu`: slot centreline.
- `squareInSquare.tessellator.json`: line-mode tessellator input.
- `SquareInSquare.tessellator.str.vtk`,
  `Discretize_Wire.tessellator.str.vtk` and
  `squareInSquare.tessellator.grid.vtk`: line-mode outputs.
- `surface/squareInSquare.embeddedSurface.tessellator.json`: embedded-surface
  tessellator input.
- `surface/*.vtk`: embedded-surface outputs used by the fdtd launch test.

Outputs were generated with `OpenSEMBA/tessellator` branch
`112-thin-slots-meshing` (PR #113):

```shell
tessellator -i squareInSquare.tessellator.json
tessellator -i surface/squareInSquare.embeddedSurface.tessellator.json
```
