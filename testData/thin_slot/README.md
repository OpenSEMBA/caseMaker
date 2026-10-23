# Thin-slot test data

- `SquareInSquare.vtu`: staircased conductive plate.
- `Discretize_Wire.vtu`: slot centreline.
- `squareInSquare.tessellator.json`: tessellator input.
- `SquareInSquare.tessellator.str.vtk`,
  `Discretize_Wire.tessellator.str.vtk` and
  `squareInSquare.tessellator.grid.vtk`: generated with

  ```shell
  tessellator -i squareInSquare.tessellator.json
  ```

  using the `OpenSEMBA/tessellator` dev commit `ce8c780`.
