import json
import os
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
from math import modf

class CaseMaker():
    # A collection of helper functions to create FDTD cases.
    def __init__(self, json_filename=""):
        self.input = {}
        if (json_filename != ""):
            self.input = json.load(open(json_filename))

    def __getitem__(self, key):
        return self.input[key]

    def _readVTK(self, path_to_vtk):
        assert os.path.isfile(path_to_vtk)

        with open(path_to_vtk, 'rb') as f:
            header = f.read(512).lstrip()

        # The input files may use a legacy VTK header even when the extension is
        # .vtu, so the format is detected from the file contents.
        if header.startswith(b'<'):
            reader = vtk.vtkXMLGenericDataObjectReader()
        else:
            reader = vtk.vtkGenericDataObjectReader()

        reader.SetFileName(path_to_vtk)
        reader.Update()
        return reader.GetOutput()

    def _checkSingleGroup(self, polyData):
        groups = polyData.GetCellData().GetArray('group')
        if groups is None:
            return
        vtkGroups = vtk_to_numpy(groups)
        if len(np.unique(vtkGroups)) > 1:
            raise ValueError("Different groups are not supported.")

    def setNumberOfTimeSteps(self, steps):
        if 'general' not in self.input:
            self.input['general'] = {}

        self.input['general']['numberOfSteps'] = steps

    def setTimeStep(self, time):
        if 'general' not in self.input:
            self.input['general'] = {}

        self.input['general']['timeStep'] = time

    def setGridFromVTK(self, path_to_grid):
        polyData = self._readVTK(path_to_grid)
        points = vtk_to_numpy(polyData.GetPoints().GetData())

        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'grid' not in self.input['mesh']:
            self.input['mesh']['grid'] = {}

        grid = self.input['mesh']['grid']
        if 'steps' not in grid:
            grid['steps'] = {}

        steps = grid['steps']
        grid['numberOfCells'] = []
        grid['origin'] = []

        index = 0
        for x in ['x', 'y', 'z']:
            steps[x] = np.diff(np.unique(points[:, index])).tolist()
            grid['numberOfCells'].append(len(steps[x]))
            grid['origin'].append(float(points[0, index]))
            index += 1

    def setAllBoundaries(self, boundary_type):
        self.input['boundary'] = {
            "all": {
                "type": boundary_type
            }
        }

    def addConformalVolumeFromVTK(self, path_to_vtk, subtype="volume", name=None):
        if subtype not in ("volume", "surface"):
            raise ValueError("subtype must be 'volume' or 'surface'.")

        polyData = self._readVTK(path_to_vtk)
        self._checkSingleGroup(polyData)

        # Stores in case input
        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'elements' not in self.input['mesh']:
            self.input['mesh']['elements'] = []
        elements = self.input['mesh']['elements']

        id = len(elements) + 1
        element = {
            "id": id,
            "type": "cell",
            "subtype": subtype,
            "intervals" : [],
            "triangles": self._getTriangles(polyData)
        }
        if name is not None:
            element["name"] = name
        elements.append(element)
        return id


    def addCellElementsFromVTK(self, path_to_vtk, name=None):
        polyData = self._readVTK(path_to_vtk)
        self._checkSingleGroup(polyData)
        return self._addCellElement(polyData, name=name)

    def addThinSlotElementsFromVTK(self, path_to_vtk, name=None):
        '''Adds line or embedded-surface cell elements for a thin slot.'''
        polyData = self._readVTK(path_to_vtk)
        self._checkSingleGroup(polyData)

        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType not in (vtk.VTK_LINE, vtk.VTK_QUAD):
                raise ValueError(
                    "Thin slot elements must be built from VTK_LINE or "
                    "VTK_QUAD cells only.")

        return self._addCellElement(polyData, name=name)

    def _addCellElement(self, polyData, name=None):
        # Stores in case input
        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'elements' not in self.input['mesh']:
            self.input['mesh']['elements'] = []
        elements = self.input['mesh']['elements']

        id = len(elements) + 1
        element = {
            "id": id,
            "type": "cell",
            "intervals": self._getCellsIntervals(polyData)
        }
        if name is not None:
            element["name"] = name
        elements.append(element)
        return id


    def addCellElementBox(self, boundingBox):
        ''' boundingBox is assumed to be in absolute coordinates.'''
        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'elements' not in self.input['mesh']:
            self.input['mesh']['elements'] = []
        elements = self.input['mesh']['elements']

        rel = self._absoluteToRelative(
            np.array([boundingBox[0], boundingBox[1]]), snap=True)

        id = len(elements) + 1
        elements.append({
            "id": id,
            "type": "cell",
            "intervals": [
                [rel[0].tolist(), rel[1, :].tolist()]
            ]
        })

        return id

    def addNodeElement(self, position):
        ''' position is assumed to be in absolute coordinates.'''
        if isinstance(position, list):
            pos = np.array(position).reshape((1, 3))
        elif isinstance(position, np.ndarray):
            pos = position.reshape((1, 3))
        else:
            raise ValueError("Invalid type for position")
        relative = self._absoluteToRelative(pos, snap=True)

        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'coordinates' not in self.input['mesh']:
            self.input['mesh']['coordinates'] = []

        coordId = len(self.input['mesh']['coordinates']) + 1
        self.input['mesh']['coordinates'].append({
            "id": coordId,
            "relativePosition": relative[0].tolist()
        })

        if 'elements' not in self.input['mesh']:
            self.input['mesh']['elements'] = []
        elementId = len(self.input['mesh']['elements']) + 1
        self.input['mesh']['elements'].append({
            "id": elementId,
            "type": "node",
            "coordinateIds": [coordId]
        })

        return elementId

    def addThinSlotMaterial(self, width, name=None):
        '''Adds a thinSlot material. Width is in metres, as used by the grid.'''
        if width <= 0:
            raise ValueError("The thin slot width must be positive.")

        if 'materials' not in self.input:
            self.input['materials'] = []

        id = len(self.input['materials']) + 1
        material = {
            "id": id,
            "type": "thinSlot",
            "width": float(width)
        }
        if name is not None:
            material["name"] = name
        self.input['materials'].append(material)

        return id

    def addPECMaterial(self):
        if 'materials' not in self.input:
            self.input['materials'] = []

        id = len(self.input['materials']) + 1
        self.input['materials'].append({
            "id": id,
            "type": "pec",
        })

        return id

    def addMaterialAssociation(self, materialId: int, elementIds: list):
        if 'materials' not in self.input:
            raise ValueError("No materials defined.")
        if 'mesh' not in self.input:
            raise ValueError("No mesh defined.")
        if 'elements' not in self.input['mesh']:
            raise ValueError("No elements defined.")
        if len(elementIds) == 0:
            raise ValueError("No elements to associate.")
        if materialId < 1 or materialId > len(self.input['materials']):
            raise ValueError("Invalid material id.")
        if any(elementId < 1 or elementId > len(self.input['mesh']['elements'])
               for elementId in elementIds):
            raise ValueError("Invalid element id.")

        if 'materialAssociations' not in self.input:
            self.input['materialAssociations'] = []

        association = {
            'materialId': materialId,
            'elementIds': elementIds
        }
        if self.input['materials'][materialId - 1]['type'] == 'thinSlot':
            association['type'] = 'thinSlot'

        self.input['materialAssociations'].append(association)

    def addPlanewaveSource(self, elementId, magnitudeFile, direction, polarization):
        assert os.path.isfile(magnitudeFile)

        if 'sources' not in self.input:
            self.input['sources'] = []

        self.input['sources'].append({
            'type': 'planewave',
            'magnitudeFile': magnitudeFile,
            'elementIds': [elementId],
            'direction': direction,
            'polarization': polarization
        })

    def addPointProbe(self, elementId, name):
        if 'probes' not in self.input:
            self.input['probes'] = []

        self.input['probes'].append({
            "name": name,
            "type": "point",
            "elementIds": [elementId]
        })

    def addFarFieldProbe(self, elementId, name, theta, phi, domain):
        if 'probes' not in self.input:
            self.input['probes'] = []

        self.input['probes'].append({
            "name": name,
            "type": "farField",
            "elementIds": [elementId],
            "theta": theta,
            "phi": phi,
            "domain": domain
        })

    def exportCase(self, case_name):
        with open(case_name + '.fdtd.json', 'w') as f:
            json.dump(self.input, f)

    def _buildGridLines(self):
        gridLines = []
        grid = self.input['mesh']['grid']
        index = 0
        for x in ['x', 'y', 'z']:
            newLines = np.cumsum(
                np.insert(grid['steps'][x], 0, grid['origin'][index]))
            index += 1
            gridLines.append(newLines)
        return gridLines

    def _relativeToAbsolute(self, relativeCoordinates):
        gridLines = self._buildGridLines()

        res = np.empty_like(relativeCoordinates, dtype=float)
        for (idx, c) in enumerate(relativeCoordinates):
            for x in range(3):
                f, i = modf(c[x])[0], int(modf(c[x])[1])
                if (i == len(gridLines[x])-1) :
                    res[idx, x] = gridLines[x][i]
                else:
                    res[idx, x] = gridLines[x][i] + f*(gridLines[x][i+1]-gridLines[x][i])
                # res[:, x] = gridLines[x][relativeCoordinates[:, x]]

        return res

    def _absoluteToRelative(self, absoluteCoordinates, arg_type = int, snap = False):
        gridLines = self._buildGridLines()
        coordinates = np.asarray(absoluteCoordinates, dtype=float)
        if coordinates.ndim == 1:
            coordinates = coordinates.reshape((1, 3))

        res = np.empty_like(coordinates, dtype=arg_type)
        for axis in range(3):
            lines = gridLines[axis]
            values = coordinates[:, axis]
            tolerance = 1e-6 * max(1.0, float(np.max(np.abs(lines))))

            if np.any(values < lines[0] - tolerance) or \
               np.any(values > lines[-1] + tolerance):
                raise ValueError(
                    "Coordinate is outside the grid in axis {}.".format(axis))

            idx = np.searchsorted(lines, values, side="right") - 1
            idx = np.clip(idx, 0, len(lines) - 2)
            lower = lines[idx]
            upper = lines[idx + 1]

            if (arg_type == int):
                nearest = np.where(
                    np.abs(values - lower) < np.abs(values - upper),
                    idx, idx + 1)
                if not snap and np.any(
                        np.abs(values - lines[nearest]) > tolerance):
                    raise ValueError(
                        "Coordinate is not aligned with a grid node in axis {}."
                        .format(axis))
                res[:, axis] = nearest
            else:
                span = upper - lower
                with np.errstate(divide="ignore", invalid="ignore"):
                    fraction = np.where(
                        span != 0, (values - lower) / span, 0.0)
                res[:, axis] = idx + fraction

        return res

    def _getTriangles(self, polyData):
        relative = self._absoluteToRelative(
            vtk_to_numpy(polyData.GetPoints().GetData()), arg_type = float)

        numberOfTriangles = 0
        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType == vtk.VTK_TRIANGLE:
                numberOfTriangles += 1
        res = [None] * numberOfTriangles
        
        c = 0
        listOfCoords = []
        vtkToCoordId = {}
        if 'coordinates' not in self.input['mesh']:
            self.input['mesh']['coordinates'] = []

        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType == vtk.VTK_TRIANGLE:
                vtk_ids = [polyData.GetCell(i).GetPointId(0),
                          polyData.GetCell(i).GetPointId(1),
                          polyData.GetCell(i).GetPointId(2)]                 
                vertices = []
                for (id,p) in enumerate(vtk_ids):
                    if p not in listOfCoords:
                        listOfCoords.append(p)

                        coordId = len(self.input['mesh']['coordinates']) + 1
                        self.input['mesh']['coordinates'].append({
                            "id": coordId,
                            "relativePosition": relative[p].tolist()
                        })
                        
                        vtkToCoordId[p] = coordId
                    vertices.append(vtkToCoordId[p])
                res[c] = vertices
                    
                c = c+1

        return res            
      

    def _getCellsIntervals(self, polyData):
        relative = self._absoluteToRelative(
            vtk_to_numpy(polyData.GetPoints().GetData()))

        res = []
        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType == vtk.VTK_QUAD:
                p1 = polyData.GetCell(i).GetPointId(0)
                p2 = polyData.GetCell(i).GetPointId(2)
                expectedVaryingDirections = 2
            elif cellType == vtk.VTK_LINE:
                p1 = polyData.GetCell(i).GetPointId(0)
                p2 = polyData.GetCell(i).GetPointId(1)
                expectedVaryingDirections = 1
            else:
                continue

            start = relative[p1]
            end = relative[p2]
            varyingDirections = np.count_nonzero(start != end)
            if varyingDirections != expectedVaryingDirections:
                raise ValueError(
                    "Cell elements must be axis-aligned with the grid. "
                    "Found a cell with {} varying directions.".format(
                        varyingDirections))

            res.append([start.tolist(), end.tolist()])

        return res
