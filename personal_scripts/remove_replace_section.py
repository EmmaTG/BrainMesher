import os
from math import ceil

from mesh.Element import HexElement, QuadElement
from mesh.PostProcessor import PostProcessor, SmoothMesh
from mesh.mesh_transformations import rotate_mesh
from mesh.refinement.Refiner import Refiner
import mesh.mesh_utils as mu
import mesh.mesh_transformations as mt
from mesh.Node import Node
from BinaryTree import TreeNode
import numpy as np
import common.helper_functions as hf
from writers.HeterogeneityConverter import MaterialsConverterFactory, Heterogeneity

from readers.Reader import Reader
from writers.Writers import Writer


def replace_duplicate_nodes(nodes_changed_map, elems_changed_maps, nodes, decimal_places=6, bounds=[-1000]):
    print("Replacing duplicated nodes")
    max_node_num = hf.calculate_max_number(nodes)
    node_location_map_x = {}
    for n, node in nodes.items():
        coords = node.getCoords()
        [xcoord, ycoord, zcoord] = [round(i, decimal_places) for i in coords]
        node_location_map_y = node_location_map_x.get(xcoord, {})
        node_location_map_z = node_location_map_y.get(ycoord, {})
        if node_location_map_z.get(zcoord,False):
            node_location_map_z[zcoord] = node
        else:
            node_location_map_z[zcoord] = node
            node_location_map_y[ycoord] = node_location_map_z
            node_location_map_x[xcoord] = node_location_map_y


    swapMap = {}
    for n, node in nodes_changed_map.items():
        orig_coords = node.getCoords()
        # if not nodeMap.__contains__(n):
        [xcoord, ycoord, zcoord] = [round(i, decimal_places) for i in orig_coords]
        node_location_map_y = node_location_map_x.get(xcoord, {})
        node_location_map_z = node_location_map_y.get(ycoord, {})
        if node_location_map_z.get(zcoord,False):
            swapMap[n] = node_location_map_z[zcoord]
        else:
            new_node = Node(n,orig_coords)
            node_location_map_z[zcoord] = new_node
            node_location_map_y[ycoord] = node_location_map_z
            node_location_map_x[xcoord] = node_location_map_y
            nodes[n] = new_node

    for e in elems_changed_maps.values():
        e_ica = e.ica
        for pos, node in enumerate(e_ica):
            n = node.number
            if swapMap.get(n, False):
                e_ica[pos] = swapMap[n]
            else:
                e_ica[pos] = node


path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/CoreFoamModels/"

# filenameIN = 'brain_base_noCSF_noSmoothing'
filenameIN = 'brain_base_noCSF_noSmoothing_converted_VTK'

reader = Reader('vtk')
reader.openReader(filenameIN, path)
mesh_base = reader.getMesh()
location = 1
diameter = 4
depth = 30
width = 8
introducer = 6

print("location: {}".format(location))
material_model = '4R'
refine_corners = False
refine_mesh = True
toBeSmoothed = True
debug = True
replacement_center_coords = [34, 0, 0]
# filename_replacement = "{}_{}mm_hole_8x8_{}mm".format(geo,diameter,depth) # Spherical hole

##############################################################################
##############################################################################

# nodeMap,elementMap,elementSetsMap = functions.readABQ(path, filenameIN)
locations = {1: [50, 18, -30], 2: [62, -24, 10], 3: [50, 16, 32], 4: [24, -30, -58], 5: [8, -20, 60],
             6: [-4, 2, -58]}

# tmp_nodeMap, tmp_elementMap, x = functions.readABQ(path, filename_replacement)
length_of_replacement = hf.round_to_interval(ceil(depth + 4 + introducer - (diameter / 2.0)), 2 / 3.)
y_direction_length = width
z_direction_length = width
safety_factor = 0.1
refinement_y = 6
refirnemnt_z = 6

# % Refine mesh
center_node = locations[location]
if center_node[2] > 0:
    side = 'R'
else:
    side = 'L'
center_node = locations[location]

if refine_mesh:
    # nodeMap, elementMap, elementSetsMap = functions.readABQ(path, filenameIN)
    bounds_to_refine = [center_node[0] - safety_factor - (length_of_replacement * 1.1),
                        center_node[0] + safety_factor + 50,  # (xmin,xmax  ...
                        center_node[1] - safety_factor - refinement_y,
                        center_node[1] + safety_factor + refinement_y,  # ymin,ymax ...
                        center_node[2] - safety_factor - refirnemnt_z,
                        center_node[2] + safety_factor + refirnemnt_z]  # zmin,zmax)
    refine_mesh = Refiner(mesh_base)
    refine_mesh.refine_within_region(bounds_to_refine)

# %% Remove and replace section in refined model with mesh with spherical hole (filename_replacement)
writer = Writer()
new_filename = "_".join([filenameIN, "refined"])
writer.openWriter('vtk', new_filename, path)
writer.writeMeshData(mesh_base)
writer.closeWriter()

# for geo, curve in [['sphere', 5], ['sphere', 6], ['sphere', 8], ['sphere', 10], ['sphere', 11]]:
for geo, curve in [['ellipse', 11]]:

    filename_replacement = "{}_in_cube_{}mm_{}mm_{}mm_curve{}".format(geo, str(diameter).replace(".", "p")[0:3],
                                                                      depth, introducer, curve)  # Spherical hole

    filenameOUT = 'Brain_{}_{}mm_{}mm_I{}mm_curve{}_Loc{}_{}'.format(geo, str(diameter).replace(".", "p")[0:3], depth,
                                                                     introducer, curve, location, material_model)
    reader = Reader('abaqus')
    reader.openReader(filename_replacement, path)
    mesh_hole = reader.getMesh()

    reader = Reader('vtk')
    new_filename = "_".join([filenameIN, "refined_VTK"])
    reader.openReader(new_filename, path)
    mesh_base = reader.getMesh()

    elementICAMap = mu.create_elements_ica_map(mesh_base.elements)
    nodeToElementMap = mu.create_node_to_elem_map(elementICAMap)
    print("Creating hole by removing elements from " + filenameIN + " and replacing with " + filename_replacement)

    bounds_to_remove = [center_node[0] - length_of_replacement - safety_factor, center_node[0] + safety_factor,
                        # (xmin,xmax,
                        center_node[1] - y_direction_length / 2 - safety_factor,
                        center_node[1] + y_direction_length / 2. + safety_factor,  # ymin,ymax
                        center_node[2] - z_direction_length / 2 - safety_factor,
                        center_node[2] + z_direction_length / 2. + safety_factor]  # zmin,zamx)

    # Create map of elements and nodes to be extracted from original mesh
    extraction_Elements = []
    extraction_Nodes = []
    for e in mesh_base.elements.values():
        centroid = e.calculate_element_centroid()
        if hf.value_in_square_bounds(centroid, bounds_to_remove):
            extraction_Elements.append(e.num)
            for n in e.ica:
                 if not extraction_Nodes.count(n.number):
                    extraction_Nodes.append(n.number)

    # extraction_NodesMap = set(list([int(x.number) for x in mesh_base.nodes.values()]))
    # Remove elements from element Map and their associated elset lists.
    # Create list of removed elements
    removed_elements = []
    # base_mesh_elements = mesh_base.elements
    for e_num in extraction_Elements:
        e = mesh_base.elements.pop(e_num)
        removed_elements.append(e)

    for n in extraction_Nodes:
        allElemsRemoved = True
        for e in nodeToElementMap[n]:
            if not extraction_Elements.count(e):
                allElemsRemoved = False
                break
        if allElemsRemoved:
            mesh_base.nodes.pop(n)

    # Create binary tree from extracted elements to assit with element material assignment later
    r = TreeNode(0, removed_elements.pop(0))
    for e in removed_elements:
        r = r.insert(r, e)

    # Increment element and node numbers based on max numbers in full mesh model
    # so as to not over accidental overwriting of nodes or elements
    replacement_nodeMap = hf.increment_numbers(mesh_hole.nodes, mesh_base.nodes)
    replacement_elementMap = hf.increment_numbers(mesh_hole.elements, mesh_base.elements)

    # Shift new mesh section to co-incide with extracted section prepared above

    for n in replacement_nodeMap.values():
        coords = n.getCoords()
        coords[0] = coords[0] + (center_node[0] - replacement_center_coords[0])
        coords[1] = coords[1] + (center_node[1] - replacement_center_coords[1])
        coords[2] = coords[2] + (center_node[2] - replacement_center_coords[2])

        # Remove overlapping/duplicate nodes within bounds, add new nodes from inserted section to nodeMap
    replace_duplicate_nodes(replacement_nodeMap, replacement_elementMap, mesh_base.nodes, 1)

    # Assign elsets in new elements based on closest element on orginal extracted mesh
    for e_num, element in replacement_elementMap.items():
        e_centroid = element.calculate_element_centroid()
        x_dim_search = r.search(r, e_centroid[0])
        y_dim_search = x_dim_search.next_dim_node.search(x_dim_search.next_dim_node, e_centroid[1])
        # if (y_dim_search.next_dim_node != None):
        z_dim_search = y_dim_search.next_dim_node.search(y_dim_search.next_dim_node, e_centroid[2])
        foundElement = z_dim_search.data.elements
        assert len(foundElement) == 1
        element.setMaterial(foundElement[0].getMaterial())

    # Merge original and new element map
    mesh_base.addElements(replacement_elementMap)

     # Defining boundary conditions
    print("Defining Boundary Conditions")
    boundary_element_map = mesh_base.locate_boundary_element_map()
    boundary_elem_to_elset_map = {}
    fixed_boundary = 4.5
    fixed = []
    prescribed = []
    boundary_elements = {}
    maxNum = max(list(mesh_base.elements.keys()))
    maxNum = int(ceil(maxNum/10)*10)
    # face_to_elem_map =
    for key, face_nodes in boundary_element_map.items():
        e_num, face_num = [int(x) for x in key.split("-")]
        volume_elem = mesh_base.elements.get(e_num)
        face_ica = []
        for n in face_nodes:
            face_ica.append(mesh_base.nodes[n])
        centroid = volume_elem.calculate_element_centroid()
        mat = volume_elem.getMaterial()
        if mat.count(16) and centroid[0] < -57:
            maxNum += 1
            boundary_elements[maxNum] = QuadElement(maxNum, face_ica, mat=2)
        elif (((side == 'L') and (centroid[2] > fixed_boundary)) or
              ((side == 'R') and (centroid[2] < -1 * fixed_boundary))):
            maxNum += 1
            boundary_elements[maxNum] = QuadElement(maxNum, face_ica, mat=2)
        elif ((centroid[0] > (center_node[0] - depth - 0.1)) and (centroid[0] < (center_node[0] + 2)) and
              (centroid[1] > (center_node[1] - diameter/2. - 1)) and (
                      centroid[1] < (center_node[1] + diameter/2. + 1)) and
              (centroid[2] > (center_node[2] - diameter/2. - 1)) and (
                      centroid[2] < (center_node[2] + diameter/2. + 1))):

            v1v2 = np.array(face_ica[0].getCoords()) - np.array(face_ica[1].getCoords())
            v1v3 = np.array(face_ica[0].getCoords()) - np.array(face_ica[2].getCoords())
            normal = np.cross(v1v2, v1v3)
            norm = normal/np.linalg.norm(normal)
            # incision conditions
            # spherical or ellipse conditins
            if round( abs(norm[0]) ) == 0:
                # print("Appplying boundary conditions to spherical and elliptical models")
                maxNum += 1
                boundary_elements[maxNum] = QuadElement(maxNum, face_ica, mat=1)

    mesh_base.addBoundaryElements(boundary_elements)

    # if (location==4) or (location==6):
    #     elem_functions.rotate_mesh(nodeMap,1,-90)
    # if location==5:
    #     elem_functions.rotate_mesh(nodeMap,1,90)

    # elementToElset_new = functions.create_element_to_elset_map(elementSetsMap)
    # %% Smoothing mesh
    if toBeSmoothed:
        bounds = [center_node[0] - depth, center_node[0] + 2,
                  center_node[1] - diameter / 2. - 0.1, center_node[1] + diameter / 2. + 0.1,
                  center_node[2] - diameter / 2. - 0.1, center_node[2] + diameter / 2. + 0.1]

        mesh_base.smooth_mesh([0.6, -0.4], 6, bounds=bounds)

    # Move mesh to be centered at (x,0,0)
    mt.translate_mesh(mesh_base.nodes, [-1*center_node[0], -1*center_node[1], -1*center_node[2]])

    mat_converter = MaterialsConverterFactory.get_converter(Heterogeneity.FOURR)
    mat_converter.convert_materials_labels(mesh_base)

    if not os.path.exists("/".join([path, filenameOUT])):
        os.mkdir("/".join([path, filenameOUT]))
    writer = Writer()
    writer.openWriter('ucd', filenameOUT, "/".join([path,filenameOUT]))
    writer.writeMeshData(mesh_base)
    writer.closeWriter()
    writer.openWriter('vtk', filenameOUT, path)
    writer.writeMeshData(mesh_base)
    writer.closeWriter()







