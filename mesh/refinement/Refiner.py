import numpy as np
from mesh import mesh_utils as mu
from mesh.Element import HexElement
from mesh.Node import Node


class Refiner:
    def __init__(self, mesh):
        self.mesh = mesh
        self.nodeMap = mu.create_node_coords_map(mesh.nodes)
        assert len(self.nodeMap) == len(mesh.nodes)
        self.elementMap = mu.create_elements_ica_map(mesh.elements)
        assert len(self.elementMap) == len(mesh.elements)

    def refine_within_region(self, bounds):
        nodes_of_refinement = self.define_refinement_area(bounds)
        sideMap = self.assign_refinement_type(nodes_of_refinement)
        self.__refine_mesh__(nodes_of_refinement, sideMap)

    def refine_around_point(self, point, radius):
        bounds = [point[0]-radius, point[0]+radius,
                  point[1]-radius, point[1]+radius,
                  point[2]-radius, point[2]+radius]
        self.refine_within_region(bounds)

    def refine_elements(self, element_numbers):
        nodes_of_refinement = []
        for e in element_numbers:
            element = self.mesh.elements[e]
            for n in element.ica:
                nodes_of_refinement.append(n.number)
        sideMap = self.assign_refinement_type(nodes_of_refinement)
        self.__refine_mesh__(nodes_of_refinement, sideMap)

    def __refine_mesh__(self, refined_nodes, sideMap):
        """
        refinement mesh in area defined by refined_nodes and categorized by sideMap

        Parameters
        ----------
        refined_nodes : TYPE
            DESCRIPTION.
        sideMap :  Map(string,array)
            Map of categories of refinement (keys) to elements assigned with that refinement

        """
        print("Refining mesh")
        self.newNodeNum = max(list(self.nodeMap.keys()))+1
        self.newElemNum = max(list(self.elementMap.keys()))+1
        nodeLocationMapX = {}
        nodeLocationMapY = {}
        nodeLocationMapZ = {}
        refined_elements = []
        swapMap = {}
        for side, side_elements in sideMap.items():
            print("Refining {} elements in category {}".format(len(side_elements), side))
            for e in side_elements:
                current_ele = self.mesh.elements.get(e)
                nodes_to_be_refined = []
                old_elem_number = e
                ica_old_elements = self.elementMap[int(e)]
                nodes = np.zeros([8, 4])
                for count, n in enumerate(ica_old_elements):
                    nodes[count, :] = np.concatenate(([n], list(self.nodeMap[n])))
                    if refined_nodes.count(n):
                        nodes_to_be_refined.append(n)
                refined_elements.append(e)
                allNodes = []
                allElements = []
                if side == "P3":
                    if len(nodes_to_be_refined) >= 3:
                        allNodes, allElements = self.__P3_refinement__(nodes_to_be_refined, nodes)

                    else:
                        print("P3 element not connected to four refined nodes")
                        print("nodes_to_be_refined length: " + str(len(nodes_to_be_refined)))

                elif side == "P2":
                    if len(nodes_to_be_refined) == 2:
                        allNodes, allElements = self.__P2_refinement__(nodes_to_be_refined, nodes)
                    else:
                        print("Error. P2 element not connected to 2 nodes")
                elif side == "P1":
                    if len(nodes_to_be_refined) == 1:
                        allNodes, allElements = self.__P1_refinement__(nodes_to_be_refined, nodes)
                    else:
                        print("Less than 1 node to be refined. Element" + str(e) + "not refined.")

                elif side == "P4":
                    allNodes, allElements = self.__P4_refinement__(nodes)
                else:
                    print("No refinement")
                    # return 0

                self.newNodeNum += np.size(allNodes, 0) + 1
                self.newElemNum += np.size(allElements, 0)

                for n in allNodes:
                    [xcoord, ycoord, zcoord] = [str(round(i, 2)) for i in n[1:]]
                    xcoord = xcoord if xcoord != '-0.0' else '0.0'
                    ycoord = ycoord if ycoord != '-0.0' else '0.0'
                    zcoord = zcoord if zcoord != '-0.0' else '0.0'
                    if nodeLocationMapX.get(xcoord, False):
                        nodeLocationMapY = nodeLocationMapX.get(xcoord)
                    else:
                        nodeLocationMapY = {}
                    if nodeLocationMapY.get(ycoord, False):
                        nodeLocationMapZ = nodeLocationMapY.get(ycoord)
                    else:
                        nodeLocationMapZ = {}
                    if nodeLocationMapZ.get(zcoord, False):
                        swapMap[str(int(n[0]))] = nodeLocationMapZ.get(zcoord)
                    else:
                        nodeLocationMapZ[zcoord] = str(int(n[0]))
                        nodeLocationMapY[ycoord] = nodeLocationMapZ
                        nodeLocationMapX[xcoord] = nodeLocationMapY
                        if not self.nodeMap.__contains__(int(n[0])):
                            newNode = Node(int(n[0]), n[1:])
                            self.mesh.nodes[int(n[0])] = newNode
                            self.nodeMap[int(n[0])] = n[1:]
                        elif ((self.nodeMap[n[0]][0] != n[1]) or
                              (self.nodeMap[n[0]][1] != n[2]) or
                              (self.nodeMap[n[0]][2] != n[3])):
                            print("Error: replicated nodes in nodeMap with different coords")
                elementNumbers = []
                for e in allElements:
                    ica = []
                    for count, nodeNum in enumerate(e[1:]):
                        num_tmp = nodeNum
                        swappedNodeNumber = int(swapMap.get(str(int(nodeNum)), -1000))
                        if swappedNodeNumber != -1000:
                            num_tmp = swappedNodeNumber
                            e[count + 1] = swappedNodeNumber
                        ica.append(self.mesh.nodes[num_tmp])
                    assert len(ica) == 8
                    newElement = HexElement(int(e[0]), ica, mat=current_ele.getMaterial())
                    self.mesh.elements[int(e[0])] = newElement
                    self.elementMap[int(e[0])] = e[1:]
                    elementNumbers.append(int(e[0]))
                del self.elementMap[old_elem_number]
                del self.mesh.elements[old_elem_number]

        print("refinement COMPLETE\n")

    def refinement_node_map(self, category, rotatedNodes, dim=-1):
        """
        Map assigning node refinement based on category provided

        Parameters
        ----------
        category : String
            Caterory of refinement to do. Options are P1|P2|P3|P4.
        rotatedNodes : Array(ints)
            Array of node numbers on current unrefined element.
        startNodeNum : int
            Starting number of new nodes.
        dim : dimension, optional
            DImension needed for edge refinement which defines dimension parrallel to edge. The default is -1.

        Returns
        -------
        newNodes : Array([int,float,float,float])
            Array of node numbers and coordinates for new refined elements.

        """
        if category == "P1":
            newNodes = self.point_refinement_nodes(rotatedNodes)
        elif category == "P2":
            newNodes = self.edge_refinement_nodes(rotatedNodes, dim)
        elif category == "P3":
            newNodes = self.face_refinement_nodes(rotatedNodes)
        elif category == "P4":
            newNodes = self.full_refinement_nodes(rotatedNodes)
        else:
            print("Error: category (" + category + ") not allowed")
            newNodes = {}
        return newNodes

    def refinement_element_map(self, category, allNodes):
        """
        Map assigning element refinement based on category provided

        Parameters
        ----------
        category : String
            Caterory of refinement to do. Options are P1|P2|P3|P4.
        allNodes : Array(Array)
            Array of new nodes created from node refinement.
        startElemNum : int
            Starting number of new elements.

        Returns
        -------
        newElements : Array(Array(ints))
            Array of node icas and element numbers for new refined elements.

        """
        if (category == "P1"):
            newElements = self.point_refinement_elements(allNodes)
        elif (category == "P2"):
            newElements = self.edge_refinement_elements(allNodes)
        elif (category == "P3"):
            newElements = self.face_refinement_elements(allNodes)
        elif (category == "P4"):
            newElements = self.full_refinement_elements(allNodes)
        else:
            print("Error: category (" + category + ") not allowed")
            newElements = {}
        return newElements

    def __P1_refinement__(self,nodes_to_be_refined, nodes):

        idx_refined_node = list(nodes[:, 0]).index(nodes_to_be_refined[0])
        coords_point_of_refinement = nodes[idx_refined_node, 1:]
        [xmin, ymin, zmin] = [min(nodes[:, 1]), min(nodes[:, 2]), min(nodes[:, 3])]
        [xmax, ymax, zmax] = [max(nodes[:, 1]), max(nodes[:, 2]), max(nodes[:, 3])]
        nodes_times_R = np.zeros([8, 4])

        R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        if coords_point_of_refinement[2] != zmin:
            R = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
            for c, n in enumerate(nodes):
                for i in range(3):
                    for j in range(3):
                        coords = n[1:]
                        nodes_times_R[c, 0] = n[0]
                        nodes_times_R[c, i + 1] += R[i][j] * coords[j]
        else:
            nodes_times_R = nodes

        [xmin, ymin, zmin] = [min(nodes_times_R[:, 1]), min(nodes_times_R[:, 2]), min(nodes_times_R[:, 3])]
        [xmax, ymax, zmax] = [max(nodes_times_R[:, 1]), max(nodes_times_R[:, 2]), max(nodes_times_R[:, 3])]
        coords_point_of_refinement = nodes_times_R[idx_refined_node, 1:]
        [x_2_node, y_2_node] = coords_point_of_refinement[0:2]
        [x_1_node, y_3_node] = [0, 0]
        x_1_node = xmax if x_2_node == xmin else xmin
        y_3_node = ymax if y_2_node == ymin else ymin

        rotatedNodes = np.zeros([8, 4])
        for n in nodes_times_R:
            if n[3] == coords_point_of_refinement[2]:
                if n[0] == nodes_to_be_refined[0]:
                    rotatedNodes[1, :] = n
                elif (n[1] == x_1_node) and (n[2] == y_2_node):
                    rotatedNodes[0, :] = n
                elif (n[1] == x_2_node) and (n[2] == y_3_node):
                    rotatedNodes[2, :] = n
                elif (n[1] == x_1_node) and (n[2] == y_3_node):
                    rotatedNodes[3, :] = n
            else:
                if (n[1] == x_1_node) and (n[2] == y_2_node):
                    rotatedNodes[4, :] = n
                elif (n[1] == x_2_node) and (n[2] == y_2_node):
                    rotatedNodes[5, :] = n
                elif (n[1] == x_2_node) and (n[2] == y_3_node):
                    rotatedNodes[6, :] = n
                elif (n[1] == x_1_node) and (n[2] == y_3_node):
                    rotatedNodes[7, :] = n

        # print("Rotated Nodes")
        # print(rotatedNodes)
        # print(rotatedNodes[0,1]-rotatedNodes[1,1])
        # print(rotatedNodes[1,2]-rotatedNodes[2,2])
        allNodes_unrot = self.refinement_node_map("P1", rotatedNodes)

        R_inv = np.linalg.inv(R)
        allNodes = np.zeros([allNodes_unrot.shape[0], 4])
        for c, n in enumerate(allNodes_unrot):
            for i in range(3):
                for j in range(3):
                    coords = n[1:]
                    allNodes[c, 0] = n[0]
                    allNodes[c, i + 1] += R_inv[i][j] * coords[j]

        flip = True
        if ((rotatedNodes[0, 1] - rotatedNodes[1, 1]) > 0 and (rotatedNodes[1, 2] - rotatedNodes[2, 2]) > 0):
            flip = False
        elif ((rotatedNodes[0, 1] - rotatedNodes[1, 1]) < 0 and (rotatedNodes[1, 2] - rotatedNodes[2, 2]) < 0):
            flip = False

        allElements = self.refinement_element_map("P1", allNodes)

        if flip:
            for i in allElements:
                e_old = list(i[1:])
                for j in range(2):
                    for pos, k in enumerate(range(3, -1, -1)):
                        section_pos = pos + (4 * j)
                        section_k = k + (4 * j)
                        i[section_pos + 1] = e_old[section_k]
        return allNodes, allElements

    def __P2_refinement__(self, nodes_to_be_refined, nodes):
        direction_of_edge = np.array(self.nodeMap[nodes_to_be_refined[0]]) - np.array(self.nodeMap[nodes_to_be_refined[1]])
        direction_of_edge = [x / np.linalg.norm(direction_of_edge) for x in direction_of_edge]
        # print("Original edge direction: " + ",".join([str(x) for x in direction_of_edge]))
        # print("nodes")
        # print(nodes)
        # Rotate nodes to match pattern
        R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        if abs(direction_of_edge[1]) == 1:
            R = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
        elif abs(direction_of_edge[2]) == 1:
            R = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]

        nodes_times_R = np.zeros([8, 4])
        for c, n in enumerate(nodes):
            for i in range(3):
                for j in range(3):
                    coords = n[1:]
                    nodes_times_R[c, 0] = n[0]
                    nodes_times_R[c, i + 1] += R[i][j] * coords[j]

        # print("nodes_times_R")
        # print(nodes_times_R)

        # Re-order to match pattern
        # [xmin,ymin,zmin] = [min(nodes[:,1]), min(nodes[:,2]), min(nodes[:,3])]
        # [xmax,ymax,zmax] = [max(nodes[:,1]), max(nodes[:,2]), max(nodes[:,3])]

        nodea_idx = list(nodes_times_R[:, 0]).index(nodes_to_be_refined[0])
        nodeb_idx = list(nodes_times_R[:, 0]).index(nodes_to_be_refined[1])
        nodea = nodes_times_R[nodea_idx][1:]
        nodeb = nodes_times_R[nodeb_idx][1:]
        node1_idx = 0
        node2_idx = 0
        dim = 0
        for k in range(3):
            if (round(nodea[k], 2) < round(nodeb[k], 2)):
                node1_idx = nodea_idx
                node2_idx = nodeb_idx
                dim = k
            elif (round(nodeb[k], 2) < round(nodea[k], 2)):
                node1_idx = nodeb_idx
                node2_idx = nodea_idx
                dim = k

        if dim == 0:
            perp_dim = 1
        elif dim == 2:
            print("Error. Line 240")
            print(direction_of_edge)
            print(nodes_to_be_refined)
        else:
            perp_dim = 0;

        rotatedNodes = np.zeros([8, 4])
        node1_coords = nodes_times_R[node1_idx, 1:]
        node2_coords = nodes_times_R[node2_idx, 1:]
        rotatedNodes[0, :] = nodes_times_R[node1_idx, :]
        rotatedNodes[1, :] = nodes_times_R[node2_idx, :]
        count = 0
        for n in nodes_times_R:
            if not nodes_to_be_refined.count(n[0]):
                if n[3] == rotatedNodes[0][3]:
                    if n[dim + 1] == node2_coords[dim]:
                        rotatedNodes[2, :] = n
                        count += 1
                    elif n[dim + 1] == node1_coords[dim]:
                        rotatedNodes[3, :] = n
                        count += 1
            if count == 2:
                break
        perp_dim_value = rotatedNodes[3, perp_dim + 1]
        for n in nodes_times_R:
            if n[3] != rotatedNodes[0][3]:
                if (n[perp_dim + 1] == node1_coords[perp_dim]) and (n[dim + 1] == node1_coords[dim]):
                    rotatedNodes[4, :] = n
                elif (n[perp_dim + 1] == node2_coords[perp_dim]) and (n[dim + 1] == node2_coords[dim]):
                    rotatedNodes[5, :] = n
                elif (n[dim + 1] == node2_coords[dim]) and (n[perp_dim + 1] == perp_dim_value):
                    rotatedNodes[6, :] = n
                elif (n[dim + 1] == node1_coords[dim]) and (n[perp_dim + 1] == perp_dim_value):
                    rotatedNodes[7, :] = n

                    # print("rotatedNodes")
        # print(rotatedNodes)
        flip = True
        if ((rotatedNodes[0, 3] - rotatedNodes[4, 3]) > 0 and (
                rotatedNodes[0, perp_dim + 1] - rotatedNodes[3, perp_dim + 1]) > 0):
            flip = False
        elif ((rotatedNodes[0, 3] - rotatedNodes[4, 3]) < 0 and (
                rotatedNodes[0, perp_dim + 1] - rotatedNodes[3, perp_dim + 1]) < 0):
            flip = False
        allNodes_unrot = self.refinement_node_map("P2", rotatedNodes, dim)
        R_inv = np.linalg.inv(R)
        allNodes = np.zeros([allNodes_unrot.shape[0], 4])
        for c, n in enumerate(allNodes_unrot):
            for i in range(3):
                for j in range(3):
                    coords = n[1:]
                    allNodes[c, 0] = n[0]
                    allNodes[c, i + 1] += R_inv[i][j] * coords[j]

        allElements = self.refinement_element_map("P2", allNodes)
        if flip:
            for i in allElements:
                e_old = list(i[1:])
                for j in range(2):
                    for pos, k in enumerate(range(3, -1, -1)):
                        section_pos = pos + (4 * j)
                        section_k = k + (4 * j)
                        i[section_pos + 1] = e_old[section_k]
        return allNodes, allElements

    def __P3_refinement__(self, nodes_to_be_refined, nodes):
        v1v2 = np.array(self.nodeMap[nodes_to_be_refined[1]]) - np.array(self.nodeMap[nodes_to_be_refined[0]])
        v1v4 = np.array(self.nodeMap[nodes_to_be_refined[2]]) - np.array(self.nodeMap[nodes_to_be_refined[0]])
        normal = np.cross(v1v2, v1v4)
        normal = [round(x) for x in (normal / np.linalg.norm(normal))]

        # check if inward or outward
        newPoint = [0, 0, 0]
        for i in range(3):
            newPoint[i] = self.nodeMap[nodes_to_be_refined[1]][i] + normal[i]

        [xmin, ymin, zmin] = [min(nodes[:, 1]), min(nodes[:, 2]), min(nodes[:, 3])]
        [xmax, ymax, zmax] = [max(nodes[:, 1]), max(nodes[:, 2]), max(nodes[:, 3])]
        if ((newPoint[0] > xmin and newPoint[0] < xmax) or
                (newPoint[1] > ymin and newPoint[1] < ymax) or
                (newPoint[2] > zmin and newPoint[2] < zmax)):
            for i in range(3):
                normal[i] = normal[i] * -1

        # Rotate nodes to match pattern
        R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        if abs(normal[0]) == 1:
            R = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
        elif abs(normal[1]) == 1:
            R = [[1, 0, 0], [0, 0, 1], [0, -1, 0]]

        nodes_times_R = np.zeros([8, 4])
        for c, n in enumerate(nodes):
            for i in range(3):
                for j in range(3):
                    coords = n[1:]
                    nodes_times_R[c, 0] = n[0]
                    nodes_times_R[c, i + 1] += R[i][j] * coords[j]

        # Re-order to match pattern
        [xmin, ymin, zmin] = [min(nodes_times_R[:, 1]), min(nodes_times_R[:, 2]), min(nodes_times_R[:, 3])]
        [xmax, ymax, zmax] = [max(nodes_times_R[:, 1]), max(nodes_times_R[:, 2]), max(nodes_times_R[:, 3])]
        node1_idx = list(nodes_times_R[:, 0]).index(nodes_to_be_refined[0])
        topZ_value = nodes_times_R[node1_idx, 3]
        rotatedNodes = np.zeros([8, 4])
        for n in nodes_times_R:
            if n[3] == topZ_value:
                if np.all(n[1:3] == [xmin, ymin]):
                    rotatedNodes[0, :] = n
                elif np.all(n[1:3] == [xmax, ymin]):
                    rotatedNodes[1, :] = n
                elif np.all(n[1:3] == [xmax, ymax]):
                    rotatedNodes[2, :] = n
                elif np.all(n[1:3] == [xmin, ymax]):
                    rotatedNodes[3, :] = n
            else:
                if np.all(n[1:3] == [xmin, ymin]):
                    rotatedNodes[4, :] = n
                elif np.all(n[1:3] == [xmax, ymin]):
                    rotatedNodes[5, :] = n
                elif np.all(n[1:3] == [xmax, ymax]):
                    rotatedNodes[6, :] = n
                elif np.all(n[1:3] == [xmin, ymax]):
                    rotatedNodes[7, :] = n

        allNodes_unrot = self.refinement_node_map("P3", rotatedNodes)
        allNodes = np.zeros([allNodes_unrot.shape[0], 4])
        # Rotate back to original locations
        R_inv = np.linalg.inv(R)
        for c, n in enumerate(allNodes_unrot):
            for i in range(3):
                for j in range(3):
                    coords = n[1:]
                    allNodes[c, 0] = n[0]
                    allNodes[c, i + 1] += R_inv[i][j] * coords[j]

                    # Element ICAs
        allElements = self.refinement_element_map("P3", allNodes)

        if (allNodes_unrot[0, 3] - allNodes_unrot[4, 3]) > 0:
            # print("Flipping node order")
            for i in allElements:
                e_old = list(i[1:])
                for j in range(2):
                    for pos, k in enumerate(range(3, -1, -1)):
                        section_pos = pos + (4 * j)
                        section_k = k + (4 * j)
                        i[section_pos + 1] = e_old[section_k]
        return allNodes, allElements

    def __P4_refinement__(self, nodes):

        [e1, e2, e3] = [1, 2, 0]
        [rev_e1, rev_e2, rev_e3] = [2, 0, 1]
        for c, n in enumerate(nodes):
            coords = n[1:]
            nodes[c, 1:] = [coords[e1], coords[e2], coords[e3]]

        # Re-order to match pattern
        [xmin, ymin, zmin] = [min(nodes[:, 1]), min(nodes[:, 2]), min(nodes[:, 3])]
        [xmax, ymax, zmax] = [max(nodes[:, 1]), max(nodes[:, 2]), max(nodes[:, 3])]
        rotatedNodes = np.zeros([8, 4])
        for n in nodes:
            if n[3] == zmin:
                if np.all(n[1:3] == [xmin, ymin]):
                    rotatedNodes[0, :] = n
                elif np.all(n[1:3] == [xmax, ymin]):
                    rotatedNodes[1, :] = n
                elif np.all(n[1:3] == [xmax, ymax]):
                    rotatedNodes[2, :] = n
                elif np.all(n[1:3] == [xmin, ymax]):
                    rotatedNodes[3, :] = n
            else:
                if np.all(n[1:3] == [xmin, ymin]):
                    rotatedNodes[4, :] = n
                elif np.all(n[1:3] == [xmax, ymin]):
                    rotatedNodes[5, :] = n
                elif np.all(n[1:3] == [xmax, ymax]):
                    rotatedNodes[6, :] = n
                elif np.all(n[1:3] == [xmin, ymax]):
                    rotatedNodes[7, :] = n

        allNodes = self.refinement_node_map("P4", rotatedNodes)

        allElements = self.refinement_element_map("P4", allNodes)

        for c, n in enumerate(allNodes):
            coords = n[1:]
            allNodes[c, 1:] = [coords[rev_e1], coords[rev_e2], coords[rev_e3]]

        return allNodes, allElements

    def define_refinement_area(self, bounds_to_refine):
        """
        Creates a list of nodes that are within the bounds given. Bounds indicate area that needs refinement

        Parameters
        ----------
        bounds_to_refine : array[6]
            [xmin, x max,ymin ,ymax, zmin, zmax] min and max cooordinates defining square within which to refine model
        nodeMap : Map(int,array)
            Map of node numbers (keys) to coordinates of nodes (values)

        Returns
        -------
        refined_nodes : Array
            List of node numbers contained in refinement area

        """
        print("Defining nodes within bounds")
        refined_nodes = []
        for node_num, n_coords in self.nodeMap.items():
            if ((n_coords[0] > bounds_to_refine[0]) and (n_coords[0] < bounds_to_refine[1]) and
                    (n_coords[1] > bounds_to_refine[2]) and (n_coords[1] < bounds_to_refine[3]) and
                    (n_coords[2] > bounds_to_refine[4]) and (n_coords[2] < bounds_to_refine[5])):
                refined_nodes.append(node_num)

        print("Number of nodes in refinement bounds: " + str(len(refined_nodes)))
        return refined_nodes

    def assign_refinement_type(self, refined_nodes):
        """
        Categorizing elements that are connected to any of the refined nodes based on four categories:
            P1: Point refinement, 1 node connected
            P2: Edge refinement, 2 nodes connected
            P3: Face refinement, 4 nodes connected
            P4: Full refined, all 8 nodes of element are within refined nodes area

        Parameters
        ----------
        refined_nodes : Array
            List of node numbers contained in refinement area

        Returns
        -------
        sideMap :  Map(string,array)
            Map of of categories of refinement (keys) to elements assigned with that refinement

        """
        elementMap = self.mesh.elements
        elementICAMap = mu.create_elements_ica_map(elementMap)
        print("Assigning refinement type based number of nodes involved")
        node_to_element_map = mu.create_node_to_elem_map(elementICAMap)
        # nodes_to_be_refined_map = {}
        new_P1 = []
        new_P2 = []
        new_P3 = []
        new_P4 = []
        element_parsed = []
        for n in refined_nodes:
            for e in node_to_element_map[n]:
                if not element_parsed.count(e):
                    element_parsed.append(e)
                    ica_old_elements = elementICAMap[e]
                    # for e, ica_old_elements in elementMap.items():
                    nodes_to_be_refined = []
                    for n in ica_old_elements:
                        if refined_nodes.count(n):
                            nodes_to_be_refined.append(n)
                            if len(nodes_to_be_refined) > 4:
                                new_P4.append(e)
                                break
                    if len(nodes_to_be_refined) == 0:
                        continue
                    elif len(nodes_to_be_refined) == 4:
                        new_P3.append(e)
                    elif len(nodes_to_be_refined) == 2:
                        new_P2.append(e)
                    elif len(nodes_to_be_refined) == 1:
                        new_P1.append(e)
            # nodes_to_be_refined_map[e] = nodes_to_be_refined
        sideMap = {}
        refine_elements = len(new_P1) + len(new_P2) + len(new_P3) + len(new_P4)
        sideMap["P1"] = new_P1
        sideMap["P2"] = new_P2
        sideMap["P3"] = new_P3
        sideMap["P4"] = new_P4
        print(str(refine_elements) + " elements to be refined")
        return sideMap

    def face_refinement_nodes(self, nodes):
        # print("face_refinement selected")
        newNodes = np.zeros([48, 4])
        for i in range(8):
            newNodes[i, :] = nodes[i, :]

        #  Upper surface nodes
        newNodes[8, 1:] = nodes[0, 1:] - (nodes[0, 1:] - nodes[1, 1:]) * (1 / 3)
        newNodes[9, 1:] = nodes[0, 1:] - (nodes[0, 1:] - nodes[1, 1:]) * (2 / 3)
        newNodes[10, 1:] = nodes[2, 1:] - (nodes[2, 1:] - nodes[1, 1:]) * (1 / 3)
        newNodes[11, 1:] = nodes[2, 1:] - (nodes[2, 1:] - nodes[1, 1:]) * (2 / 3)
        newNodes[12, 1:] = nodes[3, 1:] - (nodes[3, 1:] - nodes[2, 1:]) * (1 / 3)
        newNodes[13, 1:] = nodes[3, 1:] - (nodes[3, 1:] - nodes[2, 1:]) * (2 / 3)
        newNodes[14, 1:] = nodes[3, 1:] - (nodes[3, 1:] - nodes[0, 1:]) * (1 / 3)
        newNodes[15, 1:] = nodes[3, 1:] - (nodes[3, 1:] - nodes[0, 1:]) * (2 / 3)

        # Central uppersurface nodes
        upperSurface = nodes[1, 3];
        newNodes[16, 1:] = [newNodes[8, 1], newNodes[15, 2], upperSurface]
        newNodes[17, 1:] = [newNodes[8, 1], newNodes[14, 2], upperSurface]
        newNodes[18, 1:] = [newNodes[9, 1], newNodes[14, 2], upperSurface]
        newNodes[19, 1:] = [newNodes[9, 1], newNodes[15, 2], upperSurface]

        # Level: 1/3
        oneThird = nodes[0, 3] - (nodes[0, 3] - nodes[5, 3]) * (1 / 3)
        for i in range(20, 24):
            correspNode = i - 20
            newNodes[i, :] = [0, newNodes[correspNode, 1], newNodes[correspNode, 2], oneThird]
        for i in range(28, 40):
            correspNode = i - 20
            newNodes[i, :] = [0, newNodes[correspNode, 1], newNodes[correspNode, 2], oneThird]

        # Level: 1/2
        halfway = nodes[0, 3] - (nodes[0, 3] - nodes[5, 3]) * (1 / 2)
        newNodes[24, 1:] = [newNodes[36, 1], newNodes[36, 2], halfway]
        newNodes[25, 1:] = [newNodes[37, 1], newNodes[37, 2], halfway]
        newNodes[26, 1:] = [newNodes[38, 1], newNodes[38, 2], halfway]
        newNodes[27, 1:] = [newNodes[39, 1], newNodes[39, 2], halfway]

        # Level: 2/3
        twoThird = nodes[0, 3] - (nodes[0, 3] - nodes[5, 3]) * (2 / 3)
        newNodes[40, 1:] = [newNodes[28, 1], newNodes[28, 2], twoThird]
        newNodes[41, 1:] = [newNodes[32, 1], newNodes[32, 2], twoThird]
        newNodes[42, 1:] = [newNodes[33, 1], newNodes[33, 2], twoThird]
        newNodes[43, 1:] = [newNodes[29, 1], newNodes[29, 2], twoThird]

        newNodes[44, 1:] = [newNodes[1, 1], newNodes[11, 2], twoThird]
        newNodes[45, 1:] = [newNodes[1, 1], newNodes[10, 2], twoThird]
        newNodes[46, 1:] = [newNodes[0, 1], newNodes[15, 2], twoThird]
        newNodes[47, 1:] = [newNodes[0, 1], newNodes[14, 2], twoThird]

        for i in range(8, 48):
            newNodes[i, 0] = self.newNodeNum + i + 1

        return newNodes;

    def convertToNodeNums(self, nodes, baseNodeNums):
        """
        Converts the single element node numbers  (basenodeNumbers) to the
        global node numbers in the mesh given by (nodes) for element ica


        Parameters
        ----------
        nodes : array(Array)
            Array of nodes in mesh.
        baseNodeNums : Array(ints)
            Array of single element node numbers.

        Returns
        -------
        ica : array(ints)
            GLobal ica for element.

        """
        ica = np.zeros(8, 'int32')
        pos = 0;
        for i in baseNodeNums:
            ica[pos] = int(nodes[i - 1, 0])
            pos += 1
        return ica
    
    def face_refinement_elements(self, nodes):
        elements = np.zeros([22, 9], 'int32')
        for i in range(22):
            elements[i, 0] = self.newElemNum + i + 1

        elements[0, 1:] = self.convertToNodeNums(nodes, [1, 9, 17, 16, 21, 29, 37, 36])
        elements[1, 1:] = self.convertToNodeNums(nodes, [9, 10, 20, 17, 29, 30, 40, 37])
        elements[2, 1:] = self.convertToNodeNums(nodes, [10, 2, 12, 20, 30, 22, 32, 40])
        elements[3, 1:] = self.convertToNodeNums(nodes, [16, 17, 18, 15, 36, 37, 38, 35])
        elements[4, 1:] = self.convertToNodeNums(nodes, [17, 20, 19, 18, 37, 40, 39, 38])
        elements[5, 1:] = self.convertToNodeNums(nodes, [20, 12, 11, 19, 40, 32, 31, 39])
        elements[6, 1:] = self.convertToNodeNums(nodes, [15, 18, 13, 4, 35, 38, 33, 24])
        elements[7, 1:] = self.convertToNodeNums(nodes, [18, 19, 14, 13, 38, 39, 34, 33])
        elements[8, 1:] = self.convertToNodeNums(nodes, [19, 11, 3, 14, 39, 31, 23, 34])
        elements[9, 1:] = self.convertToNodeNums(nodes, [37, 40, 39, 38, 25, 28, 27, 26])
        elements[10, 1:] = self.convertToNodeNums(nodes, [25, 28, 27, 26, 41, 44, 43, 42])
        elements[11, 1:] = self.convertToNodeNums(nodes, [41, 44, 43, 42, 5, 6, 7, 8])
        elements[12, 1:] = self.convertToNodeNums(nodes, [29, 30, 40, 37, 41, 44, 28, 25])
        elements[13, 1:] = self.convertToNodeNums(nodes, [40, 32, 31, 39, 28, 45, 46, 27])
        elements[14, 1:] = self.convertToNodeNums(nodes, [38, 39, 34, 33, 26, 27, 43, 42])
        elements[15, 1:] = self.convertToNodeNums(nodes, [36, 37, 38, 35, 47, 25, 26, 48])
        elements[16, 1:] = self.convertToNodeNums(nodes, [21, 29, 37, 36, 5, 41, 25, 47])
        elements[17, 1:] = self.convertToNodeNums(nodes, [30, 22, 32, 40, 44, 6, 45, 28])
        elements[18, 1:] = self.convertToNodeNums(nodes, [39, 31, 23, 34, 27, 46, 7, 43])
        elements[19, 1:] = self.convertToNodeNums(nodes, [35, 38, 33, 24, 48, 26, 42, 8])
        elements[20, 1:] = self.convertToNodeNums(nodes, [47, 25, 26, 48, 5, 41, 42, 8])
        elements[21, 1:] = self.convertToNodeNums(nodes, [28, 45, 46, 27, 44, 6, 7, 43])

        return elements

    def edge_refinement_nodes(self, nodes, dim):
        newNodes = np.zeros([28, 4])
        for i in range(8, newNodes.shape[0]):
            newNodes[i, 0] = self.newNodeNum + i + 1

        for i in range(8):
            newNodes[i, :] = nodes[i, :]

        if dim == 1:
            perp_dim = 0
        else:
            perp_dim = 1

        ## Top level:
        top = nodes[0, 3]

        newNodes[8, 1:] = nodes[0, 1:] - (nodes[0, 1:] - nodes[1, 1:]) * (1 / 3)
        newNodes[9, 1:] = nodes[0, 1:] - (nodes[0, 1:] - nodes[1, 1:]) * (2 / 3)

        newNodes[14, 1:] = nodes[0, 1:] - (nodes[0, 1:] - nodes[3, 1:]) * (1 / 3)
        newNodes[15, 1:] = nodes[1, 1:] - (nodes[1, 1:] - nodes[2, 1:]) * (1 / 3)

        newNodes[10, 3] = top
        newNodes[10, dim + 1] = newNodes[8, dim + 1]
        newNodes[10, perp_dim + 1] = newNodes[14, perp_dim + 1]

        newNodes[11, 3] = top
        newNodes[11, dim + 1] = newNodes[9, dim + 1]
        newNodes[11, perp_dim + 1] = newNodes[14, perp_dim + 1]

        twoThirdY = nodes[0, :] - (nodes[0, :] - nodes[3, :]) * (2 / 3)
        newNodes[12, 3] = top
        newNodes[12, dim + 1] = newNodes[8, dim + 1]
        newNodes[12, perp_dim + 1] = twoThirdY[perp_dim + 1]
        newNodes[13, 3] = top
        newNodes[13, dim + 1] = newNodes[9, dim + 1]
        newNodes[13, perp_dim + 1] = twoThirdY[perp_dim + 1]

        ## 1/3 level
        oneThird = nodes[0, 3] - (nodes[0, 3] - nodes[4, 3]) * 1 / 3
        newNodes[16, 1:] = [newNodes[0, 1], newNodes[0, 2], oneThird]
        newNodes[17, 1:] = [newNodes[8, 1], newNodes[8, 2], oneThird]
        newNodes[18, 1:] = [newNodes[9, 1], newNodes[9, 2], oneThird]
        newNodes[19, 1:] = [newNodes[1, 1], newNodes[1, 2], oneThird]
        newNodes[20, 1:] = [newNodes[14, 1], newNodes[14, 2], oneThird]
        newNodes[21, 1:] = [newNodes[10, 1], newNodes[10, 2], oneThird]
        newNodes[22, 1:] = [newNodes[11, 1], newNodes[11, 2], oneThird]
        newNodes[23, 1:] = [newNodes[15, 1], newNodes[15, 2], oneThird]

        ## 2/3 level
        twoThird = nodes[0, 3] - (nodes[0, 3] - nodes[4, 3]) * 2 / 3
        newNodes[24, 1:] = [newNodes[8, 1], newNodes[8, 2], twoThird]
        newNodes[25, 1:] = [newNodes[9, 1], newNodes[9, 2], twoThird]
        newNodes[26, 1:] = [newNodes[12, 1], newNodes[12, 2], twoThird]
        newNodes[27, 1:] = [newNodes[13, 1], newNodes[13, 2], twoThird]

        return newNodes

    def edge_refinement_elements(self, nodes):

        elements = np.zeros([11, 9], 'int32')
        for i in range(elements.shape[0]):
            elements[i, 0] = self.newElemNum + i + 1;

        elements[0, 1:] = self.convertToNodeNums(nodes, [1, 9, 11, 15, 17, 18, 22, 21])
        elements[1, 1:] = self.convertToNodeNums(nodes, [9, 10, 12, 11, 18, 19, 23, 22])
        elements[2, 1:] = self.convertToNodeNums(nodes, [10, 2, 16, 12, 19, 20, 24, 23])
        elements[3, 1:] = self.convertToNodeNums(nodes, [15, 11, 13, 4, 21, 22, 27, 8])
        elements[4, 1:] = self.convertToNodeNums(nodes, [11, 12, 14, 13, 22, 23, 28, 27])
        elements[5, 1:] = self.convertToNodeNums(nodes, [12, 16, 3, 14, 23, 24, 7, 28])
        elements[6, 1:] = self.convertToNodeNums(nodes, [13, 14, 3, 4, 27, 28, 7, 8])
        elements[7, 1:] = self.convertToNodeNums(nodes, [17, 18, 22, 21, 5, 25, 27, 8])
        elements[8, 1:] = self.convertToNodeNums(nodes, [18, 19, 23, 22, 25, 26, 28, 27])
        elements[9, 1:] = self.convertToNodeNums(nodes, [19, 20, 24, 23, 26, 6, 7, 28])
        elements[10, 1:] = self.convertToNodeNums(nodes, [25, 26, 28, 27, 5, 6, 7, 8])
        return elements

    def point_refinement_nodes(self, nodes):

        newNodes = np.zeros([15, 4])
        for i in range(8):
            newNodes[i, :] = nodes[i, :]
        for i in range(8, newNodes.shape[0]):
            newNodes[i, 0] = self.newNodeNum + i + 1

        oneThirdX = nodes[1, 1] - (nodes[1, 1] - nodes[0, 1]) * 1 / 3
        oneThirdY = nodes[1, 2] - (nodes[1, 2] - nodes[2, 2]) * 1 / 3
        oneThirdZ = nodes[1, 3] - (nodes[1, 3] - nodes[5, 3]) * 1 / 3

        newNodes[8, 1:] = [oneThirdX, nodes[1, 2], nodes[1, 3]]
        newNodes[9, 1:] = [nodes[1, 1], oneThirdY, nodes[1, 3]]
        newNodes[10, 1:] = [oneThirdX, oneThirdY, nodes[1, 3]]
        newNodes[11, 1:] = [oneThirdX, nodes[1, 2], oneThirdZ]
        newNodes[12, 1:] = [nodes[1, 1], nodes[1, 2], oneThirdZ]
        newNodes[13, 1:] = [nodes[1, 1], oneThirdY, oneThirdZ]
        newNodes[14, 1:] = [oneThirdX, oneThirdY, oneThirdZ]

        return newNodes

    def point_refinement_elements(self, nodes):
        elements = np.zeros([4, 9], 'int32')
        for i in range(elements.shape[0]):
            elements[i, 0] = self.newElemNum + i + 1
        elements[0, 1:] = self.convertToNodeNums(nodes, [9, 2, 10, 11, 12, 13, 14, 15])
        elements[1, 1:] = self.convertToNodeNums(nodes, [1, 9, 11, 4, 5, 12, 15, 8])
        elements[2, 1:] = self.convertToNodeNums(nodes, [11, 10, 3, 4, 15, 14, 7, 8])
        elements[3, 1:] = self.convertToNodeNums(nodes, [12, 13, 14, 15, 5, 6, 7, 8])
        return elements

    def full_refinement_nodes(self, nodes):
        # print("Full refinement selected")
        newNodes = np.zeros([64, 4])
        top = nodes[0, 3]
        oneThird = nodes[0, 3] - (nodes[0, 3] - nodes[4, 3]) * 1 / 3
        twoThirds = nodes[0, 3] - (nodes[0, 3] - nodes[4, 3]) * 2 / 3
        bottom = nodes[4, 3]
        zvalues = [top, oneThird, twoThirds, bottom]

        y = nodes[0, 2]
        oneThirdy = nodes[0, 2] - (nodes[0, 2] - nodes[3, 2]) * 1 / 3
        twoThirdsy = nodes[0, 2] - (nodes[0, 2] - nodes[3, 2]) * 2 / 3
        yfull = nodes[3, 2]
        yvalues = [y, oneThirdy, twoThirdsy, yfull]

        x = nodes[0, 1]
        oneThirdx = nodes[0, 1] - (nodes[0, 1] - nodes[1, 1]) * 1 / 3
        twoThirdsx = nodes[0, 1] - (nodes[0, 1] - nodes[1, 1]) * 2 / 3
        xfull = nodes[1, 1]
        xvalues = [x, oneThirdx, twoThirdsx, xfull]

        nodeCount = 0;
        for i in range(len(zvalues)):
            for j in range(len(yvalues)):
                for k in range(len(xvalues)):
                    newNodes[nodeCount, :] = [self.newNodeNum + nodeCount + 1, xvalues[k], yvalues[j], zvalues[i]]
                    nodeCount += 1

        nodeReordering = {1: 1, 2: 4, 3: 16, 4: 13, 5: 49, 6: 52, 7: 64, 8: 61}
        for pos, nodeNum in nodeReordering.items():
            newNodes[nodeNum - 1, :] = nodes[pos - 1, :]
        return newNodes

    def full_refinement_elements(self, nodes):
        elements = np.zeros([27, 9], 'int32')
        for i in range(elements.shape[0]):
            elements[i, 0] = self.newElemNum + i + 1
        count = 0
        for k in range(3):
            for l in range(3):
                for i in range(3):
                    j = i + (k * 16) + (4 * l) + 1
                    elements[count, 1:] = self.convertToNodeNums(nodes, [j, j + 1, j + 5, j + 4, j + 16, j + 17, j + 21,
                                                                    j + 20, ])
                    count += 1
        return elements





