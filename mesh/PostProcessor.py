from abc import ABC, abstractmethod
from mesh.refinement import Refiner
from mesh.Element import QuadElement
import numpy as np


class IPostProcessor(ABC):

    def __init__(self, config, mesh):
        self.config = config
        self.mesh = mesh

    @abstractmethod
    def post_process(self):
        raise NotImplementedError

class PostProcessor(IPostProcessor):

    def post_process(self):
        return


# class PostProcessorDebug(IPostProcessor):

    # def post_process(self):

# class PostProcessorAtrophy(IPostProcessor):

    # def post_process(self):

# class PostProcessorLesion(IPostProcessor):

    # def post_process(self):


class PostProcessorDecorator(IPostProcessor):
    _post_processor: IPostProcessor = None

    def __init__(self, post_processor: IPostProcessor):
        self._post_processor = post_processor
        super().__init__(post_processor.config, post_processor.mesh)

    @property
    def post_processor(self):
        return self._post_processor

    def post_process(self):
        self._post_processor.post_process()


class SmoothMesh(PostProcessorDecorator):

    def __init__(self, post_processor: IPostProcessor, coeffs, iterations, excluded_regions=None):
        super().__init__(post_processor)
        if excluded_regions is None:
            excluded_regions = []
        self.mesh_refiner = None
        self.excluded_regions = excluded_regions
        self.coeffs = coeffs
        self.iterations = iterations

    def post_process(self):
        self.smooth_mesh()
        return super().post_process()

    def smooth_mesh(self):
        # Smooth outer surface of mesh (including CSF)
        print("########## Smoothing mesh excluding elements with material types: {} ##########"
              .format(",".join([str(x) for x in self.excluded_regions])))
        self.mesh.smooth_mesh(self.coeffs, self.iterations, elementsNotIncluded=self.excluded_regions)


class RefineMesh(PostProcessorDecorator):

    def __init__(self, post_processor: IPostProcessor, refinement_type,
                 point=None, radius=0, bounds=None, elements=None):
        super().__init__(post_processor)
        if point is None:
            point = []
        if bounds is None:
            bounds = []
        if elements is None:
            elements = []
        self.elements = elements
        self.bounds = bounds
        self.radius = radius
        self.point = point
        self.mesh_refiner = None
        self.type_ref = refinement_type

    def post_process(self):
        self.mesh_refiner = Refiner.Refiner(self.mesh)
        self.refine_mesh()
        return super().post_process()

    def refine_mesh(self):
        if self.type_ref == 'point':
            assert self.point is not None and self.radius > 0
            self.mesh_refiner.refine_around_point(self.point, self.radius)
            self.config.writeToConfig("Refined around point", ",".join([str(x) for x in self.point]))
            self.config.writeToConfig("Refined with radius", str(self.radius))
        elif self.type_ref == 'elements':
            assert self.elements is not None
            self.mesh_refiner.refine_elements(self.elements)
            self.config.writeToConfig("Refined elements", ",".join([str(x) for x in self.elements]))
        elif self.type_ref == 'bounding_box':
            assert self.bounds is not None
            self.mesh_refiner.refine_within_region(self.bounds)
            self.config.writeToConfig("Refined within bounds", ",".join([str(x) for x in self.bounds]))


class CreateBoundaryElements(PostProcessorDecorator):

    def __init__(self, post_processor: IPostProcessor, element_mat_number,
                 boundary_test_fx=None, excluded_regions=None):
        super().__init__(post_processor)
        self.mesh_refiner = None
        self.element_mat_number = element_mat_number
        self.boundary_test_fx = boundary_test_fx
        if excluded_regions is None:
            excluded_regions = []
        self.excluded_regions = excluded_regions

    def post_process(self):
        self.mesh_refiner = Refiner.Refiner(self.mesh)
        self.create_boundary()
        return super().post_process()

    def create_boundary(self):
        # boundary_elements_map = brainCreator.createCSFBoundary(mesh,200, boundaryTest=OnlyOnLabel(mesh, 24))
        # boundary_elements_map = brainCreator.createCSFBoundary(mesh,200, boundaryTest=OpenBottomCSF(mesh))
        # boundary_elements_map = brainCreator.createBoundary(mesh, 200, boundaryTest=None)
        boundary_elements_map = {}
        boundary_number = max(self.mesh.elements.keys()) if len(self.mesh.boundaryElements) == 0 else max(
            self.mesh.boundaryElements.keys())
        boundary_elements = self.mesh.locate_boundary_element_map(elementsNotIncluded=self.excluded_regions)
        print("Locating CSF boundary")
        for compoundKey, ica in boundary_elements.items():
            boundary_number += 1
            ica_nodes = [self.mesh.nodes[n] for n in ica]
            boundary_element = QuadElement(boundary_number, ica_nodes, mat=[self.element_mat_number])
            [element_num, face] = [int(x) for x in compoundKey.split("-")]
            # [xc,yc,zc,m] = e_centroids[element_num]
            boundary = True
            if not self.boundary_test_fx is None:
                boundary = self.boundary_test_fx.validElement(element_num)
            if boundary:
                boundary_elements_map[boundary_number] = boundary_element
            else:
                boundary_number -= 1
        self.mesh.addBoundaryElements(boundary_elements_map)


class ApplyAtrophyConcentration(PostProcessorDecorator):

    def post_process(self):
        self.apply_concentration()
        return super().post_process()

    def apply_concentration(self):
        print("Applying concentration field")
        # Get center of brain stem
        center_bs = self.mesh.get_center_of_region(16)
        # Get center of hippocampus
        center_h = self.mesh.get_center_of_region(17)
        bounding_box_hippo = self.mesh.getBoundingBox(regions=17)
        center = [center_bs[0], center_h[1], bounding_box_hippo[5]]
        distance = max([abs(bounding_box_hippo[0] - center_bs[0]),
                        abs(bounding_box_hippo[3] - center_bs[0]),
                        abs(bounding_box_hippo[1] - center_bs[1]),
                        abs(bounding_box_hippo[4] - center_bs[1]),
                        abs(bounding_box_hippo[2] - center_bs[2]),
                        abs(bounding_box_hippo[5] - center_bs[2])])

        mesh_bounding_box = self.mesh.getBoundingBox()
        mesh_bounding_box = [x for x in mesh_bounding_box[:3]] + [x for x in mesh_bounding_box[3:]]
        max_radius = max([abs(mesh_bounding_box[0] - center[0]),
                          abs(mesh_bounding_box[3] - center[0]),
                          abs(mesh_bounding_box[1] - center[1]),
                          abs(mesh_bounding_box[4] - center[1]),
                          abs(mesh_bounding_box[2] - center[2]),
                          abs(mesh_bounding_box[5] - center[2])])
        import random

        num = random.randrange(500, 1000)
        num /= 1000
        length = max_radius
        concentration_radius = length * num + (distance * 2)

        self.config.writeToConfig("Concentration radius", concentration_radius)
        self.config.writeToConfig("Center", ", ".join([str(x) for x in center]))
        self.mesh.dataToWrite.append("concentration")
        node_touched = []
        for n in self.mesh.nodes.values():
            centroid = n.getCoords()
            radius = np.linalg.norm(np.array(center) - np.array(centroid))
            if radius <= concentration_radius:
                c = 1 - (radius / concentration_radius)
                n.addData("concentration", [round(c, 4)])
            else:
                n.addData("concentration", [0])


class RemoveRegion(PostProcessorDecorator):

    def __init__(self, post_processor: IPostProcessor, region_label):
        super().__init__(post_processor)
        self.region_label = region_label

    def post_process(self):
        self.remove_region()
        return super().post_process()

    def remove_region(self):
        self.mesh.remove_region(self.region_label)
