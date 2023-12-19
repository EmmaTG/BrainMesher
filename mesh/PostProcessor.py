import warnings
from abc import ABC, abstractmethod

from mesh import BoundaryFunctions
from mesh.refinement import Refiner
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

    def __init__(self, post_processor: IPostProcessor, refinement_type):
        super().__init__(post_processor)
        self.mesh_refiner = None
        self.type_ref = refinement_type

    def post_process(self):
        self.mesh_refiner = Refiner.Refiner(self.mesh)
        self.refine_mesh()
        return super().post_process()

    def refine_mesh(self):
        print("########## Refining mesh using {} method ##########".format(self.type_ref))
        if self.type_ref == 'point':
            centers = self.config.get('refine.centers')
            radii = self.config.get('refine.radii')
            assert len(centers) == len(radii)
            for count in range(len(centers)):
                point = centers[count]
                radius = radii[count]
                assert point is not None and radius > 0
                self.mesh_refiner.refine_around_point(point, radius)
                self.config.write_to_config("Refined around point", ",".join([str(x) for x in point]))
                self.config.write_to_config("Refined with radius", str(radius))
        elif self.type_ref == 'elements':
            elements = self.config.get('refine.element_numbers')
            assert len(elements)>0
            self.mesh_refiner.refine_elements(elements)
            self.config.write_to_config("Refined elements", ",".join([str(x) for x in elements]))
        elif self.type_ref == 'refine.bounding_box':
            for b_box in self.config.get('bounds'):
                assert (len(b_box) == 6
                        and b_box[0] > b_box[3]
                        and b_box[1] > b_box[4]
                        and b_box[2] > b_box[5]), \
                    "Bounds must given in the format: [xmin,ymin,zmin,xmax,ymax,zmax]"
                self.mesh_refiner.refine_within_region(b_box)
                self.config.write_to_config("Refined within bounds", ",".join([str(x) for x in b_box]))


class CreateBoundaryElements(PostProcessorDecorator):

    def __init__(self, post_processor: IPostProcessor, element_mat_number,
                 boundary_test_fx='none', excluded_regions=None):
        super().__init__(post_processor)
        self.mesh_refiner = None
        self.element_mat_number = element_mat_number
        self.boundary_test = boundary_test_fx
        excluded_regions_list = []
        if excluded_regions is not None:
            all_labels = self.config.MATERIAL_LABELS.get_homogenized_labels_map()
            if excluded_regions == '':
                all_labels.clear()
            else:
                for regions in excluded_regions.split(","):
                    all_labels.pop(regions.strip().lower())
            excluded_regions_list = list(all_labels.values())
        self.excluded_regions = excluded_regions_list

    def post_process(self):
        self.create_boundary()
        return super().post_process()

    def create_boundary(self):
        boundary_test_fx = None
        if self.boundary_test.lower() != 'none':
            if 'OnlyOnLabel' in self.boundary_test:
                popped_label = self.boundary_test.split("-")[1].strip()
                labels = self.config.MATERIAL_LABELS.get_homogenized_labels_map()
                region_label = labels.pop(popped_label, False)
                if region_label:
                    boundary_test_fx = BoundaryFunctions.OnlyOnLabel(self.mesh, region_label)
                else:
                    boundary_test_fx = None
                # for e in labels.values():
                #     if not self.excluded_regions.count(e):
                #         self.excluded_regions.append(e)
            elif self.boundary_test == 'OpenBottomCSF':
                if not self.config.get('add_csf'):
                    warnings.warn("You cannot request an open open CSf boundary if CSF is not added to the model")
                    boundary_test_fx = None
                else:
                    boundary_test_fx = BoundaryFunctions.OpenBottomCSF(self.mesh)
            elif self.boundary_test == 'ExternalCSF':
                if not self.config.get('add_csf'):
                    warnings.warn("You cannot request an open open CSf boundary if CSF is not added to the model")
                    boundary_test_fx = None
                else:
                    boundary_test_fx = BoundaryFunctions.ExternalCSF(self.mesh)
            else:
                raise NotImplementedError("This boundary function has not been implemented")

        print("########## Creating boundary elements with number {} \n"
              "           excluding regions {} \n"
              "           using boundary test function {} ##########"
              .format(self.element_mat_number,
                      'None' if len(self.excluded_regions) == 0 else", ".join([str(x) for x in self.excluded_regions]),
                      'None' if self.boundary_test is None else self.boundary_test))

        self.mesh.create_boundary(boundary_test_fx, self.element_mat_number, excluded_regions=self.excluded_regions)


class ApplyAtrophyConcentration(PostProcessorDecorator):

    def post_process(self):
        self.apply_concentration()
        return super().post_process()

    def apply_concentration(self):
        print("########## Applying concentration field ##########")
        # Get center of brain stem
        center_bs = self.mesh.get_center_of_region(175)
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

        self.config.write_to_config("Concentration radius", concentration_radius)
        self.config.write_to_config("Center", ", ".join([str(x) for x in center]))
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
        print("########## Removing region {} ##########".format(self.region_label))
        self.mesh.remove_region(self.region_label)
