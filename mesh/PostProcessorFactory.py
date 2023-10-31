from mesh.PostProcessor import PostProcessor, ApplyAtrophyConcentration, SmoothMesh, RemoveRegion, \
    CreateBoundaryElements, RefineMesh


class PostProcessorFactory:

    @staticmethod
    def get_post_processor(mesh, config):
        post_processor = PostProcessor(config, mesh)
        if config.get('atrophy'):
            # Creates a concentration profile in the brain to eb used for degree fo atrophy calculations
            post_processor = ApplyAtrophyConcentration(post_processor)

        if config.get('smooth'):
            # global smoothing
            post_processor = SmoothMesh(post_processor, config.get('co_effs'), config.get('iterations'))

            if config.get('add_csf'):
                labels = config.MATERIAL_LABELS.get_homogenized_labels_map()
                label_for_csf = labels.get("csf", 24)
                post_processor = SmoothMesh(post_processor, config.get('co_effs'), config.get('iterations'),
                                            excluded_regions=[label_for_csf])

            # Smooths mesh as defined in config file
            regions = config.get("smooth_regions")
            region_coeffs = config.get("region_co_effs")
            region_iterations = config.get("region_iterations")
            count = 0
            for r in regions:
                labels = config.MATERIAL_LABELS.get_homogenized_labels_map()
                label_for_region = labels.get(r, -1)
                if label_for_region != -1:
                    excluded_materials = []
                    labels.pop(r)
                    for e in labels.values():
                        if not excluded_materials.count(e):
                            excluded_materials.append(e)
                    post_processor = SmoothMesh(post_processor, region_coeffs[count], region_iterations[count],
                                                excluded_regions=excluded_materials)
                count += 1

        # Removes regions specified at unused in the materials label description
        all_labels = config.MATERIAL_LABELS.get_homogenized_labels_map()
        label_for_unused = all_labels.get("unused")
        if label_for_unused != 0:
            post_processor = RemoveRegion(post_processor, label_for_unused)

        # Creates blank spaces in areas of the ventricles if specified
        all_labels = config.MATERIAL_LABELS.get_homogenized_labels_map()
        label_for_ventricles = all_labels.get("ventricles", all_labels.get("ventricle", -1))
        if label_for_ventricles != -1:
            post_processor = RemoveRegion(post_processor, label_for_ventricles)

        # Creates boundary elements on specified regions
        for count, e_num in enumerate(config.get('boundary_element_numbers')):
            post_processor = CreateBoundaryElements(post_processor, e_num,
                                                    boundary_test_fx=config.get('boundary_tests')[count],
                                                    excluded_regions=config.get('excluded_regions')[count])

        if config.get('refine'):
            # Refines mesh in ways defined in config file
            if config.get('refine.point'):
                post_processor = RefineMesh(post_processor, 'point')
            if config.get('refine.bounding_box'):
                post_processor = RefineMesh(post_processor, 'bounding_box', )
            if config.get('refine.elements'):
                post_processor = RefineMesh(post_processor, 'elements')

        return post_processor