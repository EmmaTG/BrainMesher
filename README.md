# BrainMesher

Create a 3D brain mesh from mri images using only hexahedral elements.

1. Import aseg.mgz from freesurfer output (after performing ,recon-all' on mri images)
2. Convert voxel image to point cloud
3. Optional Add cerebrospinal fluid
4. Convert point cloud to mesh of cube elements
5. Optional: Use laplacian smothing on the surface and/or boundaries

Run **'Brain_creation.py'** to create example brain model to be viewed using Paraview  
Edit **'model_config.ini'** file to give new input and output directory and other parameters.

Required packages:
- Nibabel
- Pyvists
- Numpy
- Scipy
- MRI (optional: only needed to view mri images)

