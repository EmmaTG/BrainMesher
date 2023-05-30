# BrainMesher

Create 3D brain mesh from mri images using only hexahedral elements.
Steps:
1. Import aseg.mgz from freesurfer output (after performing ,recon-all' on mri images)
2. Convert voxel image then point cloud -> mesh of cube elements
3. Use laplacian smothing on the surface

FUTURE WORK: 
- Addition of CSF
- Use Lapacian smoothing on the interface boundarries of these regions

Required packages:
- Nibabel
- Pyvists
- Numpy
- MRI (optional: only needed to view mri images)
