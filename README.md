# BrainMesher
Create 3D brain mesh from stl using only hexahedral elements usign MATLAB.
Steps:
1. Import stl
2. Convert stl to voxel image and then mesh of voxel/cube elements
3. Use laplacian smothing on the surface
FUTURE WORK: 
- Segment white vs grey matter and Regional areas of brain using data from FreeSurfer
- Use Lapacian smoothing on the interface boundarries of these regions
