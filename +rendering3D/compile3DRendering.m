function compile3DRendering()
    mex +rendering3D\smoothpatch_curvature_double.c -v
    mex +rendering3D\smoothpatch_inversedistance_double.c -v
    mex +rendering3D\vertex_neighbours_double.c -v

end