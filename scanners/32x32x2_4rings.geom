# Scanner file for norm-functions, in the CASToR cylindrical PET format (.geom).
# Only the fields below are read by norm-functions; other CASToR fields ('voxels number', 'field of view', ...)
# are ignored, so a reconstruction .geom describing the same geometry can be passed as well.
#
# The 59 x 59 x 10 mm module is segmented (virtually) into submodules/crystals/layers, so
# 'crystals size ...' is the voxel pitch: 1.84375 x 1.84375 x 5 mm.
# Layer-dependent fields list one value per layer (layer 0, layer 1). Lengths in mm, gaps edge to edge.

modality: PET
scanner name: 32x32x2_4rings_system
description: 4 module(s) axially, 32 submodules axially x 32 crystals transaxially, 2 depth layers
number of elements: 262144
number of layers: 2

scanner radius: 321.3, 326.3     # isocentre -> front face of each layer
number of rsectors: 32, 32
number of rsectors axial: 1, 1
number of modules transaxial: 1, 1
number of modules axial: 4, 4
number of submodules transaxial: 1, 1
number of submodules axial: 32, 32
number of crystals transaxial: 32, 32
number of crystals axial: 1, 1
crystals size trans: 1.84375, 1.84375
crystals size axial: 1.84375, 1.84375
crystals size depth: 5, 5
module gap axial: 4, 4
