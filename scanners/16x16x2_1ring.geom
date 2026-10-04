# Scanner file for norm-functions, in the CASToR cylindrical PET format (.geom).
# Only the fields below are read by norm-functions; other CASToR fields ('voxels number', 'field of view', ...)
# are ignored, so a reconstruction .geom describing the same geometry can be passed as well.
#
# The 59 x 59 x 10 mm module is segmented (virtually) into submodules/crystals/layers, so
# 'crystals size ...' is the voxel pitch: 3.6875 x 3.6875 x 5 mm.
# Layer-dependent fields list one value per layer (layer 0, layer 1). Lengths in mm, gaps edge to edge.

modality: PET
scanner name: 16x16x2_1ring_system
description: 1 module(s) axially, 16 submodules axially x 16 crystals transaxially, 2 depth layers
number of elements: 16384
number of layers: 2

scanner radius: 321.3, 326.3     # isocentre -> front face of each layer
number of rsectors: 32, 32
number of rsectors axial: 1, 1
number of modules transaxial: 1, 1
number of modules axial: 1, 1
number of submodules transaxial: 1, 1
number of submodules axial: 16, 16
number of crystals transaxial: 16, 16
number of crystals axial: 1, 1
crystals size trans: 3.6875, 3.6875
crystals size axial: 3.6875, 3.6875
crystals size depth: 5, 5
module gap axial: 4, 4
