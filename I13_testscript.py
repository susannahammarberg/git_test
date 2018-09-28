# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:23:20 2018

@author: Sanna

Copy from ptypy3d.py


3d reconstructions using ptypy.

Possibility for shifting vectors for realignment of the data according to
'real' motorpositions (aquired from STXM maps). Aligned between scans.

"""
import ptypy
from ptypy.core import Ptycho
from ptypy import utils as u
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#from ptypy.experiment.nanomax3d import NanomaxBraggJune2017 # after update need to update spec ptyScan class
from ptypy.experiment.I13_Bragg3d import I13Bragg3d # after update need to update spec ptyScan class

# TODO
# if InP
#       choose InP roi and change name

p = u.Param()
p.run = 'I13Bragg3d'   # 'XRD_InP'

sample = 'name'; scans = [155248, 155249] #range(192, 200+1)+range(205, 222+1)
#sample = 'JWX29A_NW1' #; scans =[458,459]
#scans = [458,459,460,461,462,463,464,465,466,467,468,469,470,471,518,473,474,475,476,477,478,479,480,481,482,483,484,485,486,519,488, 496,497,498, 499, 500, 501, 502, 503, 504, 505, 506,507, 508, 509, 510, 511, 512, 513, 514, 515]

p.data_type = "single"   #or "double"
# for verbose output
p.verbose_level = 5

# use special plot layout for 3d data  (but the io.home part tells ptypy where to save recons and dumps)
p.io = u.Param()
p.io.home = './'
p.io.autosave = u.Param()
p.io.autosave.interval = 1 # does not work
p.io.autoplot = u.Param()
p.io.autoplot.layout = 'bragg3d'
p.io.autoplot.dump = True
p.io.autoplot.interval = 1
# TODO make it create the plots

 
p.scans = u.Param()
p.scans.scan01 = u.Param()
p.scans.scan01.name = 'Bragg3dModel'
#p.scans.scan01.illumination = illumination
p.scans.scan01.data= u.Param()
p.scans.scan01.data.name = 'I13Bragg3d'
#base= 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/'
#p.scans.scan01.data.datapath = 'D:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/%s/' % sample
p.scans.scan01.data.datapath = 'D:/Diamond_I13_20180304_take2/testNexus_2/'
# open the nexus file for the metadata. Dont need this wor, because there is one nexus for each rotation


#print p.scans.scan01.data.datapath
p.scans.scan01.data.detfilepattern = 'scan_%04d_merlin_%04d.hdf5'
# not sure if this loads properly
p.scans.scan01.data.maskfile = 'C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/merlin_mask.h5'
p.scans.scan01.data.scans = scans
p.scans.scan01.data.center = None
p.scans.scan01.data.theta_bragg = 12.0
p.scans.scan01.data.shape = 150#512#150#60#290#128
# ptypy says: Setting center for ROI from None to [ 75.60081158  86.26238307].   bu that must be in the images that iI cut out from the detector
#p.scans.scan01.data.center = (200,270) #(512-170,245)     #(512-170,245) for 192_   #Seems like its y than x
#p.scans.scan01.data.load_parallel = 'all'
p.scans.scan01.data.psize = 55e-6
p.scans.scan01.data.energy = 9.49
p.scans.scan01.data.distance = 1

#sannas shifting parameters
##############################################################################
#S segmented NWs
##############################################################################
# the smallers values here (minus values) will include the 0:th scanning position
##TODO test with no shifting and compare result: [0]*51 #
#p.scans.scan01.data.vertical_shift =  [-1,-1,0,0,0,  0,0,2,1,0,  1,1,1,0,-1,  -1,-1,-1,-1,0,  -1,-1,0,0,1,  1,-1,0,1,0,   2,0,0,1,1,  1,0,0,1,1,  1,2,2,2,4,  3,3,3,3,3,   3]
#p.scans.scan01.data.horizontal_shift =  [3,2,0,1,2,  3,4,3,4,5,  5,6,6,5,6,  5,4,7,8,8,  8,8,10,11,12,  11,12,12,11,12,  12,11,12,13,13,  14,15,14,14,14,  13,15,16,15,14,  17,19,18,18,17,   17]
            # InGaP = [116:266,234:384]
            # InP = [116:266,80:230]
#p.scans.scan01.data.detector_roi_indices = [116,266,80,230]  #[116,266,80,230]  #    # this one should not be needed since u have shape and center...

##############################################################################
# homogenius NWs
##############################################################################
#p.scans.scan01.data.vertical_shift = [1] * len(scans) 
#p.scans.scan01.data.horizontal_shift = [1] * len(scans) 

# TODO read these in:            
p.scans.scan01.data.vertical_shift = [ 0, 0, 1, 2, 2, 2, 3, 3, 3, 1, 0, 0, 0, 0, 0, -2, -2, -2, -3, -3, -2, -2,  -3, -2, -3, -3,  -4]
p.scans.scan01.data.horizontal_shift = [ -8, -8, -8, -8, -8, -3, -2, -4, 0, -5, -3, -6, -3, -6, -3, -5, -5, -7, -4,  -7, -3, -6, -2, -7, -2, -6, -3 ] 
p.scans.scan01.data.detector_roi_indices = [275,425,150,300]  # this one should not be needed since u have shape and center...

p.scans.scan01.illumination = u.Param()
p.scans.scan01.illumination.aperture = u.Param() 
p.scans.scan01.illumination.aperture.form = 'circ'
p.scans.scan01.illumination.aperture.size = 100e-9 
p.scans.scan01.sample = u.Param()
p.scans.scan01.sample.fill = 1e-3

# to use a probe from an old reconstruction:
p.scans.scan01.illumination.model = 'recon'
p.scans.scan01.illumination.recon = u.Param()
#p.scans.scan01.illumination.recon.layer = None   #parameter valiadation failed
#p.scans.scan01.illumination.recon.ID = None      
#p.scans.scan01.illumination.recon.rfile = 'C:/Users/Sanna/Desktop/python_utilities/scan10_merlin_DM_manipulated.ptyr'#'test_hdf5_file.ptyr'#
p.scans.scan01.illumination.recon.rfile ='C:/Users/Sanna/Documents/beamtime/NanoMAX062017/Analysis_ptypy/nice_probe_ptyrfiles/scan10/scan10_pilatus_ML.ptyr'# 'C:/Users/Sanna/Desktop/ptypy_scripts/recons/scan10_pilatus/scan10_pilatus_ML.ptyr'


p.engines = u.Param()
p.engines.engine00 = u.Param()
p.engines.engine00.name = 'DM'    #Not 'DM_3dBragg' ? 
p.engines.engine00.numiter = 2
p.engines.engine00.probe_update_start = 100000
p.engines.engine00.probe_support = None
#p.engines.engine00.sample_support = None

# p.engines.engine00.sample_support = u.Param()
# p.engines.engine00.sample_support.coefficient = 0.0 
# p.engines.engine00.sample_support.type = 'rod'
# p.engines.engine00.sample_support.size = 200e-9
# p.engines.engine00.sample_support.shrinkwrap = u.Param()
# p.engines.engine00.sample_support.shrinkwrap.cutoff = .3
# p.engines.engine00.sample_support.shrinkwrap.smooth = None
# p.engines.engine00.sample_support.shrinkwrap.start = 15
# p.engines.engine00.sample_support.shrinkwrap.plot = True

# prepare and run
P = Ptycho(p,level=5)

for k, v in P.obj.S.items():
    print k, v
    
for l, m in P.model.scans.items():
    print l, m
    
for o,n in P.diff.S.items():
    print o, n
    
    
for o,n in P.exit.S.items():
    print o, n
    
P.model.scans['scan01'].ptyscan.info