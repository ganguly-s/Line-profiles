import yt
import os
import sys
import array as ar
import numpy as np
from AthenaSims import *
from yt.visualization.fixed_resolution import FixedResolutionBuffer
from yt.visualization.image_writer import write_image

# import user libraries
import AthenaInterface as AI

#---------------------------------------------------------------------
# Inputs
#---------------------------------------------------------------------  
PNG_IMAGES = False # if True, this script only generates simple png images

# set the number of pixels in x,y direction to the sim resolution
xpix = 1024    
ypix = 1024

# set the domain size
xo,xf = 0.,1.
yo,yf = 0.,1.

# specify the dumps you want to extract data from
# frames = [45,50,60,130,158,182,229]
frames1 = ar.array('i',(i for i in range(0,40,5)))
frames2 = ar.array('i',(i for i in range(40,80,2)))
frames = frames1 + frames2 

# specify which variables to extract: must be yt's syntax
if PNG_IMAGES:
	vars = ['density']
else:
	vars = ['density','velocity_x','velocity_y','pressure','temperature','vorticity_z','baroclinic_vorticity_z']
	vars = ['density','velocity_x','velocity_y','pressure','temperature']
	vars = ['pressure','temperature']
	
#---------------------------------------------------------------------
# Specify simulation info using a dictionary for each one
#---------------------------------------------------------------------
Sim2DRun1 = {}
Sim2DRun1['simpath'] = '/Volumes/Untitled/2D_LPpaper/cRFLDX_pdv5/' 
Sim2DRun1['simfolder'] = 'noPD_square'

Sim2DRun1['simtag'] = 'ti'
Sim2DRun1['timestep'] = 1.
Sim2DRun1['gamma'] = 5./3.

# path where images/data will be stored
imagepath = 'sim_data/' + Sim2DRun1['simfolder']

if PNG_IMAGES:
	imagefolder = imagepath + '/yt_png/'
else:
	imagefolder = imagepath + '/yt_dat/'

if not os.path.exists(imagefolder):
	os.makedirs(imagefolder)

#---------------------------------------------------------------------
# Figure out which dumps need to be gathered
#---------------------------------------------------------------------
dumps = []
for frame in frames:
	var = 'density'	
	if PNG_IMAGES:
		img_name = var + "%03d" % (frame,)
	else:
		img_name = var + `frame`
		
	print "Attempting to open dump:\n{}".format(imagefolder + img_name)	
	try:
		data = np.load(imagefolder + img_name)
	except:
		print "Adding {}".format(frame)
		dumps.append(frame)
		
#---------------------------------------------------------------------
# Load data by calling yt
#---------------------------------------------------------------------
for i in dumps:
	# get simulation directory
	mySim = AI.AthenaSimulation(Sim2DRun1)
	file = mySim.getFilepath(i,ext='.vtk') 
	
	# now call yt
	ds = yt.load(file) # load dataset
	gs = ds.index.select_grids(ds.index.max_level) # load all grids
	grid = gs[0] # main grid
	
	# create a frb: antialias must be True if you don't want your data altered
	prj = ds.proj('density',2)
	prj_frb = FixedResolutionBuffer(prj, (xo, xf, yo, yf),(ypix, xpix), antialias = True, periodic = True)
	
	# load all specified variables
	for var in vars:
		new_prj = prj_frb[var]
		img = np.array(new_prj) 
		if PNG_IMAGES:
			img_name = var + "%03d" % (i,)
			write_image(img,imagefolder + img_name + '.png',cmap_name='Spectral')
		else:
			img_name = var + `i`
			varfile = open(imagefolder + img_name, 'w')
			np.save(varfile,img) 
			#np.save(varfile,np.log10(img))
			varfile.close()
		print "Saving image {}...".format(img_name)
		
# save position coordinates (done just once since a fixed grid is assumed)
if not(PNG_IMAGES):
	x_prj = prj_frb['x']; x_img = np.array(x_prj)
	y_prj = prj_frb['y']; y_img = np.array(y_prj) 
	xfile = open(imagefolder + 'X', 'w'); np.save(xfile,x_img); xfile.close()
	yfile = open(imagefolder + 'Y', 'w'); np.save(yfile,y_img); yfile.close()
