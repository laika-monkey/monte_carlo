#!/usr/bin/env python
#############################################################################_80
# USAGE	:	./monte_carlo_single.py
# AUTHOR:	Walter Raymond Sessions
# DESC	:	Code for AOS 640 Project, UW-Madison, Spring 2014
# OUTPUT:	Depending on keyward options, will produce an ACSII file 
#		containing outcomes and path for all photons. Plots can
#		be enabled
#
#		This version of the code was meant to test the behavior of the 
#		photons and toy around, less to do intensive calculations. Large
#		sets of photons should be computed using monte_carlo.py instead.
#############################################################################_80

import os
from numpy import pi
d2r = pi / 180.

namelist = {
  'single'	: {
		'tau_star'	: 25.,		#_optical depth of cloud
		'ssa'		: .99,		#_single scatter albedo
		'zenith'	: 135,		#_incoming rad dir from vertical
		'azimuth'	: 0,		#_azimuthal angle
		'g'		: .85,		#_-1 to 1, 
						# desc hemi scrat dir preference
		'imgfmt'	: 'png',	#_image format, png/pdf
		'total'		: 10000,	#_total number of photons
		'group'		: 100,		#_number of photons to pass
						# each process
						# setting this to 1 and
						# drawlive True
						# allows individual photons
						# to be observed
		'rho'		: False,	#_either plot paths (False) or
						# density (True)
		'weightmin'	: .1,		#_weight at which to del photon
		'weight'	: False,	#_do we use weights and split
						# photons?
  	},
	
	'drawlive'	: True,		#_draw photons moving as it happens
					# not recommended unless doing density
					# plots or you have an intense amount
					# of memory.
					#_mostly a sanity check
	'limit_plot' : False,	#_if true, plot area will not grow
				# with photon path
}

runs = {}	#_cycle over iterations of G, weight, angle 
debug = 1	#_high for more verbosity (roughly 1-7)

#############################################################################_80
#_MAIN_######################################################################_80
#############################################################################_80

def run_main(drawlive=False,limit_plot=True,**namelist): 
	import matplotlib.pyplot as plt
	from time import sleep
	from lib import Cloud as cloud
	from lib import Photon as photon

	#_loop over runs specified in the Namelist dictionary
	for run in namelist:
	
		#_turn on interactive window
		if drawlive:
			plt.ion()
		
		#_initialize output file
		mkdir_p(run)
		fname = '/'.join((run, 'photon.dat'))
		
		#_pull out keyword args and update	
		kwargs = namelist[run]
		clouds = cloud(**kwargs)		#_initialize cloud
		kwargs.update({	'outdir'	: run,
				'figure' 	: clouds.figure,
				'fid'		: open(fname,'w'),
				'axes'		: clouds.figure.axes[0],
				'limit_plot'	: limit_plot,
				'clouds'	: clouds,
				'drawlive'	: drawlive	})
		
		#_calc how many loops need to be made
		photon_propagation(**kwargs)
	
	plt.savefig('/'.join((run,'photon.png')))
	
def photon_propagation(total=100, group=10, drawlive=False, **kwargs):
	from numpy import array, append
	
	photons = array([]) #_initialize list to carry photon objects
	figure	= kwargs.get('figure') #_get a copy to speed up drawing
	clouds	= kwargs.get('clouds') #_get a copy of clouds
	
	while photon.tot_fin < total:

		#_keep the total number of photons in list.photons == group
		need = min([total - photon.tot_fin - len(photons),
				group - len(photons)])
		photons = append(photons, [photon(**kwargs) for n 
			in range(int(need)) ])

		#_move all photons forward one step
		[ a.advance(**kwargs) for a in photons ]

		#_update number of complete
		live = array([ p.live for p in photons ]) #_find live
		photons = photons[live]	#_reduce array to living

		if drawlive: 
			clouds.figure.axes[0].set_title(str(photon.tot_fin) +
			' complete of '+str(photon.tot_gen)+' generated')
			clouds.figure.canvas.draw()
			
#############################################################################_80
#_TOOLS_#####################################################################_80
#############################################################################_80

def two_stream(ssa=1.,g=0.,tau_star=10.,**kwargs):	
	'''
	Use two stream parameterization model to calculate 
	azimuthically independent values for reflectance and transmittance
	
	ssa,	float,	[0,1], single scatter albedo
	g,		float,	[-1,1],	asymmetry parameter
	
	returns:
	refl	float,	[0,1],	reflectivity
	tran	float,	[0,1],	transmittance 

	13.63-66
	r = I^(0)/Iv(0)	= (1-g)*tau_star / [1+(1-g)*tau_star],	ssa=1
	t = Iv(tau_star)/Iv(0) = 1 / [1+(1-g)*tau_star], 				ssa=1

	r = r_ 	[exp(gamma*tau_star) - exp(-gamma*tau_star)] / 
		[exp(gamma*tau_star ) - r_**2*exp(-gamma*tau_star)] ssa<1
	t = (1-r_**2) / [exp(gamma*tau_star) - r_**2exp(-gamma*tau_star)] ssa<1
	'''
	from numpy import sqrt, exp
	
	#_calculate gamma, Petty 13.25
	G = 2*sqrt(1-ssa) * sqrt(1-ssa*g)
		
	#_calculate semi-invinit cloud reflectivity, Petty 13.45
	r8 = sqrt(1-ssa*g) - sqrt(1-ssa)
	r8 /= sqrt(1-ssa*g) + sqrt(1-ssa)
		
	#_calculate refl and trans depending on ssa value
	if ssa == 1:
		refl = (1-g)*tau_star / (1+(1-g)*tau_star) 	#_Petty 13.63
		tran = 1 / (1+(1-g)*tau_star)			#_Petty 13.64
	elif ssa < 1:
		refl = r8 * (exp(G*tau_star) - exp(-G*tau_star)) #_13.65
		refl /= (exp(G*tau_star) - r8**2*exp(-G*tau_star))
		tran = (1-r8**2) / (exp(G*tau_star) - r8**2*exp(-G*tau_star))
	else:
		raise RuntimeError, 'invalid single scatter value ' + str(ssa)

	return refl, tran

''' count photons exiting top, bottom, or absorbed '''
def count_result(phot, res):
	return sum([ p.result == res for p in phot ])

def read_photon_file(fname,**kwargs):
	dbg(fname,3)
	with open(fname,'r') as f:
		photons = [ photon(record=line,**kwargs) \
			for line in f.readlines() ]
	return photons

def kvector_init(theta,phi):
	'''
	return unit vector pointing in direction defined by theta and phi 
	In model space, this is relative to the vertical.  In photon space,
	it's relative to the propoagation direction
	theta,phi	float,	directions about axis of propagation IN RADIANS
	
	'''
	from numpy import matrix,cos,sin
	return matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),-cos(theta)]).T
		
def kvector(theta,phi):
	'''
	return unit vector in direction of theta/phi
	Initially, this will be in relationship to standard 3d cartesian coords,
	but as the propagation path shifts, it will be in relation to the last
	path

	theta,phi	float,	directions about axis of propagation IN RADIANS
	'''
	from numpy import matrix,cos,sin
	return matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)]).T

def paxes(k,phi):
	'''
	vector is horizontal in model space and orthogonal to propagation
	returns array describing coordinate system to be augmented
	'''
	from numpy import matrix, cross, dot, array, sin, cos, tile
	from numpy.linalg import norm

	z = matrix([0,0,1]).T 			#_vertical model space
	p = cross(z.T,k)			#_ortho to vert and prop
	mag = matrix([ norm(x) for x in p ]) 	#_magnitude 

	x = (p / mag.T)				#_get unit vector
	if (k[:,0]**2 + k[:,1]**2) == 0: 	#_replace nans
		x = matrix([-sin(phi),cos(phi),0.])

	y = matrix(cross(k,x))			#_calc y unit vector
	 
	return matrix([x.A.flatten(),y.A.flatten(),k.flatten()]).T

def setup_groups(objects,nproc=1,**kwargs):
	'''
	If I ever bother to parallelize photon group processing...
	
	USAGE:	GROUPS 			= setup_groups( <ITER_OBJ>, <INT> )
		[[1,2,3],[4,5,6]]	= setup_groups( range(6), 3 )

	Break iterable object into groups of NP size to be forked.
	objects	: iterable object 
	nproc	: number of items per group, int
	'''
	n_obs = len(objects)			#_total number of items
	ng = int(n_obs/nproc)	 		#_number of full groups
	extra = 1 if (n_obs%nproc) else 0	#_check for residual 
	groups = [objects[nproc*i:nproc*(i+1)] for i in range(ng)]
	if extra: groups.append(objects[nproc*ng:])
	dbg((objects,groups),l=4)
	return groups
	
def mkdir_p(path):
	'''
	USAGE: 	mkdir_p(<FULL_PATH_STRING>)
		Similar to using `mkdir -p /directory/to/make/`
	'''
	if path == '.': return
	import errno, os
	dbg(path,5)
	def mkdir(p):
		try: 
			os.makedirs(p)
			os.chmod(p,0775)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST: pass
			else: raise

	if hasattr('__iter__',path): [ mkdir(p) for p in path ]
	else: mkdir(path)

def dbg(msg,l=1,err=False):
	''' 
	if global debug is set to true, be more verbose 
	msg	: str, Message to be printed
	l	: int, Debug level of message.  Set higher for lower level 
			messages.  As debug increases, noisiness should also.
	'''
	msg = msg.__repr__()

	import inspect
	if debug >= l:
		curf 	= inspect.currentframe()
		calf 	= inspect.getouterframes(curf,2)
		file, line, method = calf[1][1:4]
		file	 = file.split('/')[-1]
		scream 	= '[%s.%s.%i] %s' % (file,method,line,msg)
		
		if not err: print scream
		else: raise RuntimeError, scream
			
if __name__ == '__main__':
	run_main(**namelist)
