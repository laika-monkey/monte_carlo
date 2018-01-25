#!/usr/bin/env python
#############################################################################_80
# USAGE	:	./monte_carlo_single.py
# AUTHOR:	Walter Raymond Sessions
# DESC	:	Code for AOS 640 Project, UW-Madison, Spring 2014
# OUTPUT:	Depending on keyward options, will produce an ACSII file 
#			containing outcomes and path for all photons. Plots can be enabled
#
#			This version of the code was meant to test the behavior of the 
#			photons and toy around, less to do intensive calculations. Large
#			sets of photons should be computed using monte_carlo.py instead.
#############################################################################_80

import os
from numpy import pi
d2r = pi / 180.

namelist = {
  'single'	: {
		'tau_star'	: 25.,		#_optical depth of cloud
		'ssa'		: .99,		#_single scatter albedo
		'zenith'	: 135,		#_incoming radiation direction from vertical 
		'azimuth'	: 0,		#_azimuthal angle
		
		'g'			: .85,		#_-1 to 1, describes hemis scrat dir preference
		
		'imgfmt'	: 'png',	#_image format, png/pdf
	
		'total'		: 10000,		#_total number of photons
		'group'		: 100,		#_number of photons to pass each process
								# setting this to 1 and drawlive True
								# allows individual photons to be observed
		
		'rho'		: False,	#_either plot paths (False) or density (True)
	
		'weightmin'	: .1,		#_weight at which to ditch photon
		'weight'	: False,	#_do we use weights and split photons?
  	},
	
	'drawlive'	: True,		#_draw photons moving as it happens
							# SUPER not recommended unless doing density
							# plots or you have an intense amount of memory.
							#_mostly a sanity check
	'limit_plot' : False,	#_if true, plot area will not grow with photon path
}

runs = {}					#_cycle over iterations of G, weight, angle 
debug = 1					#_high for more verbosity (roughly 1-7)

#############################################################################_80
#_MAIN_######################################################################_80
#############################################################################_80

def run_main(drawlive=False,limit_plot=True,**namelist): 
	import matplotlib.pyplot as plt
	from time import sleep
	for run in namelist:	
		#_turn on interactive window
		if drawlive: plt.ion()
		
		#_initialize output file
		mkdir_p(run)
		fname = '/'.join((run,'photon.dat'))
		
		#_pull out keyword args and update	
		kwargs = namelist[run]
		clouds = cloud(**kwargs)		#_initialize cloud
		kwargs.update({	'outdir'	: run,
						'figure' 	: clouds.figure, #figure,
						'fid'		: open(fname,'w'),
						'axes'		: clouds.figure.axes[0],
						'limit_plot': limit_plot,
						'clouds'	: clouds,
						'drawlive'	: drawlive	})
		
		#_calc how many loops need to be made
		photon_propagation(**kwargs)
	
	sleep(50)
	plt.savefig('/'.join((run,'photon.png')))
	
def photon_propagation(total=100,group=10,drawlive=False,**kwargs):
	from numpy import array, append
	
	photons = array([])	#_initialize list to carry photon objects
	figure	= kwargs.get('figure') #_get a copy to speed up drawing
	clouds	= kwargs.get('clouds') #_get a copy of clouds
	
	while photon.tot_fin < total:
		#_keep the total number of photons in list.photons == group
		need 	= min([total-photon.tot_fin-len(photons), group-len(photons)])
		photons = append(photons, [ photon(**kwargs) for n in range(int(need)) ])

		#_move all photons forward one step
		[ a.advance(**kwargs) for a in photons ]

		#_update number of complete
		live 	= array([ p.live for p in photons ])	#_find live
		photons = photons[live]							#_reduce array to living
		if drawlive: 
			clouds.figure.axes[0].set_title(str(photon.tot_fin) + ' complete of '
											+ str(photon.tot_gen) + ' generated')
			clouds.figure.canvas.draw()
			
#############################################################################_80
#_CLASSES_###################################################################_80
#############################################################################_80

class cloud(object):
	def __init__(self,tau_star=10.,width=40,figure=None,zenith=180,azimuth=0,
		g=0,ssa=1.,clouds=None,**kwargs):
		''' 
		main purpose is to keep track of initial cloud conditions
		and plotting axis
	
		tau_star,	float,	optical depth of cloud, constant
		width,		float,	optical width of plotting area
		zenith,		float,	incident ray zenith angle, degrees
		azimuth,	float,	incident ray azimuth angle, degrees
		'''
		import matplotlib.pyplot as plt
		from numpy import cos, tan, pi
		d2r = pi/180
		
		self.tau_star 	= tau_star
		self.zenith		= zenith
		self.azimuth	= azimuth
##		self.ssa		= ssa
##		self.g			= g

		#_initialize figure and axes
		self.figure = figure if figure != None \
			else plt.figure(figsize=[16.,4.])
		self.figure.add_subplot(111)
		
		#_initialize cloud top and base lines
		xlim = [-width/2.,width/2.]
		ylim = [-2,tau_star+2]
		self.figure.axes[0].set_ylim(ylim[::-1])
		self.figure.axes[0].set_xlim(xlim)
		
		arg = {'linewidth':.5,'color':'k'}
		a, 	= self.figure.axes[0].plot(xlim,[tau_star,tau_star],**arg)
		b, 	= self.figure.axes[0].plot(xlim,[0,0],**arg)

		#_plot incident ray path
		x 	= [0,self.incident_x(ylim[0])]
		y 	= [0,ylim[0]]
		c, 	= self.figure.axes[0].plot(x,y,linewidth=3,color='k',linestyle='--')

		#_save lines
		self.base = a
		self.top = b
		self.ray = c
		self.figure.canvas.draw()	

	def incident_x(self,o):
		''' calculate the x component of the incident ray '''
		from numpy import cos, tan, pi
		return -o*tan(self.zenith*d2r-pi)*cos(self.azimuth*d2r)

class photon(object):
	tot_gen	= 0
	tot_fin = 0
	def __init__(self,record=False,zenith=180,azimuth=0,figure=None,axes=None,
		clouds=None,**kwargs):
		'''
		basically a home brew numpy.recarray, only order is not muteable
		this should not be used as a post processing class, otherwise
		there will be giant, unwieldy arrays everywhere
		
		size,		int,	number of photons to handle at a time
		record,		str,	line of old photon.dat file
		zenith,		flt,	degrees from vertical of incidence (180==down)
							positive Z goes DOWN into CLOUD. 
		azimuth,	flt,	degrees from north of incidence (0==north)
		figure,		plt.fig	matplotlib artist object to be passed if there
							are multiple calls to this class, but only
							one desired plot area
							
		clouds,		class,	if ever get around to it, associate a photon.cloud
							attribute that points to the containing cloud class,
							which in turn holds the tau_star/g/ssa properties.
							Photon class should not contain them
		'''
		from numpy import array, zeros, ones, tile, nan, pi
		import matplotlib.pyplot as plt
		d2r = pi/180
		dbg((record,zenith,azimuth),5)
		
		#_generate object to hold #size photons
		if not record:		
			self.tau		= 0 					#_current distance from top
			self.k			= kvector_init(zenith*d2r,azimuth*d2r).T.A
													#_current direction
			self.phi		= azimuth*d2r			#_most recent phi
			self.theta		= zenith*d2r			#_most recent theta
			self.history	= zeros((3,1))			#_each coordinate traveled
			self.weight		= 1 					#_ignore
			self.result		= None					#_top, base, or absorbed
			self.scattered	= 0 					#_number of times scattered
			self.live		= True					#_still able to advance?
			self.figure		= plt.figure(figsize=([16.,4.])) if figure == None else figure
			self.axes		= self.figure.add_subplot(111) if axes==None else axes
			self.line		= None

		#_when passed a file name, returns previously run photon results
		else:
			#_break up record and put back in
			cols 			= record.split()
			shape 			= (len(cols[3:])/3,3)	#_nscat, nk
			history			= array(cols[3:],dtype='f8').reshape(shape)
			
			#_if want to be able to load and continue advancing,
			# need to add back in calculation for k vect based on hist[:,-2:]
			self.tau 		= float(cols[0])
			self.result 	= cols[1]
			self.scattered 	= int(cols[2])
			self.history	= history.T
			self.live		= False					#_don't allow advance
													# mostly for plotting
			self.figure		= plt.figure(figsize=([16.,4.])) if figure == None else figure
			self.axes		= self.figure.add_subplot(111) if axes == None else axes
			self.line		= None
			
		photon.tot_gen += 1
		
	def advance(self,ssa=1.,weight=False,force_angle=None,tau_star=10.,lag=0,
		**kwargs):
		'''
		roll to see if photon's +12 to scattering prevents grue from eating it
		
		Roll to see if absorbed (since all photons begin at cloud top, this
								begins each sequence)
		Roll to see which direction scattered
			F/B based on HG
			Azimuth uniform
		Roll to see distance
			Beer-Lambert Law
		'''
		from numpy.random import random as r
		from numpy import append, array, log, tile
		from time import sleep
		dbg((ssa,force_angle,tau_star),3)
			
		if not self.live: return
		
		#_roll to see how far until next interaction____________________________
		self.tau = -log(1-r(1)[0])	#_THIS IS NOT TAAAAUUUUU
		
		#_update current location
		point	 		= (self.history[:,-1]+self.tau*self.k).reshape(3,1)
		self.history 	= append(self.history,point,axis=1)
		self.plot(**kwargs)
		
		#_check if escaped top or bottom, update flags
		if self.history[2,-1] < 1e-8:
			self.die('top',**kwargs)
			return
		elif self.history[2,-1] > tau_star:
			self.die('base',**kwargs)
			return
					
		#_roll to see if absorbed_______________________________________________
		if r(1)[0] > ssa:
			self.die('absorbed',**kwargs)
			return
		
		#_roll to see in what direction_________________________________________
		if force_angle == None:
			phi		= 2*pi*r(1)[0]							#_randomize azimuth
			theta 	= henyey_greenstein(r(1)[0],**kwargs)	#_randomize theta
		else: #_force specific scattering angles for testing
			phi 	= force_angle[0]*d2r
			theta 	= force_angle[1]*d2r
		self.__turn__(theta,phi,**kwargs)					#_update k

		self.scattered += 1
				
		#_update number of times scattered
		dbg(self.traveled(origin=False),5)
		
		#_if we do a weighting check, add it here? above?
		#_update plot
		#self.plot(**kwargs) #_plotting here leaves out last segment
			
	def __turn__(self,theta,phi,**kwargs):
		''' augment direction of vector by theta,phi '''
		dbg((theta,phi),5)
		from numpy import ones, dot, arange
		
		#_build 3x3 array describing coordinate system orthogonal to the
		# propagating ray and with x parallel to the model horizontal plane
		aug = paxes(self.k,self.phi)

		#_get a unit vector describing new coordinate direction
		kpp 	= kvector(theta,phi)
		self.k 	= dot(aug,kpp).T.A

		#_update angle arrays
		self.phi 	= phi
		self.theta 	= theta
	
	def die(self,result,**kwargs): #_try to escape
		self.result 	= result
		self.live 		= False
		photon.tot_fin += 1
		self.dump(**kwargs)
		
	def dump(self,outdir='.',outdat='photon.dat',fid=None,**kwargs):
		'''
		remove photon from record, write 
		TAU		RESULT	SCATTERED	HISTORY
		10.5f	10s		10i 		(10.5f,10.5f,10.5f) * nscat+1
		
		index,	bool,	location of record to remove
		'''
		from numpy import array, where, empty
		dbg((outdir,outdat),3)
		mkdir_p(outdir)
			
		fopened = True if fid != None else False
		if not fopened:
			fname = '/'.join((outdir,outdat))
			fid = open(fname,'a')			#_open file handle
		
		loc = self.history.T.flatten()
		val = [self.tau,self.result,self.scattered]; val.extend(loc)
		
		#_make formatted print statement and write to file
		nk,nscat= self.history.shape
		fmt 	= '%10.5f%10s%10i ' + '%10.5f%10.5f%10.5f' * nscat + '\n'
		line 	= fmt % tuple(val)

		fid.write(line)
		
		if not fopened: fid.close()			#_close file handle
			
	def traveled(self,origin=False):
		'''
		returns geometric distance traveled by photon, total or to entry
		origin,	bool,	return only distance to origin, as the non-scattering
						photon flies
		'''
		from numpy import sqrt
		dbg(origin,5)
					
		if not origin:	#_total meanderings
			segments = self.history[:,1:] - self.history[:,:-1]
			return sqrt((segments**2).sum())
		else:			#_from cloud entrance
			return sqrt((self.history[:,-1]**2).sum())

	def plot(self,idx=None,ax=None,outdir='.',imgfmt='png',rho=False,
		limit_plot=True,clouds=None,**kwargs):
		'''
		plot photon path
		idx,	int,		if supplied, only the listed photon paths will plot
		ax,		plt.artist,	matplotlib plot space
		outdir,	str,		output directory
		imgfmt,	str,		image file type
		rho,	bool,		instead of plotting path, plot density
		'''
		from numpy import arange, array, tan, cos, pi
		mkdir_p(outdir)
		x, y, z = self.history
		if self.line == None:				#_new line
			if len(self.figure.axes) == 0:	#_new plot completely
				self.axes = self.figure.add_subplot(111)
				self.axes.set_ylim(12,-2)
				self.axes.set_xlim(-20,20)
			self.line, = self.axes.plot(x,z,linewidth=.3)
		else: 					#_update line
			self.line.set_xdata(x)
			self.line.set_ydata(z)

		#_expand plotting area if photons wander far
		if not limit_plot:
			cxmin, cxmax = self.axes.xaxis.get_view_interval()
			cymin, cymax = self.axes.yaxis.get_view_interval()
			xmax = max([cxmax,x.max()])
			xmin = min([cxmin,x.min()])
			self.axes.set_xlim(xmin,xmax)
			
			#_update cloud lines
			if clouds != None:
				x = [0,clouds.incident_x(cymax)]
				y = [0,cymax]
				clouds.ray.set_xdata(x)
				clouds.ray.set_ydata(y)
				clouds.base.set_xdata([cxmin,cxmax])
				clouds.top.set_xdata([cxmin,cxmax])	
			
		#_only call draw if this is the only line,
		# otherwise do it in a wrapper loop
		if len(self.axes.lines) == 1: 
			self.figure.canvas.draw()
			self.figure.show()
	
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
	r = I^(0)/Iv(0) 		= (1-g)*tau_star / [1+(1-g)*tau_star],	ssa=1
	t = Iv(tau_star)/Iv(0) 	= 1 / [1+(1-g)*tau_star], 				ssa=1

	r = r_ 	[exp(gamma*tau_star) - exp(-gamma*tau_star)] / 
			[exp(gamma*tau_star ) - r_**2*exp(-gamma*tau_star)]			ssa<1
	t = (1-r_**2) / [exp(gamma*tau_star) - r_**2exp(-gamma*tau_star)] 	ssa<1
	'''
	from numpy import sqrt, exp
	
	#_calculate gamma, Petty 13.25
	G = 2*sqrt(1-ssa)*sqrt(1-ssa*g)
		
	#_calculate semi-invinit cloud reflectivity, Petty 13.45
	r8 	= sqrt(1-ssa*g) - sqrt(1-ssa)
	r8 /= sqrt(1-ssa*g) + sqrt(1-ssa)
		
	#_calculate refl and trans depending on ssa value
	if ssa == 1:
		refl = (1-g)*tau_star / (1+(1-g)*tau_star)		#_Petty 13.63
		tran = 1 / (1+(1-g)*tau_star)					#_Petty 13.64
	elif ssa < 1:
		refl  = r8 * (exp(G*tau_star) - exp(-G*tau_star))	#_13.65
		refl /= (exp(G*tau_star) - r8**2*exp(-G*tau_star))
		tran  = (1-r8**2) / (exp(G*tau_star) - r8**2*exp(-G*tau_star))
	else:
		raise RuntimeError, 'invalid single scatter value ' + str(ssa)

	return refl, tran

''' count photons exiting top, bottom, or absorbed '''
def count_result(phot,res): return sum([ p.result == res for p in phot ])

def read_photon_file(fname,**kwargs):
	dbg(fname,3)
	with open(fname,'r') as f:
		photons = [ photon(record=line,**kwargs) for line in f.readlines() ]
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
	but as the propagation path shifts, it will be in relation to the last path

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

	z 	= matrix([0,0,1]).T							#_vertical model space
	p 	= cross(z.T,k)								#_ortho to vert and prop
	mag = matrix([ norm(x) for x in p ])			#_magnitude 

	x = (p / mag.T)		 							#_get unit vector
	if (k[:,0]**2 + k[:,1]**2) == 0:				#_replace nans
		x = matrix([-sin(phi),cos(phi),0.])

	y = matrix(cross(k,x))							#_calc y unit vector
	 
	return matrix([x.A.flatten(),y.A.flatten(),k.flatten()]).T

def henyey_greenstein(r,g=0,**kwargs):
	'''
	calculate scattering theta in accordance with Henyay Greenstein dist
	r,	flt,	random number, [0,1] INCLUSIVE
	g,	flt,	asymettry parameter.  0 => equal f/b hemisphere, 1, all forward

	Binzoni, T., et al. "The use of the Henyey-Greenstein phase function
	in Monte Carlo simulations in biomedical optics." Physics in medicine 
	and biology 51.17 (2006): N313.
	'''
	from numpy import arccos, isnan
	if g == 0:	#_give random isotropic direct
		cos_theta = 2.*r-1	#_this leaves out {pi,0}, 2.*(1-r)-1 leaves out pi/2
	else:		#_give hemispherically preferred direction
		cos_theta = (1+g*g-((1-g*g)/(1-g+2.*g*r))**2)/2./g
	return arccos(cos_theta)
	
def setup_groups(objects,nproc=1,**kwargs):
	'''
	If I ever bother to parallelize photon group processing...
	
	USAGE:	GROUPS 			= setup_groups( <ITER_OBJ>, <INT> )
		[[1,2,3],[4,5,6]]	= setup_groups( range(6), 3 )

	Break iterable object into groups of NP size to be forked.
	objects	: iterable object 
	nproc	: number of items per group, int
	'''
	n_obs	= len(objects)				#_total number of items
	ng		= int(n_obs/nproc)	 		#_number of full groups
	extra	= 1 if (n_obs%nproc) else 0	#_check for residual 
	groups 	= [objects[nproc*i:nproc*(i+1)] for i in range(ng)]
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
	msg = to_string(msg)
	if hasattr(msg,'__iter__'): msg = ' '.join(msg)

	import inspect
	if debug >= l:
		curf 	= inspect.currentframe()
		calf 	= inspect.getouterframes(curf,2)
		file, line, method = calf[1][1:4]
		file	 = file.split('/')[-1]
		scream 	= '[%s.%s.%i] %s' % (file,method,line,msg)
		
		if not err: print scream
		else: raise RuntimeError, scream
			
def to_string(a):
	if hasattr(a,'__iter__'): return [ str(b) for b in a ]
	else: return str(a)

if __name__ == '__main__': run_main(**namelist)