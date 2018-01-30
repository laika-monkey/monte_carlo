from tools import dbg
class Photon(object):

	def __init__(self, record=False, zenith=180, azimuth=0, figure=None,
		axes=None, clouds=None, **kwargs):
		'''
		basically a home brew numpy.recarray, only order is not muteable
		this should not be used as a post processing class, otherwise
		there will be giant, unwieldy arrays everywhere
		
		size,		int,	number of photons to handle at a time
		record,		str,	line of old photon.dat file
		zenith,		flt,	degrees from vertical of incidence
					positive Z goes DOWN into CLOUD. 
					(180 == down)
		azimuth,	flt,	degrees from north of incidence
					(0==north)
		figure,		plt.fig	matplotlib artist object to be passed if
					there are multiple calls to this class,
					but only one desired plot area
		clouds,		class,	if ever get around to it, associate a
					photon.cloud attribute that points to
					the containing cloud class,
					which in turn holds the tau_star/g/ssa
					properties.
					Photon class should not contain them
		'''
		from numpy import array, zeros, ones, tile, nan, pi
		import matplotlib.pyplot as plt
		from lib.math import kvector_init

		d2r = pi/180
		dbg((record, zenith, azimuth), 5)

		#_generate object to hold #size photons
		if not record:		
			self.tau = 0	#_current distance from top
			self.k = kvector_init(zenith*d2r, azimuth*d2r).T.A
							#_current direction
			self.phi = azimuth * d2r	#_most recent phi
			self.theta = zenith*d2r		#_most recent theta
			self.history = zeros((3, 1))	#_each coord traveled
			self.weight = 1			#_ignore
			self.result = None		#_top, base, or absorbed
			self.scattered = 0 		#_# of times scattered
			self.live = True		#_still able to advance?
			self.figure = plt.figure(figsize=([16.,4.])) \
				if figure == None else figure
			self.axes = self.figure.add_subplot(111) \
				if axes==None else axes
			self.line = None

			#_associate photon with particular cloud
			self.cloud = clouds

		#_when passed a file name, returns previously run photon results
		else:
			#_break up record and put back in
			cols = record.split()
			shape = (len(cols[3:])/3, 3)	#_nscat, nk
			history	= array(cols[3:], dtype='f8').reshape(shape)
			
			#_if want to be able to load and continue advancing,
			# need to add back in calculation for k vect based
			# on hist[:,-2:]
			self.tau = float(cols[0])
			self.result = cols[1]
			self.scattered = int(cols[2])
			self.history = history.T
			self.live = False
			#_don't allow advance
			# mostly for plotting
			self.figure = plt.figure(figsize=([16.,4.])) \
				if figure == None else figure
			self.axes = self.figure.add_subplot(111) \
				if axes == None else axes
			self.line = None
		
		self.cloud.tot_gen += 1	
		
	def advance(self, ssa=1., weight=False, force_angle=None, tau_star=10.,
		lag=0, **kwargs):
		'''
		roll to see if photon's +12 to scattering prevents grue
		from eating it. Roll to see if absorbed (since all photons
		begin at cloud top, this begins each sequence)
		Roll to see which direction scattered
			F/B based on HG
			Azimuth uniform
		Roll to see distance
			Beer-Lambert Law
		'''
		from numpy.random import random as r
		from numpy import append, array, log, tile, pi
		from time import sleep
		from lib import henyey_greenstein

		dbg((ssa, force_angle, tau_star), 5)
			
		if not self.live:
			return
		
		#_roll to see how far until next interaction____________________
		self.tau = -log(1-r(1)[0])	#_THIS IS NOT TAAAAUUUUU

		#_update current location
		point = (self.history[:,-1] + self.tau*self.k).reshape(3, 1)
		self.history = append(self.history, point, axis=1)
		self.plot(**kwargs)
		
		#_check if escaped top or bottom, update flags
		if self.history[2,-1] < 1e-8:
			self.die('top', **kwargs)
			return
		elif self.history[2,-1] > tau_star:
			self.die('base', **kwargs)
			return
					
		#_roll to see if absorbed_______________________________________
		if r(1)[0] > ssa:
			self.die('absorbed', **kwargs)
			return
		
		#_roll to see in what direction________________________________
		if force_angle == None:
			phi = 2*pi*r(1)[0]	#_randomize azimuth
			theta = henyey_greenstein(r(1)[0],**kwargs)#_now theta
		else: #_force specific scattering angles for testing
			phi = force_angle[0]*d2r
			theta = force_angle[1]*d2r
		self.__turn__(theta,phi,**kwargs) #_update k

		self.scattered += 1
				
		#_update number of times scattered
		dbg(self.traveled(origin=False), 5)
		
		#_if we do a weighting check, add it here? above?
		#_update plot
		#self.plot(**kwargs) #_plotting here leaves out last segment
			
	def __turn__(self, theta, phi, **kwargs):
		''' augment direction of vector by theta,phi '''
		dbg((theta, phi), 5)
		from numpy import ones, dot, arange
		from lib import paxes
		from lib import kvector
		
		#_build 3x3 array describing coordinate system orthogonal to the
		# propagating ray and with x parallel to the model
		# horizontal plane
		aug = paxes(self.k, self.phi)

		#_get a unit vector describing new coordinate direction
		kpp = kvector(theta,phi)
		self.k = dot(aug,kpp).T.A

		#_update angle arrays
		self.phi = phi
		self.theta = theta
	
	def die(self, result, **kwargs): #_try to escape
		self.result = result
		self.live = False
		self.cloud.tot_fin += 1
		self.dump(**kwargs)
		
	def dump(self,outdir='.', outdat='photon.dat', fid=None, **kwargs):
		'''
		remove photon from record, write 
		TAU		RESULT	SCATTERED	HISTORY
		10.5f	10s 10i	(10.5f,10.5f,10.5f) * nscat+1
		index,	bool,	location of record to remove
		'''
		from numpy import array, where, empty
		from lib import mkdir_p

		dbg((outdir, outdat), 4)
		mkdir_p(outdir)
			
		fopened = True if fid != None else False
		if not fopened:
			fname = '/'.join((outdir,outdat))
			fid = open(fname,'a')
		
		loc = self.history.T.flatten()
		val = [self.tau,self.result,self.scattered]; val.extend(loc)
		
		#_make formatted print statement and write to file
		nk,nscat= self.history.shape
		fmt = '%10.5f%10s%10i '+'%10.5f%10.5f%10.5f' * nscat + '\n'
		line = fmt % tuple(val)

		fid.write(line)
		
		if not fopened:	
			fid.close()
			
	def traveled(self, origin=False):
		'''
		returns geometric distance traveled by photon, total or to entry
		origin,	bool, 	return only distance to origin,
				as the non-scattering photon flies
		'''
		from numpy import sqrt
		dbg(origin, 5)
			
		#_total meanderings	
		if not origin:
			segments = self.history[:,1:] - self.history[:,:-1]
			return sqrt((segments**2).sum())

		#_from cloud entrance
		else:		
			return sqrt((self.history[:,-1]**2).sum())

	def plot(self, idx=None, ax=None, outdir='.', imgfmt='png', rho=False,
		limit_plot=True, clouds=None, **kwargs):
		'''
		plot photon path
		idx,	int,		if supplied, only the listed photon
					paths will plot
		ax,	plt.artist,	matplotlib plot space
		outdir,	str,		output directory
		imgfmt,	str,		image file type
		rho,	bool,		instead of plotting path, plot density
		'''
		from numpy import arange, array, tan, cos, pi
		from lib.tools import mkdir_p

		mkdir_p(outdir)
		x, y, z = self.history
		if self.line == None:	#_new line
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
				x = [0, clouds.incident_x(cymax)]
				y = [0, cymax]
				clouds.ray.set_xdata(x)
				clouds.ray.set_ydata(y)
				clouds.base.set_xdata([cxmin,cxmax])
				clouds.top.set_xdata([cxmin,cxmax])	
			
		#_only call draw if this is the only line,
		# otherwise do it in a wrapper loop
		if len(self.axes.lines) == 1: 
			self.figure.canvas.draw()
			self.figure.show()
	
