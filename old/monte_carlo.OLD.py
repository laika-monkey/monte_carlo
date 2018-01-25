#!/usr/bin/env python
#############################################################################_80
# USAGE	:	./monte_carlo.py
#			When run from the command line, the namelist dictionary controls
#			experiments. The top level of namelist.keys() represent individual
#			experiments, in t#!/usr/bin/env python
#############################################################################_80
# USAGE	:	./monte_carlo.py
#			When run from the command line, the namelist dictionary controls
#			experiments. The top level of namelist.keys() represent individual
#			experiments, in turn pointing to dictionaries that contain the
#			keyward argument settings (single scatter albedo, asymmetry param,
#			et cetera). Output will be placed in ./{runname}
# AUTHOR:	Walter Raymond Sessions
#			Code for AOS 640 Project, UW-Madison, Spring 2014
# OUTPUT:	The life of a photon after encountering a scattering medium. ASCII
#			file with each line representing a single photon. This example does
#			not produce live plots as the bookkeeping became cumbersome. For
#			prettier live updates, use monte_carlo_single.py. The class.photon
#			in that one differs from the local class.photon; this one primarily
#			focuses on matrix math of multiple photon calculations at once while
#			_single.photon advances individual photons. Results have been found
#			to be insignificantly different.
#			Really should have written photon data to binary format...
#
# DEFAULTS:
#			g			0.85
#			ssa			1.0 
#			tau_star	10.0
#			nphoton		500,000.0
#			zenith		120 (0 == upward)
#			azimuth		0
#
#############################################################################_80
from numpy import array, arange, linspace, pi
import os

args = {
	'namelist' : {		#_experiment options
		'exp_1-nphot' 	: {			#_match beer-lamber law
			'tau_star'	: 1.,
			'nphot'		: array([1e0,1e1,1e2,1e3,1e4,1e5,1e6]),
			'xscale'	: 'log',								},
####		'nphot'		: array([5e5]*100),	},	#_to test if consistent
		'exp_2-tau_star': {			#_cloud thickness
			'tau_star'	: .1*2**(arange(1,12)-1.),
			'xscale'	: 'log',
			'limit_plot': [False, True],
			'texlabel'	: r'\tau^*',							},
		'exp_3-zenith'	: {			#_solar transit
			'zenith'	: 180-arange(0,81,10),
			'texlabel'	: '\Theta',								},
		'exp_4-ssa'		: {			#_absorptivity
			'ssa'		: array([.9,.99,.999,.9999,1.]),
			'ylim'		: [0,1],
			'texlabel'	: r'\~\omega',
			'xscale'	: 'log'									},
		'exp_5-g'		: {			#_by what factors
			'g'			: linspace(0,1,11),						},
###		'noweight'		: {
###			'weightmin'	: .01,		#_weight at which to ditch photon
###			'weight'	: False, },	#_do we use weights and split photons?
	},

	#_general options
	'simulation'	: False,	#_perform monte carlo calculations?
	'postprocess'	: True,		#_perform analysis?
	
	#_photon path plotting options
	'plot_photons'	: False,	#_takes forever for large nphot
	'show'			: False,	#_show photon plots or save them?
	'pmethod'		: 'fast',	#_fast|slow, fast is all blue, paths
		 						# with bins for > 1000, slow is multi color
								#_Turning off hexbin increases speed further
	'similarity'	: True,		#_use similarity transforms
	
	#_simulation options
	'size'			: 1000,		#_how many simultaneous photons?
	'nphot'			: 5e5,		#_how many total photons. Use exp1 results
								# unless it masks behavior of interest or 
								# Unless looping over nphot, this is used
}

#############################################################################_80
#_END_NAMELIST_DO_NOT_EDIT_BELOW_############################################_80
#############################################################################_80

debug 	= 1					#_high for more verbosity (roughly 1-7),
d2r 	= pi/180 			# although I mostly ignored this since it worked

#############################################################################_80
#_MAIN_######################################################################_80
#############################################################################_80
	
def run_main(simulation=True, postprocess=False, plot_photons=False,
	figsize=[16.,4.], size=100, nphot=100, 
	namelist={'test-nphot' : {'nphot' : [1e3]}}, **args): 
	'''
	simulation,		bool,		run monte carlo simulation?
	postprocess,	bool,		run analysis and counts on data
	plot_photons,	bool,		plot individual photon paths
	figsize,		list(2),	width and height of photonplot
	size,			integer,	number of photons to process simultaneously
	nphot,			integer,	total number of photons to process
	
	namelist,		dictionary,	in the style of kwargs. This is a poor example
								of it as there is no default value set, but
								it is required, making interactive calling
								of run_main cumbersome.
	'''
	import matplotlib.pyplot as plt
	from time import time
	import matplotlib
	from numpy import tile, nan
	
	similarity = args.get('similarity')
	
	#_monte carlo simulation________________________________________
	if simulation:
		for name, experiment in namelist.iteritems():
			knob = name.split('-')[-1]
			counts 	= {}							#_dictionary for phot counts

			for value in experiment[knob]:
				dbg(('simulating', name, value))
				
				if knob == 'nphot': nphot = value	#_wonk to clean out later
						
				#_output filenames for total photon data
				fname = '.'.join(('photon',str(value),'dat'))
				fname = '/'.join((name,fname))
				mkdir_p(name)
								
				#_pull out keyword args and update in the least elegant way
				kwargs = args.copy()
				kwargs.update(experiment)
				kwargs.update({	'outdir'	: name,
								'size'		: size,
								knob		: value,	})
							
				cloud = Cloud(**kwargs)
				
				#_generate output filehandle
				fid	= open(fname, 'w')
				kwargs.update({'fid' : fid})
	
				#_calc how many loops need to be made
				nloop = int(nphot/size)
				sizes = [size] * nloop
				if nphot % size: sizes.append(int(nphot % size)) 
					
				#_initialize a set of nitem photons, advance until absorbed
				rest, scat, wght = [], [], []
				for s in sizes:
					kwargs.update({'size' : s})
					phot = Photon(cloud=cloud, **kwargs)
					while phot.live: phot.advance(**kwargs)
					
					rest.extend(phot.result)
					scat.extend(phot.scattered)
					wght.extend(phot.weight)
		
				fid.close()
		
				#_weeeeeeeee
				phot.weight 	= array(wght)
				phot.scattered	= array(scat)
				phot.result 	= array(rest)

				#_count reflected, diffuse transmitted, absorbed, dir trans
				pref = phot.count_result('top')
				pdir = phot.count_direct()
				pdif = phot.count_result('base') - pdir
				pabs = phot.count_absorbed()
				counts[value] = (pref, pdir, pdif, pabs)	#_store

			#_write table with photon counts
			write_table(counts, experiment[knob], **kwargs)

	#_for each, plot knob vs {ref,tot_tran,abs for 2str[--],mc[-]}
	if postprocess:			#_plot test ref for test param values
		for name, experiment in namelist.iteritems():
			knob 	= name.split('-')[-1]
			label	= experiment['texlabel'] \
				if 'texlabel' in experiment else knob
			mkdir_p(name)
			
			#_{mc|ts}{value}{ref|trans|absp}, store reflectance, trans, absorp
			reftrnabs 	= tile(nan, (2, len(experiment[knob]), 3))
			tname 		= '/'.join((name,'table.txt'))
			mc_rdda, v	= read_table(tname)

			if similarity:
				sm_rdda = mc_rdda.copy()
				mc_rdda, v = read_table(re.sub('^sim', 'exp', tname))

			for value in experiment[knob]:
				if knob == 'nphot': nphot = value	#_wonk to clean out later
				dbg(('postprocessing', name, value))
				
				xlabel = '$' + label + '=' + str(value) + '$'
				kwargs = args.copy()
				kwargs.update(experiment)				#_add this exp params
				kwargs.update({ 'outdir'	: name,
								'label'		: xlabel,
								knob		: value, })	#_overwrite the knob arr
				
				#_fill in values for comparison
				n 	= experiment[knob].tolist().index(value)
				key = str(value)				
				reftrnabs[0,n,:] = 								\
					mc_rdda[key][0]/nphot, 						\
					(mc_rdda[key][1]+mc_rdda[key][2])/nphot, 	\
					mc_rdda[key][3]/nphot
				
				#_generate this experiment's cloud again
				cloud = Cloud(**kwargs)
				if not similarity:	#_if similarity off
					reftrnabs[1,n,:]= cloud.two_stream()
				else:				#_read in exp data instead
					reftrnabs[1,n,:] = 								\
						sm_rdda[key][0]/nphot, 						\
						(sm_rdda[key][1]+sm_rdda[key][2])/nphot, 	\
						sm_rdda[key][3]/nphot
				
				#_read in the data we just generated for post processing
				if plot_photons:
					#_read in photon.dat file for this experiment
					fname 	= '.'.join(('photon', str(value), 'dat'))
					fname 	= '/'.join((name, fname))
					data 	= Photon(fname=fname, cloud=cloud)
					
					#_photon path plot name
					pname 	= '.'.join(('photon', str(value), 'png'))
					pname	= '/'.join((name,pname))
					
					#_time plotting of photons
					t0 = time(); data.plot(pname=pname, **kwargs); t1 = time()
					dbg(('plotting time',t1-t0), 2)
											
			#_drop values
			kwargs['label'] = '$' + label + '$'
			
			#_plot comparison between mc and two stream
			if knob == 'ssa': kwargs.update({'reverse' : True})
			plot_comparison(experiment[knob], reftrnabs, **kwargs)

#############################################################################_80
#_POST_######################################################################_80
#############################################################################_80

def plot_comparison(x, y, ylim=False, xscale='linear', label='', outdir='.',
	reverse=False, **kwargs):
	'''really pointless but i'm running otu of time'''
	import matplotlib.pyplot as plt
	import matplotlib
	
	#_initialize figure
	ax = plt.subplot(111)
	ax.grid(True)
	
	if reverse:
		xlab 	= array(x.copy())
		x 		= 1-x[::-1]
		x[0] 	= x[1]/10
	
	#_plot results
	arg = {'linewidth':.7}
	ax.plot(x, y[0,:,0], 'b-', label='reflected', **arg)
	ax.plot(x, y[1,:,0], 'b--', **arg)
	ax.plot(x, y[0,:,1], 'k-', label='transmitted', **arg)
	ax.plot(x, y[1,:,1], 'k--', **arg)
	ax.plot(x, y[0,:,2], 'r-', label='absorbed',**arg)
	ax.plot(x, y[1,:,2], 'r--', **arg)
	ax.set_ylim(0,1)
	ax.set_xlabel(label)
	
	#_add legend, reduce plot area
	lgd = ax.legend(loc=1, bbox_to_anchor=(1.35,1.0))

	#_mess up axes
	ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	if ylim: ax.set_ylim(ylim)
	shrink_ticks(ax, 8)

	ax.set_xscale(xscale)
	if reverse:
		ax.set_xticks(x)
		ax.set_xticklabels(xlab)

	tighten_plot(ax, onlyx=True)

	#_generate output name and save
	pname = '/'.join((outdir, 'comparison.png'))
	dbg(('saving', pname))
	plt.savefig(pname, bbox_extra_artists=(lgd,), bbox_inches='tight')
	plt.close()
	
def write_table(counts, values, outdir='.', label='', **kwargs):
	'''
	put output tables in latex formatting.  Essentially the same info
	as the comparison graphs...
	results		dict({'mc'|'ts'}:{'r'|'a'|'t'})
	counts, 	dict,	[value] = (nreflected, ndirect, ndiffuse, ntransmitted)
	'''
	tname = '/'.join((outdir, 'table.txt'))
	dbg(('writing',tname))
	with open(tname, 'w') as f:
		nv = len(values)
		fmt = '%30s' + ' & %10s' * nv + ' \\\\\n\\hline\n'  #_header
		fm0 = '%30s' + ' & %10i' * nv + ' \\\\\n' 		#_...?

		#_header
		head = [label]
		head.extend([ '$'+str(val)+'$' for val in values ])
		f.write(fmt % tuple(head))

		#_reflected
		refl = ['reflected']
		[ refl.append(counts[val][0]) for val in values ]
		f.write(fm0 % tuple(refl))

		#_absorbed
		absp = ['absorbed']
		[ absp.append(counts[val][3]) for val in values ]
		f.write(fm0 % tuple(absp))

		#_transmitted{dir}
		tran = ['transmitted, direct']
		[ tran.append(counts[val][1]) for val in values ]
		f.write(fm0 % tuple(tran))

		#_transmitted{dif}
		tran = ['transmitted, diffuse']
		[ tran.append(counts[val][2]) for val in values ]
		f.write(fm0 % tuple(tran))
				
def read_table(fname='table.txt', **kwargs):
	'''
	put output tables in latex formatting.  Essentially the same info
	as the comparison graphs...
	results		dict({'mc'|'ts'}:{'r'|'a'|'t'})
	counts, 	dict,	[value] = (nreflected, ndirect, ndiffuse, ntransmitted)
	'''
	dbg(('reading', fname))
	fid 	= open(fname, 'r')
	lines 	= fid.readlines()
	fid.close()

	#_remove latex markings
	lines = [ line.replace('\\','') for line in lines ]
	lines = [ line.replace('\n','') for line in lines ]
	lines = [ line.replace('$','') for line in lines ]

	#_read header
	values 	= [ str(n).strip() for n in lines[0].split('&')[1:] ]
	nval 	= len(values)

	refl = [ float(n) for n in lines[2].split('&')[1:]]
	absp = [ float(n) for n in lines[3].split('&')[1:]]
	dire = [ float(n) for n in lines[4].split('&')[1:]]
	diff = [ float(n) for n in lines[5].split('&')[1:]]

	counts = {}
	for n in xrange(nval):
		counts[values[n]] = (refl[n], dire[n], diff[n], absp[n])
		
	return counts, values
	
#############################################################################_80
#_CLASSES_###################################################################_80
#############################################################################_80
		
class Photon(object):
	def __init__(self, size=0, fname=False, zenith=120, azimuth=0, figure=None,
		figsize=[16.,4.], label=None, cloud=None, **kwargs):
		'''
		basically a  numpy.recarray, only order is not muteable
		this should not be used as a post processing class, otherwise
		there will be giant, unwieldy arrays everywhere
		
		size,		int,	number of photons to handle at a time
		fname,		str,	name of photon output file
		zenith,		flt,	degrees from vertical of incidence (0==down)
		azimuth,	flt,	degrees from north of incidence (0==north)
		figure,		plt.fig	matplotlib artist object to be passed if there
							are multiple calls to this class, but only
							one desired plot area
							
		METHODS:
		advance		controls the actual MC part of the project, calculating
					one step forward for photon.
		turn		called by advance to calculate matrix translation
		drop		writes photon to file, removes from object
		plot		tries to plot, sometimes doesn't do a great job
		'''
		from numpy import array, zeros, ones, tile, nan
		import matplotlib.pyplot as plt
		dbg((size,fname),3)
		'''
		redo i/o with pickles and see if the attributes hold
		if it doesnt, switch to simple unform binary
		'''
		
		#_generate object to hold #size photons
		if not fname:		
			self.tau		= zeros((size))
			self.k			= tile(kvector_init(zenith*d2r,azimuth*d2r).T.A,(size,1)) 														
													#_current direction
			self.phi		= tile(azimuth*d2r,size)	#_most recent phi
			self.theta		= tile(zenith*d2r,size)		#_most recent theta
			self.history	= zeros((size,3,1))		#_each coordinate traveled
			self.weight		= ones((size))			#_ignore
			self.result		= array([None]*size)	#_top, base, or absorbed
			self.scattered	= zeros((size))			#_number of times scattered
			self.live		= True					#_still able to advance?
			self.zenith		= zenith				#_next three are static 
			self.azimuth	= azimuth				# initial values 
			self.cloud		= cloud
			self.tau_star	= cloud.tau_star	
			self.tmp_result	= []					#_don't ask...
			self.tmp_weight	= []
			self.tmp_scattered = []
##			self.label		= label
			
		#_when passed a file name, returns previously run photon results
		if fname:
			res = []; sca = []; his = []; wei = []
			with open(fname,'r') as fid:		#_open file and loop over lines
				tau_star, zenith, azimuth = fid.readline().split()
				for line in fid.readlines():	#_add data back into lists
					cols = line.split(',')
					wei.append(cols[0])
					res.append(cols[1].strip())
					sca.append(cols[2])
					shape = (int(sca[-1])+2,3)	#_nscat, nk
					his.append(array(cols[3:-1],dtype='f8').reshape(shape))
			self.tau_star	= float(tau_star)
			self.zenith		= float(zenith)
			self.azimuth	= float(azimuth)
			self.weight		= array(wei, dtype='f4')
			self.result 	= array(res, dtype=str)
			self.scattered 	= array(sca, dtype='i4')
			self.history	= his
			self.axes		= None
			self.live		= False					#_don't allow advance
													# mostly for plotting
			path, data		= fname.split('/')[-2:]
##			self.label		= ' '.join((path, data[7:-4]))

	def advance(self, weight=True, force_angle=None, top=1e-8, weightmin=0.01,
		**kwargs):
		'''
		Goes through the process of advancing a photon 1 extention event 
		through a cloud. If it gets absorbed or leaves the cloud, the photon
		is marked dead.
		
		All photons initialize at {0,0,0} and must be allowed to travel
		through the cloud before an extiction event to allow for the posibility
		of direct tranmission.
		
		Advancement is done by treating beers law as a PDF of the distance
		potentially traveled by individual photons. The uniform distribution in
		python is non-inclusive, [0,1), which we subtract from 1 preventing
		the logarithm of 0 from being taken.
		
		Scattering angles are also drawn from the uniform distribution and
		applied to an integrated form of the Henyay-Greenstein distribution 
		(theta) and isotropically over 2pi (azimuth).
		
		'''
		from numpy.random import random as r
		from numpy import append, array, log, tile
		if not self.live: return
		
		ssa = self.cloud.ssa
		tau_star = self.cloud.tau_star
		
		#_get shape of single point
		nphot, nk = self.k.shape #_update before any rolls

		#_roll to see how far until interaction_________________________________
		self.tau = -log(1-r(nphot)).reshape(nphot,1) #_not k[2]
		
		#_update current location
		point	 		= (self.history[:,:,-1]+self.tau*self.k).reshape(nphot,nk,1)
		self.history 	= append(self.history,point,axis=2)
		
		#_BASE CHECK: if leaves cloud, drop_____________________________________
		bot_index = self.history[:,2,-1] > tau_star
		self.result[bot_index] = 'base'
		self.drop(bot_index, **kwargs)
		if not self.live: return
		
		#_REFLECTED CHECK: check if escaped top_________________________________
		top_index = self.history[:,2,-1] < top
		self.result[top_index] = 'top'
		self.drop(top_index, **kwargs)
		if not self.live: return
		
		nphot, nk = self.k.shape

		#_roll to see if absorbed_______________________________________________
		absorption_index = r(nphot) > ssa 
		self.result[absorption_index] = 'absorbed'
		self.drop(absorption_index, **kwargs)
		if not self.live: return
		
		#_drop weight by ssa and drop if below 0.01_____________________________
		if weight:
			self.weight *= ssa
			weight_index = self.weight < weightmin
			self.result[weight_index] = 'fractured'
			self.drop(weight_index, **kwargs)
			if not self.live: return
			
		nphot, nk = self.k.shape
			
		#_roll to see in what direction_________________________________________
		if force_angle == None:
			phi		= 2*pi*r(nphot)							#_randomize azimuth
			theta 	= self.cloud.henyey_greenstein(r(nphot))
		else: 												#_force angles
			phi 	= tile(force_angle[0]*d2r, nphot)		# for testing
			theta 	= tile(force_angle[1]*d2r, nphot)
			
		self.turn(theta,phi,**kwargs)						#_use new angles 
															# to update k
		self.scattered += 1									#_inc scat count

	def turn(self, theta, phi, **kwargs):
		''' augment direction of vector by theta,phi '''
		from numpy import ones, dot, arange
		if not self.live: return
		
		np, nk 	= self.k.shape
		shape 	= [np,nk,3]
		
		#_build 3x3 array describing coordinate system orthogonal to the
		# propagating ray and with x parallel to the model horizontal plane
		aug = paxes(self.k, self.phi)
		
		#_get a unit vector describing new coordinate direction
		kpp 	= kvector(theta, phi).A.reshape(np, nk, 1)
		self.k 	= dot(aug, kpp)[arange(np),:,arange(np),0]

		#_update angle arrays
		self.phi 	= phi
		self.theta 	= theta
	
	'''dump all photons to file'''
	def drop_all(self,**kwargs):[ self.drop(0,**kwargs) for n in self.result ]
		
	def drop(self,idx,outdir='.',outdat='photon.dat',fid=None,**kwargs):
		'''
		remove photon from record, write 
		TAU		RESULT	SCATTERED	HISTORY
		10.5f	10s		10i 		(10.5f,10.5f,10.5f) * nscat+1
		
		index,	bool,	location of record to remove
		'''
		from numpy import array, where, empty
		if not self.live: return
		mkdir_p(outdir)

		#_make an integer index of what to drop
		idx = where(idx==True)[0]
		idx.sort() #_!!!
		
		fopened = True if fid != None else False
		if not fopened:							
			fname = '/'.join((outdir, outdat))
			fid = open(fname,'a')				#_open file handle
		
		nphot, nscat, nk = self.history.shape	#_setup record format
		kidx = range(nphot)						#_make keep list

		if not fid.tell():						#_write header at start
			vals = [self.tau_star,self.zenith,self.azimuth]
			fmt	 = '%10.3f '*3 + '\n'
			line = fmt % tuple(vals)
			fid.write(line)

		#_write output ascii file
		for n in idx[::-1]:
			res = self.result[n]
			sca = self.scattered[n]
			wei = self.weight[n]
			loc = self.history[n].T.flatten(); val = [wei,res,sca]
			val.extend(loc); np = loc.shape[0] / 3
		
			#_make formatted print statement and write to file
			fmt = '%10.3f,%10s,%10i,' + '%10.5f,%10.5f,%10.5f,' * np + '\n'
			line = fmt % tuple(val)

			fid.write(line)
			del kidx[n]					#_remove offending index
		kidx = array(kidx)				#_convert to ndarray for use as index
			
		if not fopened: fid.close()		#_close file handle

		self.tmp_result.extend(self.result[idx])
		self.tmp_weight.extend(self.weight[idx])
		self.tmp_scattered.extend(self.scattered[idx])

		#_loop over relavent attributes and remove them
		attrs = ['k','history','weight','result','scattered','theta','phi']
		if self.nlive() == 0: 		#_kill when done
			[ self.__setattr__(n, empty(0)) for n in attrs ]
			self.live 		= False
			self.weight 	= array(self.tmp_weight)
			self.result 	= array(self.tmp_result)
			self.scattered	= array(self.tmp_scattered)
		else:						#_drop individual photo data
			[ self.__setattr__(n, self.__getattribute__(n)[kidx]) for n in attrs 
				if len(kidx) > 0 ]
			
	def traveled(self,origin=False):
		'''
		returns geometric distance traveled by photon, total or to entry
		origin,	bool,	return only distance to origin, as the non-scattering
						photon flies
		'''
		from numpy import sqrt
		dbg(origin,5)
		if not self.live: return
					
		if not origin:	#_total meanderings
			segments = self.history[:,:,1:] - self.history[:,:,:-1]
			return sqrt((segments**2).sum(axis=1)).sum(axis=1)
		else:			#_from cloud entrance
			return sqrt((self.history[:,:,-1]**2).sum(axis=1))

	def count_result(self,res): 
		'''
		count photons exiting top, bottom, or absorbed
		returned value is the sum of the weight - if weight
		is turned off, then it's just a count of photons
		
		'''
		from numpy import where
		idx = (self.result == res)
		return self.weight[idx].sum()

	def count_direct(self):
		'''returns integer value for photons directly transmitted through cloud'''
		from numpy import append
		size = self.result.size
		bool_base = (self.result == 'base').reshape(size,1)
		bool_scat = (self.scattered == 0).reshape(size,1)
		#_weight doesnt come into play for t_dir, no scattering, no fracturing
		return append(bool_base,bool_scat,axis=1).all(axis=1).sum()

	def count_absorbed(self):
		''' 
		count the aborbed photon weight 
		Any fractional splitting goes into this
		'''
		#_just count the number, the bits went into absorbed
		fracture = (self.result == 'fractured').sum()
		absorbed = (self.result == 'absorbed').sum()
		
		#_get the parts of the scattered out that should be in cloud
		top = (self.result == 'top').sum() - self.count_result('top')
		bot = (self.result == 'base').sum() - self.count_result('base')
		
		return fracture + absorbed + top + bot

	def nlive(self):
		'''count the number of photons still bouncing around'''
		return self.result.tolist().count(None)
		
	def plot(self, outdir='.', imgfmt='png', show=False, label='',
		limit_plot=[False, True], pmethod='fast', pname=None, 
		min_for_hexplot=1e3, exp_val='', **kwargs):
		'''
		plot photon path
		idx,		int,		if supplied, only the listed photon paths will plot
		ax,			plt.artist,	matplotlib plot space
		outdir,		str,		output directory
		imgfmt,		str,		image file type
		limit_plot,	bool(2),	don't extend plot area with data, (x, y)
		'''
		from numpy import arange, array, tan, cos, pi, vstack, tile
		import matplotlib.pyplot as plt
		from matplotlib import rcParams
		from matplotlib.cm import YlOrRd_r as YOR
		from numpy.ma import masked_outside, masked_where
		from libtools import shrink_ticks

		if self.live: return #_this version only plots the damned
		
		self.figure = plt.figure(figsize=[16., 4.])
		self.axes = self.figure.add_subplot(111)
		delt = self.tau_star/5.			
		self.axes.set_ylim(self.tau_star+delt, -delt)
		self.axes.set_xlim(-30, 30)
	
		#_if large enough, plot binned interaction sites
		nphoton	= len(self.history)
		if nphoton >= min_for_hexplot:
			#_stack all points and mask outside cloud
			zyx = vstack(self.history)						#_all points 
			zh	= masked_outside(zyx[:,2], 0, self.tau_star)
			xh 	= masked_where(zh.mask, zyx[:,0])
			hx 	= self.axes.hexbin = self.axes.hexbin(xh, zh, mincnt=10,
				bins='log', cmap=YOR, zorder=2, gridsize=200, vmax=4.)
			cb	= plt.colorbar(hx, orientation='vertical', aspect=8)
			cb.set_label('log10(N)')
		
		#_plot them
		if pmethod == 'slow':
			#_multiple colors, but slow as dirt
			[self.axes.plot(l.T[0], l.T[2], linewidth=.3) for l in self.history]
		elif pmethod == 'fast':
			#_all blue, but passed as a single line with segments separated by
			# NoneType objects
			jump	= [ nhist.shape[0] for nhist in self.history ]
			nt		= sum(jump)

			#_initialize array
			xyz	= tile([None,None,None], (nt+nphoton, 1))	#_spaced out
			last = 0
			for n in arange(nphoton):
				this 			= last + jump[n]
				xyz[last:this,:]= self.history[n].copy()
				last 			= this + 1
				xp				= xyz.T[0]
				zp 				= xyz.T[2]
				
			if nt+nphoton > 1e5:
				dbg(('modifying agg.path.chunksize! check for artifacts'), 3)
				rcParams['agg.path.chunksize'] = 1000
			
			#_plot photon paths
			self.axes.plot(xp, zp, linewidth=.3, zorder=1)
		
		cxmin, cxmax = self.axes.xaxis.get_view_interval()
		dxmin, dxmax = self.axes.xaxis.get_data_interval()
		cymin, cymax = self.axes.yaxis.get_view_interval()
		dymin, dymax = self.axes.yaxis.get_data_interval()
		
		#_expands x,y axes as needed by photon travel
		if not limit_plot[0]:	#_X
			xmax = max([cxmax,dxmax])
			xmin = min([cxmin,dxmin])
			self.axes.set_xlim(xmin,xmax)
		if not limit_plot[1]:	#_Y
			ymax = max([cymax,dymax])
			ymin = min([cymin,dymin])
			self.axes.set_ylim(ymax,ymin)
			
		#_get current size of axes
		cymin, cymax = self.axes.yaxis.get_view_interval()	
		cxmin, cxmax = self.axes.xaxis.get_view_interval()
			
		#_plot incident ray and cloud boundaries
		x = [0]; y = [0,cymax]
		x.append(-cymax*tan(self.zenith*d2r-pi)*cos(self.azimuth*d2r))
		args = { 'linewidth' : 0.9, 'color' : 'k' }
		self.axes.plot(x, y, linewidth=3, linestyle='--', color='k', zorder=3)		
		self.axes.plot([cxmin,cxmax], [0,0], zorder=3, **args)
		self.axes.plot([cxmin,cxmax], [self.tau_star,self.tau_star], 
			zorder=3, **args)
			
		#_gussy it up
		shrink_ticks(self.axes, 8)
		self.axes.set_ylabel(r'$\tau$')
		self.axes.set_xlabel(label)
		
		#_put it somewhere
		if show: 
			plt.show()
		else:
			mkdir_p(outdir)
			if pname == None:
				fname = '/'.join((outdir,'.'.join(('photon',imgfmt))))
			dbg((pname,outdir),l=2)
			plt.savefig(pname)
			
		#_the latter should release memory, but it's not. The first shouldn't
		self.figure.clear()
		plt.close()
		
class Cloud(object):
	def __init__(self, g=0.85, ssa=1., tau_star=10., similarity=False, **kwargs):
		'''
		g			float,	asymmetry parameter [0,1]
		ssa			float,	single scatter albedo, [0,1]
		tau_star	float,	optical depth of homogeneous layer, [0,inf)
		
		Having these attributed to photons was easier in the original code,
		but I prefer associating these with a class describing the medium
		'''
		self.g 			= g 		if not similarity else 0
		self.ssa 		= ssa		if not similarity else omega_prime(ssa, g)
		self.tau_star 	= tau_star	if not similarity else \
									tau_prime(ssa, g, tau_star)
		
	def two_stream(self):	
		'''
		Use two stream parameterization model to calculate 
		azimuthically independent values for reflectance and transmittance
	
		ssa,	float,	[0,1], single scatter albedo
		g,		float,	[-1,1],	asymmetry parameter
	
		returns:
		refl	float,	[0,1],	reflectivity
		tran	float,	[0,1],	transmittance 
		'''
		from numpy import sqrt, exp
	
		#_calculate gamma, Petty 13.25
		G = 2*sqrt(1-self.ssa)*sqrt(1-self.ssa*self.g)
		
		#_calculate refl and trans depending on ssa value
		if self.ssa == 1:
			#_Petty 13.63 and 64
			refl = (1-self.g)*self.tau_star / (1+(1-self.g)*self.tau_star)
			tran = 1. / (1+(1-self.g)*self.tau_star)			
##			if refl+tran != 1 or (refl+tran).mean() != 1: 
##				raise RuntimeError, 'photons magically disappearing ' +\
##					' '.join((str(refl), str(trans))) 

		elif self.ssa < 1:
			#_calculate semi-invinit cloud reflectivity, Petty 13.45
			r8 	= sqrt(1-self.ssa*self.g) - sqrt(1-self.ssa)
			r8 /= sqrt(1-self.ssa*self.g) + sqrt(1-self.ssa)
			
			#_Petty 13.65
			refl = r8 * (exp(G*self.tau_star) - exp(-G*self.tau_star))
			refl /= (exp(G*self.tau_star) - r8**2*exp(-G*self.tau_star))
			tran = (1.-r8**2) / (exp(G*self.tau_star) \
				- r8**2*exp(-G*self.tau_star))
				
		else:
			raise RuntimeError, 'invalid single scatter value ' + str(ssa)

		return refl, tran, 1-tran-refl #_r, t, a
		
	def henyey_greenstein(self, r):
		'''
		calculate scattering theta in accordance with Henyay Greenstein dist
		r,	flt,	random number, [0,1] INCLUSIVE
		g,	flt,	asymettry parameter.  0 => equal f/b hemisphere, 1, all forward

		I liked this better as an external function, but will cram it into
		a class for this project.

		Binzoni, T., et al. "The use of the Henyey-Greenstein phase function
		in Monte Carlo simulations in biomedical optics." Physics in medicine 
		and biology 51.17 (2006): N313.
		'''
		from numpy import arccos, isnan
		if self.g == 0:	#_give random isotropic direct
			cos_theta = 2.*r-1
		else:		#_give hemispherically preferred direction
			cos_theta = (1+self.g**2-((1-self.g**2)/\
				(1-self.g+2.*self.g*r))**2)/2./self.g
		return arccos(cos_theta)
		
#############################################################################_80
#_TOOLS_#####################################################################_80
#############################################################################_80

def beers(tau=1., theta=0, **kwargs): 
	from numpy import exp, cos, pi, abs
	mu = abs(cos(theta*pi/180.))
	return exp(-tau/mu)

def kvector_init(theta, phi):
	'''
	return unit vector pointing in direction defined by theta and phi 
	In model space, this is relative to the vertical.  In photon space,
	it's relative to the propoagation direction
	theta,phi	float,	directions about axis of propagation IN RADIANS
	
	'''
	from numpy import matrix,cos,sin
	return matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),-cos(theta)]).T

def paxes(k, phi):
	''' vector is horizontal in model space and orthogonal to propagation '''
	from numpy import matrix, cross, dot, array, sin, cos, tile, zeros
	from numpy.linalg import norm

	#_if k0**2+k1**2 == 0, apply different transform
	idx = array((k[:,0]**2 + k[:,1]**2) == 0)		#_indices where zero

	z 	= matrix([0,0,1]).T							#_vertical model space
	p 	= cross(z.T,k)								#_ortho to vert and prop
	mag = matrix([ norm(x) for x in p ])			#_magnitude 
	#_maybe put in a check for zero magnitude vectors (perfectly horizontal)

	x = (p / mag.T)		 							#_get unit vector
	if idx.sum() != 0: 								#_replace nans
		dbg('direction parallel to vertical', 3)
		xp = array(-sin(phi[idx]))
		yp = array( cos(phi[idx]))
		zp = tile(0, len(idx))
		x[idx] = array([xp,yp,zp]).T

	y = matrix(cross(k,x))							#_calc y unit vector

	nphot 		= k.shape[0]
	aug 		= zeros((nphot,3,3))
	aug[:,:,0] 	= x
	aug[:,:,1] 	= y 
	aug[:,:,2] 	= k # tile(z.T,(nphot,1))
	
	return aug
			
def kvector(theta, phi):
	'''
	return unit vector in direction of theta/phi
	Initially, this will be in relationship to standard 3d cartesian coords,
	but as the propagation path shifts, it will be in relation to the last path
	
	theta,phi	float,	directions about axis of propagation IN RADIANS
	'''
	from numpy import matrix,cos,sin
	return matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)]).T
	
def plot_hg(polar=True):
	'''plot what's coming out of the HG function '''
	raise RuntimeWarning, 'NOT ACTUALLY H-G'
	import matplotlib.pyplot as plt
	from numpy import cos
	theta = linspace(0,2*pi,1000)
	theta = linspace(-pi,pi,1000)
	
	for g in linspace(0,1,3):
		pT 		= (1-g*g)/(1+g*g-2*g*cos(theta)**(3/2))/4/pi
		label 	= '='.join(('g',str(g)))
		ax 		= plt.subplot(111, polar=polar)
		
		ax.plot(theta, pT, label=label, linewidth=0.4)	
		
	if not polar: 
##		ax.set_xlim(-pi,pi)
		plt.legend()
		
	plt.title('Henyey-Greenstein Phase Function')
	plt.show()

#_similarity theory
def omega_prime(ssa, g): return float(1.-g)*ssa/(1.-g*ssa)
def tau_prime(ssa, g, tau_star): return (1.-g*ssa)*tau_star
	
def setup_groups(objects,nproc=1,**kwargs):
	'''
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
	import errno
	dbg(path,5)
	def mkdir(p):
		try: 
			os.makedirs(p)
			os.chmod(p,0775)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST: pass
			else: raise

	if hasattr( '__iter__', path ): [ mkdir(p) for p in path ]
	else: mkdir(path)

def photon_test():
	''' sanity check against Petty's matrix worksheet '''
	exp = [.387,.461,-.799,.787,.388,-.479]
	
	p = Photon(size=1, zenith=37, azimuth=50, cloud=Cloud())
	got = p.k.flatten().tolist()
	
	p.advance(force_angle=[135,30], top=-10)
	got.extend(p.k.flatten().tolist())
	
	print '%5.3f '*6 % tuple(exp)
	print '%5.3f '*6 % tuple(got)
	
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
		file	 = file.split('/')[-1].split('.')[0]
		scream 	= '[%s.%s.%i] %s' % (file,method,line,msg)
		
		if not err: print scream
		else: raise RuntimeError, scream

def to_string(a):
	if hasattr(a,'__iter__'): return [ str(b) for b in a ]
	else: return str(a)

if __name__ == '__main__': 
	#_change the namelist experiment names if running similarity theory tests
	if args['similarity']:			#_change output dir name if similarity
		import re
		newargs = {}				#_new kwargs
									#_could have just put list with 
									# ['exp','sim'] to key loop over.
		for key, values in args['namelist'].iteritems():
			key = re.sub('^exp', 'sim', key)
			newargs.update({key : values})
		args['namelist'] = newargs
	run_main(**args)urn pointing to dictionaries that contain the
#			keyward argument settings (single scatter albedo, asymmetry param,
#			et cetera). Output will be placed in ./{runname}
# AUTHOR:	Walter Raymond Sessions
#			Code for AOS 640 Project, UW-Madison, Spring 2014
# OUTPUT:	The life of a photon after encountering a scattering medium. ASCII
#			file with each line representing a single photon. This example does
#			not produce live plots as the bookkeeping became cumbersome. For
#			prettier live updates, use monte_carlo_single.py. The class.photon
#			in that one differs from the local class.photon; this one primarily
#			focuses on matrix math of multiple photon calculations at once while
#			_single.photon advances individual photons. Results have been found
#			to be insignificantly different.
#			Really should have written photon data to binary format...
#
# DEFAULTS:
#			g			0.85
#			ssa			1.0 
#			tau_star	10.0
#			nphoton		500,000.0
#			zenith		120 (0 == upward)
#			azimuth		0
#
#############################################################################_80
from numpy import array, arange, linspace, pi
import os

args = {
	'namelist' : {		#_experiment options
		'exp_1-nphot' 	: {			#_match beer-lamber law
			'tau_star'	: 1.,
			'nphot'		: array([1e0,1e1,1e2,1e3,1e4,1e5,1e6]),
			'xscale'	: 'log',								},
####		'nphot'		: array([5e5]*100),	},	#_to test if consistent
		'exp_2-tau_star': {			#_cloud thickness
			'tau_star'	: .1*2**(arange(1,12)-1.),
			'xscale'	: 'log',
			'limit_plot': [False, True],
			'texlabel'	: r'\tau^*',							},
		'exp_3-zenith'	: {			#_solar transit
			'zenith'	: 180-arange(0,81,10),
			'texlabel'	: '\Theta',								},
		'exp_4-ssa'		: {			#_absorptivity
			'ssa'		: array([.9,.99,.999,.9999,1.]),
			'ylim'		: [0,1],
			'texlabel'	: r'\~\omega',
			'xscale'	: 'log'									},
		'exp_5-g'		: {			#_by what factors
			'g'			: linspace(0,1,11),						},
###		'noweight'		: {
###			'weightmin'	: .01,		#_weight at which to ditch photon
###			'weight'	: False, },	#_do we use weights and split photons?
	},

	#_general options
	'simulation'	: False,	#_perform monte carlo calculations?
	'postprocess'	: True,		#_perform analysis?
	
	#_photon path plotting options
	'plot_photons'	: False,	#_takes forever for large nphot
	'show'			: False,	#_show photon plots or save them?
	'pmethod'		: 'fast',	#_fast|slow, fast is all blue, paths
		 						# with bins for > 1000, slow is multi color
								#_Turning off hexbin increases speed further
	'similarity'	: True,		#_use similarity transforms
	
	#_simulation options
	'size'			: 1000,		#_how many simultaneous photons?
	'nphot'			: 5e5,		#_how many total photons. Use exp1 results
								# unless it masks behavior of interest or 
								# Unless looping over nphot, this is used
}

#############################################################################_80
#_END_NAMELIST_DO_NOT_EDIT_BELOW_############################################_80
#############################################################################_80

debug 	= 1					#_high for more verbosity (roughly 1-7),
d2r 	= pi/180 			# although I mostly ignored this since it worked

#############################################################################_80
#_MAIN_######################################################################_80
#############################################################################_80
	
def run_main(simulation=True, postprocess=False, plot_photons=False,
	figsize=[16.,4.], size=100, nphot=100, 
	namelist={'test-nphot' : {'nphot' : [1e3]}}, **args): 
	'''
	simulation,		bool,		run monte carlo simulation?
	postprocess,	bool,		run analysis and counts on data
	plot_photons,	bool,		plot individual photon paths
	figsize,		list(2),	width and height of photonplot
	size,			integer,	number of photons to process simultaneously
	nphot,			integer,	total number of photons to process
	
	namelist,		dictionary,	in the style of kwargs. This is a poor example
								of it as there is no default value set, but
								it is required, making interactive calling
								of run_main cumbersome.
	'''
	import matplotlib.pyplot as plt
	from time import time
	import matplotlib
	from numpy import tile, nan
	
	similarity = args.get('similarity')
	
	#_monte carlo simulation________________________________________
	if simulation:
		for name, experiment in namelist.iteritems():
			knob = name.split('-')[-1]
			counts 	= {}							#_dictionary for phot counts

			for value in experiment[knob]:
				dbg(('simulating', name, value))
				
				if knob == 'nphot': nphot = value	#_wonk to clean out later
						
				#_output filenames for total photon data
				fname = '.'.join(('photon',str(value),'dat'))
				fname = '/'.join((name,fname))
				mkdir_p(name)
								
				#_pull out keyword args and update in the least elegant way
				kwargs = args.copy()
				kwargs.update(experiment)
				kwargs.update({	'outdir'	: name,
								'size'		: size,
								knob		: value,	})
							
				cloud = Cloud(**kwargs)
				
				#_generate output filehandle
				fid	= open(fname, 'w')
				kwargs.update({'fid' : fid})
	
				#_calc how many loops need to be made
				nloop = int(nphot/size)
				sizes = [size] * nloop
				if nphot % size: sizes.append(int(nphot % size)) 
					
				#_initialize a set of nitem photons, advance until absorbed
				rest, scat, wght = [], [], []
				for s in sizes:
					kwargs.update({'size' : s})
					phot = Photon(cloud=cloud, **kwargs)
					while phot.live: phot.advance(**kwargs)
					
					rest.extend(phot.result)
					scat.extend(phot.scattered)
					wght.extend(phot.weight)
		
				fid.close()
		
				#_weeeeeeeee
				phot.weight 	= array(wght)
				phot.scattered	= array(scat)
				phot.result 	= array(rest)

				#_count reflected, diffuse transmitted, absorbed, dir trans
				pref = phot.count_result('top')
				pdir = phot.count_direct()
				pdif = phot.count_result('base') - pdir
				pabs = phot.count_absorbed()
				counts[value] = (pref, pdir, pdif, pabs)	#_store

			#_write table with photon counts
			write_table(counts, experiment[knob], **kwargs)

	#_for each, plot knob vs {ref,tot_tran,abs for 2str[--],mc[-]}
	if postprocess:			#_plot test ref for test param values
		for name, experiment in namelist.iteritems():
			knob 	= name.split('-')[-1]
			label	= experiment['texlabel'] \
				if 'texlabel' in experiment else knob
			mkdir_p(name)
			
			#_{mc|ts}{value}{ref|trans|absp}, store reflectance, trans, absorp
			reftrnabs 	= tile(nan, (2, len(experiment[knob]), 3))
			tname 		= '/'.join((name,'table.txt'))
			mc_rdda, v	= read_table(tname)

			if similarity:
				sm_rdda = mc_rdda.copy()
				mc_rdda, v = read_table(re.sub('^sim', 'exp', tname))

			for value in experiment[knob]:
				if knob == 'nphot': nphot = value	#_wonk to clean out later
				dbg(('postprocessing', name, value))
				
				xlabel = '$' + label + '=' + str(value) + '$'
				kwargs = args.copy()
				kwargs.update(experiment)				#_add this exp params
				kwargs.update({ 'outdir'	: name,
								'label'		: xlabel,
								knob		: value, })	#_overwrite the knob arr
				
				#_fill in values for comparison
				n 	= experiment[knob].tolist().index(value)
				key = str(value)				
				reftrnabs[0,n,:] = 								\
					mc_rdda[key][0]/nphot, 						\
					(mc_rdda[key][1]+mc_rdda[key][2])/nphot, 	\
					mc_rdda[key][3]/nphot
				
				#_generate this experiment's cloud again
				cloud = Cloud(**kwargs)
				if not similarity:	#_if similarity off
					reftrnabs[1,n,:]= cloud.two_stream()
				else:				#_read in exp data instead
					reftrnabs[1,n,:] = 								\
						sm_rdda[key][0]/nphot, 						\
						(sm_rdda[key][1]+sm_rdda[key][2])/nphot, 	\
						sm_rdda[key][3]/nphot
				
				#_read in the data we just generated for post processing
				if plot_photons:
					#_read in photon.dat file for this experiment
					fname 	= '.'.join(('photon', str(value), 'dat'))
					fname 	= '/'.join((name, fname))
					data 	= Photon(fname=fname, cloud=cloud)
					
					#_photon path plot name
					pname 	= '.'.join(('photon', str(value), 'png'))
					pname	= '/'.join((name,pname))
					
					#_time plotting of photons
					t0 = time(); data.plot(pname=pname, **kwargs); t1 = time()
					dbg(('plotting time',t1-t0), 2)
											
			#_drop values
			kwargs['label'] = '$' + label + '$'
			
			#_plot comparison between mc and two stream
			if knob == 'ssa': kwargs.update({'reverse' : True})
			plot_comparison(experiment[knob], reftrnabs, **kwargs)

#############################################################################_80
#_POST_######################################################################_80
#############################################################################_80

def plot_comparison(x, y, ylim=False, xscale='linear', label='', outdir='.',
	reverse=False, **kwargs):
	'''really pointless but i'm running otu of time'''
	import matplotlib.pyplot as plt
	from libtools import tighten_plot, shrink_ticks
	import matplotlib
	
	#_initialize figure
	ax = plt.subplot(111)
	ax.grid(True)
	
	if reverse:
		xlab 	= array(x.copy())
		x 		= 1-x[::-1]
		x[0] 	= x[1]/10
	
	#_plot results
	arg = {'linewidth':.7}
	ax.plot(x, y[0,:,0], 'b-', label='reflected', **arg)
	ax.plot(x, y[1,:,0], 'b--', **arg)
	ax.plot(x, y[0,:,1], 'k-', label='transmitted', **arg)
	ax.plot(x, y[1,:,1], 'k--', **arg)
	ax.plot(x, y[0,:,2], 'r-', label='absorbed',**arg)
	ax.plot(x, y[1,:,2], 'r--', **arg)
	ax.set_ylim(0,1)
	ax.set_xlabel(label)
	
	#_add legend, reduce plot area
	lgd = ax.legend(loc=1, bbox_to_anchor=(1.35,1.0))

	#_mess up axes
	ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	if ylim: ax.set_ylim(ylim)
	shrink_ticks(ax, 8)

	ax.set_xscale(xscale)
	if reverse:
		ax.set_xticks(x)
		ax.set_xticklabels(xlab)

	tighten_plot(ax, onlyx=True)

	#_generate output name and save
	pname = '/'.join((outdir, 'comparison.png'))
	dbg(('saving', pname))
	plt.savefig(pname, bbox_extra_artists=(lgd,), bbox_inches='tight')
	plt.close()
	
def write_table(counts, values, outdir='.', label='', **kwargs):
	'''
	put output tables in latex formatting.  Essentially the same info
	as the comparison graphs...
	results		dict({'mc'|'ts'}:{'r'|'a'|'t'})
	counts, 	dict,	[value] = (nreflected, ndirect, ndiffuse, ntransmitted)
	'''
	tname = '/'.join((outdir, 'table.txt'))
	dbg(('writing',tname))
	with open(tname, 'w') as f:
		nv = len(values)
		fmt = '%30s' + ' & %10s' * nv + ' \\\\\n\\hline\n'  #_header
		fm0 = '%30s' + ' & %10i' * nv + ' \\\\\n' 		#_...?

		#_header
		head = [label]
		head.extend([ '$'+str(val)+'$' for val in values ])
		f.write(fmt % tuple(head))

		#_reflected
		refl = ['reflected']
		[ refl.append(counts[val][0]) for val in values ]
		f.write(fm0 % tuple(refl))

		#_absorbed
		absp = ['absorbed']
		[ absp.append(counts[val][3]) for val in values ]
		f.write(fm0 % tuple(absp))

		#_transmitted{dir}
		tran = ['transmitted, direct']
		[ tran.append(counts[val][1]) for val in values ]
		f.write(fm0 % tuple(tran))

		#_transmitted{dif}
		tran = ['transmitted, diffuse']
		[ tran.append(counts[val][2]) for val in values ]
		f.write(fm0 % tuple(tran))
				
def read_table(fname='table.txt', **kwargs):
	'''
	put output tables in latex formatting.  Essentially the same info
	as the comparison graphs...
	results		dict({'mc'|'ts'}:{'r'|'a'|'t'})
	counts, 	dict,	[value] = (nreflected, ndirect, ndiffuse, ntransmitted)
	'''
	dbg(('reading', fname))
	fid 	= open(fname, 'r')
	lines 	= fid.readlines()
	fid.close()

	#_remove latex markings
	lines = [ line.replace('\\','') for line in lines ]
	lines = [ line.replace('\n','') for line in lines ]
	lines = [ line.replace('$','') for line in lines ]

	#_read header
	values 	= [ str(n).strip() for n in lines[0].split('&')[1:] ]
	nval 	= len(values)

	refl = [ float(n) for n in lines[2].split('&')[1:]]
	absp = [ float(n) for n in lines[3].split('&')[1:]]
	dire = [ float(n) for n in lines[4].split('&')[1:]]
	diff = [ float(n) for n in lines[5].split('&')[1:]]

	counts = {}
	for n in xrange(nval):
		counts[values[n]] = (refl[n], dire[n], diff[n], absp[n])
		
	return counts, values
	
#############################################################################_80
#_CLASSES_###################################################################_80
#############################################################################_80
		
class Photon(object):
	def __init__(self, size=0, fname=False, zenith=120, azimuth=0, figure=None,
		figsize=[16.,4.], label=None, cloud=None, **kwargs):
		'''
		basically a  numpy.recarray, only order is not muteable
		this should not be used as a post processing class, otherwise
		there will be giant, unwieldy arrays everywhere
		
		size,		int,	number of photons to handle at a time
		fname,		str,	name of photon output file
		zenith,		flt,	degrees from vertical of incidence (0==down)
		azimuth,	flt,	degrees from north of incidence (0==north)
		figure,		plt.fig	matplotlib artist object to be passed if there
							are multiple calls to this class, but only
							one desired plot area
							
		METHODS:
		advance		controls the actual MC part of the project, calculating
					one step forward for photon.
		turn		called by advance to calculate matrix translation
		drop		writes photon to file, removes from object
		plot		tries to plot, sometimes doesn't do a great job
		'''
		from numpy import array, zeros, ones, tile, nan
		import matplotlib.pyplot as plt
		dbg((size,fname),3)
		'''
		redo i/o with pickles and see if the attributes hold
		if it doesnt, switch to simple unform binary
		'''
		
		#_generate object to hold #size photons
		if not fname:		
			self.tau		= zeros((size))
			self.k			= tile(kvector_init(zenith*d2r,azimuth*d2r).T.A,(size,1)) 														
													#_current direction
			self.phi		= tile(azimuth*d2r,size)	#_most recent phi
			self.theta		= tile(zenith*d2r,size)		#_most recent theta
			self.history	= zeros((size,3,1))		#_each coordinate traveled
			self.weight		= ones((size))			#_ignore
			self.result		= array([None]*size)	#_top, base, or absorbed
			self.scattered	= zeros((size))			#_number of times scattered
			self.live		= True					#_still able to advance?
			self.zenith		= zenith				#_next three are static 
			self.azimuth	= azimuth				# initial values 
			self.cloud		= cloud
			self.tau_star	= cloud.tau_star	
			self.tmp_result	= []					#_don't ask...
			self.tmp_weight	= []
			self.tmp_scattered = []
##			self.label		= label
			
		#_when passed a file name, returns previously run photon results
		if fname:
			res = []; sca = []; his = []; wei = []
			with open(fname,'r') as fid:		#_open file and loop over lines
				tau_star, zenith, azimuth = fid.readline().split()
				for line in fid.readlines():	#_add data back into lists
					cols = line.split(',')
					wei.append(cols[0])
					res.append(cols[1].strip())
					sca.append(cols[2])
					shape = (int(sca[-1])+2,3)	#_nscat, nk
					his.append(array(cols[3:-1],dtype='f8').reshape(shape))
			self.tau_star	= float(tau_star)
			self.zenith		= float(zenith)
			self.azimuth	= float(azimuth)
			self.weight		= array(wei, dtype='f4')
			self.result 	= array(res, dtype=str)
			self.scattered 	= array(sca, dtype='i4')
			self.history	= his
			self.axes		= None
			self.live		= False					#_don't allow advance
													# mostly for plotting
			path, data		= fname.split('/')[-2:]
##			self.label		= ' '.join((path, data[7:-4]))

	def advance(self, weight=True, force_angle=None, top=1e-8, weightmin=0.01,
		**kwargs):
		'''
		Goes through the process of advancing a photon 1 extention event 
		through a cloud. If it gets absorbed or leaves the cloud, the photon
		is marked dead.
		
		All photons initialize at {0,0,0} and must be allowed to travel
		through the cloud before an extiction event to allow for the posibility
		of direct tranmission.
		
		Advancement is done by treating beers law as a PDF of the distance
		potentially traveled by individual photons. The uniform distribution in
		python is non-inclusive, [0,1), which we subtract from 1 preventing
		the logarithm of 0 from being taken.
		
		Scattering angles are also drawn from the uniform distribution and
		applied to an integrated form of the Henyay-Greenstein distribution 
		(theta) and isotropically over 2pi (azimuth).
		
		'''
		from numpy.random import random as r
		from numpy import append, array, log, tile
		if not self.live: return
		
		ssa = self.cloud.ssa
		tau_star = self.cloud.tau_star
		
		#_get shape of single point
		nphot, nk = self.k.shape #_update before any rolls

		#_roll to see how far until interaction_________________________________
		self.tau = -log(1-r(nphot)).reshape(nphot,1) #_not k[2]
		
		#_update current location
		point	 		= (self.history[:,:,-1]+self.tau*self.k).reshape(nphot,nk,1)
		self.history 	= append(self.history,point,axis=2)
		
		#_BASE CHECK: if leaves cloud, drop_____________________________________
		bot_index = self.history[:,2,-1] > tau_star
		self.result[bot_index] = 'base'
		self.drop(bot_index, **kwargs)
		if not self.live: return
		
		#_REFLECTED CHECK: check if escaped top_________________________________
		top_index = self.history[:,2,-1] < top
		self.result[top_index] = 'top'
		self.drop(top_index, **kwargs)
		if not self.live: return
		
		nphot, nk = self.k.shape

		#_roll to see if absorbed_______________________________________________
		absorption_index = r(nphot) > ssa 
		self.result[absorption_index] = 'absorbed'
		self.drop(absorption_index, **kwargs)
		if not self.live: return
		
		#_drop weight by ssa and drop if below 0.01_____________________________
		if weight:
			self.weight *= ssa
			weight_index = self.weight < weightmin
			self.result[weight_index] = 'fractured'
			self.drop(weight_index, **kwargs)
			if not self.live: return
			
		nphot, nk = self.k.shape
			
		#_roll to see in what direction_________________________________________
		if force_angle == None:
			phi		= 2*pi*r(nphot)							#_randomize azimuth
			theta 	= self.cloud.henyey_greenstein(r(nphot))
		else: 												#_force angles
			phi 	= tile(force_angle[0]*d2r, nphot)		# for testing
			theta 	= tile(force_angle[1]*d2r, nphot)
			
		self.turn(theta,phi,**kwargs)						#_use new angles 
															# to update k
		self.scattered += 1									#_inc scat count

	def turn(self, theta, phi, **kwargs):
		''' augment direction of vector by theta,phi '''
		from numpy import ones, dot, arange
		if not self.live: return
		
		np, nk 	= self.k.shape
		shape 	= [np,nk,3]
		
		#_build 3x3 array describing coordinate system orthogonal to the
		# propagating ray and with x parallel to the model horizontal plane
		aug = paxes(self.k, self.phi)
		
		#_get a unit vector describing new coordinate direction
		kpp 	= kvector(theta, phi).A.reshape(np, nk, 1)
		self.k 	= dot(aug, kpp)[arange(np),:,arange(np),0]

		#_update angle arrays
		self.phi 	= phi
		self.theta 	= theta
	
	'''dump all photons to file'''
	def drop_all(self,**kwargs):[ self.drop(0,**kwargs) for n in self.result ]
		
	def drop(self,idx,outdir='.',outdat='photon.dat',fid=None,**kwargs):
		'''
		remove photon from record, write 
		TAU		RESULT	SCATTERED	HISTORY
		10.5f	10s		10i 		(10.5f,10.5f,10.5f) * nscat+1
		
		index,	bool,	location of record to remove
		'''
		from numpy import array, where, empty
		if not self.live: return
		mkdir_p(outdir)

		#_make an integer index of what to drop
		idx = where(idx==True)[0]
		idx.sort() #_!!!
		
		fopened = True if fid != None else False
		if not fopened:							
			fname = '/'.join((outdir, outdat))
			fid = open(fname,'a')				#_open file handle
		
		nphot, nscat, nk = self.history.shape	#_setup record format
		kidx = range(nphot)						#_make keep list

		if not fid.tell():						#_write header at start
			vals = [self.tau_star,self.zenith,self.azimuth]
			fmt	 = '%10.3f '*3 + '\n'
			line = fmt % tuple(vals)
			fid.write(line)

		#_write output ascii file
		for n in idx[::-1]:
			res = self.result[n]
			sca = self.scattered[n]
			wei = self.weight[n]
			loc = self.history[n].T.flatten(); val = [wei,res,sca]
			val.extend(loc); np = loc.shape[0] / 3
		
			#_make formatted print statement and write to file
			fmt = '%10.3f,%10s,%10i,' + '%10.5f,%10.5f,%10.5f,' * np + '\n'
			line = fmt % tuple(val)

			fid.write(line)
			del kidx[n]					#_remove offending index
		kidx = array(kidx)				#_convert to ndarray for use as index
			
		if not fopened: fid.close()		#_close file handle

		self.tmp_result.extend(self.result[idx])
		self.tmp_weight.extend(self.weight[idx])
		self.tmp_scattered.extend(self.scattered[idx])

		#_loop over relavent attributes and remove them
		attrs = ['k','history','weight','result','scattered','theta','phi']
		if self.nlive() == 0: 		#_kill when done
			[ self.__setattr__(n, empty(0)) for n in attrs ]
			self.live 		= False
			self.weight 	= array(self.tmp_weight)
			self.result 	= array(self.tmp_result)
			self.scattered	= array(self.tmp_scattered)
		else:						#_drop individual photo data
			[ self.__setattr__(n, self.__getattribute__(n)[kidx]) for n in attrs 
				if len(kidx) > 0 ]
			
	def traveled(self,origin=False):
		'''
		returns geometric distance traveled by photon, total or to entry
		origin,	bool,	return only distance to origin, as the non-scattering
						photon flies
		'''
		from numpy import sqrt
		dbg(origin,5)
		if not self.live: return
					
		if not origin:	#_total meanderings
			segments = self.history[:,:,1:] - self.history[:,:,:-1]
			return sqrt((segments**2).sum(axis=1)).sum(axis=1)
		else:			#_from cloud entrance
			return sqrt((self.history[:,:,-1]**2).sum(axis=1))

	def count_result(self,res): 
		'''
		count photons exiting top, bottom, or absorbed
		returned value is the sum of the weight - if weight
		is turned off, then it's just a count of photons
		
		'''
		from numpy import where
		idx = (self.result == res)
		return self.weight[idx].sum()

	def count_direct(self):
		'''returns integer value for photons directly transmitted through cloud'''
		from numpy import append
		size = self.result.size
		bool_base = (self.result == 'base').reshape(size,1)
		bool_scat = (self.scattered == 0).reshape(size,1)
		#_weight doesnt come into play for t_dir, no scattering, no fracturing
		return append(bool_base,bool_scat,axis=1).all(axis=1).sum()

	def count_absorbed(self):
		''' 
		count the aborbed photon weight 
		Any fractional splitting goes into this
		'''
		#_just count the number, the bits went into absorbed
		fracture = (self.result == 'fractured').sum()
		absorbed = (self.result == 'absorbed').sum()
		
		#_get the parts of the scattered out that should be in cloud
		top = (self.result == 'top').sum() - self.count_result('top')
		bot = (self.result == 'base').sum() - self.count_result('base')
		
		return fracture + absorbed + top + bot

	def nlive(self):
		'''count the number of photons still bouncing around'''
		return self.result.tolist().count(None)
		
	def plot(self, outdir='.', imgfmt='png', show=False, label='',
		limit_plot=[False, True], pmethod='fast', pname=None, 
		min_for_hexplot=1e3, exp_val='', **kwargs):
		'''
		plot photon path
		idx,		int,		if supplied, only the listed photon paths will plot
		ax,			plt.artist,	matplotlib plot space
		outdir,		str,		output directory
		imgfmt,		str,		image file type
		limit_plot,	bool(2),	don't extend plot area with data, (x, y)
		'''
		from numpy import arange, array, tan, cos, pi, vstack, tile
		import matplotlib.pyplot as plt
		from matplotlib import rcParams
		from matplotlib.cm import YlOrRd_r as YOR
		from numpy.ma import masked_outside, masked_where
		from libtools import shrink_ticks

		if self.live: return #_this version only plots the damned
		
		self.figure = plt.figure(figsize=[16., 4.])
		self.axes = self.figure.add_subplot(111)
		delt = self.tau_star/5.			
		self.axes.set_ylim(self.tau_star+delt, -delt)
		self.axes.set_xlim(-30, 30)
	
		#_if large enough, plot binned interaction sites
		nphoton	= len(self.history)
		if nphoton >= min_for_hexplot:
			#_stack all points and mask outside cloud
			zyx = vstack(self.history)						#_all points 
			zh	= masked_outside(zyx[:,2], 0, self.tau_star)
			xh 	= masked_where(zh.mask, zyx[:,0])
			hx 	= self.axes.hexbin = self.axes.hexbin(xh, zh, mincnt=10,
				bins='log', cmap=YOR, zorder=2, gridsize=200, vmax=4.)
			cb	= plt.colorbar(hx, orientation='vertical', aspect=8)
			cb.set_label('log10(N)')
		
		#_plot them
		if pmethod == 'slow':
			#_multiple colors, but slow as dirt
			[self.axes.plot(l.T[0], l.T[2], linewidth=.3) for l in self.history]
		elif pmethod == 'fast':
			#_all blue, but passed as a single line with segments separated by
			# NoneType objects
			jump	= [ nhist.shape[0] for nhist in self.history ]
			nt		= sum(jump)

			#_initialize array
			xyz	= tile([None,None,None], (nt+nphoton, 1))	#_spaced out
			last = 0
			for n in arange(nphoton):
				this 			= last + jump[n]
				xyz[last:this,:]= self.history[n].copy()
				last 			= this + 1
				xp				= xyz.T[0]
				zp 				= xyz.T[2]
				
			if nt+nphoton > 1e5:
				dbg(('modifying agg.path.chunksize! check for artifacts'), 3)
				rcParams['agg.path.chunksize'] = 1000
			
			#_plot photon paths
			self.axes.plot(xp, zp, linewidth=.3, zorder=1)
		
		cxmin, cxmax = self.axes.xaxis.get_view_interval()
		dxmin, dxmax = self.axes.xaxis.get_data_interval()
		cymin, cymax = self.axes.yaxis.get_view_interval()
		dymin, dymax = self.axes.yaxis.get_data_interval()
		
		#_expands x,y axes as needed by photon travel
		if not limit_plot[0]:	#_X
			xmax = max([cxmax,dxmax])
			xmin = min([cxmin,dxmin])
			self.axes.set_xlim(xmin,xmax)
		if not limit_plot[1]:	#_Y
			ymax = max([cymax,dymax])
			ymin = min([cymin,dymin])
			self.axes.set_ylim(ymax,ymin)
			
		#_get current size of axes
		cymin, cymax = self.axes.yaxis.get_view_interval()	
		cxmin, cxmax = self.axes.xaxis.get_view_interval()
			
		#_plot incident ray and cloud boundaries
		x = [0]; y = [0,cymax]
		x.append(-cymax*tan(self.zenith*d2r-pi)*cos(self.azimuth*d2r))
		args = { 'linewidth' : 0.9, 'color' : 'k' }
		self.axes.plot(x, y, linewidth=3, linestyle='--', color='k', zorder=3)		
		self.axes.plot([cxmin,cxmax], [0,0], zorder=3, **args)
		self.axes.plot([cxmin,cxmax], [self.tau_star,self.tau_star], 
			zorder=3, **args)
			
		#_gussy it up
		shrink_ticks(self.axes, 8)
		self.axes.set_ylabel(r'$\tau$')
		self.axes.set_xlabel(label)
		
		#_put it somewhere
		if show: 
			plt.show()
		else:
			mkdir_p(outdir)
			if pname == None:
				fname = '/'.join((outdir,'.'.join(('photon',imgfmt))))
			dbg((pname,outdir),l=2)
			plt.savefig(pname)
			
		#_the latter should release memory, but it's not. The first shouldn't
		self.figure.clear()
		plt.close()
		
class Cloud(object):
	def __init__(self, g=0.85, ssa=1., tau_star=10., similarity=False, **kwargs):
		'''
		g			float,	asymmetry parameter [0,1]
		ssa			float,	single scatter albedo, [0,1]
		tau_star	float,	optical depth of homogeneous layer, [0,inf)
		
		Having these attributed to photons was easier in the original code,
		but I prefer associating these with a class describing the medium
		'''
		self.g 			= g 		if not similarity else 0
		self.ssa 		= ssa		if not similarity else omega_prime(ssa, g)
		self.tau_star 	= tau_star	if not similarity else \
									tau_prime(ssa, g, tau_star)
		
	def two_stream(self):	
		'''
		Use two stream parameterization model to calculate 
		azimuthically independent values for reflectance and transmittance
	
		ssa,	float,	[0,1], single scatter albedo
		g,		float,	[-1,1],	asymmetry parameter
	
		returns:
		refl	float,	[0,1],	reflectivity
		tran	float,	[0,1],	transmittance 
		'''
		from numpy import sqrt, exp
	
		#_calculate gamma, Petty 13.25
		G = 2*sqrt(1-self.ssa)*sqrt(1-self.ssa*self.g)
		
		#_calculate refl and trans depending on ssa value
		if self.ssa == 1:
			#_Petty 13.63 and 64
			refl = (1-self.g)*self.tau_star / (1+(1-self.g)*self.tau_star)
			tran = 1. / (1+(1-self.g)*self.tau_star)			
##			if refl+tran != 1 or (refl+tran).mean() != 1: 
##				raise RuntimeError, 'photons magically disappearing ' +\
##					' '.join((str(refl), str(trans))) 

		elif self.ssa < 1:
			#_calculate semi-invinit cloud reflectivity, Petty 13.45
			r8 	= sqrt(1-self.ssa*self.g) - sqrt(1-self.ssa)
			r8 /= sqrt(1-self.ssa*self.g) + sqrt(1-self.ssa)
			
			#_Petty 13.65
			refl = r8 * (exp(G*self.tau_star) - exp(-G*self.tau_star))
			refl /= (exp(G*self.tau_star) - r8**2*exp(-G*self.tau_star))
			tran = (1.-r8**2) / (exp(G*self.tau_star) \
				- r8**2*exp(-G*self.tau_star))
				
		else:
			raise RuntimeError, 'invalid single scatter value ' + str(ssa)

		return refl, tran, 1-tran-refl #_r, t, a
		
	def henyey_greenstein(self, r):
		'''
		calculate scattering theta in accordance with Henyay Greenstein dist
		r,	flt,	random number, [0,1] INCLUSIVE
		g,	flt,	asymettry parameter.  0 => equal f/b hemisphere, 1, all forward

		I liked this better as an external function, but will cram it into
		a class for this project.

		Binzoni, T., et al. "The use of the Henyey-Greenstein phase function
		in Monte Carlo simulations in biomedical optics." Physics in medicine 
		and biology 51.17 (2006): N313.
		'''
		from numpy import arccos, isnan
		if self.g == 0:	#_give random isotropic direct
			cos_theta = 2.*r-1
		else:			#_give hemispherically preferred direction
			cos_theta = (1+self.g**2-((1-self.g**2)/\
				(1-self.g+2.*self.g*r))**2)/2./self.g
		return arccos(cos_theta)
		
#############################################################################_80
#_TOOLS_#####################################################################_80
#############################################################################_80

def beers(tau=1., theta=0, **kwargs): 
	from numpy import exp, cos, pi, abs
	mu = abs(cos(theta*pi/180.))
	return exp(-tau/mu)

def kvector_init(theta, phi):
	'''
	return unit vector pointing in direction defined by theta and phi 
	In model space, this is relative to the vertical.  In photon space,
	it's relative to the propoagation direction
	theta,phi	float,	directions about axis of propagation IN RADIANS
	
	'''
	from numpy import matrix,cos,sin
	return matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),-cos(theta)]).T

def paxes(k, phi):
	''' vector is horizontal in model space and orthogonal to propagation '''
	from numpy import matrix, cross, dot, array, sin, cos, tile, zeros
	from numpy.linalg import norm

	#_if k0**2+k1**2 == 0, apply different transform
	idx = array((k[:,0]**2 + k[:,1]**2) == 0)		#_indices where zero

	z 	= matrix([0,0,1]).T							#_vertical model space
	p 	= cross(z.T,k)								#_ortho to vert and prop
	mag = matrix([ norm(x) for x in p ])			#_magnitude 
	#_maybe put in a check for zero magnitude vectors (perfectly horizontal)

	x = (p / mag.T)		 							#_get unit vector
	if idx.sum() != 0: 								#_replace nans
		dbg('direction parallel to vertical', 3)
		xp = array(-sin(phi[idx]))
		yp = array( cos(phi[idx]))
		zp = tile(0, len(idx))
		x[idx] = array([xp,yp,zp]).T

	y = matrix(cross(k,x))							#_calc y unit vector

	nphot 		= k.shape[0]
	aug 		= zeros((nphot,3,3))
	aug[:,:,0] 	= x
	aug[:,:,1] 	= y 
	aug[:,:,2] 	= k # tile(z.T,(nphot,1))
	
	return aug
			
def kvector(theta, phi):
	'''
	return unit vector in direction of theta/phi
	Initially, this will be in relationship to standard 3d cartesian coords,
	but as the propagation path shifts, it will be in relation to the last path
	
	theta,phi	float,	directions about axis of propagation IN RADIANS
	'''
	from numpy import matrix,cos,sin
	return matrix([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)]).T
	
def plot_hg(polar=True):
	'''plot what's coming out of the HG function '''
	raise RuntimeWarning, 'NOT ACTUALLY H-G'
	import matplotlib.pyplot as plt
	from numpy import cos
	theta = linspace(0,2*pi,1000)
	theta = linspace(-pi,pi,1000)
	
	for g in linspace(0,1,3):
		pT 		= (1-g*g)/(1+g*g-2*g*cos(theta)**(3/2))/4/pi
		label 	= '='.join(('g',str(g)))
		ax 		= plt.subplot(111, polar=polar)
		
		ax.plot(theta, pT, label=label, linewidth=0.4)	
		
	if not polar: 
##		ax.set_xlim(-pi,pi)
		plt.legend()
		
	plt.title('Henyey-Greenstein Phase Function')
	plt.show()

#_similarity theory
def omega_prime(ssa, g): return float(1.-g)*ssa/(1.-g*ssa)
def tau_prime(ssa, g, tau_star): return (1.-g*ssa)*tau_star
	
def setup_groups(objects,nproc=1,**kwargs):
	'''
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
	import errno
	dbg(path,5)
	def mkdir(p):
		try: 
			os.makedirs(p)
			os.chmod(p,0775)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST: pass
			else: raise

	if hasattr( '__iter__', path ): [ mkdir(p) for p in path ]
	else: mkdir(path)

def photon_test():
	''' sanity check against Petty's matrix worksheet '''
	exp = [.387,.461,-.799,.787,.388,-.479]
	
	p = Photon(size=1, zenith=37, azimuth=50, cloud=Cloud())
	got = p.k.flatten().tolist()
	
	p.advance(force_angle=[135,30], top=-10)
	got.extend(p.k.flatten().tolist())
	
	print '%5.3f '*6 % tuple(exp)
	print '%5.3f '*6 % tuple(got)
	
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
		file	 = file.split('/')[-1].split('.')[0]
		scream 	= '[%s.%s.%i] %s' % (file,method,line,msg)
		
		if not err: print scream
		else: raise RuntimeError, scream

def to_string(a):
	if hasattr(a,'__iter__'): return [ str(b) for b in a ]
	else: return str(a)

def shrink_ticks(ax, size=6):
     ''' reduce font size of ticks '''   
     [ t.label.set_fontsize(size) for t in ax.yaxis.get_major_ticks() ]
     [ t.label.set_fontsize(size) for t in ax.xaxis.get_major_ticks() ]

def tighten_plot(axes, rounded=False, onlyx=False, onlyy=False):
	'''
	make plot area as small as possible
	axes,			matplotlib.pyplot.axes
	rounded,		bool,	round axes up
	'''
	xmin, xmax = axes.xaxis.get_data_interval()
	ymin, ymax = axes.yaxis.get_data_interval()
	if onlyx	: axes.set_xlim(xmin,xmax)
	elif onlyy	: axes.set_ylim(ymin,ymax)
	else		: axes.axis([xmin,xmax,ymin,ymax])
	if rounded: raise ValueError, 'rounded not implemented'

if __name__ == '__main__': 
	#_change the namelist experiment names if running similarity theory tests
	if args['similarity']:			#_change output dir name if similarity
		import re
		newargs = {}				#_new kwargs
									#_could have just put list with 
									# ['exp','sim'] to key loop over.
		for key, values in args['namelist'].iteritems():
			key = re.sub('^exp', 'sim', key)
			newargs.update({key : values})
		args['namelist'] = newargs
	run_main(**args)
