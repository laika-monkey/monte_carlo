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

namelist = {
  'single'	: {
		'tau_star'	: 25.,		#_optical depth of cloud
		'ssa'		: 1.0,		#_single scatter albedo
		'zenith'	: 135,		#_incoming rad dir from vertical
		'azimuth'	: 0,		#_azimuthal angle
		'g'		: .85,		#_-1 to 1, 
						# desc hemi scrat dir preference
		'imgfmt'	: 'png',	#_image format, png/pdf
		'total'		: 10,		#_total number of photons
		'group'		: 2,		#_number of photons to pass
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
	
	'drawlive' : True,	#_draw photons moving as it happens
				# not recommended unless doing density
				# plots or you have an intense amount
				# of memory.
				#_mostly a sanity check
	'limit_plot' : False,	#_if true, plot area will not grow
				# with photon path
	}

#############################################################################_80
#_MAIN_######################################################################_80
#############################################################################_80

def run_main(drawlive=False, limit_plot=True, **namelist): 
	import os
	import matplotlib.pyplot as plt
	from lib import photon_propagation
	from lib import mkdir_p
	from lib import Cloud as cloud
	
	#_loop over runs specified in the Namelist dictionary
	for run in namelist:
	
		#_turn on interactive window
		if drawlive:
			plt.ion()
		
		#_initialize output file
		mkdir_p(run)
		fname = os.path.join(run, 'photon.dat')
		
		#_pull out keyword args and update	
		kwargs = namelist[run]
		clouds = cloud(**kwargs)		#_initialize cloud
		kwargs.update({	'outdir'	: run,
				'figure' 	: clouds.figure,
				'fid'		: open(fname, 'w'),
				'axes'		: clouds.figure.axes[0],
				'limit_plot'	: limit_plot,
				'clouds'	: clouds,
				'drawlive'	: drawlive	})
		
		#_calc how many loops need to be made
		photon_propagation(**kwargs)

	#_same image 	
	plt.savefig(os.path.join(run, 'photon.png'))
	
if __name__ == '__main__':
	run_main(**namelist)
