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

def run_main(drawlive=False, limit_plot=True, **namelist): 
	import os
	import matplotlib.pyplot as plt
	from lib import photon_propagation
	from lib import mkdir_p
	from lib import Cloud as cloud
	from time import sleep
	
	#_loop over runs specified in the Namelist dictionary
	for run in namelist:
	
		#_turn on interactive window
		if drawlive:
			plt.ion()
		
		#_initialize output file
		outdir = 'output_{}'.format(run)
		mkdir_p(outdir)
		fname = os.path.join(outdir, 'photon.dat')
		print 'Writing photon data to {}'.format(fname)
		
		#_pull out keyword args and update	
		kwargs = namelist[run]
		clouds = cloud(**kwargs)	#_initialize cloud
		kwargs.update({	'outdir'	: outdir,
				'figure' 	: clouds.figure,
				'fid'		: open(fname, 'w'),
				'axes'		: clouds.figure.axes[0],
				'limit_plot'	: limit_plot,
				'clouds'	: clouds,
				'drawlive'	: drawlive	})
		
		#_calc how many loops need to be made
		photon_propagation(**kwargs)

		#_same image 
		pname = os.path.join(outdir, 'photon.png')
		plt.savefig(pname)
		print 'Saving path image to {}'.format(pname)	

		sleep(10)
	
if __name__ == '__main__':
	import yaml
	namelist = yaml.load(open('options.yaml', 'r'))
	run_main(**namelist)
