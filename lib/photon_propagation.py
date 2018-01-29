def photon_propagation(total=100, group=10, drawlive=False, **kwargs):
	'''
	Loop over groups of photons advancing until total is met
	total		int,	How many photons to generate
	group		int,	How many photons to have alive at a time
	drawlive	bool,	Make images on the fly
	'''
	from numpy import array, append
	from photon import Photon as photon
	
	photons = array([]) #_initialize list to carry photon objects
	figure	= kwargs.get('figure') #_get a copy to speed up drawing
	clouds	= kwargs.get('clouds') #_get a copy of clouds
	
	while photon.tot_fin < total:

		#_keep the total number of photons in list.photons == group
		currently_generated = total - photon.tot_fin - len(photons)
		currently_active = group - len(photons)
		need = min([currently_generated, currently_active])

		#_generate enough photos to refill group
		photons = append(photons, [photon(**kwargs) for n \
			in range(int(need)) ])

		#_move all photons forward one step
		[ a.advance(**kwargs) for a in photons ]

		#_update number of complete
		live = array([ p.live for p in photons ]) #_find live
		photons = photons[live]	#_reduce array to living

		if drawlive:
			args = (photon.tot_fin, photon.tot_gen)
			title = '{} complete of {} generated'.format(*args) 
			clouds.figure.axes[0].set_title(title)
			clouds.figure.canvas.draw()
