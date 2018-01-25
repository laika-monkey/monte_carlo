def henyey_greenstein(r, g=0, **kwargs):
	'''
	calculate scattering theta in accordance with Henyay Greenstein dist
	r,	flt,	random number, [0,1] INCLUSIVE
	g,	flt,	asymettry parameter.  0 => equal f/b hemisphere, 1,
			all forward

	Binzoni, T., et al. "The use of the Henyey-Greenstein phase function
	in Monte Carlo simulations in biomedical optics." Physics in medicine 
	and biology 51.17 (2006): N313.
	'''
	from numpy import arccos, isnan

	#_give random isotropic direct
	if g == 0:	
		cos_theta = 2.*r-1 	#_this leaves out {pi,0},
					# 2.*(1-r)-1 leaves out pi/2
	#_give hemispherically preferred direction
	else:	
		cos_theta = (1+g*g-((1-g*g)/(1-g+2.*g*r))**2)/2./g

	return arccos(cos_theta)
