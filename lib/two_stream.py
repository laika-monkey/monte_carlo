def two_stream(ssa=1., g=0., tau_star=10., **kwargs):	
	'''
	Use two stream parameterization model to calculate 
	azimuthically independent values for reflectance and transmittance
	
	ssa,	float,	[0,1], single scatter albedo
	g,	float,	[-1,1],	asymmetry parameter
	
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
