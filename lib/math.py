def kvector_init(theta, phi):
	'''
	return unit vector pointing in direction defined by theta and phi 
	In model space, this is relative to the vertical.  In photon space,
	it's relative to the propoagation direction
	theta,phi	float,	directions about axis of propagation IN RADIANS
	
	'''
	from numpy import matrix,cos,sin
	return matrix([cos(phi)*sin(theta), sin(phi)*sin(theta), -cos(theta)]).T
		
def kvector(theta, phi):
	'''
	return unit vector in direction of theta/phi
	Initially, this will be in relationship to standard 3d cartesian coords,
	but as the propagation path shifts, it will be in relation to the last
	path

	theta,phi	float,	directions about axis of propagation IN RADIANS
	'''
	from numpy import matrix, cos, sin
	return matrix([cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]).T

def paxes(k, phi):
	'''
	vector is horizontal in model space and orthogonal to propagation
	returns array describing coordinate system to be augmented
	'''
	from numpy import matrix, cross, dot, array, sin, cos, tile
	from numpy.linalg import norm

	z = matrix([0,0, 1]).T 			#_vertical model space
	p = cross(z.T, k)			#_ortho to vert and prop
	mag = matrix([ norm(x) for x in p ]) 	#_magnitude 

	x = (p / mag.T)				#_get unit vector
	if (k[:,0]**2 + k[:,1]**2) == 0: 	#_replace nans
		x = matrix([-sin(phi), cos(phi), 0.])

	y = matrix(cross(k, x))			#_calc y unit vector
	 
	return matrix([x.A.flatten(), y.A.flatten(), k.flatten()]).T
