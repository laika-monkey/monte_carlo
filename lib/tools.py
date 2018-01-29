debug = 1

def count_result(phot, res):
	''' count photons exiting top, bottom, or absorbed '''
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

def dbg(msg, l=1, err=False):
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
