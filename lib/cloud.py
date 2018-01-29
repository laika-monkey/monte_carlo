class Cloud(object):
	def __init__(self, tau_star=10., width=40, figure=None, zenith=180, 
		azimuth=0, g=0, ssa=1., clouds=None, **kwargs):
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
		
		self.tau_star = tau_star
		self.zenith = zenith
		self.azimuth = azimuth
##		self.ssa = ssa
##		self.g = g

		#_initialize figure and axes
		self.figure = figure if figure != None \
			else plt.figure(figsize=[16., 4.])
		self.figure.add_subplot(111)
		
		#_initialize cloud top and base lines
		xlim = [-width/2., width/2.]
		ylim = [-2, tau_star+2]
		self.figure.axes[0].set_ylim(ylim[::-1])
		self.figure.axes[0].set_xlim(xlim)
		
		arg = {'linewidth':.5,'color':'k'}
		a, = self.figure.axes[0].plot(xlim, [tau_star,tau_star], **arg)
		b, = self.figure.axes[0].plot(xlim, [0,0], **arg)

		#_plot incident ray path
		x = [0, self.incident_x(ylim[0])]
		y = [0, ylim[0]]
		c,= self.figure.axes[0].plot(x, y, linewidth=3, color='k',
			linestyle='--')

		#_save lines
		self.base = a
		self.top = b
		self.ray = c
		self.figure.canvas.draw()	

	def incident_x(self,o):
		''' calculate the x component of the incident ray '''
		from numpy import cos, tan, pi
		d2r = pi/180
		return -o*tan(self.zenith*d2r-pi) * cos(self.azimuth*d2r)

