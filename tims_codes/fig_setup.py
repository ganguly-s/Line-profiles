from pylab import *
import importlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

def fix_cbar_ticks(cax):
	indx = 0
	for t in cax.get_yticklabels():
		indx += 1
		if indx%2 == 0:
			t.set_visible(False)
		else:
			t.set_fontsize(12)

# 2-panel plot -----------------------------------------------------------------------
def vel_space_plot1D(fig,line_data):

	""" Figure setup """	

	plt.subplots_adjust(left=0.12, bottom=0.08, right=0.88, top=0.95, wspace=0.1, hspace=0.)
	plt.rcParams['font.family'] = 'serif' #'STIXGeneral' #'serif'
	matplotlib.rcParams['font.size'] = '16'
	matplotlib.rcParams['ps.fonttype'] = 42 #note: fontype 42 compatible with MNRAS style file when saving figs
	matplotlib.rcParams['mathtext.fontset'] = 'stix'
	plt.rcParams['axes.labelsize'] = 16
	plt.rcParams['axes.linewidth'] = 1.
	plt.rcParams['xtick.major.size'] = 6
	plt.rcParams['xtick.minor.size'] = 3
	plt.rcParams['ytick.major.size'] = 6
	plt.rcParams['ytick.minor.size'] = 3
	plt.rcParams['xtick.labelsize'] = 14
	plt.rcParams['ytick.labelsize'] = 14
	plt.rcParams['legend.numpoints'] = 1  # uses 1 symbol instead of 2
	plt.rcParams['legend.frameon'] = False 
	plt.rcParams['legend.handletextpad'] = 0.3

	# plt.rcParams['text.usetex'] = True

	spec3 = matplotlib.gridspec.GridSpec(ncols=4, nrows=8)
	ax1 = fig.add_subplot(spec3[0:4,:])
	ax2 = fig.add_subplot(spec3[4:,:],sharex=ax1)


	# fix ticks
	for ax in [ax1,ax2]:
		ax.minorticks_on()
		ax.xaxis.set_ticks_position('both')
		ax.yaxis.set_ticks_position('both')
		ax.tick_params(axis='both', direction='in',width=1.) 
		ax.tick_params(which='minor',axis='both',direction='in',width=1.) 

	# make top axes labels invisible
	plt.setp(ax1.get_xticklabels(), visible=False)


	# labels
	ax2.set_xlabel(r'$v$' + ' [km/s]',size='large')

	#store axes in a dictionary
	axes = {'1':ax1,'2':ax2}
	
	# display ion
# 	ion_template = r'$\rm{O}$' + ' ' r'$\rm{VIII}$' + ' ' + r'$\rm{Ly\alpha}$' 
	# ion_template = r'$\rm{C}$' + ' ' r'$\rm{IV}$' + ' ' + r'$\rm{Ly\alpha}$' 
	ion_template = line_data['ion']
	ion_text = ax1.text(0.05, 0.9, ion_template, transform=ax1.transAxes, size='large')


	return fig,axes
	
# 2-panel plot -----------------------------------------------------------------------
def vel_space_plot(fig,simdata,line_data,SHOW_IMAGE=0,Y_OBSERVER=0):

	""" Figure setup """	

	plt.subplots_adjust(left=0.12, bottom=0.08, right=0.88, top=0.95, wspace=0.1, hspace=0.)
	plt.rcParams['font.family'] = 'serif' #'STIXGeneral' #'serif'
	matplotlib.rcParams['font.size'] = '16'
	matplotlib.rcParams['ps.fonttype'] = 42 #note: fontype 42 compatible with MNRAS style file when saving figs
	matplotlib.rcParams['mathtext.fontset'] = 'stix'
	plt.rcParams['axes.labelsize'] = 16
	plt.rcParams['axes.linewidth'] = 1.
	plt.rcParams['xtick.major.size'] = 6
	plt.rcParams['xtick.minor.size'] = 3
	plt.rcParams['ytick.major.size'] = 6
	plt.rcParams['ytick.minor.size'] = 3
	plt.rcParams['xtick.labelsize'] = 14
	plt.rcParams['ytick.labelsize'] = 14
	plt.rcParams['legend.numpoints'] = 1  # uses 1 symbol instead of 2
	plt.rcParams['legend.frameon'] = False 
	plt.rcParams['legend.handletextpad'] = 0.3

	# plt.rcParams['text.usetex'] = True

	spec3 = matplotlib.gridspec.GridSpec(ncols=4, nrows=8)
	ax1 = fig.add_subplot(spec3[0:4,:])
	ax2 = fig.add_subplot(spec3[4:,:],sharex=ax1)


	# fix ticks
	for ax in [ax1,ax2]:
		ax.minorticks_on()
		ax.xaxis.set_ticks_position('both')
		ax.yaxis.set_ticks_position('both')
		ax.tick_params(axis='both', direction='in',width=1.) 
		ax.tick_params(which='minor',axis='both',direction='in',width=1.) 

	# make top axes labels invisible
	plt.setp(ax1.get_xticklabels(), visible=False)


	# labels
	ax2.set_xlabel(r'$v$' + ' [km/s]',size='large')

	#store axes in a dictionary
	axes = {'1':ax1,'2':ax2}
	
	# display ion
# 	ion_template = r'$\rm{O}$' + ' ' r'$\rm{VIII}$' + ' ' + r'$\rm{Ly\alpha}$' 
	# ion_template = r'$\rm{C}$' + ' ' r'$\rm{IV}$' + ' ' + r'$\rm{Ly\alpha}$' 
	ion_template = line_data['ion']
	ion_text = ax1.text(0.05, 0.9, ion_template, transform=ax1.transAxes, size='large')


	if SHOW_IMAGE: # plot image of density field
		if Y_OBSERVER:
			imgax1 = fig.add_axes([0.71, 0.75, 0.2, 0.15], frameon=True) # right
		else:
			imgax1 = fig.add_axes([0.7, 0.75, 0.2, 0.15], frameon=True)
			
		if 1: # disable density image coordinates 
			plt.setp(imgax1.get_yticklabels(), visible=False)
			plt.setp(imgax1.get_xticklabels(), visible=False)
		else:
			xticks1 = imgax1.get_xticklabels()
			yticks1 = imgax1.get_yticklabels()
			for t in xticks1[1:-1]:
				t.set_visible(False)
			for t in yticks1[1:-1]:
				t.set_visible(False)

		imgax1.set_xlabel(r'$z$',size='large')
		imgax1.set_ylabel(r'$x$',size='large')


		divider = make_axes_locatable(imgax1)

		if Y_OBSERVER:
			img_box = [0.,1.,0,1.]
			#img_box = [0.,1.,0,0.5]
			cax = divider.append_axes("right", size="15%", pad=0.05)
			im = imgax1.imshow(simdata['density'].transpose(), cmap='Spectral',extent=img_box)
		else:
			img_box = [0.,1.,0,1.]
			cax = divider.append_axes("right", size="10%", pad=0.05)
			im = imgax1.imshow(simdata['density'], cmap='Spectral',extent=img_box)

		axes['img'] = imgax1
		axes['img_box'] = img_box
		axes['cax'] = cax


	return fig,axes


# 4-panel plot -----------------------------------------------------------------------
def pos_space_plot(fig,simdata,line_data,SHOW_IMAGE=0,Y_OBSERVER=0):
	""" Figure setup """	

	plt.subplots_adjust(left=0.16, bottom=0.08, right=0.84, top=0.96, wspace=0., hspace=0.)
	plt.rcParams['font.family'] = 'serif' #'STIXGeneral' #'serif'
	matplotlib.rcParams['font.size'] = '14'
	matplotlib.rcParams['ps.fonttype'] = 42 #note: fontype 42 compatible with MNRAS style file when saving figs
	matplotlib.rcParams['mathtext.fontset'] = 'stix'
	plt.rcParams['axes.labelsize'] = 14
	plt.rcParams['axes.linewidth'] = 1.
	plt.rcParams['xtick.major.size'] = 6
	plt.rcParams['xtick.minor.size'] = 3
	plt.rcParams['ytick.major.size'] = 6
	plt.rcParams['ytick.minor.size'] = 3
	plt.rcParams['xtick.labelsize'] = 14
	plt.rcParams['ytick.labelsize'] = 14
	plt.rcParams['legend.numpoints'] = 1  # uses 1 symbol instead of 2
	plt.rcParams['legend.frameon'] = False 
	plt.rcParams['legend.handletextpad'] = 0.3

	# plt.rcParams['text.usetex'] = True

	spec3 = matplotlib.gridspec.GridSpec(ncols=4, nrows=8)
	ax1 = fig.add_subplot(spec3[0:2,:])
	ax2 = fig.add_subplot(spec3[2:4,:],sharex=ax1)
	ax3 = fig.add_subplot(spec3[4:6,:],sharex=ax1)
	ax4 = fig.add_subplot(spec3[6:8,:],sharex=ax1)

	# fix ticks
	for ax in [ax1,ax2,ax3,ax4]:
		ax.minorticks_on()
		ax.xaxis.set_ticks_position('both')
		ax.yaxis.set_ticks_position('both')
		ax.tick_params(axis='both', direction='in',width=1.) 
		ax.tick_params(which='minor',axis='both',direction='in',width=1.)

	# make top axes labels invisible
	plt.setp(ax1.get_xticklabels(), visible=False)
	plt.setp(ax2.get_xticklabels(), visible=False)



	# twin axis showing Temperature in red
	ax2b = ax2.twinx(); #add right axis; make it red
	for tl in ax2b.get_yticklabels():
		tl.set_color('r')

	ax2b.spines['right'].set_color('red')
	ax2b.tick_params(axis='y', colors='red')


	# labels
	ax2b.set_ylabel(r'$T(x=0.5,z)$',color='r',size='large')
	ax4.set_xlabel(r'$r$',size='large')
	
	
	ax1.set_ylabel(r'$n(x=0.5,z)$',size='large')
	ax2.set_ylabel(r'$\xi(x=0.5,z)$',size='large')
	ax3.set_ylabel(r'$\alpha_{\nu_0}(x=0.5,z)$',size='large')
	ax4.set_ylabel(r'$v$' + ' [km/s]',size='large')

	#store axes in a dictionary
	axes = {'1':ax1,'2':ax2,'2b':ax2b,'3':ax3,'4':ax4}
	
	# display ion
	# ion_template = r'$\rm{O}$' + ' ' r'$\rm{VIII}$' + ' ' + r'$\rm{Ly\alpha}$' 
	# ion_text = ax4.text(0.05, 0.85, ion_template, transform=ax4.transAxes, size='large')

	if SHOW_IMAGE: # plot image of density field 
		if Y_OBSERVER:
			imgax1 = fig.add_axes([0.71, 0.8, 0.2, 0.15], frameon=True) # right
		else:
			imgax1 = fig.add_axes([0.2, 0.82, 0.2, 0.15], frameon=True)
			
		if 1: # disable density image coordinates 
			plt.setp(imgax1.get_yticklabels(), visible=False)
			plt.setp(imgax1.get_xticklabels(), visible=False)
		else:
			xticks1 = imgax1.get_xticklabels()
			yticks1 = imgax1.get_yticklabels()
			for t in xticks1[1:-1]:
				t.set_visible(False)
			for t in yticks1[1:-1]:
				t.set_visible(False)

		imgax1.set_xlabel(r'$z$',size='large')
		imgax1.set_ylabel(r'$x$',size='large')

		divider = make_axes_locatable(imgax1)

		# open or calculate line_sigmas
		# output_folder = self.out_folder + self.run_name + '/' + self.ion
		# pickled_sigmas = output_folder + '/sigmas/'
		# pickled_sigmas += 'sigmas_' + str(self.dump) + '.p'
		# line_sigmas = pickle.load(open(pickled_sigmas,"rb"))

		if Y_OBSERVER:
			img_box = [0.,1.,0,1.]
			#img_box = [0.,1.,0,0.5]
			cax = divider.append_axes("right", size="15%", pad=0.05)
			im = imgax1.imshow(simdata['density'].transpose(), cmap='Spectral',extent=img_box)
		else:
			img_box = [0.,1.,0,1.]
			cax = divider.append_axes("right", size="10%", pad=0.05)
			im = imgax1.imshow(simdata['density'], cmap='Spectral',extent=img_box)

		axes['img'] = imgax1
		axes['img_box'] = img_box
		axes['cax'] = cax


	return fig,axes

# 4-panel plot -----------------------------------------------------------------------	
def breakdown_plot(fig,simdata,line_data,SHOW_IMAGE=0,Y_OBSERVER=0):
	""" Figure setup """	

	plt.subplots_adjust(left=0.16, bottom=0.08, right=0.84, top=0.96, wspace=0., hspace=0.)
	plt.rcParams['font.family'] = 'serif' #'STIXGeneral' #'serif'
	matplotlib.rcParams['font.size'] = '14'
	matplotlib.rcParams['ps.fonttype'] = 42 #note: fontype 42 compatible with MNRAS style file when saving figs
	matplotlib.rcParams['mathtext.fontset'] = 'stix'
	plt.rcParams['axes.labelsize'] = 14
	plt.rcParams['axes.linewidth'] = 1.
	plt.rcParams['xtick.major.size'] = 6
	plt.rcParams['xtick.minor.size'] = 3
	plt.rcParams['ytick.major.size'] = 6
	plt.rcParams['ytick.minor.size'] = 3
	plt.rcParams['xtick.labelsize'] = 14
	plt.rcParams['ytick.labelsize'] = 14
	plt.rcParams['legend.numpoints'] = 1  # uses 1 symbol instead of 2
	plt.rcParams['legend.frameon'] = False 
	plt.rcParams['legend.handletextpad'] = 0.3

	# plt.rcParams['text.usetex'] = True

	spec3 = matplotlib.gridspec.GridSpec(ncols=4, nrows=9)
	ax1 = fig.add_subplot(spec3[0:2,:])
	ax2 = fig.add_subplot(spec3[2:4,:],sharex=ax1)
	ax3 = fig.add_subplot(spec3[4:6,:],sharex=ax1)
	ax4 = fig.add_subplot(spec3[7:,:])

	# fix ticks
	for ax in [ax1,ax2,ax3,ax4]:
		ax.minorticks_on()
		ax.xaxis.set_ticks_position('both')
		ax.yaxis.set_ticks_position('both')
		ax.tick_params(axis='both', direction='in',width=1.) 
		ax.tick_params(which='minor',axis='both',direction='in',width=1.)

	# make top axes labels invisible
	plt.setp(ax1.get_xticklabels(), visible=False)
	plt.setp(ax2.get_xticklabels(), visible=False)



	# twin axis showing Temperature in red
	ax2b = ax2.twinx(); #add right axis; make it red
	for tl in ax2b.get_yticklabels():
		tl.set_color('r')

	ax2b.spines['right'].set_color('red')
	ax2b.tick_params(axis='y', colors='red')


	# labels
	ax2b.set_ylabel(r'$T(x=0.5,z)$',color='r',size='large')
	ax3.set_xlabel(r'$z$',size='large')
	ax4.set_xlabel(r'$x$',size='large')
	
	
	ax1.set_ylabel(r'$n(x=0.5,z)$',size='large')
	ax2.set_ylabel(r'$\xi(x=0.5,z)$',size='large')
	ax3.set_ylabel(r'$\alpha_{\nu_0}(x=0.5,z)$',size='large')
	ax4.set_ylabel(r'$\tau_\nu(x,z=1)$',size='large')

	#store axes in a dictionary
	axes = {'1':ax1,'2':ax2,'2b':ax2b,'3':ax3,'4':ax4}
	
	# display ion
	# ion_template = r'$\rm{O}$' + ' ' r'$\rm{VIII}$' + ' ' + r'$\rm{Ly\alpha}$' 
	# ion_text = ax4.text(0.05, 0.85, ion_template, transform=ax4.transAxes, size='large')

	if SHOW_IMAGE: # plot image of density field 
		if Y_OBSERVER:
			imgax1 = fig.add_axes([0.71, 0.8, 0.2, 0.15], frameon=True) # right
		else:
			imgax1 = fig.add_axes([0.2, 0.82, 0.2, 0.15], frameon=True)
			
		if 1: # disable density image coordinates 
			plt.setp(imgax1.get_yticklabels(), visible=False)
			plt.setp(imgax1.get_xticklabels(), visible=False)
		else:
			xticks1 = imgax1.get_xticklabels()
			yticks1 = imgax1.get_yticklabels()
			for t in xticks1[1:-1]:
				t.set_visible(False)
			for t in yticks1[1:-1]:
				t.set_visible(False)

		imgax1.set_xlabel(r'$z$',size='large')
		imgax1.set_ylabel(r'$x$',size='large')

		divider = make_axes_locatable(imgax1)

		# open or calculate line_sigmas
		# output_folder = self.out_folder + self.run_name + '/' + self.ion
		# pickled_sigmas = output_folder + '/sigmas/'
		# pickled_sigmas += 'sigmas_' + str(self.dump) + '.p'
		# line_sigmas = pickle.load(open(pickled_sigmas,"rb"))

		if Y_OBSERVER:
			img_box = [0.,1.,0,1.]
			#img_box = [0.,1.,0,0.5]
			cax = divider.append_axes("right", size="15%", pad=0.05)
			im = imgax1.imshow(simdata['density'].transpose(), cmap='Spectral',extent=img_box)
		else:
			img_box = [0.,1.,0,1.]
			cax = divider.append_axes("right", size="10%", pad=0.05)
			im = imgax1.imshow(simdata['density'], cmap='Spectral',extent=img_box)

		axes['img'] = imgax1
		axes['img_box'] = img_box
		axes['cax'] = cax


	return fig,axes


