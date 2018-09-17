def configure():
	import matplotlib.pyplot as plt
	from cycler import cycler

	plt.style.use('bmh')
	font_size = 14
	# Font
	plt.rcParams['text.usetex'] = True
	plt.rcParams['text.latex.unicode'] = True
	plt.rcParams['font.size'] = font_size
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['font.sans-serif'] = "Computer Modern"

	# X-ticks
	plt.rcParams['xtick.labelsize'] = font_size
	plt.rcParams['xtick.direction'] = 'in'
	plt.rcParams['xtick.minor.visible'] = True
	plt.rcParams['xtick.major.size'] = 6
	plt.rcParams['xtick.minor.size'] = 4
	plt.rcParams['xtick.major.width'] = 2
	plt.rcParams['xtick.minor.width'] = 1

	# Y-ticks
	plt.rcParams['ytick.labelsize'] = font_size
	plt.rcParams['ytick.direction'] = 'in'
	plt.rcParams['ytick.right'] = True
	plt.rcParams['ytick.minor.visible'] = True
	plt.rcParams['ytick.major.size'] = 6
	plt.rcParams['ytick.minor.size'] = 4
	plt.rcParams['ytick.major.width'] = 2
	plt.rcParams['ytick.minor.width'] = 1

	# Axes
	plt.rcParams['axes.labelsize'] = font_size
	plt.rcParams['axes.labelweight'] = 'normal'
	plt.rcParams['axes.linewidth'] = 2
	plt.rcParams['axes.edgecolor'] = 'black'
	plt.rcParams['axes.facecolor'] = 'white'
	plt.rcParams['axes.grid'] = False

	# Color cycler
	plt.rcParams['axes.prop_cycle'] = cycler('color', ['black', '#b2182b', '#ef8a62', '#fddbc7', '#d1e5f0', '#67a9cf', '#2166ac'])

	# Legend
	plt.rcParams['legend.facecolor'] = "white"
	plt.rcParams['legend.frameon'] = False