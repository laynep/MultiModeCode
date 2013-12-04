"""
Set pyplot parameters for nice plots
"""

import matplotlib as mpl
from scipy import sqrt

fig_width_pt = 235.0  # Get this from LaTeX using \showthe\columnwidth
#fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
#fig_width_pt = 300.0  # Get this from LaTeX using \showthe\columnwidth

inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

# Explicitly set fontsizes:
hiranya=True
if hiranya:
    font_size = 13
    tick_size = 8
else:
    font_size = 10
    tick_size = 8

params = {
          #'backend': 'ps',
          'axes.labelsize': font_size,
          'xtick.labelsize': tick_size,
          'ytick.labelsize': tick_size,
          'text.fontsize': font_size,
          'legend.fontsize': font_size,
          'text.usetex': True,
          'figure.figsize': fig_size,
          }
mpl.rcParams.update(params)
