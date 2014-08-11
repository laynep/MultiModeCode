"""
Set pyplot parameters for nice plots
"""

import matplotlib as mpl
from numpy import sqrt

#fig_width_pt = 235.0  # Get this from LaTeX using \showthe\columnwidth
fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
#fig_width_pt = 300.0  # Get this from LaTeX using \showthe\columnwidth

#JCAP
#fig_width_pt = 440.0  # Get this from LaTeX using \showthe\columnwidth

inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

# Explicitly set fontsizes:
tick_size = 8
big_font=True
if big_font:
    font_size = 13
else:
    font_size = 10


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

figprops = dict(figsize=(1.0*fig_width, 1.00*fig_height))

#adjustprops = dict(left=0.13, bottom=0.12, right=0.97, top=0.96,
#    wspace=0.1, hspace=0.02)
adjustprops = dict(left=0.19, bottom=0.185, right=0.96, top=0.96,
    wspace=0.1, hspace=0.02)
