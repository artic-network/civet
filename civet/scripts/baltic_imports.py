import os
from matplotlib import font_manager as fm, rcParams
from baltic import *

from IPython.display import HTML
import re
import copy

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import numpy as np
from scipy.special import binom

import itertools
import requests

try:
    from StringIO import StringIO as sio
    from cStringIO import StringIO as csio
except ImportError:
    from io import StringIO as sio
    from io import BytesIO as csio
from Bio import Phylo
from collections import defaultdict
import matplotlib.font_manager as font_manager

#%matplotlib inline
font_dirs = ['utils/helveticaneue', ]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
font_list = font_manager.createFontList(font_files)
font_manager.fontManager.ttflist.extend(font_list)

mpl.rcParams['font.family'] = 'helveticaneue'
mpl.rcParams['font.weight']=300
mpl.rcParams['axes.labelweight']=300
mpl.rcParams['font.size']=25
