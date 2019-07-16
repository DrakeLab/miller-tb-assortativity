###### SIRS-var-sus.py simulates SIR model with variable susceptibility on assorted networks ###### 
###### This version builds on the fast_nonMarkov_SIR model 

import networkx as nx
import EoN
import random
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
from sklearn.model_selection import ParameterGrid

