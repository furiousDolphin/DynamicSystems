


import sys
import os

build_path = os.path.join(os.getcwd(), 'build')
if build_path not in sys.path:
    sys.path.append(build_path)

import numpy as np
import matplotlib.pyplot as plt
import module_uno as m  
from scripts import Plots
from typing import Callable


a_coeffs = np.array([4.0, 4.2, 1.2, 1.0]) 
b_coeffs = np.array([10.0])

system: m.System = m.System(b_coeffs, a_coeffs)
f: Callable[[float], float] = lambda t: 1.0 if np.sin(t) > 0 else -1.0  
system.set_forcing_func(f)

n: int = 1000
t_dense: np.ndarray = np.zeros(n)
y_dense: np.ndarray = np.zeros(n)

dt:float = 30/n
for i in range(n):
    (t, y) = system.do_RK4_step(dt)
    t_dense[i]=t
    y_dense[i]=y

plot = Plots.Plot((500, 500))
plot.add(t_dense, y_dense)
plot.show()