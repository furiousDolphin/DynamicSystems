


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

#-------------------------------------------------

from dataclasses import dataclass, field

FloatGetter = Callable[[], float]
FloatSetter = Callable[[float], None]

# @dataclass
# class ValueManager():
#     val: float
#     getter: FloatGetter = field(init=False)
#     setter: FloatSetter = field(init=False)

#     def __post_init__(self):
#         def _setter_func(new_val: float) -> None:
#             self.val = new_val
#         self.getter = lambda: self.val
#         self.setter = _setter_func
    
# forcing_func: ValueManager = ValueManager()
# response_func: ValueManager = ValueManager()

#-------------------------------------------------

zeta: float = 0.4
r: float = 0.1
f: float = 5.0

k1: float = zeta/(np.pi*f)
k2: float = 1/(2*np.pi*f)**2
k3: float = (zeta*r)/(2*np.pi*f)

b_coeffs = np.array([1.0, k3]) 
a_coeffs = np.array([1.0, k1, k2])

system: m.System = m.System(b_coeffs, a_coeffs)
f: Callable[[float], float] = lambda t: 1.0 if np.sin(t) > 0 else 0.0  
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