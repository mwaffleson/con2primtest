# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:30:19 2026

@author: mickl
"""

import numpy as np
import matplotlib.pyplot as plt




eps=np.finfo(float).eps


data=np.genfromtxt("data.txt")



# plt.xscale("log")
# plt.yscale("log")
steps=1000

data=data.reshape((steps,steps,3))

# %%

# plt.xlim((data[0,0],data[0,-1]))
# plt.ylim((data[1,0],data[1,-1]))




plt.pcolormesh(data[:,:,0],data[:,:,1],(data[:,:,2]),shading="nearest", norm="log")
plt.colorbar()
plt.xlabel("$log_{10}$(z)")
plt.ylabel("$log_{10}(\epsilon)$")

plt.show()
