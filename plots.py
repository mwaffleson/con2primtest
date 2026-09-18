# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:30:19 2026

@author: mickl
"""

import numpy as np
import matplotlib.pyplot as plt





data=np.genfromtxt("data.txt")




steps=1000




data=data.reshape((steps,steps,data.shape[1]))

# %%

plt.style.use('dark_background')




plt.pcolormesh(data[:,:,0],data[:,:,1],data[:,:,2],shading="nearest", norm="log")
plt.colorbar()
plt.xlabel("$log_{10}$(z)")
plt.ylabel("$log_{10}(\epsilon)$")

plt.show()

plt.xlabel("$log_{10}$(z)")
plt.ylabel("$log_{10}(\epsilon)$")
plt.pcolormesh(data[:,:,0],data[:,:,1],data[:,:,4],shading="nearest")
plt.colorbar()
plt.show()


plt.xlabel("$log_{10}$(z)")
plt.ylabel("$log_{10}(\epsilon)$")
plt.pcolormesh(data[:,:,0],data[:,:,1],data[:,:,5],shading="nearest")
plt.colorbar()
plt.show()







