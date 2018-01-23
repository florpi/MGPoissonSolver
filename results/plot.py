import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

r, acc = np.loadtxt("acceleration256.txt",unpack=True, skiprows=0)

y = np.linspace(0.00001, 1000.,100)
x =(1./256)*np.ones(100)

r_theo = np.linspace(0.002,0.6,100)


gs = gridspec.GridSpec(2, 2,
                       width_ratios=[2, 0.3],
                       height_ratios=[4, 1],hspace=0,wspace=0)

ax1 = plt.subplot(gs[0:1])
ax2 = plt.subplot(gs[2])

ax1.plot(r_theo, 2./r_theo, label='Theory', color='black', linewidth=3)
ax1.scatter(r, acc, label='Data', color='red',s=8  )
ax1.legend()
ax1.set_ylabel("Acceleration")
ax1.plot(x,y,"--",linewidth=3,color='coral')

residual = (2./r - acc)
ax2.scatter(r, residual, s=1, color='black')
ax2.set_xlabel("distance")
ax2.set_ylabel("Residual")
ax2.plot(x,y,"--",linewidth=2,color='coral')
plt.savefig("results.png")
plt.show()




