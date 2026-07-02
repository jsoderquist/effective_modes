import numpy as np
from math import pi, e
import matplotlib.pyplot as plt
from matplotlib import ticker, rcParams

# au conversions
fstoau = 41.341         # 1 fs = 41.341 a.u.
cmtoau = 4.556335e-06   # 1 cm^-1 = 4.556335e-06 a.u. - this converts a 1/lambda value to 2*pi*f and then switches to au
autoK  = 3.1577464e+05 
angtoau = 1/0.5292      # angstrom to au conversion
c = 137.0359991         # speed of light in au (137/sec)

fontSize = 12
rcParams["axes.labelsize"] = fontSize

###################### Figure 1 ############################
# system parameters
m = 1 # mode order
wwc = np.linspace(0.95,1.2,500) # omega/omega_c
alpha = np.linspace(-30,30,500)*pi/180 # incidence angle
wwc,alpha = np.meshgrid(wwc,alpha) # make it easy to do math
Q = np.array([[50,100],[200,5000]]) # quality factor
phi = pi # complex phase of rho
qzLc = m*pi*wwc*np.cos(alpha)

# plot figure 1
fig, ax = plt.subplots(2,2,figsize = (5.5,4.5),layout="constrained")
for i in range(2):
    for j in range(2):
        # calculate enhancement at a given value of Q
        rho = e**(-m*pi/2/Q[i,j])
        R = abs((1 + rho*e**(1j*(qzLc+phi)))/(1 - rho**2*e**(2*1j*(qzLc+phi))))

        # plot enhancement
        im = ax[i,j].contourf(alpha*180/pi,wwc,R,np.logspace(np.log10(1.6),np.log10(3184), 500), 
                 locator=ticker.LogLocator(subs=(1,)),cmap="inferno") 
        ax[i,j].set_box_aspect(1) # make each subplot square
        ax[i,j].text(0, 1.185, " Quality Factor: {}".format(Q[i,j]), horizontalalignment='center', verticalalignment='center',color='white') # label figures with quality factor
        # print(R.min(),R.max())
cbar = fig.colorbar(im, ax=ax,ticks=[10**0,10**1,10**2,10**3]) # log scale colorbar
plt.setp(ax[0,0].get_xticklabels(), visible=False)
plt.setp(ax[0,1].get_xticklabels(), visible=False)
plt.setp(ax[0,1].get_yticklabels(), visible=False)
plt.setp(ax[1,1].get_yticklabels(), visible=False)
ax[0,0].text(-29, 1.185, "(a)", verticalalignment='center',color='white')
ax[0,1].text(-29, 1.185, "(b)", verticalalignment='center',color='white')
ax[1,0].text(-29, 1.185, "(c)", verticalalignment='center',color='white')
ax[1,1].text(-29, 1.185, "(d)", verticalalignment='center',color='white')
cbar.set_label("|$R_q$|")
cbar.set_ticklabels(["$10^0$","$10^1$","$10^2$","$10^3$"])
ax[1,0].set_xlabel("α (deg)")
ax[1,1].set_xlabel("α (deg)")
ax[0,0].set_ylabel("$ω$/$ω_c$")
ax[1,0].set_ylabel("$ω$/$ω_c$")
plt.savefig("enhancement.png")



####################### Figure 2 ##########################
wc = 0.1 # angular frequency the cavity is tuned to (in au)
w = np.linspace(0.9,1.2,500)*wc
n = 1 # effective refractive index inside the cavity
q = n*w/c
Qs = [50, 500, 5000] # quality factors to plot 
phi = pi # complex phase of reflection coefficient

fig, ax = plt.subplots(2,figsize = (4.5,4.5),layout="constrained")
for i in range(1):
    m = i+1 # mode order
    Lc = m*pi*c/wc # length of the cavity
    xmin = m*pi+phi-pi/2+10**-8 # x is q_z*L_c+ϕ - make sure it will be evaluated from the right side
    xmax = q*Lc+phi

    for Q in Qs:
        r = e**(-m*pi/2/Q)
        prefactor = w/(8*pi**3*n*c**2*Lc*abs(r-1)*abs(r+1))*c/n
        D = prefactor*(np.arctan((r**2+1)*np.tan(xmax)/abs(r-1)/abs(r+1)) + \
            np.arctan(2*r*np.sin(xmax)/abs(r-1)/abs(r+1)) - \
            np.arctan((r**2+1)*np.tan(xmin)/abs(r-1)/abs(r+1)) - \
            np.arctan(2*r*np.sin(xmin)/abs(r-1)/abs(r+1)))
        # ax[i].semilogy(w/wc,D,label="Quality Factor: {}".format(Q))
        ax[0].semilogy(w/wc,D,label="Quality Factor: {}".format(Q),linewidth=3)
    # ax[i].text(0.91, max(D)*0.9, "m = {}".format(m), verticalalignment='center')
    # ax.text(0.91, max(D)*0.9, "m = {}".format(m), verticalalignment='center')
        
# ax.set_xlabel("$ω$/$ω_c$")
# ax.set_ylabel("AR-LDOM (au)")
plt.setp(ax[0].get_xticklabels(), visible=False)
# plt.legend(prop={'size': 13})
# ax.tick_params(axis='both', which='major', labelsize=fontSize)
# ax.tick_params(axis='both', which='minor', labelsize=8)
# plt.savefig("ARLDOM.png")



####################### Figure 3 - now included in figure 2 ##########################
wc = 0.1 # angular frequency the cavity is tuned to (in au)
w = np.linspace(0.9,1.2,1000)*wc
Δw = wc/10
n = 1 # effective refractive index inside the cavity
q = n*w/c
Δq = n*wc/c/10 # uncertainty in q - used to define Δθ
Qs = [50, 500, 5000] # quality factors to plot 
phi = pi # complex phase of reflection coefficient

# fig, ax = plt.subplots(1,figsize = (4.5,4.5),layout="constrained")
for i in range(1):
    m = i+1 # mode order
    Lc = m*pi*c/wc # length of the cavity
    xmin = m*pi+phi-pi/2+10**-8 # x is q_z*L_c+ϕ - make sure it will be evaluated from the right side
    xmax = q*Lc+phi

    for Q in Qs:
        r = e**(-m*pi/2/Q)
        prefactor = w*Δw/(8*pi**3*n*c**2*Lc*abs(r-1)*abs(r+1))
        qp = np.sqrt(abs(q**2 - (n*wc/c)**2)) # magnitude of q-parallel
        # qp[q<n*wc/c] = np.sqrt(q[q<n*wc/c]**2 - (n*wc/c/(m+1))**2) # q-parallel of peak at next mode order where no peak exists in this mode
        nm = np.floor(pi/np.arcsin(Δq/2/qp)) # number of resolvable, directional modes in x-y plane
        nm[Δq>2*qp] = 1
        Δθ = 2*pi/nm # angular extent of each mode
        # Δθ=2*pi
        lnOverV = prefactor*(np.arctan((r**2+1)*np.tan(xmax)/abs(r-1)/abs(r+1)) + \
            np.arctan(2*r*np.sin(xmax)/abs(r-1)/abs(r+1)) - \
            np.arctan((r**2+1)*np.tan(xmin)/abs(r-1)/abs(r+1)) - \
            np.arctan(2*r*np.sin(xmin)/abs(r-1)/abs(r+1)))*Δθ
        ax[1].semilogy(w/wc,lnOverV,label="Quality Factor: {}".format(Q),linewidth=3)
        
ax[1].set_xlabel("$ω$/$ω_c$")
# ax[2].set_ylabel("$N_m$")
ax[1].set_ylabel("$\ell_n$/$V$ (au)")
ax[0].set_ylabel("$AR-LDOM$ (au)")
# ax.set_ylabel("$\ell_n$/$V$ (au)")
ax[0].text(0.9, 10**-5.4, "(a)", verticalalignment='center',size=fontSize)
ax[1].text(0.9, 10**-9, "(b)", verticalalignment='center',size=fontSize)
# ax[2].text(0.9, 37, "(c)", verticalalignment='center')
# plt.setp(ax.get_xticklabels(), visible=False)
# plt.setp(ax[1].get_xticklabels(), visible=False)
# ax[2].plot(w/wc,nm)
ax[0].legend(prop={'size': 10})
# ax.tick_params(axis='both', which='major', labelsize=fontSize)
# ax.tick_params(axis='both', which='minor', labelsize=16)
# plt.savefig("lnOverV.png")
plt.savefig("figure2.svg")
