
from numpy import *
from scipy.stats import poisson
import matplotlib
matplotlib.use('gtk3agg')

from matplotlib import pylab as plt

'''
exponential for pt
dN_dpt = exp(-pt)
'''

def generate_particles_from_single( v, psi, N ):
    dN_dphi = lambda phi: 1 + sum( [ 2*v[n]*cos( (n+2)*(phi-psi[n]) ) for n in range(len(v)) ] )
    dN_dphi_max = 1 + sum( [ 2*v[n] for n in range(len(v)) ] )
    particles = []
    while True:
        phi = random.uniform(0, 2*pi)
        prob = dN_dphi(phi)/dN_dphi_max
        if random.random() < prob:
            particles.append(phi)
            if len(particles) == N:
                break
    return particles

def generate_particles_restrained( v, psi, N ):
    dN_dphi = lambda phi: 1 + sum( [ 2*v[n]*cos( (n+2)*(phi-psi[n]) ) for n in range(len(v)) ] )
    dN_dphi_max = 1 + sum( [ 2*v[n] for n in range(len(v)) ] )
    prod_prob = 1
    particles = []
    while True:
        phi = random.uniform(0, 2*pi)
        prob = dN_dphi(phi)/dN_dphi_max
        if random.random() < prob*prod_prob:
            prod_prob *= prob
            particles.append(phi)
            if len(particles) == N:
                break
    return particles

prob_func = generate_particles_from_single
N = 4*int(1e4)
Psi = [0, 0]
v = [0.05, 0.0] #v2, v3

particles_ensemble = []
reconstructed_v2 = []
reconstructed_v3 = []
Nevts = 100
for i in range(Nevts):
    particles = prob_func( v, Psi, N )
    particles_ensemble += particles
    ret = array([ sum([ exp(1j*n*phi) for phi in particles ]) for n in [2,3] ])/len(particles)
    reconstructed_v2.append( ret[0] )
    reconstructed_v3.append( ret[1] )

reconstructed_v2 = array(reconstructed_v2)
reconstructed_v3 = array(reconstructed_v3)

angular_bins = lambda N: linspace(0, 2*pi, N)

nbins=int(Nevts*N)/1000
fig = plt.figure()
ax = fig.add_subplot(111)
hist, edges = histogram( particles_ensemble, bins=angular_bins(nbins) )
x = .5*(edges[:-1]+edges[1:])
ax.step( x, hist, where='mid', label='phi' )
ax.legend()
fig.savefig('phi.png')

nbins=10
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist( abs(reconstructed_v2), bins=nbins, histtype='step', label='v2' )
ax.hist( abs(reconstructed_v3), bins=nbins, histtype='step', label='v3' )
ax.legend()
fig.savefig('vn.png')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist( angle(reconstructed_v2), bins=angular_bins(nbins), histtype='step', label='v2' )
ax.hist( angle(reconstructed_v3), bins=angular_bins(nbins), histtype='step', label='v3' )
ax.legend()
fig.savefig('angle.png')
