
from __future__ import print_function
from numpy import *
from scipy.stats import poisson
import matplotlib
matplotlib.use('gtk3agg')

from matplotlib import pylab as plt

'''
exponential for pt
dN_dpt = exp(-pt)
'''
def sign(x):
    ret = zeros_like(x)
    ret[x>0] = 1
    ret[x<0] = -1
    return ret

def angular_bins(N):
    return linspace(-pi, pi, N, endpoint=True)

def generate_particles_from_single( collective_flow, reaction_plane, N ):
    dN_dphi = lambda phi: 1 + sum( [ 2*v_n*cos( (n+2)*(phi-Psi_n) ) for n, (v_n, Psi_n) in enumerate(zip(collective_flow, reaction_plane)) ], axis=0 )
    dN_dphi_max = 1 + sum( [ 2*v_n for v_n in collective_flow ] )
    particles = []
    while True:
        phi = random.uniform(-pi, pi, N - len(particles) )
        prob = dN_dphi(phi)/dN_dphi_max
        particles = concatenate( ( particles, phi[prob > random.random( N - len(particles) )] ) )
        if len(particles) == N:
            break
    return array(particles, dtype=[('phi', float)]).view(recarray)

def compute_collective_flow( particles, number_of_harmonics = 3 ):
    harmonics = arange(0, number_of_harmonics+1, 1)
    return sum( exp( 1j*harmonics[None,:] * particles.phi[:,None]), axis=0 )/len(particles)

def make_bins( particles, number_of_bins, number_of_harmonics = 3 ):
    harmonics = arange(0, number_of_harmonics+1, 1)
    bins = angular_bins(number_of_bins)
    inds = digitize( particles.phi, bins )
    particles_in_bins = [ particles[ inds==i ] for i in unique(inds) ]
    print( 'make_bins', particles_in_bins.shape )
    return sum( exp( 1j*harmonics[None,:] * particles.phi[:,None]), axis=0 )

prob_func = generate_particles_from_single
N = 4*int(1e4)
Psi = [0, 1]
v = [0.5, 0.2] #v2, v3

all_phi = []
V_2 = []
V_3 = []
Nevts = 100
for i in range(Nevts):
    particles = prob_func( v, Psi, N )
    make_bins( particles, 20 )
    all_phi = concatenate((all_phi, particles.phi ))
    collective_flow = compute_collective_flow( particles, 3 )
    #print( collective_flow )
    V_2 = concatenate( (V_2, [collective_flow[2]]) )
    V_3 = concatenate( (V_3, [collective_flow[3]]) )

def make_pdf( X, nbins ):
    hist, edges = histogram( X, bins=nbins )
    x = .5*(edges[:-1]+edges[1:])
    return x, hist

def make_angular_pdf( X, nbins ):
    return make_pdf( X, angular_bins(nbins) )

nbins=100

def plot_histogram( fname, label_X, nbins ):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for label, Xi in label_X:
        bins = linspace(min(Xi), max(Xi), nbins)
        x, y = make_pdf( Xi, bins )
        ax.step( x, y, where='mid', label='{}\n mean={:4f}\n std={:4f}'.format( label, mean(Xi), std(Xi) ) )
    ax.legend()
    ax.grid(True)
    fig.savefig( fname )
    print( 'savefig {}'.format(fname) ) 
    
phi_bins = angular_bins(nbins)
plot_histogram( 'phi.pdf', [(r'$\phi$', all_phi)], 1000 )

Psi_2, Psi_3 = angle(V_2)/2, angle(V_3)/3
plot_histogram( 'event_plane.pdf', [(r'$\Psi_{2}$', Psi_2 ), (r'$\Psi_{3}$', Psi_3 )], 10 )
plot_histogram( 'collective_flow.pdf', [(r'$v_{2}$', abs(V_2) ), (r'$v_{3}$', abs(V_3) )], 10 )
