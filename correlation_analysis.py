
from __future__ import print_function
from numpy import *
from numpy.linalg import eig
from scipy.stats import poisson
import matplotlib
#matplotlib.use('gtk3agg')

from matplotlib import pylab as plt
from scipy.sparse import csr_matrix

from Timer import Timer
'''
exponential for pt
dN_dpt = exp(-pt)
'''

def binned_statistic(x, values, func, nbins, range):
    '''The usage is nearly the same as scipy.stats.binned_statistic''' 

    r0, r1 = range

    digitized = (float(nbins)/(r1 - r0)*(x - r0)).astype(int)
    values = values[digitized < nbins].astype(float)
    digitized = digitized[digitized < nbins]
    
    N = len(values)
    S = csr_matrix((values, [digitized, arange(N)]), shape=(nbins, N))

    return [func(group) for group in split(S.data, S.indptr[1:-1])]

class Events:
    @property
    def harmonics(self):
        try:
            return self.__harmonics
        except AttributeError:
            print('harmonics not set yet')
            raise AttributeError
    
    @harmonics.setter
    def harmonics(self, n):
        if hasattr(self, '__harmonics'):
            print('harmonics already set')
            raise AttributeError
        self.__harmonics = n
    
    @property
    def maxpt(self):
        try:
            return self.__maxpt
        except AttributeError:
            print('maxpt not set yet')
            raise AttributeError
    
    @maxpt.setter
    def maxpt(self, maxpt):
        if hasattr(self, '__maxpt'):
            print('maxpt already set')
            raise AttributeError
        self.__maxpt = maxpt
    
    @property
    def dpt(self):
        return self.__dpt
    
    @dpt.setter
    def dpt(self, dpt):
        if hasattr(self, '__dpt'):
            print('dpt already set')
            raise AttributeError
        self.__dpt = dpt
        
    @property
    def pt_bins(self):
        try:
            return self.__pt_bins
        except AttributeError:
            self.__pt_bins = arange( 0, self.maxpt+self.dpt, self.dpt )
        return self.__pt_bins
        
    @property
    def pt_center(self):
        try:
            return self.__pt_center
        except AttributeError:
            self.__pt_center = (self.pt_bins[1:] + self.pt_bins[:-1])/2.
        return self.__pt_center
    
    def compute_Q_n( self ):
        print('computing Q_n')
        inds = [ digitize( event.pt, self.pt_bins ) for event in self.events ]
        bin_inds = unique(inds)[:-1]
        M = array( [ [ count_nonzero( ind == i ) for i in bin_inds ] for ind in inds ] )
        Q = lambda phi: sum( exp( (1j*self.harmonics[None,:] * phi[:,None]).astype(complex) ), axis=0 )
        Q_n = array( [ [ Q( event.phi[ind == i] ) for i in bin_inds ] for ind, event in zip(inds, self.events) ] )
        return M, Q_n/(2*pi*self.dpt)

    def compute_Q_n__csr( self ):
        print('computing Q_n csr')
        Q = lambda phi: sum( exp( (1j*self.harmonics[None,:] * phi[:,None]).astype(complex) ), axis=0 )
        Q_n = array( [ binned_statistic(event.pt, event.phi, Q, len(self.pt_bins), [0, self.maxpt+self.dpt] ) for event in self.events ] )
        return Q_n/(2*pi*self.dpt)

    @property
    def M_i_p(self):
        try:
            return self.__M_i_p
        except AttributeError:
            self.__M_i_p, self.__Q_i_p_n = self.compute_Q_n__csr()
        return self.__M_i_p

    @property
    def std_M_p(self):
        try:
            return self.__std_M_p
        except AttributeError:
            self.__std_M_p = std(self.M_i_p, axis=0)
            return self.__std_M_p

    @property
    def mean_M_p(self):
        try:
            return self.__mean_M_p
        except AttributeError:
            self.__mean_M_p = mean(self.M_i_p, axis=0)
            return self.__mean_M_p

    @property
    def Q_i_p_n(self):
        try:
            return self.__Q_i_p_n
        except AttributeError:
            self.__M_i_p, self.__Q_i_p_n = self.compute_Q_n()
        return self.__Q_i_p_n

    @property
    def mean_Q_p_n(self):
        try:
            return self.__mean_Q_p_n
        except AttributeError:
            self.__mean_Q_p_n = mean(self.__Q_i_p_n, axis=0)
        return self.__mean_Q_p_n
    

    @property
    def V_i_p_p_nDelta(self):
        try:
            return self.__V_i_p_p_nDelta
        except AttributeError:
            self.__V_i_p_p_nDelta = self.Q_i_p_n[:,:,None,:]*self.Q_i_p_n.conj()[:,None,:,:]
        return self.__V_i_p_p_nDelta
    
    @property
    def V_p_p_nDelta(self):
        try:
            return self.__V_p_p_nDelta
        except AttributeError:
            mean_prod = mean( self.Q_i_p_n[:,:,None,:]*self.Q_i_p_n.conj()[:,None,:,:], axis=0 )
            prod_mean = mean( self.Q_i_p_n[:,:,None,:], axis=0) * mean(self.Q_i_p_n.conj()[:,None,:,:], axis=0)
            self.__V_p_p_nDelta = mean_prod - prod_mean
            for ipt in range(len(self.pt_bins)-1):
                self.__V_p_p_nDelta[ipt,ipt,:] -= self.mean_M_p[ipt]/(2*pi*self.dpt)**2
        return self.__V_p_p_nDelta
    
    def V_p_p(self, n):
        return self.V_p_p_nDelta[:,:,n]
    
    def __compute_pca(self):
        print('computing PCA')
        return zip( *[ eig( self.V_p_p_nDelta[:,:,n] ) for n in self.harmonics ] )
        
    @property
    def eigvec_n_p_alpha(self):
        try:
            return self.__eigvec_n_p_alpha
        except AttributeError:
            self.__eigval_n_alpha, self.__eigvec_n_p_alpha = self.__compute_pca()
        return self.__eigvec_n_p_alpha
    
    @property
    def eigval_n_alpha(self):
        try:
            return self.__eigval_n_alpha
        except AttributeError:
            self.__eigval_n_alpha, self.__eigvec_n_p_alpha = self.__compute_pca()
        return self.__eigval_n_alpha
    
    @property
    def V_n_p_alpha(self):
        try:
            return self.__V_n_p_alpha
        except AttributeError:
            item = lambda n: (self.eigval_n_alpha[n]/sqrt(abs(self.eigval_n_alpha[n])))[None,:] * self.eigvec_n_p_alpha[n]
            self.__V_n_p_alpha = array([ item(n) for n in self.harmonics ])
        return self.__V_n_p_alpha
    
    def V_p(self, n, alpha):
        return self.V_n_p_alpha[n,:,alpha]
    
    @property
    def v_n_p_alpha(self):
        try:
            return self.__v_n_p_alpha
        except AttributeError:
            self.__v_n_p_alpha = self.V_n_p_alpha/(self.mean_M_p[None,:,None]/(2*pi*self.dpt))
        return self.__v_n_p_alpha
    
    def v_p(self, n, alpha):
        return self.v_n_p_alpha[n,:,alpha]



class plot:
    def __init__(self, fname):
        self.fname = fname
        
    def __enter__(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        return self.ax
    
    def __exit__(self, type, value, traceback):
        self.fig.savefig(self.fname)
        print( 'saved {}'.format(self.fname) )

def sign(x):
    ret = zeros_like(x)
    ret[x>0] = 1
    ret[x<0] = -1
    return ret

def angular_bins(N, endpoint=True):
    return linspace(-pi, pi, N, endpoint=endpoint)

def random_func( f, f_max, a, b, N ):
    x = []
    while True:
        _x = random.uniform(a, b, N - len(x) )
        prob = f(_x)/f_max
        x = concatenate( ( x, _x[prob > random.random( N - len(x) )] ) )
        if len(x) == N:
            break
    return array(x)

def random_func2( y, f, f_max, a, b, N ):
    x = zeros_like(y)*nan
    unassigned = isnan(x)
    while True:
        _x = random.uniform(a, b, sum(unassigned) )
        _y = y[unassigned]
        prob = f(_y,_x)/f_max(_y)
        assign = prob > random.random( sum(unassigned) )
        to_assign = array(unassigned)
        to_assign[to_assign] = assign
        x[to_assign] = _x[assign]
        unassigned = isnan(x)
        if sum(unassigned) == 0:
            break
    return array(x, dtype=float)
    
    
def generate_particles_from_single( beta, collective_flow_func, N ):
    Psi = random.uniform( -pi, pi, len(collective_flow_func) ) 
    def dN_dptdphi( pt, phi ):
        return 1 + sum( [ 2*cos( n*(phi-Psi[n]) ) * v[n](pt) for n in range(2,len(Psi)) ], axis=0 )
    dN_dptdphi_max = lambda pt: 1 + sum( [ 2*abs(v[n](pt)) for n in range(2,len(Psi)) ], axis=0 )
    pt = random.exponential( beta, N )
    phi = random_func2( pt, dN_dptdphi, dN_dptdphi_max, -pi, pi, N )
    return array( zip(pt, phi), dtype=[('pt', float), ('phi', float)]).view(recarray)


N = int(2e5)
Nevts = 200

p2max = 3
v2max = 0.1

p3max = 2.
v3max = 0.05

v = [
    lambda p: ones_like(p),
    lambda p: 0,
    lambda p, p2max=p2max, v2max=v2max: (p2max**2 - (p-p2max)**2 )*v2max/p2max**2,
    lambda p, p3max=p3max, v3max=v3max: (p3max**2 - (p-p3max)**2 )*v3max/p3max**2,
]

data = Events()
number_of_harmonics = 3
harmonics = arange(0, number_of_harmonics+1, 1)
data.harmonics = harmonics

beta = .75
data.dpt = .2
data.maxpt = 5.

events = array( [ generate_particles_from_single( beta, v, N ) for i in range(Nevts) ], dtype=[('pt', object), ('phi', object)] ).view(recarray)
data.events = events

with Timer('csr'):
    Q__scr = data.compute_Q_n__csr()
    print( Q__scr.shape )
with Timer('Q'):
    _, Q = data.compute_Q_n()
    print( Q.shape )


with plot('M_p') as p:
    p.set_title('M(p)')
    p.errorbar( data.pt_center, data.mean_M_p, yerr = data.std_M_p, fmt='.' )
    p.errorbar( data.pt_center, data.mean_Q_p_n[:,0].real*(2*pi*data.dpt), yerr = std(data.Q_i_p_n[:,:,0]*(2*pi*data.dpt).real, axis=0), fmt='.' )
    p.errorbar( data.pt_center, data.mean_Q_p_n[:,1].real*(2*pi*data.dpt), yerr = std(data.Q_i_p_n[:,:,1]*(2*pi*data.dpt).real, axis=0), fmt='.' )
    p.set_xlabel(r'$p_{\rm T}$')
    p.set_yscale('log')

print( 'M(p)', std(data.M_i_p, axis=0) )


for n in [2,3]:
    with plot('v{}_p_0_.png'.format(n)) as p:
        p.plot( data.pt_center, v[n](data.pt_center), label=r'$v_{}(p_t)$'.format(n) )
        p.errorbar( data.pt_center, data.v_p(n,0).real, fmt='.', label=r'$v_{}(p_t)^{{ (0) }}$'.format(n), xerr=data.dpt/2 )
        p.errorbar( data.pt_center, data.v_p(n,1).real, fmt='.', label=r'$v_{}(p_t)^{{ (1) }}$'.format(n), xerr=data.dpt/2 )
        p.errorbar( data.pt_center, data.v_p(n,2).real, fmt='.', label=r'$v_{}(p_t)^{{ (2) }}$'.format(n), xerr=data.dpt/2 )
        p.legend()
        p.grid(True)

exit()
