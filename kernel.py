import numpy as np

def cspline_unormalized( q ):
    o = np.zeros_like(q)
    less_than_1 = q < 1
    q_1 = q[less_than_1]
    between_1_and_2 = np.all( [q >= 1, q < 2], axis=0 )
    q_1_2 = q[between_1_and_2]
    o[less_than_1] = 1 - 1.5*q_1**2+.75*q_1**3
    o[between_1_and_2] = .25*(2-q_1_2)**3
    return o

def cspline( dist, dim, h ):
    return ( 2./3, 10*np.pi/7, 1./np.pi )[dim-1] / h**dim * cspline_unormalized( dist/h )

#def cspline(h,dim):
    #K = (2./3, 10*np.pi/7, 1./np.pi)[dim-1]/h**dim
    #return (
        #np.vectorize( lambda q: K* ( 0 if q > 2 else ( .25*(2-q)**3 if q > 1 else 1 - 1.5*q**2+.75*q**3 ) ) ), 
        #np.vectorize( lambda q: K* ( 0 if q > 2 else ( -3*.25*(2-q)**2 if q > 1 else - 3*q+ 3*.75*q**2 ) ) ),
        #np.vectorize( lambda q: K* ( 0 if q > 2 else ( 6*.25*(2-q) if q > 1 else -3 + 6*.75*q ) ) )
        #)

class Kernel:
    def __init__():
        self.h = None
        self.mode = None
        self.dimension = None
        self.r = None
        self.W = None
        self.normalization = None
    
    def set_dimension( d ):
        self.dimension = d
        return
    
    def set_positions( r ):
        self.r = r
        return
    
    def set_length( h ):
        self.h = h
        self.normalization = self.compute_normalization()
        return
    
    def set_mode( mode ):
        if mode is 'cspline':
            self.mode = mode
            self.W = self.__cspline__
            return
        raise Exception('mode "%s" is not defined' % mode)
            
    
    def compute_normalization(self):
        if self.mode is 'cspline':
            return ( 2./3, 10*np.pi/7, 1./np.pi )[self.dimension-1] / self.h**self.dimension
        return None
    
    def __cspline__( q ):
        o = np.zeros_like(q)
        less_than_1 = q < 1
        between_1_and_2 = np.all( [q >= 1, q < 2], axis=0 )
        o[less_than_1] = 1 - 1.5*q[less_than_1]**2+.75*q[less_than_1]**3
        o[between_1_and_2] = .25*(2-q[between_1_and_2])**3
        return normalization*o
        

def generate_kernel_function(mode='cspline', length=.1, dimension=1):
    inv_length = 1./length
    if mode is 'cspline':
        normalization = ( 2./3, 10*np.pi/7, 1./np.pi )[dimension-1] * inv_length**dimension
        def func( x, r ):
            q = (r - x)*inv_length
            q = np.sqrt( np.einsum('ij,ij->i', q, q) )
            o = np.zeros_like(q)
            less_than_1 = q < 1
            between_1_and_2 = np.all( [q >= 1, q < 2], axis=0 )
            o[less_than_1] = 1 - 1.5*q[less_than_1]**2+.75*q[less_than_1]**3
            o[between_1_and_2] = .25*(2-q[between_1_and_2])**3
            return normalization*o
        return func
    raise Exception('mode "%s" is not defined' % mode)

def cspline_opt(h,dim):
    K = (2./3, 10*np.pi/7, 1./np.pi)[dim-1]/h**dim
    def d0( q ):
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        res[less_than_1] = 1 - 1.5*q[less_than_1]**2+.75*q[less_than_1]**3
        res[between_1_and_2] = .25*(2-q[between_1_and_2])**3
        return K*res
    def d1( q ):
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        res[less_than_1] = -3*q[less_than_1]+ 3*.75*q[less_than_1]**2
        res[between_1_and_2] = -3*.25*(2-q[between_1_and_2])**2
        return K*res
    def d2( q ):
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        res[less_than_1] = -3 + 6*.75*q[less_than_1]
        res[between_1_and_2] = 6*.25*(2-q[between_1_and_2])
        return K*res
    return (d0,d1,d2)

def cspline_opth(dim):
    K = (2./3, 10*np.pi/7, 1./np.pi)[dim-1]
    def d0( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        res[less_than_1] = 1 - 1.5*q[less_than_1]**2+.75*q[less_than_1]**3
        res[between_1_and_2] = .25*(2-q[between_1_and_2])**3
        return K*res/h**dim
    def d1( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        res[less_than_1] = -3*q[less_than_1]+ 3*.75*q[less_than_1]**2
        res[between_1_and_2] = -3*.25*(2-q[between_1_and_2])**2
        return K*res/h**dim
    def d2( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        res[less_than_1] = -3 + 6*.75*q[less_than_1]
        res[between_1_and_2] = 6*.25*(2-q[between_1_and_2])
        return K*res/h**dim
    def d3( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        res[less_than_1] = 6*.75
        res[between_1_and_2] = -6*.25
        return K*res/h**dim
    def dint( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_2 = q<2
        res[less_than_2] = 1
        return res
    return (d0,d1,d2,d3)

def qspline_opth(dim):
    K = (1./120., 7/(478*np.pi), 3/(359*np.pi))[dim-1]
    def d0( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        between_2_and_3 = np.all( [q>=2,q<3], axis=0 )
        res[less_than_1] = (3 - q[less_than_1])**5 - 6*(2 - q[less_than_1])**5 + 15*(1 - q[less_than_1])**5
        res[between_1_and_2] = (3 - q[between_1_and_2])**5 - 6*(2 - q[between_1_and_2])**5
        res[between_2_and_3] = (3 - q[between_2_and_3])**5
        return K*res/h**dim
    def d1( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        between_2_and_3 = np.all( [q>=2,q<3], axis=0 )
        res[less_than_1] = -5*( (3 - q[less_than_1])**4 - 6*(2 - q[less_than_1])**4 + 15*(1 - q[less_than_1])**4 )
        res[between_1_and_2] = -5*( (3 - q[between_1_and_2])**4 - 6*(2 - q[between_1_and_2])**4 )
        res[between_2_and_3] = -5*(3 - q[between_2_and_3])**4
        return K*res/h**dim
    def d2( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        between_2_and_3 = np.all( [q>=2,q<3], axis=0 )
        res[less_than_1] = 4*5*( (3 - q[less_than_1])**3 - 6*(2 - q[less_than_1])**3 + 15*(1 - q[less_than_1])**3 )
        res[between_1_and_2] = 4*5*( (3 - q[between_1_and_2])**3 - 6*(2 - q[between_1_and_2])**3 )
        res[between_2_and_3] = 4*5*(3 - q[between_2_and_3])**3
        return K*res/h**dim
    def d3( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_1 = q<1
        between_1_and_2 = np.all( [q>=1,q<2], axis=0 )
        between_2_and_3 = np.all( [q>=2,q<3], axis=0 )
        res[less_than_1] = -3*4*5*( (3 - q[less_than_1])**2 - 6*(2 - q[less_than_1])**2 + 15*(1 - q[less_than_1])**2 )
        res[between_1_and_2] = -3*4*5*( (3 - q[between_1_and_2])**2 - 6*(2 - q[between_1_and_2])**2 )
        res[between_2_and_3] = -3*4*5*(3 - q[between_2_and_3])**2
        return K*res/h**dim
    def dint( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_2 = q<2
        res[less_than_2] = 1
        return res
    return (d0,d1,d2,d3)

def exp_opth(dim):
    #K = (2./3, 10*np.pi/7, 1./np.pi)[dim-1]
    K = 1./(2*np.pi)
    def d0( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_2 = q<2
        res[less_than_2] = np.exp(-q[less_than_2]**2)
        return K*res/h**dim
    return (d0,0,0,0)
    
def liq(h,dim):
    A = -1.458
    B = 3.790
    C = -2.624
    D = -0.2915
    E = 0.5831
    F = 0.6500
    W = lambda q: 0 if q > 2 else ( A*(q/2)**4 + B*(q/2)**3 + C*(q/2)**2 + D*(q/2) + E if q > .6 else F-q/2 )
    Wp = lambda q: 0 if q > 2 else ( 2*A*(q/2)**3 + 3*B/2*(q/2)**2 + C*(q/2) + D/2 if q > .6 else -.5 )
    Wpp = lambda q: 0 if q > 2 else ( 3*A*(q/2)**2 + 3*B/2*(q/2) + C/2 if q > .6 else 0 )
    #N = 2*scipy.integrate.quad( lambda q: W(q), 0, 3 )[0]
    N = 1
    return (
        np.vectorize( lambda q: W(q)/N/h**dim ),
        np.vectorize( lambda q: Wp(q)/N/h**dim ),
        np.vectorize( lambda q: Wpp(q)/N/h**dim )
        )

def liq_opth(dim):
    A = -1.458
    B = 3.790
    C = -2.624
    D = -0.2915
    E = 0.5831
    F = 0.6500
    #W = lambda q: 0 if q > 2 else ( A*(q/2)**4 + B*(q/2)**3 + C*(q/2)**2 + D*(q/2) + E if q > .6 else F-q/2 )
    #Wp = lambda q: 0 if q > 2 else ( 2*A*(q/2)**3 + 3*B/2*(q/2)**2 + C*(q/2) + D/2 if q > .6 else -.5 )
    #Wpp = lambda q: 0 if q > 2 else ( 3*A*(q/2)**2 + 3*B/2*(q/2) + C/2 if q > .6 else 0 )
    #N = 2*scipy.integrate.quad( lambda q: W(q), 0, 3 )[0]
    N = 1
    def _0( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_06 = q<.6
        between_06_and_2 = np.all( [q>=.6,q<2], axis=0 )
        res[less_than_06] = F-q[less_than_06]/2
        res[between_06_and_2] = A*(q[between_06_and_2]/2)**4 + B*(q[between_06_and_2]/2)**3 + C*(q[between_06_and_2]/2)**2 + D*(q[between_06_and_2]/2) + E
        return res/N/h**dim
    def _1( d, h ):
        q = d/h
        res = np.zeros_like(q)
        less_than_06 = q<.6
        between_06_and_2 = np.all( [q>=.6,q<2], axis=0 )
        res[less_than_06] = -1./2
        res[between_06_and_2] = 2*A*(q[between_06_and_2]/2)**3 + 3*B/2*(q[between_06_and_2]/2)**2 + C*(q[between_06_and_2]/2) + D/2
        return res/N/h**dim
    return (
        _0,
        _1,
        #np.vectorize( lambda q: Wpp(q)/N/h**dim )
        )

def derivatives(method, D=1):
    if method == 'cspline':
        method = cspline_opth(D)
    if method == 'qspline':
        method = qspline_opth(D)
    if method == 'liq':
        method = liq_opth(D)
    
    return [lambda x, r, h: method[0]( np.abs(x-r), h ), 
        lambda x, r, h: method[1]( np.abs(x-r), h )*np.sign(x-r)/h, 
        lambda x, r, h: method[2]( np.abs(x-r), h )*np.sign(x-r)**2/h**2, 
        lambda x, r, h: method[3]( np.abs(x-r), h )*np.sign(x-r)**3/h**3, #not complete
        ]

