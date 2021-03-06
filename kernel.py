import numpy as np
import weave

pi = np.pi

def cspline_unormalized( q ):
    o = np.zeros_like(q)
    less_than_1 = q < 1
    q_1 = q[less_than_1]
    between_1_and_2 = np.all( [q >= 1, q < 2], axis=0 )
    q_1_2 = q[between_1_and_2]
    o[less_than_1] = 1 - 1.5*q_1**2+.75*q_1**3
    o[between_1_and_2] = .25*(2-q_1_2)**3
    return o

def cspline_unormalized_scalar( q ):
    if q > 2: return 0
    if q > 1: return .25*( 2 - q )**3
    return 1 - 1.5 * q**2 + .75 * q**3

def cspline_vectorized( dist, dim, h ):
    return ( 2./3, 10*pi/7, 1./pi )[dim-1] / h**dim * np.vectorize( cspline_unormalized_scalar )( dist/h )

def cspline_vectorized_k( dist, dim, h ):
    k = 1./h**dim
    if dim == 1: k *= 2./3
    elif dim == 2: k *= 10*pi/7
    elif dim == 3: k *= 1./pi
    return k * np.vectorize( cspline_unormalized_scalar )( dist/h )

def cspline( dist, dim, h ):
    return ( 2./3, 10*np.pi/7, 1./np.pi )[dim-1] / h**dim * cspline_unormalized( dist/h )

def generate_cspline( dim, h = None ):
    k = ( 2./3, 10*np.pi/7, 1./np.pi )[dim-1] / h**dim
    return lambda dist, h = h, k = k: k*cspline_unormalized_weave( dist/h )

def generate_cspline_prime( dim, h = None ):
    k = ( 2./3, 10*np.pi/7, 1./np.pi )[dim-1] / h**(dim+1)
    return lambda dist, h = h, k = k: k*cspline_prime_unormalized_weave( dist/h )

def generate_qspline( dim, h = None ):
    k = (1./120., 7./(478*np.pi), 3./(359*np.pi))[dim-1] / h**dim
    return lambda dist, h = h, k = k: k*qspline_unormalized_weave( dist/h )

def generate_qspline_prime( dim, h = None ):
    k = (1./120., 7./(478*np.pi), 3./(359*np.pi))[dim-1] / h**(dim+1)
    return lambda dist, h = h, k = k: k*qspline_prime_unormalized_weave( dist/h )


def cspline_weave( dist, dim, h ):
    code = '''
    double k = 0;
    if( dim == 1 ){ k = 2./3; }
    if( dim == 2 ){ k = 10*pi/7; }
    if( dim == 3 ){ k = 1./pi; }
    for( int i=0; i < shape0; ++i){
        for( int j=0; j < shape1; ++j){
            int index = i + j*shape0;
            double Q = q[index];
            if( Q > 2 ){ 
                ret[index] = 0; 
            }
            else { 
                if( Q > 1 ){ 
                    ret[index] = k*.25*(2-Q)*(2-Q)*(2-Q); 
                    }
                else{
                    ret[index] = k*(1 - 1.5*Q*Q + .75*Q*Q*Q);
                }
            }
        }
    }
    '''
    q = dist/h
    ret = np.zeros_like( q )
    shape0 = ret.shape[0]
    shape1 = ret.shape[1]
    weave.inline( code, ['q', 'dim', 'ret', 'pi', 'shape0', 'shape1'], verbose=1 )
    return ret.reshape((shape0,shape1))/h**dim

def cspline_unormalized_weave( q ):
    code = '''
    for( int i=0; i < shape0; ++i ){
        for( int j=0; j < shape1; ++j ){
            int index = i + j*shape0;
            double Q = q[index];
            if( Q > 2 ){ ret[index] = 0; }
            else { 
                if( Q > 1 ){ ret[index] = .25*(2-Q)*(2-Q)*(2-Q); }
                else{ ret[index] = 1 - 1.5*Q*Q + .75*Q*Q*Q; }
            }
        }
    }
    '''
    ret = np.zeros_like( q )
    shape0 = ret.shape[0]
    shape1 = ret.shape[1]
    weave.inline( code, ['ret', 'q', 'shape0', 'shape1'], verbose=1, compiler = 'gcc', extra_compile_args=['-O3'] )
    return ret.reshape((shape0,shape1))

def cspline_prime_unormalized_weave( q ):
    code = '''
    for( int i=0; i < shape0; ++i ){
        for( int j=0; j < shape1; ++j ){
            int index = i + j*shape0;
            double Q = q[index];
            if( Q > 2 ){ ret[index] = 0; }
            else { 
                if( Q > 1 ){ ret[index] = - 3*.25*(2-Q)*(2-Q); }
                else{ ret[index] = - 3.*Q + 3*.75*Q*Q; }
            }
        }
    }
    '''
    ret = np.zeros_like( q )
    shape0 = ret.shape[0]
    shape1 = ret.shape[1]
    weave.inline( code, ['ret', 'q', 'shape0', 'shape1'], verbose=1, compiler = 'gcc', extra_compile_args=['-O3'] )
    return ret.reshape((shape0,shape1))

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

def qspline_unormalized_weave( q ):
    code = '''
    for( int i=0; i < shape0; ++i ){
        for( int j=0; j < shape1; ++j ){
            int index = i + j*shape0;
            double Q = q[index];
            if( Q > 3 ){ ret[index] = 0; }
            else {
                if( Q > 2 ){ ret[index] = pow(3-Q, 5); }
                else { 
                    if( Q > 1 ){ ret[index] = pow(3-Q, 5) - 6*pow(2-Q, 5); }
                    else{ ret[index] = pow(3-Q, 5) - 6*pow(2-Q, 5) + 15*(1-Q, 5); }
                }
            }
        }
    }
    '''
    ret = np.zeros_like( q )
    shape0 = ret.shape[0]
    shape1 = ret.shape[1]
    weave.inline( code, ['ret', 'q', 'shape0', 'shape1'], verbose=1, compiler = 'gcc', extra_compile_args=['-O3'] )
    return ret.reshape((shape0,shape1))

def qspline_prime_unormalized_weave( q ):
    code = '''
    for( int i=0; i < shape0; ++i ){
        for( int j=0; j < shape1; ++j ){
            int index = i + j*shape0;
            double Q = q[index];
            if( Q > 3 ){ ret[index] = 0; }
            else {
                if( Q > 2 ){ ret[index] = -5*pow(3-Q, 4); }
                else { 
                    if( Q > 1 ){ ret[index] = -5*pow(3-Q, 5) + 5*6*pow(2-Q, 4); }
                    else{ ret[index] = -5*pow(3-Q, 4) + 5*6*pow(2-Q, 4) - 5*15*(1-Q, 4); }
                }
            }
        }
    }
    '''
    ret = np.zeros_like( q )
    shape0 = ret.shape[0]
    shape1 = ret.shape[1]
    weave.inline( code, ['ret', 'q', 'shape0', 'shape1'], verbose=1, compiler = 'gcc', extra_compile_args=['-O3'] )
    return ret.reshape((shape0,shape1))

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

def liq_unormalized_weave( q ):
    code = '''
    double A = -1.458;
    double B = 3.790;
    double C = -2.624;
    double D = -0.2915;
    double E = 0.5831; //changed here
    double F = 0.6500;   
    for( int i=0; i < shape0; ++i ){
        for( int j=0; j < shape1; ++j ){
            int index = i + j*shape0;
            double u = q[index]/2.;
            if( u >= 1 ){ ret[index] = 0; }
            else {
                if( u >= 0.3 ){ ret[index] = A*pow(u,4) + B*pow(u,3) + C*pow(u,2) + D*u + E; }
                else { 
                    ret[index] = F - u;
                }
            }
        }
    }
    '''
    ret = np.zeros_like( q )
    shape0 = ret.shape[0]
    shape1 = ret.shape[1]
    weave.inline( code, ['ret', 'q', 'shape0', 'shape1'], verbose=1, compiler = 'gcc', extra_compile_args=['-O3'] )
    return ret.reshape((shape0,shape1))

def generate_liq( dim, h = None ):
    k = (1., 2.962, 3.947)[dim-1] / h**dim
    return lambda dist, h = h, k = k: k*liq_unormalized_weave( dist/h )

def liq_prime_unormalized_weave( q ):
    code = '''
    double A = -1.458;
    double B = 3.790;
    double C = -2.624;
    double D = -0.2915;
    double E = 0.5831;
    double F = 0.6500;   
    for( int i=0; i < shape0; ++i ){
        for( int j=0; j < shape1; ++j ){
            int index = i + j*shape0;
            double u = q[index]/2;
            if( u > 1 ){ ret[index] = 0; }
            else {
                if( u > 0.3 ){ ret[index] = 2*A*pow(u,3) + 1.5*B*pow(u,2) + C*u + D/2; }
                else { 
                    ret[index] = - .5;
                }
            }
        }
    }
    '''
    ret = np.zeros_like( q )
    shape0 = ret.shape[0]
    shape1 = ret.shape[1]
    weave.inline( code, ['ret', 'q', 'shape0', 'shape1'], verbose=1, compiler = 'gcc', extra_compile_args=['-O3'] )
    return ret.reshape((shape0,shape1))

def generate_liq_prime( dim, h = None ):
    k = (1., 2.962, 3.947)[dim-1] / h**(dim+1)
    return lambda dist, h = h, k = k: k*liq_prime_unormalized_weave( dist/h )

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

