# coding: utf-8
from __future__ import print_function

import numpy as np
import scipy.integrate
from math import sqrt, factorial
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time
import os, sys
import traceback

import summation
import kernel

#def SL2(fp,fv):
    #return lambda p, q, dt: ( lambda dp1: (
                             #lambda dq1: ( dp1/2 + dt/2*fp( p+dp1/2, q+dq1 ), dq1 ) 
                             #) ( dt * fv( p+dp1/2, q ) ) #dq1
                             #) ( dt * fp(p,q) ) #dp1

#def SEuler(fp,fv):
    #return lambda p, q, dt: ( lambda dp1: (
                             #lambda dq1: ( dp1, dq1 ) 
                             #) ( dt * fv( p+dp1, q ) ) #dq1
                             #) ( dt * fp(p,q) ) #dp1

class TimeIntegrators:
    @staticmethod
    def euler(f):
        return lambda t, y, dt: y + dt * f(t,y)

    @staticmethod
    def midpoint(f):
        return lambda t, y, dt: y + (lambda dy1: dt * f( t + dt/2, y + dy1/2 ) )( dt * f( t, y ) )

    @staticmethod
    def RK4(f):
        def integrator( t, y, dt ):
            dy1 = dt * f( t, y )
            dy2 = dt * f( t + dt/2., y + dy1/2. )
            dy3 = dt * f( t + dt/2., y + dy2/2. )
            dy4 = dt * f( t + dt, y + dy3)
            return t + dt, y + (dy1 + 2.*dy2 + 2.*dy3 + dy4)/6.
        return integrator

    class Sympletic:
        @staticmethod
        def Euler( Q, P ):
            def integrator( t, q, p, dt ):
                dp = dt * P( t, q, p )
                dq = dt * Q( t, q, p + dp )
                return t + dt, q + dq, p + dp
            return integrator

        @staticmethod
        def Euler2( fp, fv ):
            def integrator( p, q, dt ):
                p_1 = p + dt * fp( p, q )
                q_1 = q + dt * fv( p_1, q )
                return p_1, q_1
            return integrator

        @staticmethod
        def LeapFrog( fp, fv ):
            def integrator( p, q, dt ):
                dp = dt/2 * fp(p,q)
                dq = dt * fv( p + dp, q )
                dp = dt/2 * fp( p + dp, q + dq )
                return dp, dq
            return integrator
            #return lambda p, q, dt: ( lambda p1: (
                             #lambda q2: ( p1 + dt/2*fp(p1,q2) - p, q2 - q ) 
                             #) ( q + dt * fv(p1,q) ) #q2
                             #) ( p + dt/2 * fp(p,q) ) #p1
        @staticmethod
        def LeapFrog2( fp, fv ):
            def integrator( p, q, dt ):
                p_1 = p + dt/2 * fp(p,q)
                q_2 = q + dt * fv( p_1, q )
                p_2 = p_1 + dt/2 * fp( p_1, q_2 )
                return p_2, q_2
            return integrator
        
        @staticmethod
        def RungeKutta4( dot_q, dot_p ):
            def integrator( t, q, p, dt ):
                P1 = dot_p( t, q, p )
                Q1 = dot_q( t, q, p )
                
                P2 = dot_p( t + dt/2., q + dt/2.*Q1, p + dt/2.*P1 )
                Q2 = dot_q( t + dt/2., q + dt/2.*Q1, p + dt/2.*P1 )
                
                P3 = dot_p( t + dt/2., q + dt/2.*Q2, p + dt/2.*P2 )
                Q3 = dot_q( t + dt/2., q + dt/2.*Q2, p + dt/2.*P2 )
                
                P4 = dot_p( t + dt, q + dt*Q3, p + dt*P3 )
                Q4 = dot_q( t + dt, q + dt*Q3, p + dt*P3 )
                
                dp = (P1 + 2*P2 + 2*P3 + P4)*dt/6.
                dq = (Q1 + 2*Q2 + 2*Q3 + Q4)*dt/6.
                
                return t + dt, q + dq, p + dp
            return integrator

        @staticmethod
        def RungeKutta4a( Q, P ):
            def integrator( t, q, p, dt ):
                P1 = P( t, q, p )
                Q1 = Q( t, q, p ) + dt/2.*P1
                
                P2 = P( t + dt/2., q + dt/2.*Q1, p + dt/2.*P1 )
                Q2 = Q( t + dt/2., q + dt/2.*Q1, p + dt/2.*P1 ) + dt/4.*P1
                
                P3 = P( t + dt/2., q + dt/2.*Q2, p + dt/2.*P2 )
                Q3 = Q( t + dt/2., q + dt/2.*Q2, p + dt/2.*P2 ) + dt/4.*P2
                
                P4 = P( t + dt, q + dt*Q3, p + dt*P3 )
                Q4 = Q( t + dt, q + dt*Q3, p + dt*P3 ) + dt/2.*P3
                
                dp = (P1 + 2*P2 + 2*P3 + P4)*dt/6.
                dq = (Q1 + 2*Q2 + 2*Q3 + Q4)*dt/6.
                
                return t + dt, q + dq, p + dp+p
            return integrator

        @staticmethod
        def RungeKutta4Single( fv ):
            def integrator( q, dt ):
                dq_1 = dt * fv( q )
                dq_2 = dt * fv( q + dq_1/2 )
                dq_3 = dt * fv( q + dq_2/2 )
                dq_4 = dt * fv( q + dq_3 )
                q_1 = q + (dq_1 + 2.*dq_2 + 2.*dq_3 + dq_4)/6.
                return q_1
            return integrator

        @staticmethod
        def Yoshida4( fp, fv ):
            def integrator( p, q, dt ):
                d = [1.351207191959657, -1.702414383919315]
                p_1, q_1 = TimeIntegrators.Sympletic.LeapFrog2( fp, fv )( p, q, dt*d[0] )
                p_2, q_2 = TimeIntegrators.Sympletic.LeapFrog2( fp, fv )( p_1, q_1, dt*d[1] )
                p_3, q_3 = TimeIntegrators.Sympletic.LeapFrog2( fp, fv )( p_2, q_2, dt*d[0] )
                return p_3, q_3
            return integrator

        @staticmethod
        def Yoshida8( fp, fv ):
            def integrator( p, q, dt ):
                d = [0.104242620869991e1, 0.182020630970714e1, 0.157739928123617e0, 0.244002732616735e1, -0.716989419708120e-2, -0.244699182370524e1, -0.161582374150097e1, -0.17808286265894516e1]
                for i in range(6):
                    p, q = TimeIntegrators.Sympletic.LeapFrog2( fp, fv )( p, q, dt*d[i] )
                p, q = TimeIntegrators.Sympletic.LeapFrog2( fp, fv )( p, q, dt*d[7] )
                for i in range(6):
                    p, q = TimeIntegrators.Sympletic.LeapFrog2( fp, fv )( p, q, dt*d[6-i] )
                return p, q
            return integrator

_i=np.s_[:,None]
_j=np.s_[None,:]

class AnalyticalSolution:
    @staticmethod
    def sum_converge_array( f, n, eps=1e-40 ):
        F = f(n)
        while 1:
            n += 1
            fn = f(n)
            F += fn
            if np.all( np.abs(fn) <= np.abs(eps*F) ):
                return F

    def __init__(self, visc, u_0 = 1, L = 1 ):
        self.visc = visc
        self.L = L
        self.u_0 = u_0
        self.k = np.pi / self.L
        self.K = 4 * self.k * self.visc
        self.I0 = scipy.special.i0( u_0*L/(2*np.pi*visc) )
        
    def Iv(self,n): return scipy.special.iv( n, self.u_0*self.L/(2*np.pi*self.visc) )
    
    def t_term(self, n, t ): return np.exp( -n**2 * self.k**2 * self.visc*t )
    
    def x_N(self, n, x, diff = 0 ): return np.imag( np.exp( 1.j * n*self.k*x ) * (1.j * n*self.k)**diff )

    def x_D(self, n, x, diff = 0 ): return np.real( np.exp( 1.j * n*self.k*x ) * (1.j * n*self.k)**diff )
    
    def N(self, x, t, diff_x = 0 ):
        return AnalyticalSolution.sum_converge_array(
                lambda n, x=x, t=t, diff_x=diff_x: n * self.Iv(n) * self.t_term( n, t ) * self.x_N( n, x, diff_x ), 
                n=1
            )

    def D(self, x, t, diff_x = 0 ):
        return (self.I0 if diff_x == 0 else 0) +\
            2*AnalyticalSolution.sum_converge_array(
                lambda n, x=x, t=t, diff_x=diff_x: self.Iv(n) * self.t_term( n, t ) * self.x_D( n, x, diff_x ), 
                n=1
            )

    def u(self, t, x ):
        return self.K * self.N(x,t)/self.D(x,t)
    
    def du_dx(self, t, x ):
        '''
        d(N/D) = dN/D - N/D**2*dD 
               = (dN*D - N*dD)/D**2
        '''
        return self.K * ( 
            self.N(x,t,1)*self.D(x,t) 
            - self.N(x,t)*self.D(x,t,1) 
            )/self.D(x,t)**2

    def d2u_dx2(self, x, t ):
        '''
        d2(N/D) = d(dN*D - N*dD)/D**2 + 2(dN*D - N*dD)/D**3*dD
                = (d2N*D + dN*dD - dN*dD - N*d2D)/D**2 + 2(dN*D - N*dD)/D**3*dD
                = (d2N*D**2 + dN*dD*D - dN*dD*D - N*d2D*D + 2dN*D*dD - 2N*dD*dD)/D**3
                = (d2N*D**2 + dN*dD*D - dN*dD*D - N*d2D*D + 2dN*D*dD - 2N*dD*dD)/D**3
        '''
        return self.K* (
            self.N(x,t,2) * self.D(x,t)**2
            - 2*self.D(x,t) * self.N(x,t,1)*self.D(x,t,1) 
            - self.N(x,t) * self.D(x,t) * self.D(x,t,2)
            + 2*self.N(x,t)*self.D(x,t,1)**2
            )/self.D(x,t)**3
        
    def d3u_dx3(self, x, t ):
        return self.K* (
            self.N(x,t,3) * self.D(x,t)**3
            - 3*self.D(x,t)**2*self.N(x,t,2)*self.D(x,t,1) 
            - 3*self.D(x,t)**2*self.N(x,t,1)*self.D(x,t,2) 
            + 6*self.D(x,t)*self.N(x,t,1)*self.D(x,t,1)**2
            - self.N(x,t)*self.D(x,t,3)*self.D(x,t)**2
            - 6*self.N(x,t)*self.D(x,t,1)**3
            + 6*self.N(x,t)*self.D(x,t)*self.D(x,t,1)*self.D(x,t,2)
            )/self.D(x,t)**4

            
def diff_array( a, x ):
    return (a[2:] - a[:-2])/(x[2:] - x[:-2])

def calculate_h( h, q, N, status, dim = 1, eps = 1e-1 ):
    hnew = np.array(h)
    NN = np.ones_like(q)*N
    while 1:
        prevh = np.array(hnew)
        hnew = h0*NN/rho(q,hnew,0,q)**(1./dim)
        dist_h = np.abs( q[_i] - q[_j] )/hnew[_i]
        dist_h[ dist_h <= 2 ] = 1
        dist_h[ dist_h > 2 ] = 0
        count = np.sum( dist_h, axis = 0 )
        #print( 'count', count[ status_border == 0 ], len(count), np.all( count[ status_border == 0 ] == N ) )
        if np.all( count[ status == 0 ] == N ): break
        maxdiff = np.max( np.abs(prevh-hnew) )
        print( 'diff', np.max( np.abs(prevh-hnew) ) )
        prevh = np.array(hnew)
        if maxdiff/h0 < 1e-2: break
        break
    return hnew

def calculate_h_dist( h, q, p, dt, N, dim = 1, eps = 1e-1 ):
    fraction = .1
    q1 = q+p*dt
    dist_h = np.abs( q1[_i]-q1[_j] )/h[_i]
    dist_h[dist_h < .1*fraction] = np.max(dist_h)
    #dist_h_min = np.min( dist_h, axis = 1 )
    dist_h_min = np.min( dist_h, axis = 0 )
    
    hnew = np.array(h)
    #hnew[dist_h_min < fraction] *= dist_h_min[dist_h_min < fraction]/fraction
    hnew *= dist_h_min/fraction
    return hnew
    
import timeit

colors = ['r','g','b','m','y','c']

def generate_border( t, r, s, ds, h, method_smooth = 'self', method_sum = 'np_sum' ):
    
    dists = None
    if method_smooth == 'self':
        H = h[_j]
    elif method_smooth == 'other':
        H = h[_i]
    elif method_smooth == 'mean':
        H = .5*(h[_i]+h[_j])

    dists = (r[_i] - r[_j])/H
    dists[np.abs(dists)>3] = 0
    #inds = border_indices(r, h, interpolator )
    numer = summation.sum_2( dists, axis=0, method = method_sum )
    denom = summation.sum_2( np.abs(dists), axis=0, method = method_sum )
    #print(numer)
    #print(denom)
    uni = numer/denom
    inds = np.argwhere( np.abs(uni) == 1 )
    print('border inds',inds.flatten(), uni[inds].flatten() )
    newt = []
    newr = []
    news = []
    newds = []
    newh = []

    for ind in inds:
        neighbors = np.argwhere( np.abs(r-r[ind])<3*h[ind] )
        dist = r[neighbors] - r[ind]
        mask = np.all( [np.abs(dist)>0, dist*uni[ind] > 0 ], axis=0 )
        newr = np.append( newr, r[ind] - dist[mask] )
        
        delta = s[neighbors] - s[ind]
        news = np.append( news, s[ind] - delta[mask] )

        delta = ds[neighbors] - ds[ind]
        newds = np.append( news, ds[ind] + delta[mask] )

        deltah = h[neighbors] - h[ind]
        newh = np.append( newh, h[ind] + deltah[mask] )

        deltat = t[neighbors] - t[ind]
        newt = np.append( newt, t[ind] + deltat[mask] )

    #print('newr',newr)
    return newt, newr, news, newds, newh, -1*np.ones_like(newr)
    

def merge( r, p, h, status, fuse=False ):
    fraction = .2
    inds = np.argwhere( h[status==0] < fraction*h0 )
    print( 'h<h0', inds.flatten() , end='; ')
    return r, p, h, status
    
    for ind in inds:
        print( 'h/h0', h[ind]/h0, 'id', ind )
    
    neighbors = np.argwhere( np.abs(r-r[ind]) < .5*fraction*h0 )
    print( 'neighbors', neighbors.flatten() )
    
    if h[ind] > fraction*h0:
        return r, p, h, status

    if fuse:
        print( 'fuse!!!' )
        r[neighbors] = np.mean( r[neighbors] )
        p[neighbors] = np.mean( p[neighbors] )
        h[neighbors] = np.mean( h[neighbors] )
        return r, p, h, status
    
    print( 'before', r[ind], p[ind], h[ind] )
    r[ind] = np.mean( r[neighbors] )
    #p[ind] = np.mean( p[neighbors] )
    #h[ind] = np.mean( h[neighbors] )
    print( 'after', r[ind], p[ind], h[ind] )
    
    print( 'merge!!!' )
    print( 'before', r.shape )
    r = np.delete(r, remove)
    p = np.delete(p, remove)
    h = np.delete(h, remove)

    status = np.delete(status, remove)
    print( 'after', r.shape )

    return r, p, h, status

class BaseSPH:
    def __init__( self, h0, W, density, phi='ref' ):
        self.W = kernel.derivatives( W )
        self.h0 = h0
        self.density = density
        if phi == 'ref':
            self.phi = self.ref
        if phi == 'id':
            self.phi = self.id
    
    def interpolation( self, f, r, D=0, function=False ):
        W = self.W[D]
        def inner( x ):
            if type(f) is float:
                val = W( x[_j], r[_i], self.h0 )*f
            elif len(f.shape) == 1:
                val = W( x[_j], r[_i], self.h0 )*f[_i]
            else:
                val = W( x[_j], r[_i], self.h0 )*f
            
            valp = np.zeros_like(val)
            valp[val > 0] = val[val>0]
            valn = np.zeros_like(val)
            valn[val < 0] = val[val<0]
            return np.sum( valp, axis=0 ) + np.sum( valn, axis=0 )
        if function:
            return inner
        else:
            return inner(r)

    def id( self, x, hx = None ):
        return x
    
    def delta( self, r, D=0, power=1, function=False ):
        def inner( x ):
            if power == 0: diff = 1
            else: diff = r[_i] - x[_j]
            if power < 0: diff[ diff == 0 ] = self.h0
            return self.interpolation( diff**power/self.phi(r)[_i], x, D=D, function=function )
        if function:
            return inner
        else:
            return inner(r)

    def deltaAbs( self, x, D=0, power=1 ):
        if power == 0: diff = 1
        else: diff = x[_i] - x[_j]
        return self.interpolation( np.abs(diff)**power/self.phi(x)[_i], x, D=D )

    def fdelta( self, f, x, D=0, power=1 ):
        if f is None: f = np.ones_like(x)
        if power == 0: diff = 1
        else: diff = x[_i] - x[_j]
        if power < 0: diff[ diff == 0] = self.h0
        return self.interpolation( diff**power * f[_i]/self.phi(x)[_i], x, D=D )
        
    def ref( self, x, r=None, D=0, function=False ):
        return self.interpolation( 1., x, D=D, function=function )
    
    def sph( self, f, x, D=0, function=False ):
        return self.interpolation( f/self.phi(x), x, D=D, function=function )
    
    def uni( self, x ):
        return self.delta( x )/self.deltaAbs( x )
    
class StandardSPH:
    def __init__( self, sph ):
        self.sph = sph
    
    def diff( self, f, r ):
        return self.sph.sph( f, r, D=1 )

    def diff2( self, f, r ):
        return self.sph.sph( f, r, D=2 )
        
    def interp( self, f, r ):
        return self.sph.sph( f, r, D=0 )
    
class KernelGradientFree:
    def compute_invM( self, r ):
        m = np.array( [[ self.sph.delta( r, power=i+j+self.basepower, D=0 ) for j in range(self.order)] for i in range(self.order) ] ).T
        self.det = np.abs( np.linalg.det( m ) )
        self.invM = np.zeros_like( m )
        self.invM[self.det > 0] = np.linalg.inv( m[self.det > 0] )
        return self.invM
    
    def compute_V( self, f, r ):
        self.V = np.array([ self.sph.fdelta( f, r, power=i+self.basepower, D=0 ) for i in range(self.order) ] ).T
        return self.V

    def compute_F( self, f, r ):
        self.F = np.einsum( '...ij,...i->...j', self.compute_invM( r ), self.compute_V( f, r ) ).T
        return self.F

    def compute_F2( self, f, r ):
        m = np.array( [[ self.sph.delta( r, power=i+j+self.basepower, D=0 ) for j in range(self.order)] for i in range(self.order) ] ).T
        self.det = np.abs( np.linalg.det( m ) )
        self.invM = np.zeros_like( m )
        self.invM[self.det > 0] = np.linalg.inv( m[self.det > 0] )

        self.V = np.array([ self.sph.fdelta( f, r, power=i+self.basepower, D=0 ) for i in range(self.order) ]).T
        self.F = np.zeros_like( self.V )
        self.F[ self.det>0 ] = np.einsum( '...ij,...i->...j', self.invM[ self.det>0 ], self.V[ self.det>0 ] )
        
        if np.sum( (self.det==0).astype(int) ) > 0:
            print( 'det==0', np.sum( (self.det==0).astype(int) ) )
            zeros = np.zeros_like( f[ self.det==0 ] )
            self.F[ self.det==0 ] = np.array([ f[ self.det==0 ] if i==0 else zeros for i in range(self.order) ]).T

        return self.F.T
    
    def __init__( self, order = 4, basepower = 0, sph = None ):
        self.order = order
        self.basepower = basepower
        self.sph = sph
        self.r = 0
        self.f = 0
        self.invM = None
        self.V = None
        self.F = None
    
    def diff( self, f, r ):
        return self.compute_F2( f, r )[1]
    
    def interp( self, f, r, x ):
         return self.compute_F2( f, r )[0]

class Output:
    def makeName( self, s ):
        return 'plots/sph_%svst_%s'%(s, self.name)
        
    def __init__( self, name ):
        self.name = name
        self.doneHeader = False
        
    def header( self, N ):
        for s in [ 'r', 'u' ]:
            f = open( self.makeName(s), 'w' )
            f.write("# t " + ' '.join( map(lambda x: '%s[%d]'%(s,x), range(N) ) ) + '\n' )
            f.close()
        self.doneHeader = True
    
    def write( self, t, r, u ):
        if not self.doneHeader:
            self.header( len(r) )
        for s, arr in [ ('r', r), ('u', u) ]:
            f = open( self.makeName(s), 'a' )
            f.write("%f "%t + ' '.join( map(lambda x: '%f'%x, arr ) ) + '\n' )
            f.close()

    def read( self ):
        t = []
        r = []
        u = []
        for s in [ 'r', 'u' ]:
            data = open( self.makeName(s), 'r' ).readlines()
            for line in data:
                if line[0] == '#': continue
                linesplit = line.replace('  ', ' ').strip(' \n').split(' ')
                if s == 'r':
                    r += [ map( float, linesplit[1:] ) ]
                    t += [ float(linesplit[0]) ]
                elif s == 'u':
                    u += [ map( float, linesplit[1:] ) ]
        t = np.array(t)
        r = np.array(r)
        u = np.array(u)
        return t, r, u
        
        for s, arr in [ ('r', r), ('u', p) ]:
            fname = 'plots/sph_%svst_%s'%(s, self.name)
            f = open( fname, 'a' )
            f.write("%f "%t + ' '.join( map(lambda x: '%f'%x, arr ) ) + '\n' )
            f.close()

class Burger:
    def __init__(self, viscosityCoefficient = None, diff = None ):
        self.viscosityCoefficient = viscosityCoefficient
        self.diff = diff
    
    def set_viscosityCoefficient( self, v ):
        self.viscosityCoefficient = v
        return self
    
    def dot_r( self ):
        return lambda t, r, u: u #self.interp( u, r ) #u
    
    def dot_u( self ):
        return lambda t, r, u: self.viscosityCoefficient * self.diff( self.diff(u,r), r )
    
    def interp( self ):
        return lambda t, u, r, x: self.interpolator( u, r, x )
            
    def set_diff( self, method ):
        self.interpolator = method.interp
        self.diff = method.diff
        return self
    
class LagrangianSystem:
    def __init__( self ):
        self.dot_r = None
        self.dot_u = None
        
    def set_timeStep( self, dt0 ):
        self.timeStep = dt0
        return self
    
    def set_endTime( self, tf ):
        self.endTime = tf
        return self
    
    def set_outputFile( self, out ):
        self.outputFile = Output
        return self
    
    def set_equations( self, eqs ):
        self.dot_r = eqs.dot_r()
        self.dot_u = eqs.dot_u()
        return self

    def set_interpolator( self, interpolator ):
        self.interp = interpolator.interp
        return self
    
    def set_integrator( self, integrator ):
        self.timeIntegrator = integrator
        return self

    def set_initialCondition( self, r0, u0 ):
        self.t = 0
        self.r = r0
        self.u = u0( self.r )
        return self
        
    def boundaryCondition( self, t, r, u ):
        u[0] = self.u[0]
        u[-1] = self.u[-1]
        r[0] = self.r[0]
        r[-1] = self.r[-1]
        return r, u
    
    def addBoundary( self, t, r, u ):
        br1 = 2*r[0] - r[:5]
        br2 = 2*r[-1] - r[-5:]
        
        return r, u
    
    def removeBoundary( self, t, r, u ):
        return r, u

    def integrate( self, fname ):
        
        t, r, u = float(self.t), self.r.copy(), self.u.copy()

        print( 'r', len(r) )
        print( 'u', len(u) )
        
        stream = Output( fname )
        stream.write(t,r,u)

        loopnumber = 0
        
        while 1:
            loopnumber += 1
            print( '%d'%loopnumber, end=' ' )
            t0 = time.time()
            dt = self.timeStep
            t, r, u = self.timeIntegrator( self.dot_r, self.dot_u )( t, r, u, dt )
            r, u = self.boundaryCondition( t, r, u )
            r, u = self.fixCrossing( r, u )

            stream.write( t, r, u )
            elapsed = time.time() - t0
            print( 'clock %.2e (%.2e/N) t=%s'%( elapsed, elapsed/len(r), t ), end='; ' )

            crossed = np.argwhere( r[1:] < r[:-1] )
            if len(crossed) > 0:
                print( '(crossed %d)'%len(crossed), end='; ' )
                print( crossed )
                exit(0)

            if np.min(t) >= self.endTime:
                print( 'reached tf; simulation ended' )
                break
            plotOutputs( fname, fname.replace( 'dat', 'pdf' ) )
            print('')
        print( 'finished' )
        plotOutputs( fname, fname.replace( 'dat', 'pdf' ) )
        return
    
    def fixCrossing( self, r, u ):
        ufix = u.copy()
        n = 0
        while 1:
            n+=1
            rnext = r + ufix*self.timeStep*2
            frozenCluster = np.zeros_like(r)
            
            crossed = np.argwhere( rnext[1:] < rnext[:-1] ).T[0]
            if len(crossed) == 0:
                return r, ufix
            print( 'crossed', len(crossed) )
            for crossPair in crossed:
                for cross in [crossPair, crossPair+1]:
                    newCluster = int(np.max( frozenCluster ))+1
                    if frozenCluster[cross-1] == 0 and frozenCluster[cross+1] == 0:
                        frozenCluster[cross] = newCluster
                        
                    elif frozenCluster[cross-1] == 0:
                        frozenCluster[cross] = frozenCluster[cross+1]
                        
                    elif frozenCluster[cross+1] == 0:
                        frozenCluster[cross] = frozenCluster[cross-1]
                    
                    elif frozenCluster[cross-1] < frozenCluster[cross+1]:
                        frozenCluster[cross] = frozenCluster[cross-1]
                    
                    else:
                        frozenCluster[cross] = frozenCluster[cross+1]
            
            for i in range(1, int(np.max( frozenCluster ))+1 ):
                args = np.argwhere( frozenCluster == i )
                if len(args) == 0: continue
                print('freeze', i, args.T )
                uCM = np.sum( u[args] )/len(args)
                ufix[ args ] = uCM + ( ufix[args] - uCM )*.9**n

def append( d, e, v ):
    if e in d:
        d[e].append(v)
    else:
        d[e] = [v]
    return 

def readFile( fname ):
    data = open( fname, 'r' ).readlines()
    t = []
    r = []
    for line in data:
        if line[0] == '#': continue
        linesplit = line.replace('  ', ' ').strip(' \n').split(' ')
        t += [ float(linesplit[0]) ]
        r += [ map( float, linesplit[1:] ) ]
    t = np.array(t)
    r = np.array(r)
    print(r.shape)
    return t, r
    
def plotOutputs( fname, outname ):
    plt.rcParams['font.size'] = 12.
    plt.rcParams['font.family'] = "serif"
    plt.rcParams["xtick.labelsize"] = 'xx-small'
    plt.rcParams["ytick.labelsize"] = 'xx-small'
    width = plt.rcParams['figure.figsize'][0]/2

    fig = plt.figure( figsize = (width,width) )
    axr = fig.add_subplot(211)
    axu = fig.add_subplot(212)

    stream = Output( fname )
    t, r, u = stream.read()

    indices = np.argwhere( np.abs(r[0,:] - 1) < 1.1 )
    tmax = 1.
    for i in indices:
        axr.plot( r[:,i], t, '-', lw=.5 )
        axu.plot( u[:,i], t, '-', lw=.5 )
        
    axr.set_ylabel('t')
    axr.set_xlabel( r'r(t)' )
    axu.set_ylabel('t')
    axu.set_xlabel( r'u(t)' )

    fig.tight_layout()
    fig.savefig( outname )
    plt.close()

if __name__ == '__main__':
    print(sys.argv)
    if len(sys.argv) > 1:
        plotOutputs( sys.argv[-1] )
        exit(0)
    
    L = 1.
    u_0 = 1.
    N = 251

    def u0( x, u_0, L ):
        return u_0*np.sin(np.pi*x/L)

    v = .05
    system = LagrangianSystem()
    system.set_timeStep( .01 )
    system.set_endTime( 3. )
    system.set_initialCondition( r0 = np.linspace( 0., 2.*L, N ), u0 = lambda x: u0( x, u_0, L ) )
    
    if False:
        anaSol = AnalyticalSolution( v )
        def AnalyticalIntegrator( dot_r, dot_u ):
            def integrator( t, r, u, dt ):
                t, r = TimeIntegrators.RK4( anaSol.u )( t, r, dt )
                return t, r, anaSol.u( t, r )
            return integrator
        system.set_integrator( AnalyticalIntegrator )
        system.integrate( 'analyticalv%s.dat'%v )

    h0 = .05
    burger = Burger()
    burger.set_viscosityCoefficient( v )
    burger.set_diff( KernelGradientFree( sph = BaseSPH( h0 = h0, W = 'qspline', density = N/L*h0 ), order = 3, basepower = 0 ) )
    system.set_equations( burger )
    system.set_integrator( TimeIntegrators.Sympletic.RungeKutta4 )
    system.integrate( 'simulationv%s.dat'%v )
    
    exit(0)

exit(0)
