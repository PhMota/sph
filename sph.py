# coding: utf-8

from __future__ import print_function

from numpy import *
import scipy.integrate
import scipy.stats as stats

from math import factorial
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time
import os, sys
import traceback
import glob

import summation
import kernel
import threading

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib

from Timer import Timer
import weave

from scipy.sparse import csr_matrix

def binned_statistic(x, values, func, nbins, range):
    '''The usage is nearly the same as scipy.stats.binned_statistic'''

    r0, r1 = range
    digitized = (float(nbins)/(r1 - r0)*(x - r0)).astype(int)
    values = values[digitized < nbins].astype(float)
    digitized = digitized[digitized < nbins]

    N = len(values)
    S = csr_matrix((values, [digitized, arange(N)]), shape=(nbins, N))

    return [func(group) for group in split(S.data, S.indptr[1:-1])]

# def binned_statistic_fast(x, values, func, nbins, range):
#     '''The usage is nearly the same as scipy.stats.binned_statistic'''
#
#     N = len(values)
#     r0, r1 = range
#
#     digitized = (float(nbins)/(r1 - r0)*(x - r0)).astype(int)
#     S = csr_matrix((values, [digitized, arange(N)]), shape=(nbins, N))
#
#     return [func(group) for group in split(S.data, S.indptr[1:-1])]

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

_i=s_[:,None]
_j=s_[None,:]

class AnalyticalSolution:
    @staticmethod
    def sum_converge_array( f, n, eps=1e-40 ):
        F = f(n)
        while 1:
            n += 1
            fn = f(n)
            F += fn
            if all( abs(fn) <= abs(eps*F) ):
                return F

    def __init__(self, visc, u_0 = 1, L = 1 ):
        self.visc = visc
        self.L = L
        self.u_0 = u_0
        self.k = pi / self.L
        self.K = 4 * self.k * self.visc
        self.I0 = scipy.special.i0( u_0*L/(2*pi*visc) )

    def Iv(self,n): return scipy.special.iv( n, self.u_0*self.L/(2*pi*self.visc) )

    def t_term(self, n, t ): return exp( -n**2 * self.k**2 * self.visc*t )

    def x_N(self, n, x, diff = 0 ): return imag( exp( 1.j * n*self.k*x ) * (1.j * n*self.k)**diff )

    def x_D(self, n, x, diff = 0 ): return real( exp( 1.j * n*self.k*x ) * (1.j * n*self.k)**diff )

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
    hnew = array(h)
    NN = ones_like(q)*N
    while 1:
        prevh = array(hnew)
        hnew = h0*NN/rho(q,hnew,0,q)**(1./dim)
        dist_h = abs( q[_i] - q[_j] )/hnew[_i]
        dist_h[ dist_h <= 2 ] = 1
        dist_h[ dist_h > 2 ] = 0
        count = sum( dist_h, axis = 0 )
        #print( 'count', count[ status_border == 0 ], len(count), all( count[ status_border == 0 ] == N ) )
        if all( count[ status == 0 ] == N ): break
        maxdiff = max( abs(prevh-hnew) )
        print( 'diff', max( abs(prevh-hnew) ) )
        prevh = array(hnew)
        if maxdiff/h0 < 1e-2: break
        break
    return hnew

def calculate_h_dist( h, q, p, dt, N, dim = 1, eps = 1e-1 ):
    fraction = .1
    q1 = q+p*dt
    dist_h = abs( q1[_i]-q1[_j] )/h[_i]
    dist_h[dist_h < .1*fraction] = max(dist_h)
    #dist_h_min = min( dist_h, axis = 1 )
    dist_h_min = min( dist_h, axis = 0 )

    hnew = array(h)
    #hnew[dist_h_min < fraction] *= dist_h_min[dist_h_min < fraction]/fraction
    hnew *= dist_h_min/fraction
    return hnew

import timeit

colors = ['r','g','b','m','y','c']

def generate_border( t, r, s, ds, h, method_smooth = 'self', method_sum = 'sum' ):

    dists = None
    if method_smooth == 'self':
        H = h[_j]
    elif method_smooth == 'other':
        H = h[_i]
    elif method_smooth == 'mean':
        H = .5*(h[_i]+h[_j])

    dists = (r[_i] - r[_j])/H
    dists[abs(dists)>3] = 0
    #inds = border_indices(r, h, interpolator )
    numer = summation.sum_2( dists, axis=0, method = method_sum )
    denom = summation.sum_2( abs(dists), axis=0, method = method_sum )
    #print(numer)
    #print(denom)
    uni = numer/denom
    inds = argwhere( abs(uni) == 1 )
    print('border inds',inds.flatten(), uni[inds].flatten() )
    newt = []
    newr = []
    news = []
    newds = []
    newh = []

    for ind in inds:
        neighbors = argwhere( abs(r-r[ind])<3*h[ind] )
        dist = r[neighbors] - r[ind]
        mask = all( [abs(dist)>0, dist*uni[ind] > 0 ], axis=0 )
        newr = append( newr, r[ind] - dist[mask] )

        delta = s[neighbors] - s[ind]
        news = append( news, s[ind] - delta[mask] )

        delta = ds[neighbors] - ds[ind]
        newds = append( news, ds[ind] + delta[mask] )

        deltah = h[neighbors] - h[ind]
        newh = append( newh, h[ind] + deltah[mask] )

        deltat = t[neighbors] - t[ind]
        newt = append( newt, t[ind] + deltat[mask] )

    #print('newr',newr)
    return newt, newr, news, newds, newh, -1*ones_like(newr)


def merge( r, p, h, status, fuse=False ):
    fraction = .2
    inds = argwhere( h[status==0] < fraction*h0 )
    print( 'h<h0', inds.flatten() , end='; ')
    return r, p, h, status

    for ind in inds:
        print( 'h/h0', h[ind]/h0, 'id', ind )

    neighbors = argwhere( abs(r-r[ind]) < .5*fraction*h0 )
    print( 'neighbors', neighbors.flatten() )

    if h[ind] > fraction*h0:
        return r, p, h, status

    if fuse:
        print( 'fuse!!!' )
        r[neighbors] = mean( r[neighbors] )
        p[neighbors] = mean( p[neighbors] )
        h[neighbors] = mean( h[neighbors] )
        return r, p, h, status

    print( 'before', r[ind], p[ind], h[ind] )
    r[ind] = mean( r[neighbors] )
    #p[ind] = mean( p[neighbors] )
    #h[ind] = mean( h[neighbors] )
    print( 'after', r[ind], p[ind], h[ind] )

    print( 'merge!!!' )

    print( 'before', r.shape )
    r = delete(r, remove)
    p = delete(p, remove)
    h = delete(h, remove)

    status = delete(status, remove)
    print( 'after', r.shape )

    return r, p, h, status

class BaseInterpolator:
    def set_positions( self, r ):
        self.r = r
        self.i = ((r - np.amin(r, axis=1))/self.h).astype(int)
        print( 'max', np.amax(self.i, axis=1) )
        #self.neighbors = np.zeros( () )


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

            valp = zeros_like(val)
            valp[val > 0] = val[val>0]
            valn = zeros_like(val)
            valn[val < 0] = val[val<0]
            return sum( valp, axis=0 ) + sum( valn, axis=0 )
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
        return self.interpolation( abs(diff)**power/self.phi(x)[_i], x, D=D )

    def fdelta( self, f, x, D=0, power=1 ):
        if f is None: f = ones_like(x)
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
        m = array( [[ self.sph.delta( r, power=i+j+self.basepower, D=0 ) for j in range(self.order)] for i in range(self.order) ] ).T
        self.det = abs( linalg.det( m ) )
        self.invM = zeros_like( m )
        self.invM[self.det > 0] = linalg.inv( m[self.det > 0] )
        return self.invM

    def compute_V( self, f, r ):
        self.V = array([ self.sph.fdelta( f, r, power=i+self.basepower, D=0 ) for i in range(self.order) ] ).T
        return self.V

    def compute_F( self, f, r ):
        self.F = einsum( '...ij,...i->...j', self.compute_invM( r ), self.compute_V( f, r ) ).T
        return self.F

    def compute_F2( self, f, r ):
        m = array( [[ self.sph.delta( r, power=i+j+self.basepower, D=0 ) for j in range(self.order)] for i in range(self.order) ] ).T
        self.det = abs( linalg.det( m ) )
        self.invM = zeros_like( m )
        self.invM[self.det > 0] = linalg.inv( m[self.det > 0] )

        self.V = array([ self.sph.fdelta( f, r, power=i+self.basepower, D=0 ) for i in range(self.order) ]).T
        self.F = zeros_like( self.V )
        self.F[ self.det>0 ] = einsum( '...ij,...i->...j', self.invM[ self.det>0 ], self.V[ self.det>0 ] )

        if sum( (self.det==0).astype(int) ) > 0:
            print( 'det==0', sum( (self.det==0).astype(int) ) )
            zeros = zeros_like( f[ self.det==0 ] )
            self.F[ self.det==0 ] = array([ f[ self.det==0 ] if i==0 else zeros for i in range(self.order) ]).T

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

    def interp( self, f, r ):
        return self.compute_F2( f, r )[0]

class Output:
    def makeName( self, s ):
        return self.name.replace('.dat', '.%s.dat'%s )

    def figname( self ):
        return self.name.replace('.dat', '.png' )

    def __del__( self ):
        Output.cleanLocks( os.path.dirname(self.name) )

    def __init__( self, name ):
        self.name = name
        self.doneHeader = False
        if self.islocked():
            print('file is locked')
            exit(0)
        self.lock()

    def lock(self):
        open( self.name+'.lock', 'w')

    def islocked(self):
        return os.path.exists( self.name+'.lock' )

    @staticmethod
    def cleanLocks( dir ):
        print( glob.glob( dir+'/*.lock' ) )
        for lock in glob.glob( dir+'/*.lock' ):
            print('removing lock', lock)
            os.remove( lock )

    def header( self, N ):
        for s in [ 'r', 'u', 'u_smooth' ]:
            f = open( self.makeName(s), 'w' )
            f.write("# t " + ' '.join( map(lambda x: '%s[%d]'%(s,x), range(N) ) ) + '\n' )
            f.close()
        self.doneHeader = True

    def write( self, t, r, u, u_smooth ):
        if not self.doneHeader:
            self.header( len(r) )
        for s, arr in [ ('r', r), ('u', u), ('u_smooth', u_smooth) ]:
            f = open( self.makeName(s), 'a' )
            f.write("%f "%t + ' '.join( map(lambda x: '%s'%x, arr ) ) + '\n' )
            f.close()

    def read( self ):
        t = []
        r = []
        u = []
        u_smooth = []
        for s in [ 'r', 'u', 'u_smooth' ]:
            data = open( self.makeName(s), 'r' ).readlines()
            for line in data:
                if line[0] == '#': continue
                linesplit = line.replace('  ', ' ').strip(' \n').split(' ')
                if s == 'r':
                    r += [ map( float, linesplit[1:] ) ]
                    t += [ float(linesplit[0]) ]
                elif s == 'u':
                    u += [ map( float, linesplit[1:] ) ]
                elif s == 'u_smooth':
                    u_smooth += [ map( float, linesplit[1:] ) ]
        t = array(t)
        r = array(r)
        u = array(u)
        u_smooth = array(u_smooth)
        return t, r, u, u_smooth

        #for s, arr in [ ('r', r), ('u', p) ]:
            #fname = 'plots/sph_%svst_%s'%(s, self.name)
            #f = open( fname, 'a' )
            #f.write("%f "%t + ' '.join( map(lambda x: '%f'%x, arr ) ) + '\n' )
            #f.close()

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
        return lambda t, u, r: self.interpolator( u, r )

    def set_diff( self, method ):
        self.interpolator = method.interp
        self.diff = method.diff
        return self

class Hydro:
    def dot_r(self):
        return lambda t, r, u: u

    def dot_u(self):
        return lambda t, r, u: - self.gradP(u,r)/self.s(u,r)

    def gradP(self):
        return -( self.diff( self.P()/self.s(), r ) + self.P()/self.s()**2 * self.diff( self.s(), r ) )

class LagrangianSystem:
    def __init__( self ):
        self.dot_r = None
        self.dot_u = None
        self.loopnumber = 0

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
        print('r', r0)
        print('r0shape', r0.shape)
        print('u0shape', self.u.shape)
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

    def integrate( self, fname, signal=None ):

        t, r, u = float(self.t), self.r.copy(), self.u.copy()

        self.linklist(r,)
        print( 'r', len(r) )
        print( 'u', len(u) )

        stream = Output( fname )
        stream.write(t,r,u,u)

        loopnumber = 0

        while 1:
            loopnumber += 1
            print( '%d'%loopnumber, end=' ' )
            t0 = time.time()
            dt = self.timeStep
            t, r, u = self.timeIntegrator( self.dot_r, self.dot_u )( t, r, u, dt )
            r, u = self.boundaryCondition( t, r, u )
            r, u = self.fixCrossingPair( r, u )
            u_smooth = self.interp( u, r )

            stream.write( t, r, u, u_smooth )
            elapsed = time.time() - t0
            print( 'clock %.2e (%.2e/N) t=%s'%( elapsed, elapsed/len(r), t ), end='; ' )

            crossed = argwhere( r[1:] < r[:-1] )
            if len(crossed) > 0:
                print( '(crossed %d)'%len(crossed), end='; ' )
                print( 'c', crossed.T )
                #exit(0)

            if np.min(t) >= self.endTime:
                print( 'reached tf; simulation ended' )
                break
            plotOutputs( stream )
            if signal is not None: signal()
            print('')
        print( 'finished' )
        plotOutputs( stream )
        return

    def fixCrossingPair( self, r, u ):
        ufix = u[:]
        n = 0
        while 1:
            n+=1
            rnext = r + ufix*self.timeStep*2
            frozenCluster = zeros_like(r)

            crossed = argwhere( rnext[1:] < rnext[:-1] ).T[0]
            if len(crossed) == 0:
                return r, ufix
            print( 'will cross', len(crossed) )
            for cross in crossed:
                newCluster = int(max( frozenCluster ))+1
                if frozenCluster[cross] == 0:
                    frozenCluster[cross] = newCluster
                    frozenCluster[cross+1] = newCluster

                else:
                    frozenCluster[cross+1] = frozenCluster[cross]

            u_interp = self.interp( ufix, r )
            ufixnew = ufix[:]
            for i in range(1, int(max( frozenCluster ))+1 ):
                args = argwhere( frozenCluster == i )
                if len(args) == 0: continue
                print('freeze', i, args.T )
                ufixnew[ args ] = u_interp[ args ]
            ufix = ufixnew[:]
        return r, ufix

    def fixCrossing( self, r, u ):
        ufix = u.copy()
        n = 0
        while 1:
            n+=1
            rnext = r + ufix*self.timeStep*2
            frozenCluster = zeros_like(r)

            crossed = argwhere( rnext[1:] < rnext[:-1] ).T[0]
            if len(crossed) == 0:
                return r, ufix
            print( 'crossed', len(crossed) )
            for crossPair in crossed:
                for cross in [crossPair, crossPair+1]:
                    newCluster = int(max( frozenCluster ))+1
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

            for i in range(1, int(max( frozenCluster ))+1 ):
                args = argwhere( frozenCluster == i )
                if len(args) == 0: continue
                print('freeze', i, args.T )
                uCM = sum( u[args] )/len(args)
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
    t = array(t)
    r = array(r)
    print('rshape', r.shape)
    return t, r

def plotOutputs( stream ):
    plt.rcParams['font.size'] = 12.
    plt.rcParams['font.family'] = "serif"
    plt.rcParams["xtick.labelsize"] = 'xx-small'
    plt.rcParams["ytick.labelsize"] = 'xx-small'
    width = plt.rcParams['figure.figsize'][0]/2

    fig = plt.figure( figsize = (3*width,width*.75) )
    grid = plt.GridSpec(2, 2, hspace=0, wspace=0)
    axr = fig.add_subplot(grid[0,0])
    axu = fig.add_subplot(grid[1,1])
    axprofile = fig.add_subplot(grid[1,0], sharex=axr, sharey=axu)
    fig.subplots_adjust(wspace=0)
    axr.grid(True)
    axu.grid(True)
    axprofile.grid(True)

    t, r, u, u_smooth = stream.read()

    indices = argwhere( abs(r[0,:] - 1) < 1.1 )
    tmax = 1.
    for i in indices:
        axr.plot( r[:,i], t, '-', lw=.5 )
        #print( 't.x.u', t, u[:,i] )
        axu.plot( t, u[:,i], '-', lw=.5 )
    axprofile.plot( r[-1,:], u_smooth[-1,:], '-', label='t=%s'%t[-1] )
    axprofile.plot( r[0,:], u_smooth[0,:], ':', label='t=%s'%t[0] )
    axr.set_ylabel('t')
    axu.set_xlabel('t')
    axprofile.set_xlabel( r'r(t)' )
    axprofile.set_ylabel( r'u[smooth](t)' )
    #axu.set_xlabel( r'u(t)' )
    axprofile.legend( fancybox=True, framealpha=0 )

    plt.setp(axu.get_yticklabels(), visible=False)
    plt.setp(axr.get_xticklabels(), visible=False)
    #plt.setp(axu.get_yticklabels(), visible=False)
    # remove last tick label for the second subplot
    #yticks = axu.yaxis.get_major_ticks()
    #yticks[-1].label1.set_visible(False)

    fig.tight_layout()
    fig.savefig( stream.figname() )
    print('savefig', stream.figname() )
    plt.close()

class SimulationApp(Gtk.Window):

    def __init__(self):
        Gtk.Window.__init__( self, title="SPH Simulation" )
        Output.cleanLocks( 'simulations/' )
        self.connect( 'destroy', Gtk.main_quit )
        self.set_border_width(3)
        self.maximize()

        self.vbox = Gtk.VBox()

        self.simulationPanels = []
        self.paramsList = ['dt', 'h', 'N', 'IC', 'v', 'SPH' ]
        self.paramsEntry = []
        self.images = []
        self.pixbufs = []
        self.threads = []
        self.createSimulationPanel( self.simulationPanels )
        self.createSimulationPanel( self.simulationPanels )
        self.createSimulationPanel( self.simulationPanels, empty = True )

        scrolled = Gtk.ScrolledWindow()
        scrolled.add_with_viewport( self.vbox )
        scrolled.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        self.add( scrolled )

        self.show_all()

    def onaddclicked( self, widget ):
        self.simulationPanels[-1].destroy()
        del self.simulationPanels[-1]
        self.createSimulationPanel( self.simulationPanels )
        self.createSimulationPanel( self.simulationPanels, True )
        self.show_all()

    def onminusclicked( self, widget, box ):
        box.destroy()
        ind = self.simulationPanels.index(box)
        self.simulationPanels[ ind ].destroy()
        del self.simulationPanels[ ind ]
        del self.paramsEntry[ ind ]
        self.show_all()

    def createSimulationPanel(self, panel, empty = False):
        box = Gtk.HBox()
        panel.append( box )
        this = panel[-1]
        self.vbox.pack_start(box, False, False, 3 )

        if empty:
            button = Gtk.Button()
            box.pack_start( button,  False, False, 3 )
            button.set_label(' + ')
            button.connect( 'clicked', self.onaddclicked )
            return

        options = Gtk.VBox()
        box.pack_start( options, False, False, 3 )

        firstLine = Gtk.HBox()
        title = Gtk.Label()
        destroyButton = Gtk.Button()
        options.pack_start(firstLine, False, False, 3 )

        firstLine.pack_start(title, True, True, 3)
        firstLine.pack_end(destroyButton, False, False, 3)

        title.set_label('simulation')
        destroyButton.set_label(' x ')
        destroyButton.connect( 'clicked', self.onminusclicked, box )

        secondLine = Gtk.HBox()
        thirdLine = Gtk.HBox()
        fourthLine = Gtk.HBox()
        options.pack_start(secondLine, False, False, 3 )
        options.pack_start(thirdLine, False, False, 3 )
        options.pack_start(fourthLine, False, False, 3 )
        label = {}
        self.paramsEntry.append( {} )
        for key in self.paramsList:
            label[key] = Gtk.Label()
            label[key].set_label(key)
            self.paramsEntry[-1][key] = Gtk.Entry()
            self.paramsEntry[-1][key].set_width_chars(5)
            if key == 'IC':
                thirdLine.pack_start( label[key], False, False, 1)
                thirdLine.pack_start( self.paramsEntry[-1][key], True, True, 1)
            elif key == 'SPH':
                fourthLine.pack_start( label[key], False, False, 1)
                fourthLine.pack_start( self.paramsEntry[-1][key], True, True, 1)
            else:
                secondLine.pack_start( label[key], False, False, 1)
                secondLine.pack_start( self.paramsEntry[-1][key], True, True, 1)
        if len( self.paramsEntry ) == 1:
            self.paramsEntry[-1]['dt'].set_text( '0.01' )
            self.paramsEntry[-1]['h'].set_text( '0.05' )
            self.paramsEntry[-1]['N'].set_text( '125' )
            self.paramsEntry[-1]['v'].set_text( '0.1' )
            self.paramsEntry[-1]['IC'].set_text( 'sin(pi*x)' )
            self.paramsEntry[-1]['SPH'].set_text( 'kgf3' )
        elif len( self.paramsEntry ) == 2:
            self.paramsEntry[-1]['dt'].set_text( '0.01' )
            self.paramsEntry[-1]['h'].set_text( '0.05' )
            self.paramsEntry[-1]['N'].set_text( '125' )
            self.paramsEntry[-1]['v'].set_text( '0.01' )
            self.paramsEntry[-1]['IC'].set_text( 'norm(x,.5,.2)' )
            self.paramsEntry[-1]['SPH'].set_text( 'kgf3' )
        else:
            self.paramsEntry[-1]['dt'].set_text( self.paramsEntry[-2]['dt'].get_text() )
            self.paramsEntry[-1]['h'].set_text( self.paramsEntry[-2]['h'].get_text() )
            self.paramsEntry[-1]['N'].set_text( self.paramsEntry[-2]['N'].get_text() )
            self.paramsEntry[-1]['v'].set_text( self.paramsEntry[-2]['v'].get_text() )
            self.paramsEntry[-1]['IC'].set_text( self.paramsEntry[-2]['IC'].get_text() )
            self.paramsEntry[-1]['SPH'].set_text( self.paramsEntry[-2]['SPH'].get_text() )

        runButton = Gtk.Button()
        options.pack_end(runButton, False, False, 3 )
        runButton.set_label('run')
        runButton.connect( 'clicked', self.onclickrun, len(panel)-1 )

        self.images.append( Gtk.Image() )
        self.pixbufs.append( None )
        #scrolled = Gtk.ScrolledWindow()
        #scrolled.add_with_viewport(self.images[-1])
        #scrolled.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)

        box.pack_start( self.images[-1], True, True, 3 )
        return

    def refreshImage( self, panelIndex, figname ):
        self.images[panelIndex].set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file(figname) )
        print('signal', figname )
        self.images[panelIndex].show()

    def onclickrun( self, widget, panelIndex ):
        print( 'run panelIndex', panelIndex )
        L = 1.
        u_0 = 1.
        N = int(self.paramsEntry[panelIndex]['N'].get_text())
        dt = float(self.paramsEntry[panelIndex]['dt'].get_text())
        h0 = float(self.paramsEntry[panelIndex]['h'].get_text())
        v = float(self.paramsEntry[panelIndex]['v'].get_text())
        ic_str = self.paramsEntry[panelIndex]['IC'].get_text()
        sph_str = self.paramsEntry[panelIndex]['SPH'].get_text()

        outname = 'simulations/dt%sh%sN%sv%s_ic=%s_sph%s.dat'%(dt,h0,N,v,ic_str,sph_str)
        figname = outname.replace('.dat', '.png')
        print( 'outname', outname, figname )
        if not os.path.exists('simulations/'): os.makedirs('simulations/')

        ic = lambda x: eval( ic_str )
        #def u0( x, u_0, L ):
            #return u_0*sin(pi*x/L)

        system = LagrangianSystem()
        system.set_timeStep( dt )
        system.set_endTime( 1. )
        r0 = np.array( [ [2*L/N*i, 2*L/N*j] for i in range(N) for j in range(N)] )
        system.set_initialCondition( r0 = r0, u0 = ic )

        burger = Burger()
        burger.set_viscosityCoefficient( v )
        if sph_str.startswith('standard'):
            interpolator = StandardSPH( sph = BaseSPH( h0 = h0, W = 'qspline', density = N/L*h0 ) )
        elif sph_str.startswith('kgf'):
            if sph_str[-1].isdigit(): order=int(sph_str[-1])
            else: order=3
            interpolator = KernelGradientFree( sph = BaseSPH( h0 = h0, W = 'qspline', density = N/L*h0 ), order = order, basepower = 0 )
        else:
            print('unpredicted SPH method')
            return
        burger.set_diff( interpolator )
        system.set_equations( burger )
        system.set_interpolator( interpolator )
        system.set_integrator( TimeIntegrators.Sympletic.RungeKutta4 )

        def run():
            return system.integrate( outname, signal=lambda: GLib.idle_add( lambda panelIndex=panelIndex, figname=figname: self.refreshImage(panelIndex, figname) ) )
        self.threads.append( threading.Thread(target=run) )
        self.threads[-1].daemon = True
        self.threads[-1].start()
        return

norm = lambda x, loc=0, scale=1.: exp(-.5*(x-loc)**2/scale**2)/sqrt(2*pi)/scale

import itertools
def array_to_int( r, h=1 ):
    return np.rint(r/h).astype(int)

def array_to_int_tuple( r, h=1):
    return tuple(array_to_int(r,h))

def get_relative_bin_shifts( d, radius ):
    shifts = range(-radius,radius+1)
    list_of_shifts = [ shifts for i in range(d) ]
    relative_bin_shifts = itertools.product(*list_of_shifts)
    relative_bin_shifts = [ relative_bin_shift for relative_bin_shift in relative_bin_shifts if np.sum(np.power(relative_bin_shift,2)) < (radius + np.sqrt(d))**2 ]
    return tuple( relative_bin_shifts )

def build_sum2( positions, h, d=1, radius=2 ):
    particle_bins = array_to_int( positions, h )
    unique = np.unique( particle_bins, axis=0 )
    print( 'particles/bins', float(len(positions))/len(unique) )

    relative_bin_shifts = get_relative_bin_shifts( d, radius )

    indexes_bin = {}
    number_of_occupied_bins = 0
    for particle_index, particle_bin in enumerate( particle_bins ):
        in_bin = tuple( particle_bin )
        try:
            indexes_bin[in_bin]['in'].append( particle_index )
        except KeyError:
            indexes_bin[in_bin] = { 'in': [ particle_index ], 'around': [] }

        for relative_bin_shift in relative_bin_shifts:
            around_bin = tuple( particle_bin + relative_bin_shift )
            try:
                indexes_bin[around_bin]['around'].append( particle_index )
            except KeyError:
                indexes_bin[around_bin] = { 'in': [], 'around': [ particle_index ] }

    lengths = [ len(v['around']) for k, v in indexes_bin.items() ]
    print( 'link list', min(lengths), max(lengths), np.median(lengths), len(lengths) )
    indexes_bin = { key: index_bin for key, index_bin in indexes_bin.items() if len(index_bin['in']) > 0 }
    lengths = [ len(v['around']) for k, v in indexes_bin.items() ]
    print( 'link list', min(lengths), max(lengths), np.median(lengths), len(lengths) )

    def __sum__( W, v = None ):
        if v is None:
            v = np.ones( len(positions) )
        ret = np.zeros( len(positions) )
        for bin_indexes in indexes_bin.values():
            dists = dist_vectorized( positions[ bin_indexes['in'] ], positions[ bin_indexes['around'] ] )
            ret[ bin_indexes['in'] ] = np.sum( v[ bin_indexes['around'] ][:,None] * W( dists ), axis=0 )
        return ret

    return __sum__


def build_sum( positions, h, d, dist, radius=2 ):
    particle_bins = array_to_int( positions, h)

    particle_bins_tuple = map( tuple, particle_bins )
    occupied_bins, particle_indexes_in_bin = np.unique( particle_bins_tuple, axis=0, return_inverse = True )
    print( 'particles/bins', float(len(positions))/len(occupied_bins) )

    relative_bin_shifts = get_relative_bin_shifts( d, radius )

    indexes_bin = { tuple(occupied_bin): {'in':[], 'around':[]} for occupied_bin in occupied_bins }
    first = True
    for bin_index, occupied_bin in enumerate(occupied_bins):
        particle_indexes = np.nonzero( particle_indexes_in_bin == bin_index )[0].tolist()
        in_bin = tuple(occupied_bin)
        indexes_bin[in_bin]['in'].extend( particle_indexes )

        for relative_bin_shift in relative_bin_shifts:
            around_bin = tuple( occupied_bin + relative_bin_shift )
            try:
                indexes_bin[around_bin]['around'].extend( particle_indexes )
            except KeyError:
                pass

    lengths = [ len(v['around']) for k, v in indexes_bin.items() ]
    print( 'link list statistics', np.median(lengths), len(lengths) )

    def __sum__( v = None, W = None, prime = False ):
        if v is None:
            v = np.ones( len(positions) )[:,None]
        elif type(v) is float:
            v = v*np.ones( len(positions) )[:,None]
        elif len(v.shape) == 1:
            v = v[:,None]
        else:
            raise Exception( 'error %s' % type(v) )
        if not prime:
            ret = np.zeros( len(positions) )
            for bin_indexes in indexes_bin.values():
                dists = dist( positions[ bin_indexes['in'] ], positions[ bin_indexes['around'] ] )
                ret[ bin_indexes['in'] ] = np.sum( v[ bin_indexes['around'] ] * W( dists ), axis=0 )
        else:
            ret = np.zeros( positions.shape )
            for bin_indexes in indexes_bin.values():
                dists = dist( positions[ bin_indexes['in'] ], positions[ bin_indexes['around'] ] )
                directions = positions[ bin_indexes['in'] ][None,:,:] - positions[ bin_indexes['around'] ][:,None,:]
                directions[ dists>0 ] /= dists[:,:,None][ dists>0 ]*h
                ret[ bin_indexes['in'] ] = np.sum( v[ bin_indexes['around'] ][:,:,None] * directions * W( dists )[:,:,None], axis=0 )
        return ret
    return __sum__


def interpolation( r, h, W, x, f ):
    return np.sum( f*W( x, r, h ) )

def generate_IC(mode='random', dimension=1, number_of_fluid_elements=100):
    modes = {'random': 'generate random r and u'}
    if mode is 'random':
        r = np.random.uniform(-1,1, dimension * number_of_fluid_elements)
        u = np.random.uniform(-1,1, dimension * number_of_fluid_elements)
        r = np.reshape(r,(-1,dimension))
        u = np.reshape(u,(-1,dimension))
        return r, u
    raise Exception( '''IC mode "%s" not defined\ndefined modes are %s''' % (mode, modes) )

def run_app():
    app = SimulationApp()
    Gtk.main()

def run_simulation():
    h = .05
    d = 3
    r, u = generate_IC(mode='random', dimension = d, number_of_fluid_elements = 1000)
    kernel_function = kernel.generate_kernel_function(mode='cspline', length=h, dimension=d)
    print( kernel_function(r[0], r[:10]) )
    print( r[0], r, 'r' )
    link_list = make_link_list( r, h, d )
    #print( 'list index 0', link_list[ array_to_int_tuple(r[0],h) ] )
    W = kernel.derivatives('qspline', d)
    print(interpolation( r, h, W[0], r[0], 1 ))

    L = 1.
    u_0 = 1.
    N = 251

    def u0( x, u_0, L ):
        return u_0*sin(pi*x/L)

    v = .05
    system = LagrangianSystem()
    system.set_timeStep( .01 )
    system.set_endTime( 3. )
    system.set_initialCondition( r0 = np.reshape( linspace( 0., 2.*L, N ), (-1,1) ), u0 = lambda x: u0( x, u_0, L ) )

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

    interpolator = KernelGradientFree( sph = BaseSPH( h0 = h0, W = 'qspline', density = N/L*h0 ), order = order, basepower = 0 )
    system.set_interpolator( interpolator )

    system.set_integrator( TimeIntegrators.Sympletic.RungeKutta4 )
    system.integrate( 'simulationv%s.dat'%v )

def dist( a, b ):
    return np.sqrt( np.sum( (a-b)**2, axis=-1 ) )

def dist_vectorized( a, b ):
    return np.sqrt( np.sum( (a[None,:,:] - b[:,None,:])**2, axis=-1 ) )

def dist_weave( x, y ):
    code = '''
    for(int a=0; a < Na; ++a){
        for(int b=0; b < Nb; ++b){
            float sum = 0;
            for(int i=0; i<d; ++i){
                int index_a = i + a*d;
                int index_b = i + b*d;
                sum += pow(x[index_a] - y[index_b], 2);
            }
            int index = a + b*Na;
            ret[index] = pow( sum, .5 );
        }
    }
    '''
    Na = x.shape[0]
    Nb = y.shape[0]
    d = x.shape[1]
    ret = np.zeros((Nb,Na))
    weave.inline( code, ['ret', 'x', 'Na', 'y', 'Nb', 'd'], compiler = 'gcc', extra_compile_args=['-O3'] )
    return ret

class Base:
    def __init__( self, array, name, signature):
        if not isinstance(name, str): raise Exception('%s cannot be used as name' % self.name )
        self.name = name
        self.array = array
        self.signature = signature

    def __getitem__( self, axes ):
        if not type(axes) is tuple: axes = (axes,)
        indices = []
        slices = {}
        if len(axes) is not len(self.signature):
            raise Exception('dimension mismatch %s and %s' %(len(self.signature), len(axes)) )
        for axis, t in zip(axes, self.signature):
            if type(axis) is slice:
                ind = axis.step
                indices.append(ind)
                slices[ind] = slice(axis.start,axis.stop)
            else:
                ind = axis
                indices.append(ind)
                slices[ind] = slice(None)
            if not ind.signature is t:
                raise Exception('signature mismatch %s and %s' %(t, ind.signature) )

        #print( 'getitem', slices )
        return Indexed( self, slices, indices )

    def __setitem__( self, axes, expr ):
        print( 'setitem {}'.format(self.name), axes, expr )
        self.array = expr.eval( axes )

    def __str__(self):
        print('str')
        return str(self.array)
    def __repr__(self):
        print('repr')
        return repr(self.array)
        #return '%s' % ( self.name )

class Indexed:
    def __init__(self, base, slices, indices ):
        self.base = base
        self.slices = slices
        self.indices = indices

    def __repr__(self):
        return '%s[%s]' % ( self.base.name, ','.join(map(repr,self.indices)) )
        #return '%s' % ( self.base.name )

    def __getitem__(self, s):
        print( 'Indexed getitem {}'.format(self.base.name), s )
        return self.base.array[s]

    def __setitem__(self, s, val):
        print( 'Indexed setitem {}'.format(self.base.name), s, val )
        return self.base.array[s]

    def __mul__(self, other):
        return Expr( Mul, (self, other) )

    def __add__(self, other):
        return Expr( Add, (self, other) )

class Function:
    def __init__(self, name, func):
        self.name = name
        self.func = func

    def __call__(self, *args):
        return Expr( self, args )

Add = Function('add', lambda x,y: x+y )
Mul = Function('mul', lambda x,y: x*y )

class Expr:
    def __init__( self, func, args ):
        #print('constructed expr', func, args )
        self.func = func
        self.args = args

    def eval( self, axes, level=0 ):
        if type(axes) is not tuple: axes = (axes,)
        if level is 0: print('eval')
        print( '  '*level + repr(self) )
        new_args = []
        for arg in self.args:
            if isinstance(arg, Base):
                thisslice = [ arg.slices[axis] if axis in arg.slices.keys() else None for axis in axes ]
                new_args.append( arg[thisslice] )
            elif isinstance(arg, Expr):
                #print( 'arg', arg )
                new_args.append( arg.eval(axes, level=level+1) )
        ret = self.func( *new_args )
        #print( 'ret', ret.shape )
        return ret

    def __repr__( self ):
        #print( self.args )
        return '%s(%s)' %(self.func.name, ', '.join( map( repr, self.args )) )

    def __mul__( self, other ):
        return Expr( Mul, (self, other) )

    def __add__( self, other ):
        return Expr( Add, (self, other) )

class Index:
    def __init__(self, name, signature):
        self.name = name
        self.signature = signature

    def __repr__(self):
        return '%s' % (self.name)

class BaseSum:
    def __init__(self):
        pass

    def __getitem__(self, axis):
        return lambda expr: Expr( Sum(axis), expr )


def get_neighbor_positions( particle_position, all_positions, link_list, h ):
    neighbor_indexes = link_list[ array_to_int_tuple( particle_position, h) ]
    return all_positions[ neighbor_indexes ]

def get_h( number_of_fluid_elements, fluid_elements_per_bin, volume, d ):
    number_of_bins = float(number_of_fluid_elements)/(fluid_elements_per_bin/5.**d)
    h = np.power( volume/number_of_bins, 1./d )
    return h

class StateBase:
    def __init__( self, r, u ):
        self.__r = r
        self.__u = u

class SystemBase:
    def __init__( self, energy_distribution, external_potential ):
        self.__e = energy_distribution
        self.__V = external_potential

    def evolve( self, dt, exit_condition ):
        pass

def generateIC( e, N, h ):
    r = zeros(N)*nan
    unassigned = isnan(r)
    e_max = max(e)
    while True:
        new_r = random.uniform( e.shape )
        test = random.uniform( new_r.size )



def test_link_list():
    L = np.array([1,1,1])
    number_of_fluid_elements = int(1e3)
    d = 1
    volume = np.prod(L)
    h = get_h( number_of_fluid_elements, 50., volume, d )
    number_of_bins = volume/h**d
    #h = .02
    print('h',h)
    N = 10
    mask = slice( None,N )

    positions, velocities = generate_IC( mode='random', dimension=d, number_of_fluid_elements = number_of_fluid_elements )
    print( 'shape', positions.shape, np.prod(L), number_of_bins, number_of_fluid_elements/number_of_bins )

    with Timer('brute force') as timer:
        dists =  dist( positions[None,:N,:], positions[:,None,:] )
        ret = np.sum( kernel.cspline_weave( dists, d, h ), axis = 0 )
        print( ret, ret.shape, timer.time( number_of_fluid_elements/N ) )
    print()

    W = kernel.generate_cspline( d, h )
    W_prime = kernel.generate_cspline_prime( d, h )

    print()
    dt = h
    velocities = np.zeros_like(velocities)
    for i in range(10):
        a = ( slice(None,None), None )
        b = ( None, slice(None,None) )
        axis_a = 0
        axis_b = 1
        mask = positions[:,0] > .99
        with Timer('build sum weave'):
            print( 'r', positions[mask] )
            print( 'v', velocities[mask] )
            with Timer('build sum weave'):
                _sum_ = build_sum( positions, h, d, dist = dist_vectorized )
            with Timer('rho'):
                rho = _sum_( None, W )
                print( 'rho', rho[mask], rho.shape )
            #with Timer('vol'):
                #vol = _sum_( 1./rho, W )
                ##print( 'vol', vol[mask], vol.shape )
            with Timer('P'):
                P = rho**(1./3)
                F = _sum_( (P/rho**2), W_prime, prime=True ) + (P/rho**2)[a] * _sum_( 1., W_prime, prime=True )
                print( 'F', F[mask], F.shape )
            dv_dt = -F
            dr_dt = velocities

            positions += dr_dt * dt
            velocities += dv_dt * dt
    print()

def test_IC():
    Lx = 10
    Npoints = 10
    L = array([[ -Lx, Lx, Npoints ], [ -Lx, Lx, Npoints ]])
    ix = [ linspace(*Li) for Li in L ]
    xv = meshgrid( *ix )
    e = scipy.stats.norm.pdf(xv[0])*scipy.stats.norm.pdf(xv[1])
    print( e.shape )

    e_max = amax(e)
    print( e )
    print( e[Lx/2,Lx/2]/e_max )
    Nelements = 10
    print( len(L) )
    r = random.random( size=Nelements*len(L) ).reshape( (len(L),-1) )*(L[:,1] - L[:,0])[:,None] + L[:,0][:,None]
    print( r.shape )

class System:
    def __init__( self, positions, canonical_momenta, kernel, kernel_prime ):
        if len(positions) != len(canonical_momenta):
            print( 'different sizes')
            exit(0)
        self.__positions = positions
        self.__canonical_momenta = canonical_momenta
        self.size = len(self.__positions)
        self.kernel = kernel
        self.kernel_prime = kernel_prime
        return

    def build_sum( self, support, dimension, distance_func, radius=2 ):
        digitized_positions = map( tuple, rint(self.positions/support).astype(int) )

        occupied_bins, inverse_indices = np.unique( digitized_positions, axis=0, return_inverse = True )
        print( 'particles/bins', float( self.size )/len(occupied_bins) )

        relative_bin_shifts = get_relative_bin_shifts( dimension, radius )

        # particles_in_bin = { tuple(occupied_bin): ([],[]) for occupied_bin in occupied_bins }
        particles_in_bin = { tuple(occupied_bin): [] for occupied_bin in occupied_bins }
        neighbor_particles_of_bin = { tuple(occupied_bin): [] for occupied_bin in occupied_bins }

        for bin_index, occupied_bin in enumerate( occupied_bins ):
            particles = np.nonzero( inverse_indices == bin_index )[0].tolist()

            particles_in_bin[occupied_bin].extend( particles )

            for relative_bin_shift in relative_bin_shifts:
                around_bin = tuple( occupied_bin + relative_bin_shift )
                try:
                    # particles_in_bin[around_bin][1].extend( particles )
                    neighbors_particles_of_bin[around_bin].extend( particles )
                except KeyError:
                    pass

        lengths = map(len, neighbors_particles_of_bin.values() )
        print( 'link list statistics',
                np.mean(lenghts),
                np.median(lengths),
                len(lengths)
                )

        def __sum_W__( v = None, W = None, prime = False ):
            ret = np.zeros_like( self.positions )
            for occupied_bin in occupied_bins:
                a = particles_in_bin[occupied_bin]
                b = neighbor_particles_of_bin[occupied_bin]
                r_a = self.positions[a]
                r_b = self.positions[b]
                distances = distance_function( r_a, r_b )
                ret[a] = sum( v[b] * W( distances ), axis=0 )
            return ret

        def __sum_W_prime__():
            ret = np.zeros_like( self.positions )
            for occupied_bin in occupied_bins:
                a = particles_in_bin[occupied_bin]
                b = neighbor_particles_of_bin[occupied_bin]
                r_a = self.positions[a]
                r_b = self.positions[b]
                distances = distance_function( r_a, r_b )
                directions = direction_function( r_a, r_b )
                directions[ distances>0 ] /= distances[:,:,None][ distances>0 ]
                ret[a] = sum( v[b] * directions * W_prime( dists ), axis=0 )
            return ret

        return __sum_W__, __sum_W_prime__

    @property
    def positions(self):
        return self.__positions

    @property
    def canonical_momenta(self):
        return self.__canonical_momenta

    def sum(self, a):
        return sum( a, axis=1 )

    def system_summation(self, a):
        return sum( a, axis=1 )

    @property
    def dr_dt(self):
        try:
            return self.__positions_time_differentials
        except AttributeError:
            self.__positions_time_differential = self.canonical_momenta#/self.rho[:,None]**2
            return self.__positions_time_differential

    @property
    def dp_dt(self):
        try:
            return self.__canonical_momenta_time_differential
        except AttributeError:
            self.__canonical_momenta_time_differential = -self.sum( ( self.forces[None,:,None] + self.forces[:,None,None] ) * self.gradW(self.positions, self.positions) )
            return self.__canonical_momenta_time_differential

    @property
    def forces(self):
        try:
            return self.__forces
        except AttributeError:
            self.__forces = self.pressures/self.neighbor_densities**2# - self.pSqr/self.rho**3
            return self.__forces

    @property
    def pressures(self):
        try:
            return self.__pressures
        except AttributeError:
            self.__pressures = self.neighbor_densities**(4./3)
            return self.__pressures

    @property
    def energy_densities(self):
        return 3*self.pressures

    @property
    def pSqr(self):
        return sum(self.canonical_momenta**2, axis=-1)

    @property
    def neighbor_densities(self):
        try:
            return self.__neighbor_densities
        except AttributeError:
            self.__neighbor_densities = self.sum( self.W(self.positions, self.positions) )
            #print( 'rho.shape', self.__rho.shape )
            return self.__neighbor_densities

    def W(self, x, x2):
        diffs = x[:,None,:] - x2[None,:,:]
        dists = sqrt( sum(diffs**2, axis=-1) )
        return self.kernel( dists )

    #@property
    #def W(self):
        #try:
            #return self.__W
        #except AttributeError:
            #self.__W = self.kernel(self.dists)
            #return self.__W

    def gradW(self, x, x2):
        diffs = x[:,None,:] - x2[None,:,:]
        dists = sqrt( sum(diffs**2, axis=-1) )
        dirs = diffs
        dirs[ dists!=0 ] /= dists[ dists!=0 ][:,None]
        return self.kernel_prime(dists)[:,:,None] * dirs


    @property
    def Hamiltonian(self):
        try:
            return self.__Hamiltonian
        except AttributeError:
            self.__Hamiltonian = .5*sum( self.pSqr ) + sum( self.energy_densities/self.neighbor_densities )
            return self.__Hamiltonian

    @property
    def J(self):
        try:
            return self.__J
        except AttributeError:
            self.__J = sum( self.canonical_momenta, axis=0 )
            return self.__J

    @property
    def velocities(self):
        try:
            return self.__velocities
        except AttributeError:
            self.__velocities = self.canonical_momenta/self.neighbor_densities**2

    def ascii_e(self, a, b, l=100):
        x = linspace(a, b, l)
        W = self.kernel( abs(x[:,None] - self.r[None,:,0]) )
        rho = sum( W, axis=1 )
        #print( W )
        e = 3*rho**(4./3)
        e_max = amax(e)
        s = (' ',' ̣','.','·','˙','^')
        line1 = ''.join( [ s[int((len(s)-1)*ie/e_max)] for ie in e ] )
        line2 = [' ']*l
        chars = range(ord('0'), ord('9')) + range(ord('a'), ord('z')) + range(ord('A'), ord('Z'))
        for i, ri in enumerate(self.__r):
            r__ = int( (ri[0] - a)/(b - a) * l )
            if r__ < len(line2):
                line2[r__] = '%c' % chars[i % len(chars)]
        line2 = ''.join(line2)
        return '\n'.join((line1,line2))

    def ascii_e_2d(self, a, b, length = 100, maximum_energy_density = None):
        if maximum_energy_density is None:
            maximum_energy_density = amax(self.energy_densities)
        bg = lambda text, color: "\33[48;5;" + str(color) + "m" + text + "\33[0m"
        factor = .4
        x_ = linspace(a, b, length)
        y_ = linspace(a, b, int(length*factor))
        r__ = array( [ [xi, yi, 0] for yi in y_ for xi in x_ ] )
        pos_min = array([a,a,a])[None,:]
        pos_max = array([b,b,b])[None,:]
        pos = (self.positions - pos_min)/(pos_max - pos_min)
        pos[:,0] *= length
        pos[:,1] = (pos[:,1]*length*factor).astype(int)
        pos = pos.astype(int)
        pos = [ (ipos[0], ipos[1]) for ipos in pos ]
        dists = sqrt(sum( (r__[:,None,:] - self.positions[None,:,:])**2, axis=-1 ))
        W = self.kernel( dists )
        neighbor_densities = sum( W, axis=1 )
        energy_densities = 3 * neighbor_densities**(4./3)
        s = ''
        #exit(0)
        for y in range(len(y_)):
            for x in range(len(x_)):
                #print( (e/e_max)/(195-160) + 160 )
                #e_ = log( e[x+y*l]+1e-3 )/log( e_max+1e-3 )
                scaled_energy_density = energy_densities[x+y*length]/maximum_energy_density
                color_code = int( scaled_energy_density*(255-232) ) + 232
                #print( int(e__), end=' ' )
                count = pos.count((x,y))
                marker = (str(count) if count < 10 else '9') if count>0 else ' '
                s += bg(marker, int(color_code))
            s += '\n'
        os.system('cls || clear')
        print( s )
        print('maximum_energy_density', maximum_energy_density, end=' ')
        return


def collision_test():
    number_of_fluid_elements = 1000
    ri = scipy.stats.norm.rvs(0,.1, number_of_fluid_elements)
    rj = scipy.stats.norm.rvs(0,.1, number_of_fluid_elements)
    r1 = array([[ri_-.8,rj_,0] for ri_, rj_ in zip(ri,rj) ]).astype(float)
    p1 = array([[20,0,0] for ri_,rj_ in zip(ri,rj) ]).astype(float)

    ri = scipy.stats.norm.rvs(0,.1, number_of_fluid_elements)
    rj = scipy.stats.norm.rvs(0,.1, number_of_fluid_elements)
    r2 = array([[ri_+.8,rj_,0] for ri_, rj_ in zip(ri,rj) ]).astype(float)
    p2 = array([[-20,0,0] for ri_,rj_ in zip(ri,rj) ]).astype(float)

    positions = concatenate( (r1,r2) )
    canonical_momenta = concatenate( (p1,p2) )
    print( 'r.shape', positions.shape )
    print( 'p.shape', canonical_momenta.shape )
    dimension = positions.shape[1]
    support = .05
    time_step = support/10.

    #W = kernel.generate_liq( dim, h )
    #W_prime = kernel.generate_liq_prime( dim, h )

    kernel_function = kernel.generate_cspline( dimension, support )
    kernel_function_derivative = kernel.generate_cspline_prime( dimension, support )

    next_canonical_momenta = canonical_momenta
    next_positions = positions

    #c = [1, -2./3, 2./3]
    #d = [-1./24, 3./4, 7./24]

    c = [
        .5/(2-2**(1./3)),
        .5*(1-2**(1./3))/(2-2**(1./3)),
        .5*(1-2**(1./3))/(2-2**(1./3)),
        .5/(2-2**(1./3)),
        ]
    d = [
        1./(2-2**(1./3)),
        -2**(1./3)/(2-2**(1./3)),
        1./(2-2**(1./3)),
        0
        ]

    system = System( positions, canonical_momenta, kernel_function, kernel_function_derivative )
    initial_Hamiltonian = system.Hamiltonian
    maximal_energy_density = amax( system.energy_densities )
    time = 0
    L = 1
    while True:
        for i, (ci, di) in enumerate(zip(c, d)):
            print( i )
            system = System( next_positions, next_canonical_momenta, kernel_function, kernel_function_derivative )
            next_positions = system.positions + ci * time_step * system.dr_dt
            system = System( next_positions, next_canonical_momenta, kernel_function, kernel_function_derivative )
            next_canonical_momenta = system.canonical_momenta + di * time_step * system.dp_dt
        time += time_step

        positions, canonical_momenta = system.positions, system.canonical_momenta
        if amax( abs(positions[:,0]) ) > L or amax( abs(positions[:,1]) ) > L:
            L *= 2
        #for a in range(len(r)):
            #for i in range(dim):
                #if r[a][i] < -1 and p[a][i] < 0:
                    #r[a][i] = -1
                    #p[a] = abs(p[a])
                #elif r[a][i] > 1 and p[a][i] > 0:
                    #r[a][i] = 1
                    #p[a] = -abs(p[a])
        next_canonical_momenta = canonical_momenta
        next_positions = positions
        system.ascii_e_2d( -L, L, 70, None)
        print( 'time=', time,
              '[{}]'.format( (system.Hamiltonian - initial_Hamiltonian)/initial_Hamiltonian ),
              'pmax', amax( sqrt(system.pSqr) ),
              'vmax', amax( sqrt(sum(system.dr_dt**2, axis=-1)) )
              )

if __name__ == '__main__':
    #test_link_list()
    #test_IC()
    collision_test()
    exit(0)

exit(0)
