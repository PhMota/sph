import os
from sympy import *
from sympy.core.singleton import S

init_printing()
opening_string = r'''
\documentclass{article}
\usepackage{amsmath}
%\usepackage{nath}
%\delimgrowth=1
\usepackage{breqn}
\delimitershortfall=-1pt
\def\operatorname{}
\begin{document}
'''
open('output.tex','w').write( opening_string )
def print2file( str ):
    output = open('output.tex','a')
    output.write( str + '\n' )
    return

def section(title):
    print2file( r'\section{%s}' % title )

def par(s):
    print2file( s )

def eq( lhs, rhs = None ):
    #s = r'$$'
    s = r'\begin{dmath}'
    if rhs is None:
        s += '\n' + latex( simplify(lhs.doit()).expand() )
    else:
        s += '\n' + latex( lhs )
        s += '\n=\n' + latex( simplify(rhs.doit() ).expand() )
    s += '\n$$'
    s += '\n'+r'\end{dmath}'
    print2file( s )

def convert2pdf():
    os.system('pdflatex output.tex')

def const_Diff(self, var, **kwargs):
    return S.Zero
Integer.Diff = const_Diff
Rational.Diff = const_Diff

def symbol_Diff( self, var, **kwargs ):
    if self == var:
        return S.One
    return S.Zero
Symbol.Diff = symbol_Diff

def mul_Diff( self, var, **kwargs ):
    df = 0
    prod = 1
    for i, arg in enumerate(self.args):
        prod *= arg
    for i, arg in enumerate(self.args):
        df += prod/arg*arg.Diff( var, **kwargs )
    return df
Mul.Diff = mul_Diff

def pow_Diff( self, var, **kwargs ):
    a, b = self.args
    return a**b * ( a.Diff(var, **kwargs)*b/a + b.Diff(var, **kwargs)*log(a) )
Pow.Diff = pow_Diff

def add_Diff( self, var, **kwargs ):
    df = 0
    for i, arg in enumerate(self.args):
        df += arg.Diff( var, **kwargs )
    return df
Add.Diff = add_Diff

def sum_Diff( self, var, **kwargs ):
    return Sum( self.args[0].Diff( var, **kwargs ), self.args[1] )
Sum.Diff = sum_Diff

def integral_Diff( self, var, **kwargs ):
    return Integral( self.args[0].Diff( var, **kwargs ), self.args[1] )
Integral.Diff = integral_Diff

def func_diff( self, var, **kwargs ):
    return self.diff(var, **kwargs)
Function.Diff = func_diff

def Diff( self, var, **kwargs ):
    return self.Diff(var, **kwargs)

Idx.name = 'idx'

N = Symbol('N', integer=True)
i = Idx('i', N+1)
j = Idx('j', N+1)
k = Idx('k', N+1)

x = Symbol('x', real=True)
y = Symbol('y')
t = Symbol('t')
h = Symbol('h')

class nu(Function):
    is_real = True
    def Diff( self, var, **kwargs ):
        return S.Zero
    def diff( self, var ):
        return S.Zero
    def _latex(self, *args, **kwargs):
        return r'\nu_{%s}' % (self.args[0])

class W(Function):
    def fdiff( self, argindex ):
        return Function('W_%s' % argindex )(*self.args)
    
    def _eval_Integral(self, arg):
        return S.One

    def diff( self, var ):
        #print 'W.diff', 'args=', self.args, 'self=', self, 'var=', var
        if isinstance( var, Dummy ):
            #print 'Dummy found in W.diff'
            return Derivative( self.func( *self.args ), var )
        #print 'W.diff', self.args[0], var 
        return self.fdiff(1) * self.args[0].diff( var )
    
    def Diff( self, var, **kwargs ):
        return self.fdiff(1) * self.args[0].Diff( var )

class v(Function):
    def Diff( self, var, **kwargs ):
        if isinstance(var, r):
            return S.Zero
        if isinstance(var, Idx):
            return S.Zero
        if isinstance(var, v):
            #print 'v', self.args, var.args
            if self.args[1] == var.args[1]:
                return S.One
            return KroneckerDelta( self.args[1], var.args[1] )
        return self.fdiff(1) * self.args[0].Diff( var )
    diff = Diff

    def _latex(self, *args, **kwargs):
        if len(self.args) == 1:
            return 'v(%s)' % (self.args[0])
        return 'v_{%s}(%s)' % (self.args[1], self.args[0])

class r(Function):    
    is_real = True
    _diff_wrt = True
    def fdiff(self, argindex=1):
        if argindex == 1:
            return v(*self.args)
        raiseExepction( 'cannot derivate arg %s' % argindex )
    
    def Diff( self, var, **kwargs ):
        if isinstance(var, v):
            return S.Zero
        if isinstance(var, r):
            if self.args[1] == var.args[1]:
                return S.One
            return KroneckerDelta( self.args[1], var.args[1] )
        return self.fdiff(1) * self.args[0].Diff( var )
    diff = Diff

    def _latex(self, a):
        return 'r_{%s}(%s)' % (self.args[1], self.args[0])

class sph_interpolation(Function):
    @classmethod
    def eval( cls, arg, var, index ):
        i = symbols('i', cls=Idx)
        
        kernel = nu(index)*W( (var-r(t, index))/h )
        range = (index, 0, N)
        if arg is None:
            return Sum( kernel, range )
        
        elif isinstance(arg, FunctionClass):
            return Sum( arg.func(t, index)*kernel, range )
        
        return Sum( arg*kernel, range )
        print 'unpredicted instance of arg', type(arg)

class rho(Function):
    _diff_wrt = True
    def explicit( self ):
        arg = self.args[0]
        index = Idx( '%s\'' % arg.args[-1], N )
        return sph_interpolation( None, arg, index)
    
    def fdiff(self, argindex ):
        return Function('rho_%s' % argindex )(*self.args)
    
    def diff(self, var):
        if isinstance(var, Dummy):
            return Derivative( self.func(*self.args), var )
        return self.fdiff(1) * self.args[0].diff( var )

    def Diff(self, var ):
        if isinstance(self.args[0], r):
            arg = self.args[0]
            index = Idx( '%s_1' % arg.args[-1], N+1 )
            return sph_interpolation( None, arg, index ).Diff( var )
        return self.fdiff(1) * self.args[0].Diff( var )

class Monaghan(Function):
    @classmethod
    def eval(cls, arg, var, index):
        innerIndex = symbols('%s\'' % index.label, cls=Idx)
        if arg is None:
            return sph_interpolation( 1/rho(r(t,index)), var, index )
        #pprint( arg )
        arg = arg.replace( x, r(t,index) )
        #pprint( arg )
        return sph_interpolation( arg/rho(r(t,index)), var, index )

class FirstOrder(Function):
    @classmethod
    def eval(cls, arg, var, index):
        innerIndex = symbols('%s\'' % index.label, cls=Idx)
        if arg is None:
            return sph_interpolation( 1/rho(r(t,index)), var, index )
        arg = arg.replace( x, r(t,index) )
        return sph_interpolation( arg, var, index )/rho(var)
    
def explicit( expr, level=0 ):
    if hasattr( expr, 'explicit' ):
        newexpr = expr.explicit()
        return newexpr
    if len(expr.args) == 0:
        return expr
        
    newargs = tuple( explicit(arg, level=level+1) for arg in expr.args )
    return expr.func(*newargs)

class Conservation_Eq(Function):
    @classmethod
    def eval(cls, arg):
        eq = Dmat( arg )
        return solve( eq, v(t, i) )

class Sum(Sum):
    def _eval_Integral(self, arg):
        return Sum( Integral( self.function, arg ).doit(), (i,0,N) )
    
    def diff(self, var):
        #pprint( self.args[0], var )
        #pprint( self.args[0].diff(var) )
        return Sum( self.args[0].diff(var), self.args[1] )
    
class EulerLagrange_Eq(Function):
    @classmethod
    def eval(cls, arg, var, dvar, **kwargs ):
        eq = arg.Diff( var ) - arg.Diff( dvar ).Diff(t)
        return eq
        #return simplify( eq.doit() )

class epsilon(Function):

    def fdiff( self, argindex ):
        return Function('epsilon_%s' % argindex )(*self.args)
    
    def Diff( self, var, implicit=True, **kwargs ):
        fdiff = None
        arg = self.args[0]
        if isinstance(arg, Function):
            fdiff = Function('epsilon_%s' % arg.func )(*self.args)
        else:
            fdiff = self.fdiff(1)
        return fdiff * self.args[0].Diff( var )

#eq( epsilon( r(t,i) ).Diff( r(t,j) ) )
##eq( Derivative( epsilon( rho( rho( r(t,i) ) ) ), r(t,j) ), Diff( epsilon( rho( rho( r(t,i) ) ) ), r(t,j) ) )
#epsilon_over_rho = lambda _: epsilon(rho(_))/rho(_)
#eq( Derivative( epsilon_over_rho( x ), x ), ( epsilon( rho( x ))/rho( x ) ).Diff( x ) )
#eq( Derivative( epsilon( r(t,i) )/rho( r(t,i) ), r(t,j) ), ( epsilon( r(t,i) )/rho( r(t,i) ) ).Diff( r(t,j) ) )

#eq( epsilon( r(t,i) ).Diff(t) )
#eq( (epsilon( r(t,i) )/rho(r(t,i))).Diff(t) )

#eq( Diff( epsilon( r(t,i) ), t ) )

par('All calculations here are performed automatically SymPy')
section( 'kernel' )
par('The Kernel function is defined by the following properties. Its integral is')
eq( integrate( W((x - r(t,i)/h)), x ) )
par('Its derivative in respect to $x$')
eq( Derivative( W((x - r(t,i))/h), x ), W( (x - r(t,i))/h ).Diff( x ) ) 
par('Its derivative in respect to $r(t,j)$')
eq( Derivative( W(Abs(x - r(t,i))/h), r(t,j) ), W( (x - r(t,i))/h ).Diff( r(t,j) ) ) 

section('Interpolation')
par('the interpolation reference is')
eq( rho(x), sph_interpolation( None, x, i) )
par('and its derivative is')
eq( Derivative( rho(x), x ), simplify( sph_interpolation( None, x, i).Diff( x ).doit() ) )
par('on the other hand, when taking the variational approach one gets')
eq( Derivative( rho( r(t,k) ), r(t,j) ), simplify( sph_interpolation( None, r(t,k), i).Diff( r(t,j) ).doit() ) )

_L = Function('\mathcal{L}')
L = Function('L')

section('Monaghan recipe')
par('The Lagrangean density is')
eq( _L(x,t), S.Half*rho(x)*v(x)**2 - rho(x)*epsilon(rho(x)) )

LagrangeanDensity = Monaghan( S.Half*rho(x)*v(t,i)**2 - rho(x)*epsilon(rho(x)), x, i )
par('Using the Monaghan recipe it transcribes into with the transformation of $v(x)$ into $v(t,i)$')
eq( _L(x,t), simplify(LagrangeanDensity) )
Lagrangean = Integral( LagrangeanDensity, x )
par('The integrated Lagrangean density is')
eq( L(x,t), simplify(Lagrangean) )
par('The variation in respect to r')
eq( Derivative( L(x, t), r(t,j) ), simplify( Lagrangean.Diff( r(t,j) ).doit() ) )
par('The variation in respect to v')
eq( Derivative( L(x, t), v(t,j) ), simplify( Lagrangean.Diff( v(t,j) ).doit() ) )
par('and its time derivative')
eq( simplify( Lagrangean.Diff( v(t,j) ).Diff(t).doit() ) )
par('Finally the Euler-Lagrange relation')
eq( 0, simplify( EulerLagrange_Eq( Lagrangean, r(t,j), v(t,j) ).doit() ) )

open('output.tex','a').write( r'\end{document}' + '\n' )
convert2pdf()
