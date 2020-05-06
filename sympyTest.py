import subprocess
from sympy import *
from sympy.core.singleton import S

class PDF:
    def writelines( self, s, mode = 'a' ):
        open( '%s.tex' % self.fname, mode ).writelines( s )
    
    def __init__(self, fname='calculations'):
        init_printing()
        self.fname = fname
        preamble = [
            r'\documentclass{article}',
            r'\usepackage{amsmath}',
            r'\usepackage{breqn}',
            r'\delimitershortfall=-1pt',
            r'\def\operatorname{}',
            r'\begin{document}'
            ]
        self.writelines( preamble, mode = 'w' )
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.writelines( [r'\end{document}'] )
        self.pdflatex()
    
    def par( self, s ):
        self.writelines( [s] )

    equation_evironment = 'equation'
    def pdflatex( self ):
        subprocess.check_output(['pdflatex', '%s.tex' % self.fname])

    def section( self, text ):
        self.writelines( [ r'\section{%s}' % text] )
    
    def environment( self, s, env = 'equation' ):
        self.writelines( [r'\begin{%s}' % env] )
        self.writelines( s )
        self.writelines( [r'\end{%s}' % env] )

    def eq( self, lhs, rhs = None ):
        print 'eq', lhs
        r = [r'\begin{dmath}', latex(lhs)]
        if not rhs is None:
            r += ['=', latex(rhs)]
        r += [r'\end{dmath}']
        self.writelines( r )

    def math( self, lhs, rhs = None, label = None ):
        print 'eq', label
        r = [r'\begin{equation}', latex(lhs)]
        if not rhs is None:
            r += ['=', latex(rhs)]
        if not label is None:
            r += [r'\label{%s}'%label]
        r += [r'\end{equation}']
        self.writelines( r )

    def matheval( self, s ):
        print 'eq', label
        self.writelines( [r'\begin{equation}', latex(eval(m)), r'\label{%s}'%m, r'\end{equation}'] )


def openPDF( fname ): return PDF( fname )

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

#def func_diff( self, var, **kwargs ):
    #return self.Diff(var, **kwargs)
#Function.Diff = func_diff

KroneckerDelta.Diff = lambda *args: S.Zero

def Diff( self, var, **kwargs ):
    return self.Diff(var, **kwargs)

def deriv_Diff(self, var):
    return self.diff(var)
Derivative.Diff = deriv_Diff

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
        return Function( r'W^{\prime}' )(*self.args)
    
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
    #def fdiff(self, argindex):
        ##return Derivative( self, 
        #return Function('v_{%s}' % self.args[argindex-1])(*self.args)
    
    def Diff( self, var, **kwargs ):
        if isinstance(var, r):
            return S.Zero
        if isinstance(var, Idx):
            return S.Zero
        if isinstance(var, v):
            if len(self.args) > 1:
                if self.args[1] == var.args[1]:
                    return S.One
                return KroneckerDelta( self.args[1], var.args[1] )
        return Derivative(self, self.args[0]) * self.args[0].Diff( var )
    diff = Diff

    #def _latex(self, *args, **kwargs):
        ##print self.base
        #if len(self.args) == 1:
            #return 'v(%s)' % (self.args[0])
        #return 'v_{%s}(%s)' % (self.args[1], self.args[0])

class dot_r(Function):
    def Diff( self, var, **kwargs ):
        if isinstance(var, r):
            return S.Zero
        if isinstance(var, dot_r):
            if self.args[1] == var.args[1]:
                return S.One
            return KroneckerDelta( self.args[1], var.args[1] )
        if var == t:
            return Function(r'\ddot{r}')(*self.args)
        return S.Zero
    diff = Diff
    
class r(Function):
    is_real = True
    _diff_wrt = True
    
    def Diff( self, var, **kwargs ):
        if isinstance(var, dot_r):
            return S.Zero
        if isinstance(var, r):
            if self.args[1] == var.args[1]:
                return S.One
            return KroneckerDelta( self.args[1], var.args[1] )
        if var == t:
            return dot_r(*self.args)
        return S.Zero
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
        raise Exception('called diff in rho')

    def Diff(self, var ):
        #if isinstance(self.args[0], r):
            #arg = self.args[0]
            #return sph_interpolation( None, arg, index ).Diff( var )
        index = None
        if isinstance(self.args[0], r):
            index = Idx( '%s_1' % self.args[0].args[-1], N+1 )
        else:
            index = i
        return sph_interpolation( None, self.args[0], index ).Diff( var )
        
        #return self.fdiff(1) * self.args[0].Diff( var )

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
    def Diff( self, var, implicit=True, **kwargs ):
        fdiff = None
        arg = self.args[0]
        if isinstance( arg, rho ):
            return epsilon_rho(*self.args) * self.args[0].Diff( var )
        raise Exception( 'problem in Diff arg=%s, var=%s, argtype=%s' % (arg, var, type(arg) ) )

class epsilon_rho(Function):
    def Diff( self, var ):
        return Function('epsilon_2' )(*self.args) * self.args[0].Diff( var )
    
#eq( epsilon( r(t,i) ).Diff( r(t,j) ) )
##eq( Derivative( epsilon( rho( rho( r(t,i) ) ) ), r(t,j) ), Diff( epsilon( rho( rho( r(t,i) ) ) ), r(t,j) ) )
#epsilon_over_rho = lambda _: epsilon(rho(_))/rho(_)
#eq( Derivative( epsilon_over_rho( x ), x ), ( epsilon( rho( x ))/rho( x ) ).Diff( x ) )
#eq( Derivative( epsilon( r(t,i) )/rho( r(t,i) ), r(t,j) ), ( epsilon( r(t,i) )/rho( r(t,i) ) ).Diff( r(t,j) ) )

#eq( epsilon( r(t,i) ).Diff(t) )
#eq( (epsilon( r(t,i) )/rho(r(t,i))).Diff(t) )

#eq( Diff( epsilon( r(t,i) ), t ) )

def write_kernel( document ):
    document.section( 'kernel' )

    document.par('The Kernel is determined by the desired property of the density integral form,')
    f = Function('f')
    i = Symbol('i', real=True, positive=True)
    N = Symbol('N', real=True, positive=True)
    F = IndexedBase('F')
    V = IndexedBase('V')
    sum = lambda func, ind: Sum( func, (ind, 0, N) )
    document.math( Integral(f(x), x), sum( F[i]*V[i], i ), label='int.f' )
    
    r = IndexedBase('r')

    document.par('this is needed so that the expression of the integral of the density does not depend on x explicitly.')
    document.par('The quantities $F$ and $V$ are not defined yet, but we suppose that $F$ is related to the function and $V$ is not')
    document.par('The general form of the interpolation using a general kernel function $G$ is')
    G = Function('G')
    interpG = lambda func, x, ind: sum( func * G(x - r[ind]), ind )
    interpG_f = lambda x, i: interpG( F[i], x, i )
    document.math( f(x), interpG_f(x, i), label='f' )

    document.par(r'The relations \eqref{f} and \eqref{int.f} bind the integral of the kernel function,')
    document.eq( Integral( Function('G')(x-Indexed('r', i) ), x), Indexed('V', i) )

    document.par('Therefore we can write the general kernel as a normalized kernel $W$ times the fluid element volume $V_i$')
    W = Function('W')
    document.eq( G(x-r[i]), V[i]*W(x-r[i]) )

    interp = lambda func, x, ind: sum( func * W(x - r[ind]), ind )
    interp_f = lambda x, ind: interp( F[ind]* V[ind], x, ind )
    document.par(r'Taking the interpolation computed on $r_j$,')
    j = Symbol('j', real=True, positive=True)
    document.eq( f(r[i]), interp_f( r[i], j ) )
    document.par(r'and summing it with the kernel gives')
    document.eq( interp( f(r[i]), x, i ), interp( interp_f( r[j], i ), x, j).doit().factor() )

    document.par(r'which is satisfied if')
    document.eq( f(r[i])*W(x-r[i]), F[i]*V[i]*interp( W(r[j]-r[i]), x, j) )

    document.par(r'the integral of the previous relation is')
    rho = lambda i,j: interp( 1, r[i], j)
    document.eq( f(r[i]), F[i]*V[i]*rho(i,j) )
    
    #document.par(r'let us define,')
    #v = Function('v')
    #document.eq( v(r[i]), interp( V[j], r[i], j) )
    
    document.par(r'substituting back to \eqref{f},')
    document.eq( f(x), interp( f(r[i])/rho(i,j), x, i ) )

    document.par(r'for the special case that $f(x)=1$ the integral of the rhs is the total system volume, which can be written as the sum of all the volumes of the fluid elements')
    document.eq( sum( V[i], i ), sum( 1/rho(i,j), i ) )


    document.par('It can be shown that the only way to compute the element volume with a non-recurrent formula')
    document.eq( Sum( Indexed('V', i)**(-1)*Function('G')(x-Indexed('r', i) ), (i,0,N) ), Sum( Function('W')(x-Indexed('r', i) ), (i,0,N) ) )

    document.par(r'This is the fluid element density $\rho=V^{-1}$')
    document.eq( Function('rho')(x), Sum( Function('W')(x-Indexed('r', i) ), (i,0,N) ) )
    document.par(r'Therefore the expression of the fluid element density ultimately arises from the requirement that the integral of the density does not depend on the kernel (and that is the only explicit definition)')

with openPDF( 'analytical' ) as document:
    document.par('All calculations here are performed automatically SymPy')
    
    write_kernel(document)
    document.section('Kernel')
    document.par('The Kernel function is defined by the following properties. Its integral is')
    document.eq( integrate( W((x - r(t,i)/h)), x ) )
    document.par('Its derivative in respect to $x$')
    document.eq( Derivative( W((x - r(t,i))/h), x ), W( (x - r(t,i))/h ).Diff( x ) ) 
    document.par('Its derivative in respect to $r(t,j)$')
    document.eq( Derivative( W(Abs(x - r(t,i))/h), r(t,j) ), W( (x - r(t,i))/h ).Diff( r(t,j) ) ) 

    document.section('Interpolation')
    document.par('the interpolation reference is')
    document.eq( rho(x), sph_interpolation( None, x, i) )
    document.par('and its derivative is')
    document.eq( Derivative( rho(x), x ), simplify( sph_interpolation( None, x, i).Diff( x ).doit() ) )
    document.par('on the other hand, when taking the variational approach one gets')
    document.eq( Derivative( rho( r(t,k) ), r(t,j) ), simplify( sph_interpolation( None, r(t,k), i).Diff( r(t,j) ).doit() ) )
    document.par('and its derivative is')
    document.eq( Derivative( rho(x), t ), simplify( sph_interpolation( None, x, i ).Diff( t ).doit() ) )
    document.par('and its derivative is')
    document.eq( Derivative( rho(r(t,k)), t ), simplify( sph_interpolation( None, r(t,k), i ).Diff( t ).doit() ) )

    _L = Function('\mathcal{L}')
    L = Function('L')

    document.section('Monaghan recipe')
    document.par('The Lagrangean density is')
    document.eq( _L(x,t), S.Half*rho(x)*v(x)*v(x) - rho(x)*epsilon(rho(x)) )

    LagrangeanDensity = Monaghan( S.Half*rho(x)*dot_r(t,i)**2 - rho(x)*epsilon(rho(x)), x, i )
    document.par('Using the Monaghan recipe it transcribes into with the transformation of $v(x)$ into $v(t,i)$')
    document.eq( _L(x,t), simplify(LagrangeanDensity) )
    Lagrangean = Integral( LagrangeanDensity, x )
    document.par('The integrated Lagrangean density is')
    document.eq( L(x,t), Lagrangean )
    document.par('The variation in respect to r')
    document.eq( Derivative( L(x, t), r(t,j) ), Lagrangean.Diff( r(t,j) ) )
    document.par('The variation in respect to v')
    document.eq( Derivative( L(x, t), dot_r(t,j) ), simplify( Lagrangean.Diff( dot_r(t,j) ).doit() ) )
    document.par('and its time derivative')
    document.eq( simplify( Lagrangean.Diff( dot_r(t,j) ).Diff(t).doit() ) )
    document.par('Finally the Euler-Lagrange relation')
    EL = EulerLagrange_Eq( Lagrangean, r(t,j), dot_r(t,j) )
    document.eq( 0, simplify( EL.doit() ) )
    document.eq( 0, EL.expand().doit() )

    document.section('Beyond Monaghan')
    document.par(r"In Monaghan's recipe, the fluid element velocity is set to be the energy velocity by hand. We will attempt to derive the relation between these two velocities out of the variational principle directly")
    lamb_ = Function('lambda')
    LagrangeanDensity = S.Half*rho(x)*v(x)**2 - rho(x)*epsilon(rho(x)) + lamb_(x)*( rho(x).Diff(t) + (rho(x)*v(x)).Diff(x) )
    document.eq( _L(x,t), simplify(LagrangeanDensity) )
    LagrangeanDensity = S.Half*rho(x)*v(x)**2 - rho(x)*epsilon(rho(x)) + lamb_(x)*( rho(x).Diff(t) + (rho(x)*v(x)).Diff(x) )
    document.par(r"Differential with respect to $v(x)$")
    document.eq( _L(x,t), LagrangeanDensity.Diff( v(x) ) )

    document.par(r"Differential with respect to $\lambda(x)$")
    document.eq( _L(x,t), LagrangeanDensity.Diff( lamb_(x) ) )
    document.par(r"Differential with respect to $\rho(x)$")
    document.eq( _L(x,t), LagrangeanDensity.Diff( rho(x) ) )

    document.par(r"Interpolation of the lagrangean")
    document.par(r"Differential of $\rho$ with respect to $t$")
    rhoDiff = rho( r(t,i) ).Diff(t)
    document.eq( rhoDiff )
    LagrangeanDensityMon = Monaghan( S.Half*rho(x)*v(x)*v(x) - rho(x)*epsilon(rho(x)) + lamb_(x)*(rho(x)*v(x)).Diff(x), x, i )
    document.eq( rhoDiff )
    #print Monaghan( lamb_( x )*, x, i )
    document.eq( _L(x,t), simplify(LagrangeanDensityMon) )
    LagrangeanMon = Integral( LagrangeanDensityMon, x )
    document.par(r"Differential with respect to $r_j(t)$")
    document.eq( _L(x,t), simplify( LagrangeanMon.Diff( r(t,j) ).doit() ) )
