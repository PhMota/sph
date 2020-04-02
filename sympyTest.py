from sympy import *
from sympy.core.singleton import S

init_printing()
N = Symbol('N', integer=True)
i = Idx('i', N+1)
j = Idx('j', N+1)
k = Idx('k', N+1)

print type(i)

x = Symbol('x')
y = Symbol('y')
t = Symbol('t')
h = Symbol('h')
nu = Function('nu')
u = Function('u')

class W(Function):
    def fdiff(self, argindex=1):
        if argindex==1:
            return Function('W\'')(*self.args)
        else:
            return 0
    
    def _eval_Integral(self, arg):
        print 'W._eval_Integral', arg
        return S.One

r_t = Function('r_t')
class r(Function):
    def fdiff(self, argindex=1):
        if argindex==1:
            return r_t(*self.args)
        else:
            return 0
    def diff(self, var):
        print 'r.diff', type(var)
        if isinstance(var, r):
            print 'is r'
            if self.args[1] == var.args[1]:
                return 1
            return KroneckerDelta( self.args[1], var.args[1] )
        return Function.diff(self, var)

class v(Function):
    #@classmethod
    def as_base_exp(self):
        #print 'v.as_base_exp'
        return self, S.One

class sph_interpolation(Function):
    @classmethod
    def eval( cls, arg, var, index ):
        print 'sph_interpolation', arg, type(arg)
        i = symbols('i', cls=Idx)
        
        kernel = nu(index)*W( (var-r(t, index))/h )
        range = (index, 0, N)
        if arg is None:
            return Sum( kernel, range )
        
        elif isinstance(arg, FunctionClass):
            return Sum( arg.func(t, index)*kernel, range)
        
        return Sum( arg*kernel, range)
        print 'unpredicted instance of arg', type(arg)
    
class rho(Function):
    def explicit( self ):
        print 'rho.expand'
        arg = self.args[0]
        index = Idx( '%s\'' % arg.args[-1], N )
        
        print arg, index
        return sph_interpolation( None, arg, index)
    
    def fdiff(self, var):
        return Function('rho_x')(*self.args)
    #def diff(self, var):
        #return Function('rho_x')(var) * self.args[0].diff( var )

pprint( rho(x) )
pprint( rho(r(t,j)) )
#exit(0)

class Monaghan(Function):
    @classmethod
    def eval(cls, arg, var, index):
        innerIndex = symbols('%s\'' % index.label, cls=Idx)
        if arg is None:
            return sph_interpolation( 1/rho(r(t,index)), var, index )
        pprint( arg )
        arg = arg.replace( x, r(t,index) )
        pprint( arg )
        return sph_interpolation( arg/rho(r(t,index)), var, index )
    
def explicit( expr, level=0 ):
    if hasattr( expr, 'explicit' ):
        newexpr = expr.explicit()
        return newexpr
    if len(expr.args) == 0:
        return expr
        
    newargs = tuple( explicit(arg, level=level+1) for arg in expr.args )
    return expr.func(*newargs)

#pprint( Monaghan(None,x,i) )
#pprint( Monaghan(rho(x),x,i) )
#pprint( Monaghan(v(t,x),x,i) )
#pprint( explicit( Monaghan(v(t,x),x,i) ) )
#exit(0)

#class D4(Function):
    #@classmethod
    #def eval(cls, arg):
        #print 'SPH.eval', cls, arg
        #return diff( arg, t ) - diff( arg, x )

#class Dmat(Function):
    #@classmethod
    #def eval(cls, arg):
        #print 'eval', cls, arg
        #return diff( arg, t ) + diff( v*arg, x )
    #@classmethod
    #def as_base_exp(cls, arg=None):
        #print 'as_base_exp', cls, arg

#class Dmat4(Function):
    #@classmethod
    #def eval(cls, arg):
        #print 'SPH.eval', cls, arg
        #return diff( gamma*arg, t ) + diff( v*arg, x )
    
class Conservation_Eq(Function):
    @classmethod
    def eval(cls, arg):
        print 'SPH.eval', cls, arg
        eq = Dmat( arg )
        return solve( eq, v(t, i) )

class Sum(Sum):
    def _eval_Integral(self, arg):
        print 'sum _eval_Integral', arg, self.function
        return Sum( Integral( self.function, arg ).doit(), (i,0,N) )
    
    #def diff(self, var):
        #print 'Sum.diff'
        #pprint( self.args[0], var )
        #pprint( self.args[0].diff(var) )
        #return Sum( diff( self.args[0], var), self.args[1] )
    
class EulerLagrange_Eq(Function):
    @classmethod
    def eval(cls, arg, var, dvar):
        print 'EL.eval', cls, arg, var, dvar
        eq = arg.diff( var ) #- arg.diff( dvar ).Sum( self.args[0].diff(), self.args[1] )diff( t )
        #eq = Derivative( arg, var )
        
        return eq

term = Sum( rho(r(t,i)).diff(r(t,j)), (i,0,N)).doit()

class epsilon(Function):
    #pass
    #def _eval_derivative(self, var):
        ##print 'epsilon._eval_derivative', self.args, self, var
        #return Function('epsilon_t')(var, self.args[1]) * self.args[0].diff( var ) \
            #+ Function('epsilon_x')(self.args[0],var) * self.args[1].diff( var )

    def fdiff(self, var):
        print 'epsilon.fdiff', self.args
        return Function('epsilon_x')(*self.args)
    #def diff(self, var):
        #print 'epsilon.diff', self.args, self, var
        #return Function('epsilon_t')(var, self.args[1]) * self.args[0].diff( var ) \
            #+ Function('epsilon_x')(self.args[0],var) * self.args[1].diff( var )

print 'epsilon.diff'
pprint( epsilon( r(t,i) ).diff(r(t,j)) )
pprint( diff( epsilon( r(t,i) ), r(t,j)) )
pprint( epsilon( rho(x) ).diff(x) )

print 'de/dt'
pprint( epsilon( r(t,i) ).diff(t) )
print 'de/dt'
pprint( diff( epsilon( r(t,i) ), t ) )
print 'done'

exit(0)

LagrangeanDensity = Monaghan( epsilon(x), x, i )
print 'LagrangeanDensity'
pprint( LagrangeanDensity )
Lagrangean = Integral( LagrangeanDensity, x )
print 'Lagrangean'
pprint( Lagrangean )
print 'EulerLagrange'
pprint( EulerLagrange_Eq( Lagrangean, r(t,j), r_t(t,j) ) )
exit(0)
EulerLagrange_Eq( LagrangeanDensity, r, a )
