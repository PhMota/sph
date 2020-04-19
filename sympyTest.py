from sympy import *
from sympy.core.singleton import S

init_printing()

Dummy.init = Dummy.__init__
def new__init__( self, a='', dummy_index='', commutative=True ):
    if a == 'xi_0':
        print 'new Dummy', a, dummy_index, commutative
        #raise Exception('here')
    Dummy.init(self, a='', dummy_index=dummy_index, commutative=commutative)
    
Dummy.__init__ = new__init__

def new_fdiff( self, argindex ):
    """
    Returns the first derivative of the function.
    """
    if not (1 <= argindex <= len(self.args)):
        raise ArgumentIndexError(self, argindex)
    ix = argindex - 1
    A = self.args[ix]
    if A._diff_wrt:
        if len(self.args) == 1:
            return Derivative(self, A)
        if A.is_Symbol:
            for i, v in enumerate(self.args):
                if i != ix and A in v.free_symbols:
                    # it can't be in any other argument's free symbols
                    # issue 8510
                    break
            else:
                return Derivative(self, A)
        else:
            free = A.free_symbols
            for i, a in enumerate(self.args):
                if ix != i and a.free_symbols & free:
                    break
            else:
                # there is no possible interaction bewtween args
                return Derivative(self, A)
    # See issue 4624 and issue 4719, 5600 and 8510
    #D = Dummy('xi_%i' % argindex, dummy_index=hash(A))
    #args = self.args[:ix] + (D,) + self.args[ix + 1:]
    #return Subs(Derivative(self.func(*args), D), D, A)
    print 'no xi here'
    return Derivative( self.func(*self.args), A )


def new_eval_subs(self, old, new):
    # The substitution (old, new) cannot be done inside
    # Derivative(expr, vars) for a variety of reasons
    # as handled below.
    if old in self._wrt_variables:
        # first handle the counts
        expr = self.func(self.expr, *[(v, c.subs(old, new))
            for v, c in self.variable_count])
        if expr != self:
            return expr._eval_subs(old, new)
        # quick exit case
        if not getattr(new, '_diff_wrt', False):
            # case (0): new is not a valid variable of
            # differentiation
            if isinstance(old, Symbol):
                # don't introduce a new symbol if the old will do
                return Subs(self, old, new)
            else:
                xi = Dummy('xi')
                return Subs(self.xreplace({old: xi}), xi, new)

    # If both are Derivatives with the same expr, check if old is
    # equivalent to self or if old is a subderivative of self.
    if old.is_Derivative and old.expr == self.expr:
        if self.canonical == old.canonical:
            return new

        # collections.Counter doesn't have __le__
        def _subset(a, b):
            return all((a[i] <= b[i]) == True for i in a)

        old_vars = Counter(dict(reversed(old.variable_count)))
        self_vars = Counter(dict(reversed(self.variable_count)))
        if _subset(old_vars, self_vars):
            return Derivative(new, *(self_vars - old_vars).items()).canonical

    args = list(self.args)
    newargs = list(x._subs(old, new) for x in args)
    if args[0] == old:
        # complete replacement of self.expr
        # we already checked that the new is valid so we know
        # it won't be a problem should it appear in variables
        return Derivative(*newargs)

    if newargs[0] != args[0]:
        # case (1) can't change expr by introducing something that is in
        # the _wrt_variables if it was already in the expr
        # e.g.
        # for Derivative(f(x, g(y)), y), x cannot be replaced with
        # anything that has y in it; for f(g(x), g(y)).diff(g(y))
        # g(x) cannot be replaced with anything that has g(y)
        syms = {vi: Dummy() for vi in self._wrt_variables
            if not vi.is_Symbol}
        wrt = set(syms.get(vi, vi) for vi in self._wrt_variables)
        forbidden = args[0].xreplace(syms).free_symbols & wrt
        nfree = new.xreplace(syms).free_symbols
        ofree = old.xreplace(syms).free_symbols
        if (nfree - ofree) & forbidden:
            return Subs(self, old, new)

    viter = ((i, j) for ((i, _), (j, _)) in zip(newargs[1:], args[1:]))
    if any(i != j for i, j in viter):  # a wrt-variable change
        # case (2) can't change vars by introducing a variable
        # that is contained in expr, e.g.
        # for Derivative(f(z, g(h(x), y)), y), y cannot be changed to
        # x, h(x), or g(h(x), y)
        for a in _atomic(self.expr, recursive=True):
            for i in range(1, len(newargs)):
                vi, _ = newargs[i]
                if a == vi and vi != args[i][0]:
                    return Subs(self, old, new)
        # more arg-wise checks
        vc = newargs[1:]
        oldv = self._wrt_variables
        newe = self.expr
        subs = []
        for i, (vi, ci) in enumerate(vc):
            if not vi._diff_wrt:
                # case (3) invalid differentiation expression so
                # create a replacement dummy
                print 'created a xi here!'
                xi = Dummy('xi_%i' % i)
                # replace the old valid variable with the dummy
                # in the expression
                newe = newe.xreplace({oldv[i]: xi})
                # and replace the bad variable with the dummy
                vc[i] = (xi, ci)
                # and record the dummy with the new (invalid)
                # differentiation expression
                subs.append((xi, vi))

        if subs:
            # handle any residual substitution in the expression
            newe = newe._subs(old, new)
            # return the Subs-wrapped derivative
            return Subs(Derivative(newe, *vc), *zip(*subs))

    # everything was ok
    return Derivative(*newargs)

Function.fdiff = new_fdiff

Function.__eval_subs = Function._eval_subs
def new_eval_subs(self, old, new):
    print 'new_eval_subs'
    return Function.__eval_subs(self, old, new)

Function._eval_subs = new_eval_subs

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

U = Function('U')
rho = Function('rho')
v = Function('v')
r = Function('r')

class r(Function):
    _diff_wrt = True

class U(Function):
    #def __init__(self, arg):
        #for key, func in Function.__dict__.items():
            #def _(self, **kwargs):
                #print key, kwargs
                #return Function.__dict__[key](self,**kwargs)
            #self.__dict__[key] = _
        #Function.__init__(self, arg )
    def fdiff(self, argindex=1):
        print 'U.fdiff', argindex
        ix = argindex - 1
        A = self.args[ix]
        return Derivative( self.func(*self.args), A )

class W(Function):
    _diff_wrt = True
    def fdiff( self, argindex ):
        return Function('W_%s' % argindex )(*self.args)
    
    def _eval_Integral(self, arg):
        print 'W._eval_Integral', arg
        return S.One

    def diff( self, var ):
        print 'W.diff', 'args=', self.args, 'self=', self, 'var=', var
        if isinstance( var, Dummy ):
            print 'Dummy found in W.diff'
            return Derivative( self.func( *self.args ), var )
        pprint( self.args[0], var )
        return self.fdiff(1) * self.args[0].diff( var )

class r(Function):
    _diff_wrt = True
    def fdiff(self, argindex=1):
        if argindex == 1:
            return Function('v')
        raiseExepction( 'cannot derivate arg %s' % argindex )
    
    def diff( self, var ):
        #print 'r.diff', var, type(var)
        if isinstance(var, r):
            #print 'is r'
            if self.args[1] == var.args[1]:
                return S.One
            return KroneckerDelta( self.args[1], var.args[1] )
        #print 'r.diff last return'
        return self.fdiff(1)(*self.args) * self.args[0].diff(var)

print '*** kernel derivatives ***'
pprint( r(t,i).diff( r(t,j) ) )
pprint( ((x - r(t,i))/h).diff(x) )
pprint( ((x - r(t,i))/h).diff( r(t,j) ) )

pprint( W( (x - r(t,i))/h ) )
pprint( W( (x - r(t,i))/h ).diff(x) )
pprint( W( (x - r(t,i))/h ).diff( r(t,j) ) )
print '*** done kernel ***'
exit(0)

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
    _diff_wrt = True
    def explicit( self ):
        print 'rho.expand'
        arg = self.args[0]
        index = Idx( '%s\'' % arg.args[-1], N )
        
        print arg, index
        return sph_interpolation( None, arg, index)
    
    def fdiff(self, argindex ):
        return Function('rho_%s' % argindex )(*self.args)
    
    def diff(self, var):
        if isinstance(var, Dummy):
            return Derivative( self.func(*self.args), var )
        return self.fdiff(1) * self.args[0].diff( var )

pprint( rho(x) )
pprint( rho(r(t,j)) )
pprint( sph_interpolation( None, r(t,j), i) )
print 'done here'
#pprint( rho(r(t,j)).explicit() )
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
    
    def diff(self, var):
        print 'Sum.diff', var
        pprint( self.args[0], var )
        pprint( self.args[0].diff(var) )
        return Sum( self.args[0].diff(var), self.args[1] )
    
class EulerLagrange_Eq(Function):
    @classmethod
    def eval(cls, arg, var, dvar):
        print 'EL.eval', cls, arg, var, dvar
        eq = arg.diff( var )
        
        return eq

#term = Sum( rho(r(t,i)).diff(r(t,j)), (i,0,N)).doit()

class epsilon(Function):

    def fdiff( self, argindex ):
        return Function('epsilon_%s' % argindex )(*self.args)
    
    def diff( self, var ):
        print 'epsilon.diff', 'args=', self.args, 'self=', self, 'var=', var
        if isinstance(var, Dummy):
            return Derivative( self.func(*self.args), var )
        return self.fdiff(1) * self.args[0].diff(var)

print
print '*** we are here!!! ***'
pprint( epsilon( r(t,i) ).diff( r(t,j) ) )
pprint( epsilon( rho( rho( r(t,i) ) ) ).diff( r(t,j) ) )
pprint( ( epsilon( x )/rho( x ) ).diff( x ) )
pprint( ( epsilon( r(t,i) )/rho( r(t,i) ) ).diff( r(t,j) ) )

print type( epsilon( x )/rho( x ) )
#exit(0)

print 'de/dt'
pprint( epsilon( r(t,i) ).diff(t) )
pprint( (epsilon( r(t,i) )/rho(r(t,i))).diff(t) )

print 'de/dt'
pprint( diff( epsilon( r(t,i) ), t ) )
print 'done'

LagrangeanDensity = Monaghan( epsilon(x), x, i )
print '####LagrangeanDensity'
pprint( LagrangeanDensity )
Lagrangean = Integral( LagrangeanDensity, x )
print '###########Lagrangean'
pprint( Lagrangean )
pprint( Lagrangean.diff( r(t,j) ) )
print '########EulerLagrange'
pprint( EulerLagrange_Eq( Lagrangean, r(t,j), v(t,j) ) )

#print '########d(r[i], r[j])'
#pprint( r(t,j).diff( r(t,i) ) )
exit(0)
EulerLagrange_Eq( LagrangeanDensity, r, a )
