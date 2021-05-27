from cas.utils import delim
from cas.rule import Rule, RuleList
from cas.expr import Expr
from cas.symbol import Symbol
from cas.number import Number
from cas.indexed import Indexed
from cas.neg import Neg
from cas.add import Add
from cas.mul import Mul
from cas.pow import Pow
from cas.appliedfunction import AppliedFunction
from cas.appliedindexedfunction import AppliedIndexedFunction
from cas.summation import Summation
from cas.kroneckerdelta import KroneckerDelta
from cas.diracdelta import DiracDelta

def rules_derivative(obj, v):
    if isinstance(obj.arg, Indexed):
        if isinstance(obj.wrt, Indexed):
            if obj.arg.base == obj.wrt.base:
                if v:
                    print('kroneckerDelta', obj.arg.indices, obj.wrt.indices)
                return Mul( *[ KroneckerDelta(ind0, ind1) for ind0, ind1 in zip(obj.arg.indices, obj.wrt.indices)] )

    if isinstance(obj.arg, AppliedIndexedFunction):
        if isinstance(obj.wrt, AppliedIndexedFunction):
            if obj.arg.base == obj.wrt.base:
                if v:
                    print('kroneckerDelta+dirac', obj.arg.indices, obj.wrt.indices)
                return Mul( *[ KroneckerDelta(ind0, ind1) for ind0, ind1 in zip(obj.arg.indices, obj.wrt.indices)], *[ DiracDelta(arg0 - arg1) for arg0, arg1 in zip(obj.arg.args, obj.wrt.args)] )
    
#     if isinstance(obj.arg, Neg):
#         if v: print('neg deriv')
#         return Neg( Derivative(obj.arg.arg, obj.wrt) )
#     if isinstance(obj.arg, Summation):
#         if v: print('sum deriv')
#         return Summation( Derivative(obj.arg.arg, obj.wrt), *obj.arg.lims )
#     if isinstance(obj.arg, Add):
#         if v: print('add deriv')
#         return Add( *[Derivative(arg, obj.wrt) for arg in obj.arg.args] )
    if isinstance(obj.arg, Mul):
        if v: print('mul deriv')        
        return Add( *[ 
            Mul(*[ 
                jarg if i!=j else Derivative(jarg, obj.wrt) for j, jarg in enumerate(obj.arg.args) 
            ] )
            for i, iarg in enumerate(obj.arg.args)
        ] )
    if isinstance(obj.arg, Pow):
        if v: print('pow deriv')
        return obj.arg.power * obj.arg.base**(obj.arg.power-1) * Derivative(obj.arg.base, obj.wrt)

    if isinstance(obj.arg, AppliedFunction):
        if v: print('func deriv')
        if isinstance(obj.wrt, AppliedFunction):
            if v: print('func deriv')
            if obj.arg.operator == obj.wrt.operator:
                if v: print('diracDelta deriv')
                return Mul( *[ DiracDelta( arg0 - arg1 ) for arg0, arg1 in zip(obj.arg.args, obj.wrt.args) ])
        return Add(*[ PartialDerivative(obj.arg.operator, i)(*obj.arg.args) * Derivative(arg, obj.wrt) for i, arg in enumerate(obj.arg.args) ])
    
    if v: print('nothing deriv')
    return obj

class Derivative(Expr):
    @staticmethod
    def add_order():
        return Mul.add_order()-1
    @staticmethod
    def mul_order():
        return AppliedFunction.mul_order()+1
    
    rules = RuleList(
        Rule(
            'd(-...) -> -d(...)',
            lambda obj: isinstance(obj.arg, Neg),
            lambda obj, v: Neg( Derivative(obj.arg.arg, obj.wrt) )
        ),
        Rule(
            'd(sum(...)) -> sum(d(...))',
            lambda obj: isinstance(obj.arg, Summation),
            lambda obj, v: Summation( Derivative(obj.arg.arg, obj.wrt), *obj.arg.lims )
        ),
        Rule(
            'd(add(...)) -> add(d(...))',
            lambda obj: isinstance(obj.arg, Add),
            lambda obj, v: Add( *[Derivative(arg, obj.wrt) for arg in obj.arg.args] )
        ),
        Rule(
            'derivative rule',
            lambda obj: True,
            lambda obj, v: rules_derivative(obj, v)
        )
    )
    
    def __new__(cls, *args, **kwargs):
        arg, wrt, = args
        if isinstance(wrt, Number):
            raise Exception('trying to derivate with respect to a number')
        if isinstance(arg, Number):
            return Number(0)
        if arg == wrt:
            return Number(1)
        obj = super().__new__(cls, *args)
        if not hasattr(obj, 'wrt'):
            obj.arg = arg
            obj.wrt = wrt
        return obj

    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        wrt = self.wrt
        if not isinstance(wrt, Symbol):
            wrt = delim(wrt)
        if isinstance(self.arg, Expr):
            return fr'\frac{{ {{\rm d}} }}{{ {{\rm d}}{self.wrt} }}{delim(self.arg)}'
        return fr'\frac{{ {{\rm d}}{self.arg} }}{{ {{\rm d}}{self.wrt} }}'

    
class PartialDerivative(Expr):
    def __new__(cls, symbol, index, **kwargs):
        name = fr'\partial_{{{index}}}{symbol}'
        obj = super().__new__(cls, name, index)
        if not hasattr(obj, 'index'):
            obj.arg = symbol
            obj.index = index
            obj.name = name
        return obj

    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        return delim(self.name)
    