from cas.utils import delim
from cas.rule import Rule, RuleList
from cas.operable import translate
from cas.number import Number
from cas.expr import Expr
from cas.neg import Neg
from cas.add import Add
from cas.mul import Mul
from cas.pow import Pow
from cas.symbol import Symbol
from cas.indexed import Indexed
from cas.appliedfunction import AppliedFunction
from cas.appliedindexedfunction import AppliedIndexedFunction
from cas.kroneckerdelta import KroneckerDelta

def hasindex(obj, index):
    if isinstance(obj, (Indexed, AppliedIndexedFunction)):
        if index in obj.indices:
            return True
    if isinstance(obj, Symbol):
        if obj == index:
            return True
    if not hasattr(obj,'args'):
        return False
    for arg in obj.args:
        if hasindex(arg, index):
            return True
    return False

def update_indices(obj, ind_old, index, verbose=False):
    if isinstance(obj, Indexed):
        if verbose: print( '<update index', index, obj )
        indices = [ index if ind == ind_old else ind for ind in obj.indices ]
        obj = obj.base[indices]
        if verbose: print( '>update index', indices, obj )
        return obj
    if isinstance(obj, (Indexed, AppliedIndexedFunction)):
        if verbose: print( '<update index', index, obj )
        indices = [ index if ind == ind_old else ind for ind in obj.indices ]
        obj = obj.base[indices](*obj.args)
        if verbose: print( '<update index', indices, obj )
        return obj
    if hasattr(obj, 'args'):
        new_args = []
        for arg in obj.args:
            new_args.append( update_indices(arg, ind_old, index) )
        obj = obj.operator(*new_args)
    return obj
    
def rules_summation(obj, v):
    if isinstance(obj.arg, Add):
        if v: print('add summation')
        inside_args = [ arg for arg in obj.arg.args if hasindex(arg, obj.lims[0]) ]
        outside_args = [ arg for arg in obj.arg.args if not arg in inside_args ]
        if len(inside_args) == 0:
            inside_args = [Number(1)]        
        return Add( *[ arg*Summation(Number(1), *obj.lims) for arg in outside_args], *[ Summation( arg, *obj.lims ) for arg in inside_args] )
    
    if isinstance(obj.arg, Neg):
        if v: print('neg summation')
        return Neg(Summation(obj.arg.arg, *obj.lims))
    
    if isinstance(obj.arg, Mul):
        if v: print('mul summation')
        inside_args = [ arg for arg in obj.arg.args if hasindex(arg, obj.lims[0]) ]
        outside_args = [ arg for arg in obj.arg.args if not arg in inside_args ]
        if len(inside_args) == 0:
            inside_args = [Number(1)]
        if KroneckerDelta in list(map(lambda x: x.__class__, inside_args)):
            new_ind = [ ind for delta in inside_args if isinstance(delta, KroneckerDelta) for ind in delta.indices if not ind == obj.lims[0] ][0]
            return Mul( *outside_args, *[ update_indices(arg, obj.lims[0], new_ind, v) for arg in inside_args if not isinstance(arg, KroneckerDelta)] )
        return Mul( *outside_args, Summation(Mul(*[ arg for arg in inside_args]), *obj.lims ) )
    
    if not hasindex(obj.arg, obj.lims[0]):
        if v: print('no index summation')        
        return obj.arg*(obj.lims[2] - obj.lims[1])
    
    if v: print('nothing summation')
    return obj

class Summation(Expr):
    @staticmethod
    def add_order():
        return Number.add_order()+1
    @staticmethod
    def mul_order():
        return AppliedFunction.mul_order()+1
    
    rules = RuleList(
        Rule(
            'summation rule',
            lambda obj: True,
            lambda obj, v: rules_summation(obj, v)
        )
    )
    
    def __new__(cls, arg, *lims, **kwargs):
        lims = list(map(translate, lims))
        if len(lims) < 2:
            lims = (lims[0], -Infinity, Infinity)
        if len(lims) < 3:
            lims = (lims[0], lims[1], Infinity)
        if arg == Number(1):
            return lims[2]-lims[1]
        obj = super().__new__(cls, arg, *lims)
        if not hasattr(obj,'name'):
            obj.arg = arg
            obj.lims = lims
        return obj

    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        lims = self.lims
        arg = self.arg
        if isinstance(arg, Add):
            arg = delim(arg)
        return fr'\sum_{{{lims[0]}={lims[1]}}}^{{{lims[2]}}}{arg}'
