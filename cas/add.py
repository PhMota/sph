from cas.rule import RuleList, Rule
from cas.expr import Expr
from cas.number import Number, getnumbervalues, isnumber, Zero
from cas.utils import delim, flatten_args, areinstance
from collections import Counter


def rules_add(obj, verbose=False):
    args = flatten_args(obj)
    terms = Counter()
    numeric = 0
    for arg in args:
        if isinstance(arg, Neg):
            arg0 = arg.args[0]
            if isinstance(arg0, Number):
                numeric += -arg0.value
            else:
                terms.update({Mul(arg0):-1})
        elif isinstance(arg, Number):
            numeric += arg.value
        elif isinstance(arg, Mul):
            arg0 = arg.args[0]
            if isinstance(arg0, Neg):
                if isinstance(arg0.args[0], Number):
                    terms.update({Mul(*arg.args[1:]): -arg0.args[0].value})
            elif isinstance(arg0, Number):
                terms.update({Mul(*arg.args[1:]): arg0.value})
            else:
                terms.update({arg:1})
        else:
            terms.update({arg:1})
    new_args = [ v*s if v!=1 else s for s, v in terms.items() if v != 0]
    if len(new_args) == 0:
        return Number(numeric)
    if len(new_args) == 1 and numeric == 0:
        return new_args[0]    
    if numeric == 0:
        return Add(*new_args)
    return Add( *new_args, Number(numeric) )

class Add(Expr):
    @staticmethod
    def add_order():
        return Symbol.add_order()-1
    @staticmethod
    def mul_order():
        return Symbol.mul_order()+2
    
    rules = RuleList(
        Rule(
            r'flatten args',
            lambda obj: any(areinstance(obj.args, Add)),
            lambda obj, v: Add( *flatten_args(obj) )
        ),
        Rule(
            r'add numbers',
            lambda obj: any(areinstance(obj.args, Number)),
            lambda obj, v: Add(
                Number(sum(getnumbervalues(obj))),
                *[ arg for arg in obj.args if not isnumber(arg) ]
            )
        ),
        Rule(
            r'add rules',
            lambda obj: True,
            lambda obj, v: rules_add(obj, v)
        ),
    )

    def __new__(cls, *args, **kwargs):
        if len(args) == 0:
            return Zero
        if len(args) == 1:
            return args[0]
        if Zero in args:
            return Add(*[arg for arg in args if not arg == Zero])
        args = sorted( args, key=lambda x: x.add_order() )
        obj = super().__new__(cls, *args)
        return obj
    
    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        s = '+'.join( [ 
            str(arg) if not isinstance(arg, Add) else delim(str(arg)) 
            for arg in self.args 
        ] )
        return s.replace('+-', '-')

from cas.symbol import Symbol
from cas.neg import Neg
from cas.mul import Mul