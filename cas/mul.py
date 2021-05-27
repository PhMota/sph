from cas.rule import RuleList, Rule
from cas.expr import Expr
from cas.number import Number, getnumbervalues, Zero, One
from cas.utils import flatten_args, areinstance, delim
from collections import Counter


def rules_mul(obj, verbose=False):
    args = flatten_args(obj)    
    terms = Counter()
    numeric = 1
    for arg in args:
        if isinstance(arg, Add):
            return Add(*[ Mul(add_arg, *[ other for other in args if not other == arg ]) for add_arg in arg.args ])
        if isinstance(arg, Neg):
            iarg = arg.args[0]
            numeric *= -1
            if isinstance(iarg, Number):
                numeric *= iarg.value
            else:
                terms.update({iarg:1})
        elif isinstance(arg, Number):
            numeric *= arg.value
        elif isinstance(arg, Pow):
            if isinstance(arg.power, Neg):
                if isinstance(arg.power.arg, Number):
                    terms.update({arg.base: -arg.power.arg.value})
            elif isinstance(arg.power, Number):
                terms.update({arg.base: arg.power.value})
            else:
                terms.update({arg.base: arg.power})
        else:
            terms.update({arg:1})
    new_args = [ s**v if v != 1 else s for s, v in terms.items() if v != 0]
    if numeric == 0:
        return Zero
    if len(new_args) == 0:
        return Number(numeric)
    if len(args) == 1 and numeric == 1:
        return new_args[0]
    if len(args) == 1 and numeric == -1:
        return Neg(new_args[0])
    if numeric == 1:
        return Mul( *new_args )
    if numeric == -1:
        return Neg(Mul( *new_args ))
    ret = Mul( Number(numeric), *new_args )
    return ret

class Mul(Expr):
    @staticmethod
    def add_order():
        return Symbol.add_order()+2
    @staticmethod
    def mul_order():
        return Symbol.mul_order()+1
    
    rules = RuleList(
        Rule(
            r'mul rules',
            lambda obj: True,
            lambda obj, v: rules_mul(obj, v)
        ),
    )

    def __new__(cls, *args, **kwargs):
        if len(args) == 0:
            return Zero
        if len(args) == 1:
            return args[0]
        if Zero in args:
            return Zero
        if One in args:
            return Mul(*[arg for arg in args if not arg == One])
        args = sorted( args, key=lambda x: x.mul_order() )
        obj = super().__new__(cls, *args )
        return obj
    
    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        s = ''
        for arg in self.args:
            if isinstance(arg, (Add, Neg, Mul)):
                s += delim(str(arg))
            else:
                s += str(arg)
        return s

from cas.symbol import Symbol
from cas.number import Number
from cas.neg import Neg
from cas.add import Add
from cas.pow import Pow