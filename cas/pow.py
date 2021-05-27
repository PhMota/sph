from cas.rule import Rule, RuleList
from cas.utils import delim
from cas.symbol import Symbol
from cas.expr import Expr
from cas.number import Number
from cas.neg import Neg
from cas.mul import Mul
from cas.add import Add

class Pow(Expr):
    @staticmethod
    def add_order():
        return Mul.add_order()-1
    @staticmethod
    def mul_order():
        return Add.mul_order()-1
    
    rules = RuleList(
        Rule(
            '(x!=0)^0 -> 1',
            lambda obj: obj.base == Number(0) and not obj.power == Number(0),
            lambda obj, v: Number(1)
        ),
        Rule(
            '0^{x!=0} -> 0',
            lambda obj: obj.power == Number(0) and not obj.base == Number(0),
            lambda obj, v: Number(0)
        ),
        Rule(
            'x^1 -> x',
            lambda obj: obj.power == Number(1),
            lambda obj, v: obj.base
        ),
        Rule(
            '1^x -> 1',
            lambda obj: obj.base == Number(1),
            lambda obj, v: Number(1)
        ),
        Rule(
            '(x^a)^b -> x^{a*b}',
            lambda obj: isinstance(obj.base, Pow),
            lambda obj, v: obj.base.base**( obj.base.power*obj.power )
        ),
    )

    def __new__(cls, *args, **kwargs):
        obj = super().__new__(cls, *args)
        if not hasattr(obj, 'base'):
            obj.base = args[0]
            obj.power = args[1]
        return obj
    
    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        base = self.base
        if not isinstance(base, (Symbol, Indexed)):
            base = delim(base)
        power = self.power
        if isinstance(power, Neg):
            if power == -Number(1):
                return fr'\frac{{1}}{{{base}}}'
            return fr'\frac{{1}}{{{base}^{{{power.arg}}} }}'
        return fr'{base}^{{ {power} }}'
    
from cas.indexed import Indexed
