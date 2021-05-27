from cas.rule import Rule, RuleList
from cas.symbol import Symbol
from cas.number import Number
from cas.indexed import Indexed
from cas.mul import Mul
from cas.appliedfunction import AppliedFunction

class KroneckerDelta(Indexed):
    @staticmethod
    def add_order():
        return Mul.add_order()-1
    @staticmethod
    def mul_order():
        return AppliedFunction.mul_order()+2
    
    rules = RuleList(
        Rule(
            'KroneckerDelta rule',
            lambda obj: obj.indices[0] == obj.indices[1],
            lambda obj, v: Number(1)
        )
    )
    
    def __new__(cls, ind0, ind1, **kwargs):
        if ind0 == ind1:
            return Number(1)
        symbol = Symbol(r'{\delta}')
        obj = super().__new__(cls, symbol, [ind0, ind1])
        return obj
