from cas.rule import RuleList, Rule
from cas.number import Number
from cas.symbol import Symbol
from cas.appliedfunction import AppliedFunction

def rules_diracdelta(obj, v):
    if obj.arg == Number(0):
        return Number(1)
    return obj

class DiracDelta(AppliedFunction):
    rules = RuleList(
        Rule(
            'delta dirac zero',
            lambda obj: obj.arg == Number(0),
            lambda obj, v: Number(1)
        )
    )
    
    def __new__(cls, arg, **kwargs):
        if arg == Number(0):
            return Number(1)
        symbol = Symbol(r'{\delta}')
        obj = super().__new__(cls, symbol, arg)
        if not hasattr(obj, 'arg'):
            obj.arg = arg
        return obj

