from cas.rule import RuleList, Rule
from cas.expr import Expr
from cas.utils import delim

class Neg(Expr):
    def add_order(self):
        return self.arg.add_order()+.1
    def mul_order(self):
        return self.arg.mul_order()
    
    rules = RuleList(
#         Rule(
#             r'negative number',
#             lambda obj: isinstance(obj.arg, Number),
#             lambda obj, v: NegNumber(obj.arg.value)
#         ),
        Rule(
            r'double negative',
            lambda obj: isinstance(obj.arg, Neg),
            lambda obj, v: obj.arg.arg
        )
    )

    def __new__(cls, arg):
        obj = super().__new__(cls, arg)
        if not hasattr(obj, 'arg'):
            obj.arg = arg
        return obj
    
    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        arg = self.arg
        s = str(arg)
        if not isinstance(arg, Symbol):
            if not hasattr(arg, 'base'):
                s = delim(s)
        return fr'-{s}'

from cas.symbol import Symbol