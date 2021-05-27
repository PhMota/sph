from cas.rule import RuleList
from cas.operable import Operable

class Expr(Operable):
    def __new__(cls, *args):
        label = fr'{cls.__name__}({", ".join(map(repr, args))})'
        obj = super().__new__(cls, label)
        if not hasattr(obj, 'operator'):
            obj.operator = cls
            obj.args = args
            obj.rules = RuleList()
        return obj
    
    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        return self.label

