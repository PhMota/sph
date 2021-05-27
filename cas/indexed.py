from cas.rule import RuleList
from cas.operable import Operable, translate
from cas.symbol import Symbol

class Indexed(Operable):
    @staticmethod
    def add_order():
        return Symbol.add_order()+.5
    @staticmethod
    def mul_order():
        return Symbol.mul_order()+.5

    rules = RuleList()
    def __new__(cls, base, indices, **kwargs):
        if not hasattr(indices, '__iter__'):
            indices = [indices]
        label = fr'{base}_{{{" ".join(map(str,indices))}}}'
        name = fr'{base.name}_{{{" ".join(map(str,indices))}}}'
        new_obj = super().__new__( cls, label )
        if not hasattr(new_obj, 'name'):
            new_obj.base = base
            new_obj.name = name
            new_obj.label = label
            new_obj.indices = indices
            new_obj.rules = RuleList()
        return new_obj
    
    def __init__(self, *args, **kwargs):
        pass
    
    def __call__(self, *args):
        return AppliedIndexedFunction( self, *map(translate, args) )
    
from cas.appliedindexedfunction import AppliedIndexedFunction