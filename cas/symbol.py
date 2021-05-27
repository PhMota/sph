from cas.rule import RuleList
from cas.operable import Operable

class Symbol(Operable):
    @staticmethod
    def add_order():
        return 0
    @staticmethod
    def mul_order():
        return 0
    
    rules = RuleList()
    def __new__(cls, label):
        name = label
        label = fr'{cls.__name__}({repr(label)})'        
        new_obj = super().__new__(cls, label)
        if not hasattr(new_obj, 'name'):
            new_obj.name = name
            new_obj.rules = RuleList()
        return new_obj
    
    def __init__(self, label):
        pass

Infinity = Symbol(r'{\infty}')
