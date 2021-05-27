from cas.rule import Rule, RuleList
from cas.number import Number
from cas.expr import Expr
from cas.symbol import Symbol
from cas.add import Add
from cas.mul import Mul
from cas.pow import Pow

class AppliedFunction(Expr):
    @staticmethod
    def add_order():
        return Mul.add_order()+1
    @staticmethod
    def mul_order():
        return Pow.mul_order()+1
    
    rules = RuleList()
    
    def __new__(cls, symbol, *args, **kwargs):
        name = symbol.name
        obj = super().__new__(cls, symbol.name, *args)
        if not hasattr(obj,'name'):
            obj.operator = symbol
            obj.args = args
            obj.name = name
        return obj

    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        return fr'{self.operator}({", ".join(map(str, self.args))})'
    
    def __lshift__(self, arg):
        def _(obj, v):
            print( arg.label )
            print( self.args[0].label )
            print( obj.args[0].label )
            string = arg.label.replace(self.args[0].label, obj.args[0].label)
            print( 'string', string )
            return eval( arg.label.replace(self.args[0].label, obj.args[0].label) )
        
        rule = Rule(
            fr'{self} := {arg}',
            lambda obj: obj.operator == self.operator,
            lambda obj, v: _(obj, v)
        )
        self.__class__.rules.append(rule)
        return
    
