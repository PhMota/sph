from cas.symbol import Symbol

class Number(Symbol):
    @staticmethod
    def add_order():
        return Symbol.add_order()+1
    @staticmethod
    def mul_order():
        return Symbol.mul_order()-1
    
    def __new__(cls, label):
        if isinstance(label, str):
            label = int(label)
        obj = super().__new__(cls, repr(label))
        if not hasattr(obj, 'value'):
            obj.value = label
        return obj

    def __init__(self, label):
        pass

Zero = Number(0)
One = Number(1)

def isnumber(obj):
    return isinstance(obj, Number) or (isinstance(obj, Neg) and isinstance(obj.arg, Number))

def getnumbervalues(obj):
    positive = [ arg.value for arg in obj.args if isinstance(arg, Number) ]
    negative = [ -arg.arg.value for arg in obj.args if isinstance(arg, Neg) and isinstance(arg.arg, Number) ]
    return [*positive, *negative]

from cas.neg import Neg
