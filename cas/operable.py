from cas.rule import Rule, RuleList

def translate( arg ):
    if isinstance(arg, Operable):
        return arg
    if isinstance(arg, (int, float)):
        if arg < 0:
            return Neg(Number(abs(arg)))
        return Number(arg)
    if isinstance(arg, str):
        return Symbol(arg)
    raise Exception('translate not implemented for', arg.__class__)


_instances = {}
class Operable:
    '''
    this is the basic class that defines the interface with Python
    mainly, it converts Python operators into Functions and inputs
    '''
    def __hash__(self):
        return hash(self.label)
    
    def __new__(cls, label):
        _hash = hash(label)
        if not _hash in _instances.keys():
            new_obj = object().__new__(cls)
            new_obj.label = label
            new_obj.rules = RuleList()
            _instances[_hash] = new_obj
            return new_obj
        return _instances[_hash]

    def __init__(self, label):
        pass
    
    @classmethod
    def instance(cls, obj):
        return isinstance(obj, cls)
    
    @classmethod
    def list(cls):
        print( cls._instances )
    
    def __del__(self):
        print('destroyed obj', str(self), hash(self.label) )
        del _instances[hash(self.label)]
        del self
        
    def __repr__(self):
        return str(self.label)

    def __str__(self):
        return str(self.name)
        
    def __eq__(self, arg):
        return hash(self.label) == hash(translate(arg).label)
    
    def __neg__(self):
        return Neg( self )
    
    def __add__(self, arg):
        return Add( self, translate(arg) )

    def __radd__(self, arg):
        return Add( translate(arg), self)

    def __sub__(self, arg):
        return Add( self, Neg(translate(arg)) ) 

    def __rsub__(self, arg):
        return Add( translate(arg), Neg(self) )
    
    def __mul__(self, arg):
        return Mul( self, translate(arg) )

    def __rmul__(self, arg):
        return Mul( translate(arg), self )

    def __truediv__(self, arg):
        return Mul( self, translate(arg)**(translate(-1)) )

    def __rtruediv__(self, arg):
        return Mul( translate(arg), self**(translate(-1)) )

    def __pow__(self, arg):
        return Pow( self, translate(arg))
    
    def __rpow__(self, arg):
        return Pow( translate(arg), self )
    
    def __call__(self, *args):
        if hasattr(args[0], '__iter__'):
            args = args[0]
        return AppliedFunction( self, *map(translate, args) )
    
    def __getitem__(self, indices):
        return Indexed(self, indices)

    def eval(self, tab=0, verbose=False):
        obj = self
        while True:
            new_obj = Operable.recursive_eval(obj, tab=0, verbose=verbose)
            if new_obj == obj:
                break
            obj = new_obj
        return new_obj
    
    @classmethod
    def recursive_eval(cls, obj, tab=0, verbose=False):
        new_obj = obj.rules( obj, verbose=verbose )
        new_obj = new_obj.__class__.rules( new_obj, verbose=verbose )
        if hasattr(new_obj, 'args'):
            new_args = []
            for arg in new_obj.args:
                arg = Operable.recursive_eval(arg, tab=tab+1, verbose=verbose)
                new_args.append( arg )
            new_obj = new_obj.operator(*new_args)
        return new_obj

from cas.number import Number
from cas.neg import Neg
from cas.add import Add
from cas.mul import Mul
from cas.pow import Pow
from cas.indexed import Indexed
from cas.appliedfunction import AppliedFunction