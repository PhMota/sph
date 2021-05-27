import cas.utils as utils

class Rule:
    def __init__(self, name, condition, func):
        self.name = name
        self.condition = condition
        self.func = func
        
    def __call__(self, arg, verbose=False):
        if self.condition(arg):
            return self.func(arg, verbose)
        return arg
    
    def __str__(self):
        return name

class RuleList:
    def __init__(self, *rules):
        self.rules = list(rules)
    
    def empty(self):
        return len(self.rules) == 0
    
    def append(self, *args):
        rule = None
        if not isinstance(args[0], Rule):
            rule = Rule(*args)
        else:
            rule = args[0]
        self.rules.append( rule )
        return self
    
    def __call__(self, arg, verbose=False):
        for rule in self.rules:
            new_arg = rule( arg, verbose=verbose )
            if not new_arg == arg:
                if verbose: 
                    print( rule.name, ':')
                    print( repr(arg), '->', repr(new_arg) )
                    utils.pprint( '\quad', arg, r'\rightarrow', new_arg )                    
                return new_arg
            arg = new_arg
        return arg