from IPython.display import display, Markdown, Latex, Math

def pprint(*s):
    _s = '\quad '.join( map(str,s) )
    return display(Math(_s))

def eprint(s, verbose=False):
    pprint( f'{s} = {s.eval(verbose=verbose)}' )

def delim(s):
    return fr'\left( {s}\right)'

def areinstance(objs, classinfo):
    return [ isinstance(obj, classinfo) for obj in objs ]
    
def flatten_args( obj ):
    new_args = []
    for arg in obj.args:
        if isinstance(arg, obj.__class__):
            new_args.extend( flatten_args(arg) )
        else:
            new_args.append( arg )
    new_args = sorted( new_args, key=lambda x: hash(x.label) )
    new_args = sorted( new_args, key=lambda x: str(x) )
    return new_args

