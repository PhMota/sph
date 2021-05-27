from cas.appliedfunction import AppliedFunction

class AppliedIndexedFunction(AppliedFunction):
    def __new__(cls, indexed, *args, **kwargs):
        name = indexed.name
        obj = super().__new__(cls, indexed, *args)
        if not hasattr(obj, 'indices'):
            obj.base = indexed.base
            obj.indices = indexed.indices
        return obj

