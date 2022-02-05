import types,string
from re import escape

def SADFormat(x):
    typ=type(x)
    if typ == types.TupleType or typ == types.ListType:
        return "{%s}"%string.join(map(SADFormat, x),",\n")
    elif typ == types.StringType:
        return '"%s"'%escape(x)
    elif typ == types.FloatType or \
         typ == types.IntType or typ == types.IntType:
        return str(x)
    else:
        return str(x)
    
UTCEPOCHSAD=2208988800L

