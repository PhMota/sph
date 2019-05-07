import numpy as np
#cimport numpy as np

DTYPE = np.double

#ctypedef np.double_t DTYPE_t

def Kahan( a, axis=0 ):
    '''Kahan summation of the numpy array along an axis.
    '''
    #cdef list shape = a.shape[:axis]
    #shape += a.shape[axis+1:]
    s = np.zeros( a.shape[:axis] + a.shape[axis+1:], dtype = DTYPE )
    #cdef np.ndarray s = np.zeros( shape, dtype = DTYPE )
    c = np.zeros( s.shape, dtype = DTYPE )
    for i in range(a.shape[axis]):
        # https://stackoverflow.com/a/42817610/353337
        y = a[(slice(None),) * axis + (i,)] - c
        t = s + y
        c = (t - s) - y
        s = t.copy()
    return s

def Neumaier( a, axis=0 ):
    '''
    Neumaier summation of the numpy array along an axis.
    '''
    s = np.zeros(a.shape[:axis] + a.shape[axis+1:])
    c = np.zeros(s.shape)
    for i in range(a.shape[axis]):
        y = a[(slice(None),) * axis + (i,)]
        t = s + y
        c[ np.abs(s) >= np.abs(y) ] += ( (s - t) + y )[ np.abs(s) >= np.abs(y) ]
        c[ np.abs(s) < np.abs(y) ] += ( (y - t) + s )[ np.abs(s) < np.abs(y) ]
        s = t.copy()
    return s + c

def sum_2(a, axis=0, method = None ):
    ap = a.copy()
    ap[ ap < 0 ] = 0
    an = a.copy()
    an[ an > 0 ] = 0
    if method == 'np_sum':
        method = np.sum
    if method == 'Neumaier':
        method = Neumaier
    return method(ap, axis=axis) + method(an, axis=axis)
