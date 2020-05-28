from __future__ import print_function
import datetime

class Timer:
    def __init__(self, msg='elapsed'):
        self.msg = msg
        
    def __enter__(self):
        self.start = datetime.datetime.now()
        return self
    
    def __exit__(self, type, value, traceback):
        s = [ self.msg ]
        s += [ self.time() ]
        if value:
            s += [ '%s'%type ]
            s += [ '%s'%value ]
        print( ' '.join(s) )
    
    def seconds(self, N=1):
        return (datetime.datetime.now() - self.start).total_seconds()*N
    
    def time(self, N=1):
        t = self.seconds(N)
        if t > 60*60:
            s = '%sh%sm%ss%sms' %( int(t)/60/60, int(t/60)%60, int(t)%60, int(t*1e3)%1000 )
        else: 
            if t > 60: s = '%sm%ss%sms' %( int(t)/60, int(t)%60, int(t*1e3)%1000 )
            else: s = '%ss%sms' % (int(t), int(t*1e3)%1000 )
        return s

