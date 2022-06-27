import sys,time,inspect,thread,threading,os,errno


from  geom_filter import DoubleArray


##def handle_stackframe_without_leak():
##  frame = inspect.currentframe()
##  print "stack:"
##  try:      
##    while frame!=None:
##      print inspect.getframeinfo(frame)
##      frame=frame.f_back
##  finally:      
##    del frame

        
class VerbLock:  
  """ to debug deadlocks..."""
  
  def __init__(self,lock):      
    self.lock=lock
  def acquire(self,*x):      
    sys.stdout.write("acquire...")
    handle_stackframe_without_leak()
    sys.stdout.write("<")    
    sys.stdout.flush()
    self.lock.acquire(*x)
    sys.stdout.write(">\n")
    
  def release(self):
    sys.stdout.write("release<")
    sys.stdout.flush()
    self.lock.release()
    sys.stdout.write(">\n")

class State:
  """ state of a sequence of operations. Default implementation just displays
  stage on stdout"""

  def next(self,name,**kwargs):
    " step to a stage called name "
    print " "*70+"\rstage "+name,
    for k,v in kwargs.iteritems():
      print "%s=%s"%(k,v),
    print

  def set_frac(self,frac):
    " set fraction completed for current stage "
    n=int(frac*70)
    sys.stdout.write("+"*n+"-"*(70-n)+"\r")
    sys.stdout.flush()

class DummyState:
  def next(self,name,**kwargs):
    pass
  set_frac=None

    
class ThreadsafeState(State):
  """ Protected against multiple access """
  
  def __init__(self):
    self.n=0
    self.name='init'
    self.frac=0.0
    self.t0=time.time()
    self.times=[]
    self._lock=thread.allocate_lock()

  def next(self,name,**kwargs):
    self._lock.acquire()
    self.times.append(time.time()-self.t0)
    self.n+=1
    self.name=name
    self.frac=0.0
    for (k,v) in kwargs.iteritems():
      setattr(self,k,v)
    self._lock.release()
      
  def set_frac(self,frac):
    self._lock.acquire()    
    self.frac=frac
    self._lock.release()

  # serialization

  def get_as_dict(self):
    self._lock.acquire()    
    d=self.__dict__.copy()
    self._lock.release()
    for k in d.keys():
      if k[0]=='_': del d[k]
    return d

  def set_from_dict(self,d):
    self._lock.acquire()    
    for (k,v) in d.iteritems():
      setattr(self,k,v)
    self._lock.release()
    


class WritePriorityLock:
  """ Protects an object against concurrent read/writes:
  - reads can be concurrent
  - writes must be atomic
  - writes have priority over reads
  """
  
  def __init__(self):
    self.lock=thread.allocate_lock()    
    self.r_cond=threading.Condition(self.lock)
    self.w_cond=threading.Condition(self.lock) # waiting readers
    self.nr=0   # nb of readers (-1 = currently writing)
    self.nww=0  # nb of waiting writers

  def get(self):
    self.lock.acquire()    
    while self.nr==-1 or self.nww>0:
      self.r_cond.wait()
    self.nr+=1
    self.lock.release()

  def get_w(self):
    self.lock.acquire()
    while self.nr!=0:
      self.nww+=1
      self.w_cond.wait()
      self.nww-=1
    self.nr=-1
    self.lock.release()

  def release(self):
    self.lock.acquire()
    assert self.nr>0
    self.nr-=1
    if self.nr==0 and self.nww>0:
      self.w_cond.notify()
    self.lock.release()

  def release_w(self):
    self.lock.acquire()
    assert self.nr==-1
    self.nr=0
    if self.nww>0:
      self.w_cond.notify()
    else:
      self.r_cond.notifyAll()
    self.lock.release()



def try_mkdir(dirname):
  """ returns False when the directory already existed """
  try:
    os.mkdir(dirname)
  except OSError,e:
    if e.errno==errno.EEXIST:
      return False
    else:
      raise e
  return True


def copy_file(src,dst):
  # print "copy %s -> %s"%(src,dst)
  open(dst,"w").write(open(src,"r").read())



def aff_to_py(sd):
  da=DoubleArray.frompointer(sd)
  return [da[i] for i in range(6)]

def gen_thumbnail(infile,outfile):
  """ uses ImageMagick """
  os.system('convert -strip -quality 50  -geometry "150x150>" "%s" "%s"'%
            (infile,outfile))



class Pool:
  """ pool = array of (object,lock) where each object can be used by one
  thread only. If a new object is needed, cons() is called """

  def __init__(self,cons):
    self.cons=cons
    self.plock=thread.allocate_lock()
    self.pool=[]


  def get_one(self):
    """ get a non-locked object from the pool.
    If there is none, allocate a new one with cons() """
    self.plock.acquire()
    for o,lock in self.pool:
      if lock.acquire(0):
        break
    else:
      o=self.cons()
      lock=thread.allocate_lock()
      lock.acquire()
      self.pool.append((o,lock))
    self.plock.release()
    return o                  

  def release(self,o_in):
    """ release a locked object from the pool """
    self.plock.acquire()
    for o,lock in self.pool:
      if o_in==o:
        lock.release()
        break
    else:
      print "!!! lock to release not found o=%s pool=%s"%(o_in,pool)
    self.plock.release()

  def clear(self):
    "removes all elements in the pool"
    self.plock.acquire()
    for o,lock in self.pool:
      assert lock.acquire(0)
    self.pool=[]
    self.plock.release()
    
  def size(self):
    return len(self.pool)
