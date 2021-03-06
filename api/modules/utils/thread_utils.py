
import sys,time,inspect,thread,threading,os,errno,types,traceback
from yael import yael
#import utils

def handle_stackframe_without_leak():
  frame = inspect.currentframe()
  print "stack:"
  try:      
    while frame!=None:
      print inspect.getframeinfo(frame)
      frame=frame.f_back
  finally:      
    del frame
  sys.stdout.flush()

        
class VerbLock:  
  """ to debug deadlocks..."""

  id=0
  pool=[]
  
  def __init__(self,lock):      
    self.lock=lock
    self.id=VerbLock.id
    VerbLock.id+=1    
    
    
  def acquire(self,*x):      
    sys.stdout.write("acquire %d pool=%s..."%(self.id,VerbLock.pool))
    handle_stackframe_without_leak()
    sys.stdout.write("<")    
    sys.stdout.flush()
    self.lock.acquire(*x)
    sys.stdout.write(">\n")
    sys.stdout.flush()
    VerbLock.pool.append(self.id)
    
  def release(self):
    sys.stdout.write("release %d pool=%s<"%(self.id,VerbLock.pool))
    sys.stdout.flush()
    self.lock.release()
    sys.stdout.write(">\n")
    sys.stdout.flush()
    VerbLock.pool.remove(self.id)
    

the_allocate_lock=thread.allocate_lock

def alt_allocate_lock():
  return VerbLock(the_allocate_lock())

# thread.allocate_lock=alt_allocate_lock

class DummyState:
  def next(self,name,**kwargs):
    pass
  set_frac=None


class ThreadsafeState:
  """ Protected against multiple access """
  
  def __init__(self):
    self.n=0
    self.name='init'
    self.frac=0.0
    self.t0=time.time()
    self.times=[]
    self._lock=thread.allocate_lock()
    self._logf=None

  def next(self,name,**kwargs):
    self._lock.acquire()
    self.times.append(time.time()-self.t0)
    self.n+=1
    self.name=name
    self.frac=0.0
    logs="state %d: %s "%(self.n,name)
    for (k,v) in kwargs.iteritems():
      setattr(self,k,v)
      logs+="%s:%s "%(k,v)
    if self._logf: self._logf.log(logs)
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


class RunOnSet:
  """ run function f(i) on all i in l (l is an iterable) with n threads """
  
  def __init__(self,n,l,f):
    if n==1:
      for i in l:
        f(i)
    else:
      (self.n,self.l,self.f)=(n,iter(l),f)
      self.exception=None
      self.lock=thread.allocate_lock()
      self.lock2=thread.allocate_lock()
      self.lock2.acquire()      
      for i in range(n):
        thread.start_new_thread(self.loop,())
      self.lock2.acquire()
      if self.exception:
        raise self.exception
      
  def loop(self):
    while True:
      self.lock.acquire()
      try:
        i=self.l.next()
      except StopIteration:
        self.n-=1
        if self.n==0:
          self.lock2.release()
        self.lock.release()
        return
      self.lock.release()
      
      try:
        self.f(i)
      except Exception,e:
        traceback.print_exc(50,sys.stderr)        
        self.lock.acquire()
        self.l=[] # stop all
        self.exception=e
        self.lock2.release() # avoid deadlock
        self.lock.release()

class RunOnSetWithPeek(RunOnSet):
  """ same as RunOnSet, with a peek function (l must have a length)"""

  def __init__(self,n,l,f,peek_fun=None):
    self.tot_l=len(l)
    self.l_no=0
    self.peek_fun=peek_fun
    self.orig_f=f
    self.lock3=thread.allocate_lock()
    RunOnSet.__init__(self,n,l,self.fx)
    del self.f # avoid circular ref

  def fx(self,i):
    self.orig_f(i)
    if self.peek_fun:
      self.lock3.acquire()
      self.l_no+=1
      self.peek_fun(float(self.l_no)/self.tot_l)        
      self.lock3.release()

class ParallelMap:

   def __init__(self,n,l,f):
     self.orig_f=f
     self.result=[]
     self.reslock=thread.allocate_lock()
     RunOnSet(n,enumerate(l),self.f)

   def f(self,(i,x)):
     y=self.orig_f(x)
     self.reslock.acquire()
     if i>=len(self.result):
       self.result+=[None]*(i-len(self.result)+1)
     self.result[i]=y
     self.reslock.release()

def parallel_map(n,l,f):
  return ParallelMap(n,l,f).result

class ParallelMapWithPeek:

   def __init__(self,n,l,f,peek_fun):
     self.orig_f2=f
     self.l2=l
     self.result=[None]*len(l)
     RunOnSetWithPeek(n,range(len(l)),self.f2,peek_fun)

   def f2(self,i):
     self.result[i]=self.orig_f2(self.l2[i])


# clear to end of line, nice in combination with \r
clr_eol="\033[K"
      


class ProducerConsumer:
  """ Producer-consumer with a maximum buffer size and
  a merchanism to stop consumers when there are no producers left"""  

  def __init__(self,maxbuf,nprod):
    self.maxbuf=maxbuf
    self.nprod=nprod
    self.file=[]
    self.l=thread.allocate_lock()
    self.cv_full=threading.Condition(self.l)
    self.cv_empty=threading.Condition(self.l)
    self.cv_end=threading.Condition(self.l)
    
  def prod(self,x):
    " Produce something "
    self.l.acquire()
    while len(self.file)==self.maxbuf:
      self.cv_full.wait()
    self.file.append(x)
    if len(self.file)==1:
      self.cv_empty.notify()
    self.l.release()

  def cons(self):
    " get something to consume. Return None if there are no producers left "
    self.l.acquire()
    while self.file==[] and self.nprod>0:
      self.cv_empty.wait()
    if self.file!=[]:
      x=self.file.pop(0)
    else:
      x=None
    if len(self.file)==self.maxbuf-1:
      self.cv_full.notify()
    self.l.release()
    return x

  def prod_end(self):
    " Notify end of one producer "
    self.l.acquire()
    self.nprod-=1
    if self.nprod==0:
      print "prod_end"
      self.cv_empty.notifyAll()
      self.cv_end.notifyAll()
    self.l.release()
    
  def wait_end(self):
    " Wait for end of producers "
    self.l.acquire()
    while self.nprod>0:
      self.cv_end.wait()
    self.l.release()


def parse_nt(nt):
  """ Parse a -nt argument, that can be a nb of threads or "all" = as
  many threads as cpus"""  
  if nt=="all":
    nt=yael.count_cpu()
    print "using %d threads"%nt
  else:
    nt=int(nt)
  return nt



class Cache:
  """ Cache with upper bound on size, LRU policy. Use like a dictionary
  TODO: replace the LRU queue with a linked list
  """

  def __init__(self,maxsz=32768):
    self.cache={}
    self.queue=[]
    self.maxsz=maxsz
    self.nget=0
    self.nfail=0
    self.ndel=0
    
  def __setitem__(self,i,v):
    # print "cache add %s "%i
    if i not in self.cache:
      if len(self.queue)>=self.maxsz:
        # print "cache del ",self.queue[-1]
        self.ndel+=1
        del self.cache[self.queue[-1]]
        del self.queue[-1]        
    else:
      self.queue.remove(i)
    self.cache[i]=v
    self.queue.insert(0,i)
    
  def __getitem__(self,i):
    self.nget+=1
    # print "cache get %s "%i
    if i not in self:
      self.nfail+=1
      raise KeyError("%s not in cache"%i)
    self.queue.remove(i)
    self.queue.insert(0,i)
    return self.cache[i]

  def __contains__(self,i):
    return i in self.cache

  def stats(self):
    print "size: %d nget: %d nfail: %d ndel: %d"%(len(self.cache),self.nget,self.nfail,self.ndel)


class ThreadsafeCache(Cache):
  """ LRU cache with upper bound on size
  TODO: replace the LRU queue with a linked list
  """

  def __init__(self,maxsz=32768):
    Cache.__init__(self,maxsz)
    self.lock=thread.allocate_lock()
    
  def __setitem__(self,*args):
    self.lock.acquire()
    try:
      Cache.__setitem__(self,*args)
    finally:
      self.lock.release()
    
  def __getitem__(self,*args):
    self.lock.acquire()
    try:
      return Cache.__getitem__(self,*args)
    finally:
      self.lock.release()



class ThreadsafeCache2(Cache):
  "give a function to compute/load the cached value"

  def __init__(self,f,maxsz=32768):
    Cache.__init__(self,maxsz=maxsz)
    self.lock=thread.allocate_lock()
    self.f=f
    
  def __getitem__(self,k):
    self.lock.acquire()
    try:
      try:
        return Cache.__getitem__(self,k)
      except KeyError:
        v=self.f(k)
        Cache.__setitem__(self,k,v)
        return v
    finally:
      self.lock.release()

  def __setitem__(self,*args):
    assert False
    
