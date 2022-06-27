import os,pdb,cPickle,time,errno,sys,thread,traceback,socket,threading,gc

# default
PORT=12032


#########################################################################
# simple I/O functions 



def inline_send_handle(f,conn):
  st=os.fstat(f.fileno())
  size=st.st_size
  cPickle.dump(size,conn)
  conn.write(f.read(size))
     

def inline_send(filename,conn):
  inline_send_handle(open(filename,'r'),conn)

def inline_recv_handle(f,conn):
  size=cPickle.load(conn)
  rd=0
  while rd<size:
    sz=size-rd
    if sz>4096: sz=4096
    buf=conn.read(sz)
    f.write(buf)
    rd+=len(buf)    

def inline_recv(filename,conn):
  inline_recv_handle(open(filename,"w"),conn)

class FileSock:
  " wraps a socket so that it is usable by pickle/cPickle "
  
  def __init__(self,sock):
    assert type(sock)!=type('')
    self.sock=sock
    
  def write(self,buf):
    assert type(self.sock)!=type('')
    self.sock.sendall(buf)
    
  def read(self,bs=1024):
    return self.sock.recv(bs)

  def readline(self):
    """may be optimized..."""
    s=''
    while True:
      c=self.read(1)
      s+=c
      if c=='\n' or len(c)==0:
        return s

class SocketTransferred(Exception):
  " the controlling socket of this RID has been given to another one "
  def __init__(self,to_rid):
    self.to_rid=to_rid

class CannotChangeRID(Exception):
  pass

class ClientExit(Exception):
  pass

class EndOfAsyncCall(Exception):
  pass

class ServerException(Exception):
  pass


class Server:
  """
  server protocol. Methods from classes that subclass Server can be called
  transparently from a client
  """

  instance_pool={}
  important_rids=[] 
  instance_pool_lock=thread.allocate_lock()
  
  def __init__(self,s,rid=0):
    self.rid=rid
    self.datadir="data/commondata/"

    # register instance
    Server.instance_pool_lock.acquire()
    Server.instance_pool[self.rid]=self
    Server.instance_pool_lock.release()

    # connection

    self.conn=s
    self.fs=FileSock(s)

    # if data should be transmitted after the return message,
    # register a callback that does it here
    self.call_after_return=None

    # should the next function call be done in detached mode?
    # 0: no
    # 1,5: yes (at next function call)
    # 2,6: yes (at this function call)
    # if next_detach & 4: at end of function, we'll wait to
    #    return the result to an attach() call
    self.next_detach=0

    self.logf=sys.stderr

  def log(self,s):
    self.logf.write("rid %d %s\n"%(self.rid,s))

  def get_by_rid(self,rid):
    Server.instance_pool_lock.acquire()
    # print "instance pool: ",Server.instance_pool.keys()
    try:
      other=Server.instance_pool[rid]
    except KeyError:
      other=None
    Server.instance_pool_lock.release()
    return other

  def resume_rid(self,rid):
    """ resumes a finished RID
    does NOT set the rid """
    Server.instance_pool_lock.acquire()
    if Server.instance_pool.has_key(rid):
      Server.instance_pool_lock.release()
      raise CannotChangeRID('RID busy')
    # silently move to new id
    del Server.instance_pool[self.rid]
    Server.instance_pool[rid]=self    
    Server.instance_pool_lock.release()
    
  def detach_at_next_call(self,join=True):
    """ next call will close the socket """
    self.next_detach=join and 5 or 1
    
  def detach(self):
    self.det_lock=thread.allocate_lock()
    self.det_lock.acquire()
    self.waiting_cond=None
    self.fs=None
    self.conn.close(); self.conn=None    
    self.det_lock.release()
    
  def attach(self,rid):
    other=self.get_by_rid(rid)
    if not other:
      raise CannotChangeRID("cannot find other")

    self.log("found other RID")

    other.det_lock.acquire()
    if not other.waiting_cond:
      other.det_lock.release()
      raise CannotChangeRID("other not in wait_attach")

    # transfer connection 
    other.conn=self.conn;    self.conn=None
    other.fs=self.fs;        self.fs=None
    
    other.waiting_cond.notify()
    other.det_lock.release()

    # exit gracefully
    raise SocketTransferred(rid)

  def wait_attach(self):
    self.det_lock.acquire()
    self.waiting_cond=threading.Condition(self.det_lock)
    self.log("wait_attach wait")
    self.waiting_cond.wait()
    # shoud set timeout for dead clients
    self.log("wait_attach end of wait")
    self.det_lock.release()
    del self.waiting_cond
    del self.det_lock

  def one_function(self):
    """
    Executes a single function with associated I/O.
    Protocol:
    - the arguments and results are serialized with the pickle protocol
    - client sends : (fname,args)
        fname = method name to call
        args = tuple of arguments
    - server sends result: (rid,st,ret)
        rid = request id
        st = None, or exception if there was during execution
        ret = return value or None if st!=None
    - if data must be sent raw over the socket after the function returns,
      call_after_return must be set to a function that does this
    - if the client wants to close a connection during a function, it
      should call detach_at_next_call(join)
      If join then
        the client should reatach from another rid with attach(rid), which
        will return the function's result
      else
        the server exits at the end of the function call and drops the result
    """
    
    try:
      (fname,args)=cPickle.load(self.fs)
    except EOFError:
      raise ClientExit("read args")
    self.log("executing method %s with args %s"%(fname,args))
    st=None
    ret=None
    try:
      f=getattr(self,fname)
    except AttributeError:
      st=AttributeError("unknown method "+fname)
      self.log("unknown method ")
    else:

      if self.next_detach in [2,6]:
        cPickle.dump((self.rid,None,None),self.fs)
        self.detach()

      try:
        ret=f(*args)
      except SocketTransferred,e:
        raise e
      except Exception,e:
        st=ServerException(
          "".join(traceback.format_tb(sys.exc_info()[2]))+
          str(e))
        self.log("exception in method")
        traceback.print_exc(50,self.logf)
        self.logf.flush()

      if self.next_detach==6:
        self.wait_attach()
        self.next_detach=0
      if self.next_detach==2:
        raise EndOfAsyncCall()
      elif self.next_detach in [1,5]:
        self.next_detach+=1

    # print "return",ret
    try:
      cPickle.dump((self.rid,st,ret),self.fs)
    except EOFError:
      raise ClientExit("function return")
      
    # pdb.set_trace()
    if self.call_after_return!=None:
      self.log("doing call_after_return")
      try:
        self.call_after_return()
      except Exception,e:
        # this is a problem: we can't propagate the exception to the client
        # so he will probably crash...
        self.log("!!! exception in call_after_return")
        traceback.print_exc(50,self.logf)
        traceback.print_exc(50,sys.stderr)
        self.logf.flush()          
      self.call_after_return=None
    
  def exec_loop(self):
    """ main execution loop. Loops and handles exit states
    """

    self.log("in exec_loop")
    try:
      while True:
        self.one_function()
    except ClientExit,e:
      self.log("ClientExit %s"%e)
    except SocketTransferred,e:
      self.log("socket transferred to RID thread %d"%e.to_rid)
    except EndOfAsyncCall,e:
      self.log("EndOfAsyncCall %s"%e)      
    except socket.error,e:
      self.log("socket error %s"%e)
      traceback.print_exc(50,self.logf)
    except EOFError:
      self.log("EOF during communication")
      traceback.print_exc(50,self.logf)  
    except:
      # unexpected
      traceback.print_exc(50,sys.stderr)  
      sys.exit(1)

    Server.instance_pool_lock.acquire()
    del Server.instance_pool[self.rid]
    Server.instance_pool_lock.release()
    try:
      self.exec_loop_cleanup()
    except:
      traceback.print_exc(50,sys.stderr)  
          
    print "exit rid %d "%self.rid

  def exec_loop_cleanup(self):
    pass

  def get_datafile(self,fname):
    " sends a file from the datadir to the client "
    fh=open(fname,"r")
    self.call_after_return=lambda *x: inline_send_handle(fh,self.fs)    
    return None

  def put_datafile(self,fname):
    " puts a file coming from the client in the datadir "
    fname=self.datadir+fname
    fh=open(fname,"w")
    self.call_after_return=lambda *x: inline_recv_handle(fh,self.fs)
    return fname

  ###################################################################
  # spying stuff

  def get_ps_stats(self):
    ret=''
    f=os.popen("echo ============ uptime:; uptime;"+
               "echo ============ self:; "+
               "ps -p %d -o pid,vsize,rss,%%cpu,nlwp,psr; "%os.getpid()+
               "echo ============ run queue:;"+
               "ps ar -o user,pid,%cpu,%mem,ni,nlwp,psr,vsz,rss,cputime,command")
    for l in f:
      ret+=l
    return ret
 
  def get_server_pool(self):
    Server.instance_pool_lock.acquire()
    rids=Server.instance_pool.keys()
    Server.instance_pool_lock.release()
    return rids

  def iam_important(self,key):
    Server.instance_pool_lock.acquire()
    Server.important_rids.append((self.rid,key))
    Server.instance_pool_lock.release()

  def get_importants(self):
    Server.instance_pool_lock.acquire()
    i=Server.important_rids[:]
    Server.instance_pool_lock.release()
    return i
  
class Client:
  """
  Methods of the server object can be called transparently. Exceptions are
  re-raised.
  """
  def __init__(self,HOST,port=PORT):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    print "connecting",(HOST, port)
    sock.connect((HOST, port))
    
    self.sock=sock
    self.fs=FileSock(sock)
    self.rid=-1

    self.async_state=0
    
  def generic_fun(self,fname,args):
    # int "gen fun",fname
    if self.async_state==2:
      raise RuntimeError("async call in progress")
    
    cPickle.dump((fname,args),self.fs)
    if self.async_state==1:
      self.async_state=2
    else:
      return self.get_result()

  def get_result(self):
    (rid,st,ret)=cPickle.load(self.fs)
    self.async_state=0
    self.rid=rid
    if st!=None:
      raise st
    else:
      return ret

  def async_at_next_call(self):
    self.async_state=1
  
  def __getattr__(self,name):
    return lambda *x: self.generic_fun(name,x)

  def get_datafile_handle(self,dist_name,local_fh):
    self.generic_fun('get_datafile',(dist_name,))
    inline_recv_handle(local_fh,self.fs)

  def put_datafile_handle(self,local_fh,dist_name):
    dist_name_2=self.generic_fun('put_datafile',(dist_name,))
    inline_send_handle(local_fh,self.fs)
    return dist_name_2

  def get_datafile(self,dist_name,local_name):
    return self.get_datafile_handle(dist_name,open(local_name,"w"))

  def put_datafile(self,local_name,dist_name):
    return self.put_datafile_handle(open(local_name,"r"),dist_name)


def run_server(new_handler,port=PORT):

  HOST = ''                 # Symbolic name meaning the local host
  s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
  s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

  print "bind ",(HOST, port)
  s.bind((HOST, port))
  s.listen(5)

  print "accepting connections"

  while True:
    conn, addr = s.accept()
    print 'Connected by', addr,

    ibs=new_handler(conn)

    print "handled by rid ",ibs.rid,

    tid=thread.start_new_thread(ibs.exec_loop,())

    print "tid",tid

