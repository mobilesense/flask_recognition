import sys,socket,thread,cPickle,errno,time,copy,inspect,types

from protocol import *
from imbase import *
from config import select_config
from utils import *




class ThreadsafeVWAssignement(VWAssignement):

  def __init__(self,config):
    VWAssignement.__init__(self,config)
    self.scores_pool=Pool(lambda : ihash_clone(self.scores))
    
  def get_ann(self):
    return (self.clusteri,self.ann,self.scores_pool.get_one())

  def release_scores(self,scores_in):
    self.scores_pool.release(scores_in)



class ThreadsafeBigBase(BigBase):
  
  def __init__(self,config):
    self.ivf_lock=WritePriorityLock()
    self.ivfq_pool=Pool(
      lambda: invfile_query_new(self.config.q_nbcells))
    BigBase.__init__(self,config)
    
  def get_ivf(self,w=False):
    if w:
      self.ivf_lock.get_w()
    else:
      self.ivf_lock.get()
    (ivf,ivfm)=BigBase.get_ivf(self)    
    return (ivf,ivfm)

  def release_ivf(self,w=False):
    if w:
      self.ivf_lock.release_w()
    else:
      self.ivf_lock.release()
      
  def prepare_ivfm(self):
    self.ivfq_pool.clear()
    BigBase.prepare_ivfm(self)

  def get_ivf_query(self):
    return self.ivfq_pool.get_one()
    
  def release_ivf_query(self,ivfq):
    self.ivfq_pool.release(ivfq)


class ThreadsafeCachedData(CachedData):
  """ data that is common between request handler (and costly to load)
  access must be protected against multiple thread access.
  """

  def __init__(self,config):
    CachedData.__init__(self,config)
    self.load_lock=thread.allocate_lock()
    self.file_lock=thread.allocate_lock()
    self.file_pipe=None
        
  def get_FeaturesComputer(self):
    self.load_lock.acquire()
    fc=CachedData.get_FeaturesComputer(self)
    self.load_lock.release()
    return fc

  def get_VWAssignement(self):
    self.load_lock.acquire()    
    if self.the_VWAssignement==None:
      self.the_VWAssignement=ThreadsafeVWAssignement(self.config)
    self.load_lock.release()
    return self.the_VWAssignement

  def get_BigBase(self):
    self.load_lock.acquire()    
    if self.the_BigBase==None:
      self.the_BigBase=ThreadsafeBigBase(self.config)
    self.load_lock.release()
    return self.the_BigBase

  def get_GeometricVerification(self):
    self.load_lock.acquire()
    gv=CachedData.get_GeometricVerification(self)
    self.load_lock.release()
    return gv


  def get_mimetype(self,filename):
    self.file_lock.acquire()
    if self.file_pipe==None:
      self.file_pipe=os.popen2("file -b -n -i -L -f -")
    (file_in,file_out)=self.file_pipe
    file_in.write("%s\n"%filename)
    file_in.flush()
    line=file_out.readline()
    self.file_lock.release()   
    return line.strip()

  def alloc_nos(self,nf):
    self.load_lock.acquire()
    a=CachedData.alloc_nos(self,nf)
    self.load_lock.release()
    return a


class ImBaseServer(ImBase,Server):
  """ The server can be used to call the individual functions or, with
  process_search, in a detached remote function call """

  def __init__(self,cachedData,s):
    ImBase.__init__(self,cachedData)
    # those 2 are overridden by Server.__init__
    # so save + restore
    logf=self.logf
    datadir=self.datadir
    Server.__init__(self,s,self.rid)
    self.logf=logf
    self.datadir=datadir
    self.state=ThreadsafeState()
    self.move_datadir=None
    
  def process_search(self,image_file):
    return self.search(image_file,self.state)

  def process_add(self,imnames):
    return self.add_images(imnames,True,self.state)

  def process_add_from_dir(self,dirname):
    # self.state.next("scan dir")
    imnames=[]
    for f in os.listdir(dirname):
      imname=dirname+"/"+f
      if self.cachedData.get_mimetype(imname) in [
        "image/jpeg", "image/png", "image/gif",
        "image/x-portable-pixmap","image/tiff"]:
        imnames.append(imname)
    self.log("adding %s"%imnames)
    return self.process_add(imnames)
          
  def get_state(self):
    return self.state.get_as_dict()
  
  def get_other_state(self,rid):
    other=self.get_by_rid(rid)
    if other:
      return other.get_state()
    else:
      return None
    
  def get_mimetype(self,filename):
    return self.cachedData.get_mimetype(filename)

  def atexit_move_datadir(self,other_rid):
    if other_rid==self.rid:  return
    self.move_datadir=RequestHandler.datadir_pattern%other_rid

  def get_config_dict(self):
    return self.cachedData.config.__dict__
        
  def get_cache_state(self):
    s=""
    if self.cachedData.the_VWAssignement!=None:
      vwa=self.cachedData.get_VWAssignement()
      s+="clusters(%d scores pool)"%vwa.scores_pool.size()
    if self.cachedData.the_BigBase!=None:
      bb=self.cachedData.get_BigBase()
      s+="ivf(%d pool) "%bb.ivfq_pool.size()
    if self.cachedData.the_GeometricVerification!=None:
      s+="PCA "
    if self.cachedData.file_pipe!=None:
      s+="file_pipe "      
    return s


  ############################################################
  # persistent variable management (can be used after resume_rid)
  # stored in the state object

  # parameters that should be preserved => copied to state variable
  # when dumped
  params=['n_short','geom_thr','geom_fthr']
  
  def copy_params(self,src,dst):
    for param in ImBaseServer.params:
      try:
        setattr(dst,param,getattr(src,param))
      except AttributeError:
        pass

  def resume_rid(self,rid):
    Server.resume_rid(self,rid)
    ImBase.resume_rid(self,rid)
    try:
      f=open(self.datadir+"state.raw","r")
    except IOError:
      self.log("no variables to load")
    else:
      self.log("loading vars")
      self.state=ThreadsafeState()
      self.state.set_from_dict(cPickle.load(f))
      self.copy_params(self.state,self)
            
  def exec_loop_cleanup(self):
    if True:
      self.log("saving vars")
      self.copy_params(self,self.state)
      cPickle.dump(self.state.get_as_dict(),
                   open(self.datadir+"state.raw","w"))
    if self.move_datadir!=None:
      dd2=self.move_datadir+'joined'
      try_mkdir(dd2)
      jn=0
      while True:
        dd3=dd2+"/%d"%jn
        if try_mkdir(dd3): break
        jn+=1
      # print "ren %s => %s"%(self.datadir,dd3+"/%05d"%self.rid)
      os.rename(self.datadir,dd3+"/%05d"%self.rid)
    
  def set_var(self,name,val):
    setattr(self.state,name,val)
  
  def get_var(self,name):
    return getattr(self.state,name)
  

def usage():
  sys.stderr.write("%s [-port port] dbname\n"%sys.argv[0])
  sys.exit(1)

    
# if we're in the main program
if __name__=="__main__":

  dbname=None
  port=12032
  args=sys.argv[1:]
  while args:
    if args[0] in ('-h','--help'):
      usage()
    elif args[0]=='-port':
      del args[0]
      port=int(args[0])
    else:
      if not dbname: dbname=args[0]
      else: usage()
    del args[0]

  if not dbname: usage()
 
  cd=ThreadsafeCachedData(select_config(dbname))
  
  run_server(lambda conn: ImBaseServer(cd,conn),port=port)
  

