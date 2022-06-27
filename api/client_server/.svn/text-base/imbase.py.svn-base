

import sys,socket,errno,time,os,pdb

from invfile import *

import geom_filter

import obsidian

from utils import *



class RequestHandler:
  """
  A request has a request id (rid), and a corresponding data directory for
  a log file and intermediate results.
  
  """

  id_counter=0

  datadir_pattern="data/%05d/"

  def __init__(self):    
    """ s is the r/w communication socket with the client
    may be None if exec_loop is not used
    """
    # request id & datadir
    
    base_rid=RequestHandler.id_counter

    # find a directory that is not used
    while True:
      datadir=RequestHandler.datadir_pattern%base_rid      
      if try_mkdir(datadir): break
      base_rid+=1

    RequestHandler.id_counter=base_rid+1

    self.rid=base_rid
    self.datadir=datadir

    # log file
    
    self.verbose=0
    self.logf=open(self.datadir+"log.txt","w")
    self.log("rid=%d"%self.rid)


  def log(self,s):
    " write s to the log file "
    self.logf.write("%s: %s\n"%(time.strftime("%Y-%m-%d %H:%M:%S"),s))
    if self.verbose:
      print "log: %s"%s
    self.logf.flush()

  def resume_rid(self,rid):
    """ changes the RID to the one of a *finished* RequestHandler """
    self.log("quitting current rid... new one=%d"%rid)
    self.logf.close(); self.logf=None
    old_rid=self.rid
    self.rid=rid
    self.datadir=RequestHandler.datadir_pattern%rid
    self.logf=open(self.datadir+"log.txt","a")
    self.log("====================== Resumed from rid %d"%old_rid)

  def dump_log_file(self,rid):
    datadir=RequestHandler.datadir_pattern%rid
    return open(datadir+"log.txt","r").read()




#############################################################################
# the four processing stages for a search


class ObsidianStage:
  " stage of processing"
  pass

class ObsidianOptimStage(ObsidianStage):
  " stage that has several levels of optimization"
  def __init__(self):
    self.optim='cpu'
    # Alternatives('lava|cpu')

class ObsidianFeaturesComputer:
  """ Computes descriptors with Obsidian code """

  def __init__(self,config):
    " default values "
    self.ins=ObsidianStage()
    self.ins.max_pixel=-1

    # pyramid generation
    self.pyr=ObsidianOptimStage()
    self.pyr.scale_fact=1.2
    self.pyr.min_size=20
    self.pyr.smooth=True

    # detector
    self.det=ObsidianOptimStage()
    self.det.type='harris_laplace'
    #   Alternatives('harris_laplace|harris|hessian|log')
    self.det.threshold=1000.0
    self.det.max_desc=-1

    # normalizations
    self.aff=ObsidianOptimStage()
    self.aff.skip=False

    self.ori=ObsidianOptimStage()
    self.ori.skip=False

    # descriptor
    self.desc=ObsidianOptimStage()
    self.desc.type='cs_lbp'
    # Alternatives('cs_lbp|sift')
    self.desc.window_size=64
    self.desc.scale_one=16.0
    self.desc.index_size=4
    self.desc.T=0.015
    self.desc.ori_size=8
    self.desc.norm=0
    self.desc.dim=128

    # output file
    self.out=ObsidianStage()
    # self.out.file='a.desc'
    self.out.format='geo'
    # Alternatives('geo|kristian')

    # customize through provided configuration
    self.config=config
    self.parse_args(self.config.det_obsidian_args)

    self.image_pyramid_build=self.obsidian_func(
      'image_pyramid_build_%s'%self.pyr.optim)
    
    self.detector=self.obsidian_func(
      'detector_%s_%s'%(self.det.type,self.det.optim))

    if not self.aff.skip:
      self.normalizer_second_moment=self.obsidian_func(
        'normalizer_second_moment_%s'%self.aff.optim)
    else:
      self.normalizer_second_moment=lambda x: None

    if not self.ori.skip:
      self.normalizer_grad_ori=self.obsidian_func(
        'normalizer_grad_ori_%s'%self.aff.optim)
    else:
      self.normalizer_grad_ori=lambda x: None

    self.descriptor=self.obsidian_func(
      'descriptor_%s_%s'%(self.desc.type,self.det.optim))

  def obsidian_func(self,funcname):
    return obsidian.__dict__[funcname]
      
  def parse_args(self,args):
    """ parse arguments, passes as a string in the form
    "pyr-min_size=45 ori-skip=1"
    NB that booleans must be passed as 0 or 1
    """
    argv=args.split()
    while argv:
      opt=argv.pop(0)
      [k,v]=opt.split('=')
      [stname,parname]=k.split('-')
      stage=getattr(self,stname)
      param=getattr(stage,parname)
      setattr(stage,parname,(type(param))(v))

  def compute_descriptors(self,img,peek_fun=None):
    """ returns initial nb of descriptors & descriptor list """

    pstep=1/4.0
    if peek_fun: peek_fun(0)
    
    pyr=self.image_pyramid_build(
      img,self.pyr.scale_fact,self.pyr.min_size,self.pyr.smooth)

    if peek_fun: peek_fun(pstep)
    
    descs=self.detector(pyr,self.det.threshold)
    n_interest=descs.size

    obsidian.local_desc_selector(descs,self.det.max_desc)

    if peek_fun: peek_fun(2*pstep)

    if not self.aff.skip:
      self.normalizer_second_moment(descs,pyr)

    if not self.ori.skip:
      self.normalizer_grad_ori(descs,pyr)

    if peek_fun: peek_fun(3*pstep)

    if str(self.desc.type)=='cs_lbp':
      self.descriptor(descs,pyr,self.desc.window_size,self.desc.scale_one,
                 self.desc.index_size,self.desc.T,self.desc.norm,
                 self.desc.dim)
    elif str(self.desc.type)=='sift':
      self.descriptor(descs,pyr,self.desc.window_size,self.desc.scale_one,
                 self.desc.index_size)

    if peek_fun: peek_fun(4*pstep)

    return (n_interest,descs.size,
            obsidian.local_desc_list_as_string(descs))
    
  def load_image(self,fname):
    # return obsidian.image_load_max(fname,self.ins.max_pixel)
    return obsidian.image_load(fname)

  def store_descriptors(self,descs,filename):
    open(filename,'w').write(descs)
    

class MikoFeaturesComputer:
  """ Computes features with Krystian Mikolajczyk's programs.
  Images and descriptors are passed with files """

  def __init__(self,config):
    self.config=config
  
  def exec_miko(self,cmd,peek_fun=None):
    """ execute one of Krystian's description programs and read n_desc
    and n_interest from the stdout"""
    
    # self.log("execing %s"%cmd)
    
    f=os.popen(cmd,'r')

    n_interest=None
    n_desc=None
    line_no=0
    for l in f:
      # self.log("stdout %s"%l[:-1])
      l=l.split()
      if len(l)>=3 and l[0]=="interest" and l[1]=="points":
        n_interest=int(l[2])
      elif len(l)>=3 and l[0]=="saving" and l[2]=="feaures":
        n_desc=int(l[1])
      # peek does not work well because of buffering
      if peek_fun!=None:
        peek_fun(line_no/9.0)
      line_no+=1

    ret=f.close()
    if ret:
      raise RuntimeError("error in exec_miko: %s "%ret)

    return (n_interest,n_desc)
  
  def compute_descriptors_files(self,imagefile,descfile,peek_fun=None):
    """ computes the descriptors for the image stored in filename
    returns the name of the bsiftgeo descriptor file """
    
    bsiftgeo=descfile

    if self.config.det_prog=='compute_descriptors':
      cmd="./compute_descriptors -hesaff -thres %g -i2 %s -sift  -max %d -o4 %s"%(
        self.config.det_thresh,imagefile,self.config.det_max_sift,bsiftgeo)
      (n_interest,n_desc)=self.exec_miko(cmd,peek_fun)

    else:
      assert False
     
    return (n_interest,n_desc)
  
  def load_image(self,fname):
    " dummy load" 
    return fname

  def compute_descriptors(self,img,peek_fun=None):
    """ returns initial nb of descriptors & descriptor list """
    descfile=os.tempnam()
    (n_interest,n_desc)=self.compute_descriptors_files(img,descfile,peek_fun)
    descs=open(descfile,"r").read()
    os.unlink(descfile)
    return (n_interest,n_desc,descs)
    
  def store_descriptors(self,descs,filename):
    open(filename,'w').write(descs)


# unfortunately, interest point lists appear in 3 memory formats:
# - obsidian local_desc_list_t
# - imat when they are to be assigned to clusters
# - PhotoMole ImageDescriptor
# the common format to store them is bsiftgeo
# Therefore, I decided that the common interchange format would be
# a python string with the contents of a bsiftgeo. Caevat: if the nb
# of descirptors is not a multiple of 4, floats may get misaligned!

class VWAssignement:
  """ assigns local descriptors to clusters (visual words). NB that all
  encapsulated objects must be freed explicitly"""

  def __init__(self,config):
    self.config=config
    self.load_clusteri()
    
  def default_ann(self,clusteri,lattice):
    """ translated from ann_construct in siftgeo2vw.c"""
    w=15
    w=(w<<SIFTGEO_PRECISION)/256
    vs=6
    nh=200
    nbcells=10000
    l=200
    d=imat_width(clusteri)
    seed=666
    if not lattice:
      ann=ann_new_projection (d, nh, nbcells, w, l, vs, 1, seed)
    else:
      ann=ann_new_latticeE8 (d, nbcells, w, l, vs, 1, seed)
    for i in range(imat_height(clusteri)):
      ann_add(ann,imat_getitem(clusteri,i),i)
    return ann
    
  def load_clusteri(self):
    print "reading clusters from %s"%self.config.clusterfile
    if self.config.cluster_format=='float':        
      self.clusteri=load_clusters(self.config.clusterfile)
    elif self.config.cluster_format=='imat':
      f=open(self.config.clusterfile,"r")
      self.clusteri=imat_new_fread(f)
    elif self.config.cluster_format=='ivecs':
      self.clusteri=load_clusters_ivecs(self.config.clusterfile)
    self.ann=self.default_ann(self.clusteri,False)        
    self.scores=new_basic_ihash(10000);

  def get_ann(self):
    return (self.clusteri,self.ann,self.scores)

  def release_scores(self,scores):
    """ useless in single thread context"""    
    pass

  def load_descs(self,filename):
    return open(filename,"r").read()

  def compute_vw(self,descs,peek_f=None):
    """ computes the visual words corresponding to the bsiftgeo file
    returns the visual word file""" 
    
    (clusteri,ann,scores)=self.get_ann()
    nbvw=imat_height(clusteri)
    if descs: # coords may be empty (None)
      vw=ivec_new(nbvw)

      coords=read_bsiftgeo_string(descs)

      if self.config.nn_method=='simple':
        assign_vw_peek(clusteri,coords,vw,peek_f)        
      elif self.config.nn_method=='ann':
        assign_vw_ann_peek(ann,scores,clusteri,coords,vw,100,peek_f)

      imat_delete(coords)    
    else:      
      vw=ivec_new_zeros(nbvw)

    self.release_scores(scores)

    vws=ivec_to_sparse(vw)
    ivec_delete(vw)
        
    return vws
   
  def store_vw(self,(svi,sv),outfname):
    if self.config.sparse_vw:
      f=open(outfname,"w")
      ivec_fwrite(f,svi)
      ivec_fwrite(f,sv)
    else:
      vw=ivec_from_sparse(svi,sv)
      ivec_fwrite(open(outfname,"w"),vw)
      ivec_delete(vw)

  def dealloc_vw(self,(svi,sv)):
    ivec_delete(svi)
    ivec_delete(sv)
  

class BigBase:
  """ does searches of sparse vectors. Encapsulates an
  inverted file"""

  def __init__(self,config):
    self.config=config
    print "reading inverted file %s"%self.config.invfilename
    self.ivf=invfile_read(open(self.config.invfilename,"r"))
    self.ivfm=invfile_meta_new(self.config.q_nbcells)
    self.ivfq=invfile_query_new(self.config.q_nbcells)
    self.prepare_ivfm()
    self.load_cdm(self.config.cdminvfilename, self.config.cdmfactorsname)
  
  def prepare_ivfm(self):
    """ no ivfq should be active """
    print "preparing inverted file..."
    invfile_compute_vec_norm_L1_tfidf(self.ivf,self.ivfm)
    
  def load_cdm(self, ivffile, hashfile):
    if ivffile:
      fd=open(ivffile,"r")
      invfile_cdm_support_read(fd, self.ivfm)
      print "load file %s"%ivffile
    if hashfile:
      fd=open(hashfile,"r")
      invfile_cdm_factors_read(fd, self.ivf, self.ivfm)
      print "load file %s"%hashfile

  def get_ivf(self,for_write=False):
    return (self.ivf,self.ivfm)
    
  def get_ivf_query(self):
    return self.ivfq
    
  def release_ivf_query(self,ivfq_in):
    """ useless in single thread context"""    
    pass

  def release_ivf(self,for_write=False):
    """ useless in single thread context"""    
    pass

  def load_vw(self,vwf):
    if not self.config.sparse_vw:
      v=ivec_new_fread(open(vwf,"r"))
      (svi,sv)=ivec_to_sparse(v); ivec_delete(v)
    else:
      f=open(vwf,"r")
      svi=ivec_new_fread(f)
      sv=ivec_new_fread(f)
    return (svi,sv)

  def dealloc_vw(self,(svi,sv)):
    ivec_delete(svi)
    ivec_delete(sv)
  
  def query(self,vw,n_short=20,peek_f=None):
    """ queries a visual word file and returns the corresponding
    image id's + distances """

    (svi,sv)=vw
    
    (ivf,ivfm)=self.get_ivf()
    ivfq=self.get_ivf_query()
    
    (res_ivec,dis_vec)=invfile_query_L1_norm_peek(
        ivf,ivfm,ivfq,svi,sv,n_short,peek_f)

    self.release_ivf_query(ivfq)    
    self.release_ivf()
    
    res=ivec_to_python(res_ivec); ivec_delete(res_ivec)
    dis=vec_to_python(dis_vec); vec_delete(dis_vec)

    res1=zip(res,dis)
  
    return res1

  def cdm_load(self, ivffile, hashfile):    
    (ivf,ivfm)=self.get_ivf(True)
    self.load_cdm(ivffile, hashfile)
    self.release_ivf(True) 

  def cdm_reset(self,alpha):
    (ivf,ivfm)=self.get_ivf(True)
    invfile_cdm_factors_reset(ivfm, alpha)
    self.release_ivf(True) 

  def cdm_iterate(self,alpha):
    self.cdm_reset(alpha)
    (ivf,ivfm)=self.get_ivf(True)
    ivfq=self.get_ivf_query()

    invfile_cdm_compute_one_iteration(ivf, ivfm, ivfq, self.cdm_ngb)

    self.release_ivf_query(ivfq)    
    self.release_ivf(True) 


class GeometricVerification:
  """ does point-basesd image matching for a short list of images """

  def __init__(self,config):
    self.config=config
    print "reading PCA transform file %s"%self.config.pcafile
    self.pcat=geom_filter.load_pcatransform(self.config.pcafile)
    
  def load_imdesc(self,filename):
    return open(filename,'r').read()

  def filter(self,short_list,query,peek_f):
    """ filters a result from query() with a geometric filter
    returns an ordered list of tuples (imno,npt,nmatch,score), where:
      imno   = image number
      npt    = nb of descriptors in the image
      nmatch = # of matching points
      score  = computed from matching points and distances
    """

    pcat=self.pcat

    imdesc_l=[geom_filter.loadImdesc(self.config.desc_filename(i),geom_filter.DescType_bsiftgeo)
              for i in short_list]
    imdescs=geom_filter.ImageDescriptorVector() 
    for imdesc in imdesc_l:      
      imdescs.push_back(imdesc)
    # imdescs_l is needed so that the descriptors are not erased before
    # function output

    results=geom_filter.ImageMatchVector()

    sp=geom_filter.SearchParams()
    sp.st=self.config.geom_thr
    sp.stf=self.config.geom_fthr
    sp.pcat=pcat

    qdesc=geom_filter.imdescFromString(query,geom_filter.DescType_bsiftgeo)

    geom_filter.search_peek(sp,imdescs,qdesc,results,peek_f)

    res1=[results[i] for i in range(results.size())]
    
    res2=[(r.final_votes,r.imno,r.npt,r.stage3_nmatch,
           aff_to_py(r.affine)) for r in res1 if r.final_votes>0]

    res2.sort()
    res2.reverse()

    return [(short_list[imno],npt,nmatch,score,aff) for
            (score,imno,npt,nmatch,aff) in res2]


    
############################################################################
# Common data between instances

class CachedData:
  """ Manages the processing stage objects. One instance of each is created,
  lazily.
  """

  def __init__(self,config):
    self.config=config

    self.the_FeaturesComputer=None
    self.the_VWAssignement=None
    self.the_BigBase=None
    self.the_GeometricVerification=None

  def get_FeaturesComputer(self):
    if self.the_FeaturesComputer==None:
      if self.config.det_prog=='obsidian':        
        self.the_FeaturesComputer=ObsidianFeaturesComputer(self.config)
      elif self.config.det_prog=='compute_descriptors':        
        self.the_FeaturesComputer=MikoFeaturesComputer(self.config)
    return self.the_FeaturesComputer

  def get_VWAssignement(self):
    if self.the_VWAssignement==None:
      self.the_VWAssignement=VWAssignement(self.config)
    return self.the_VWAssignement

  def get_BigBase(self):
    if self.the_BigBase==None:
      self.the_BigBase=BigBase(self.config)
    return self.the_BigBase

  def get_GeometricVerification(self):
    if self.the_GeometricVerification==None:
      self.the_GeometricVerification=GeometricVerification(self.config)
    return self.the_GeometricVerification

  def alloc_nos(self,nf):
    """ allocates nf new image numbers """
    n0=self.config.nfixed+self.config.nalloc
    self.config.nalloc+=nf
    n1=self.config.nfixed+self.config.nalloc    
    return range(n0,n1)


  

class ImBase(RequestHandler):
  """
  Implements all the image database functions. The inputs & outputs are
  all stored on disk, so they are given as filenames. This may be suboptimal
  but allows to trace all intermediate results
  """
 
  def __init__(self,cachedData):
    RequestHandler.__init__(self)
    self.cachedData=cachedData
    self.config=cachedData.config
    # size of short list
    
    self.n_short=20

    # CDM parameters

    self.cdm_ngb = 10
    
  def set_n_short(self,n_short):
    self.n_short=n_short

  def get_nimg(self):
    return self.config.nimg

  def get_nfixed(self):
    return self.config.nfixed

  def db_im_name(self,imno):
    return self.config.img_filename(imno)

  def db_vw_name(self,imno):
    return self.config.vw_filename(imno)

  def db_thumb_name(self,imno):
    return self.config.thumb_filename(imno)

  def db_desc_name(self,imno):
    return self.config.desc_filename(imno)

  #################################################
  # Decomposed query with I/O through files

  def compute_descriptors(self,filename,peek_fun=None):
    """ computes the descriptors for the image stored in filename
    returns the name of the bsiftgeo descriptor file """
    
    bsiftgeo=self.datadir+"descs.bsiftgeo"

    fc=self.cachedData.get_FeaturesComputer()

    img=fc.load_image(filename)     
    (n_interest,n_desc,descs)=fc.compute_descriptors(img,peek_fun)
    fc.store_descriptors(descs,bsiftgeo)
     
    return (n_interest,n_desc,bsiftgeo)

  def compute_vw(self,bsiftgeo,peek_f=None):
    """ computes the visual words corresponding to the bsiftgeo file
    returns the visual word file"""
    vwa=self.cachedData.get_VWAssignement()
    outfname=self.datadir+"descs.vw"
    coords=vwa.load_descs(bsiftgeo)
    vw=vwa.compute_vw(coords,peek_f)
    vwa.store_vw(vw,outfname)
    vwa.dealloc_vw(vw)
    return outfname

  def query(self,vwf,peek_f=None):
    bb=self.cachedData.get_BigBase()
    vw=bb.load_vw(vwf)
    res=bb.query(vw,self.n_short,peek_f)
    bb.dealloc_vw(vw)
    return res

  def filter(self,bsiftgeo,res,peek_f=None):
    gv=self.cachedData.get_GeometricVerification()
    imdesc=gv.load_imdesc(bsiftgeo)
    return gv.filter(res,imdesc,peek_f)

  ########################################################
  # Show intermediate results
  
  def draw_pts(self,image_name,bsiftgeo,out_name,max_size=-1):
    """ return an image with interest points drawn on it """
    im=geom_filter.CRGBImage(image_name)
    sz=max(im.getxsize(),im.getysize())
    if max_size>0 and sz>max_size:
      scale=max_size/float(sz)
      im=geom_filter.downscale(im,scale); im.thisown=True
    else:
      scale=1.0
    imdesc=geom_filter.loadImdesc(bsiftgeo,True); imdesc.thisown=True
    geom_filter.drawInterestPoints(
      im,imdesc,scale,geom_filter.CRGBPixel(255,0,0),True)
    im.store(out_name)
  
  def draw_query_descs(self,image_name,bsiftgeo,max_size=-1):
    out_name=self.datadir+"query_descs.png"
    self.draw_pts(image_name,bsiftgeo,out_name,max_size)
    return out_name

  def draw_db_descs(self,dbno,max_size=-1):
    image_name=self.db_im_name(dbno)
    bsiftgeo=self.config.desc_filename(dbno)
    out_name=self.datadir+"db_descs_%d.png"%dbno
    self.draw_pts(image_name,bsiftgeo,out_name,max_size)
    return out_name
      
  def draw_superpose(self,query_image_name,dbno,aff):
    """ draw images superposed with the found affine transform """
    im0=geom_filter.CRGBImage(query_image_name)
    im1=geom_filter.CRGBImage(self.db_im_name(dbno))
    da=geom_filter.DoubleArray(6)
    for i in range(6): da[i]=aff[i]
    geom_filter.drawContrasted(im0,im1,da.cast())
    outname=self.datadir+"superposed_%d.png"%dbno
    im0.store(outname)
    return outname


  ########################################################
  # Chained processing

  def search(self,image_file,state=None):
    """ calls all the functions to do a query
    if image_file is a number, it is assumed to be a db image number
    """
    if state==None: state=DummyStateObj()

    if type(image_file)!=type(1):

      state.next('computing descriptors',image_file=image_file)
      fc=self.cachedData.get_FeaturesComputer()
      image=fc.load_image(image_file)
      (n_int,n_desc,descs)=fc.compute_descriptors(image,state.set_frac)

      state.next('assign VW',n_int=n_int,n_desc=n_desc)
      vwa=self.cachedData.get_VWAssignement()   
      vw=vwa.compute_vw(descs,state.set_frac)

    else: # this is a db image number
      imno=image_name
      
      vwa=self.cachedData.get_VWAssignement()   
      descs=vwa.load_descs(self.db_desc_name(imno))
      
      bb=self.cachedData.get_BigBase()
      vw=bb.load_vw(self.db_vw_name(imno))
    
    state.next('query')
    bb=self.cachedData.get_BigBase()
    res=bb.query(vw,self.n_short,state.set_frac)
    bb.dealloc_vw(vw)

    state.next('filter',res=res)
    gv=self.cachedData.get_GeometricVerification()
    res2=gv.filter([imno for (imno,dist) in res],descs,state.set_frac)
    
    state.next('end',res2=res2)

    return res2

    
  def add_images(self,imlist,gen_thumb=False,state_obj=None):
    """ adds images from a list of files """

    if state_obj==None:
      state_obj=DummyState()
    
    config=self.config
    vws=[]
    nf=len(imlist)
    n_map=self.cachedData.alloc_nos(nf)

    self.log("allocated n_map=%s"%n_map)

    state_obj.next("begin",n_map=n_map,imlist=imlist,thumbs_ok=0)
    
    self.log("begin enumerate")

    for no,im in enumerate(imlist):
      n=n_map[no]
      
      self.log("handle %s %s"%(no,im))

      imname=config.img_filename(n)
      copy_file(im,imname)
      
      state_obj.next("compute descs",imname=im,imno=no)
      (npt,ndesc,desc)=self.compute_descriptors(
        imname,peek_fun=state_obj.set_frac)      
      descname=config.desc_filename(n)
      copy_file(desc,descname)

      state_obj.next("assign VW") 
      vw=self.compute_vw(desc,state_obj.set_frac)
      vwname=config.vw_filename(n)
      copy_file(vw,vwname)
      vws.append(vwname)

      if gen_thumb:
        state_obj.next("compute thumbnail")
        thumbname=config.thumb_filename(n)
        gen_thumbnail(im,thumbname)
    
    state_obj.next("lock invfile",thumbs_ok=1)

    bb=self.cachedData.get_BigBase()
            
    (ivf,ivfm)=bb.get_ivf(True)
    # ivf now blocked
    
    for no,vwf in enumerate(vws):        
      n=n_map[no]
      (svi,sv)=bb.load_vw(vwf)
      invfile_add(ivf,n,svi,sv)
      #invfile_cdm_add(ivfm, n, 1.0)
      
      state_obj.next("add image",imno=no)

      ivfq=bb.get_ivf_query()
      #invfile_cdm_compute_one_iteration(ivf, ivfm, ivfq, self.cdm_ngb)
      factor = invfile_cdm_compute_one_factor(ivfm, ivfq, svi, sv, self.cdm_ngb)
      invfile_cdm_add(ivfm, n, factor)
      bb.release_ivf_query(ivfq)
      bb.dealloc_vw((svi,sv))

    config.nimg+=len(vws)

    # ivfm must be updated
    state_obj.next("update invfile metadata")
    
    bb.prepare_ivfm()
    bb.release_ivf(True) 

    state_obj.next("end")      

    return n_map

  def store_imbase(self,basedir):
    """ store the database in a directory.
    Generates python ocde that can be added to config.py to use the
    stored data """
    
    try_mkdir(basedir)

    bb=self.cachedData.get_BigBase()
   
    (ivf,ivfm)=bb.get_ivf()

    config=self.config

    nim=config.nimg

    print "storing %d-%d"%(config.nfixed,nim)

    datadir=basedir+"/data"
    try_mkdir(datadir)
   
    for i in range(config.nfixed,nim):
      copy_file(config.img_filename(i),  datadir+"/%09d.img"%i)
      copy_file(config.desc_filename(i), datadir+"/%09d.bsifgeo"%i)
      copy_file(config.vw_filename(i),   datadir+"/%09d.vw"%i)
      copy_file(config.thumb_filename(i),datadir+"/%09d.thumb.jpg"%i)

    parentclassname=config.__class__.__name__

    invfile_write(open(basedir+"/extended.ivf", "w"),ivf)
    
    invfile_cdm_factors_write(open(basedir+"/extended.cdmfac.hash","w"),ivfm)
    invfile_cdm_support_write(open(basedir+"/extended.cdmsup.ivf","w"),ivfm)
    
    f=open(basedir+"/extended.py","w")
    f.write("""
    
# add this to your config.py
class Extended%s(%s):
  def __init__(self):
    %s.__init__(self)
    self.extend("%s",%d,%d)
    
"""%(parentclassname,parentclassname,parentclassname,
     basedir,config.nfixed,nim))
    f.close()
    
    bb.release_ivf() 

  def cdm_op(self,opname,args):
    """
    do imbase.cdm_op("reset",(0.2,))
    to call cdm_reset(0.2)
    """
    bb=self.cachedData.get_BigBase()
    f=getattr(bb,'cdm_'+opname)
    f(*args)
    
