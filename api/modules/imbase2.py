

import sys,socket,errno,time,os,pdb

from invfile import *

import geom_filter

import obsidian

from utils import *

#############################################################################
# the four processing stages for a search

class MikoFeaturesComputer:
  """ Computes features with Krystian Mikolajczyk's programs.
  Images and descriptors are passed with files """

  def __init__(self,config):
    self.config=config
  
  def exec_miko(self,cmd,peek_fun=None):
    """ execute one of Krystian's description programs and read n_desc
    and n_interest from the stdout"""

    while True:
      try:
        f=os.popen(cmd,'r')
      except OSError,e:
        # ERESTARTNOINTR raised randomly by FC6 kernel
        if e.errno==513:
          print >>sys.stderr,"exec_miko got ERESTARTNOINTR, retrying..."
          continue
        raise e
      break
    
    n_interest=None
    n_desc=None
    line_no=0
    for l in f:
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
      print >>sys.stderr,"error in exec_miko: %s "%ret
      
      if (ret&0xff)==0 and (ret>>8)==2:
        print >>sys.stderr,"Assuming image is too small"
        return (-1,-1)
               
      raise RuntimeError("error in exec_miko: %s "%ret)

    return (n_interest,n_desc)
  
  def compute_descriptors_files(self,fname,descfile,peek_fun=None):
    """ computes the descriptors for the image stored in filename
    returns the name of the bsiftgeo descriptor file """
    if self.config.det_pre_filter:
      tmpname=os.tempnam()
      cmd=self.config.det_pre_filter%(fname,tmpname)
      # print "execing",cmd
      os.system(cmd)
      fname=tmpname
    else:
      tmpname=None      
    try:
      cmd="%s %s -i2 %s -o4 %s"%(self.config.miko_exec,self.config.det_args,fname,descfile)
      (n_interest,n_desc)=self.exec_miko(cmd,peek_fun)      
    finally:
      if tmpname and os.access(tmpname,os.R_OK):
        os.unlink(tmpname)      
    return (n_interest,n_desc)
  
  def load_image(self,fname):
    " dummy load" 
    return fname

  def save_image(self,img,fname):
    copy_file(img, fname)

  def compute_descriptors(self,img,peek_fun=None):
    """ returns initial nb of descriptors & descriptor list """
    descfile=os.tempnam()
    try:
      (n_interest,n_desc)=self.compute_descriptors_files(img,descfile,peek_fun)
      if n_interest==-1:
        pointset=pointset_new()
      else:
        stringdescs=open(descfile,"r").read()
        pointset=self.string_to_pointset(stringdescs)      
        assert n_desc==pointset.n
    finally:
      if os.access(descfile,os.R_OK):
        os.unlink(descfile)
      
    return (n_interest,n_desc,pointset)


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
    
