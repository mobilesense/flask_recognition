
import re,sys,socket,errno,time,os,pickle,pdb,math
import basics
# mandatory ones

from siftgeo.siftgeo import *
from basics import *
from yael import yael
from utils import utils

import ivfgeo.ivfgeo as ivfgeo

import geom_filter.geom_filter as geom_filter

import os

import pdb

#############################################################################
# Descriptor computation
# classes with load_image() and compute_descriptor()


class FeaturesComputer:
  """ functions to extract features from images Features are put in a
  pointset_t. Subclasses should define load_image() and
  compute_descriptors()  
  """

  def store_descriptors(self,descs,filename):
    write_points(filename,descs.pts,descs.n,0)

  def string_to_pointset(self,stringdescs):
    pointset=pointset_t()
    pointset.name=None
    (pointset.pts,pointset.n)=read_points_string(stringdescs,0)
    return pointset
  
  def load_image(self,fname):
    pass

  def compute_descriptors(self,img,peek_fun=None):
    """ returns (nb_interest,nb_descs,desc_list) """
    pass
    
class ObsidianStage:
  " stage of processing"
  pass

class ObsidianOptimStage(ObsidianStage):
  " stage that has several levels of optimization"
  def __init__(self):
    self.optim='cpu'
    # Alternatives('lava|cpu')

class ObsidianFeaturesComputer(FeaturesComputer):
  """ Computes descriptors with Obsidian code """

  def __init__(self,config):
    " default values "
    self.ins=ObsidianStage()
    self.ins.max_pixel=-1

    # pyramid generation
    self.pyr=ObsidianOptimStage()
    self.pyr.scale_fact=1.2
    self.pyr.min_size=12
    self.pyr.smooth=True

    # detector
    self.det=ObsidianOptimStage()
    self.det.type='hessian'
    self.det.fname=''
    #   Alternatives('harris_laplace|harris|hessian|log')
    self.det.threshold=500.0
    self.det.max_desc=-1

    self.det.mser_min_size=30
    self.det.mser_max_size=0.02
    self.det.mser_delta=10
    self.det.mser_qmax=0.5
    self.det.mser_exfac=1.3

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
    self.desc.scale_one=15.0
    self.desc.index_size=4
    self.desc.T=0.0001
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
    self.parse_args(self.config.det_args)

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
        'normalizer_grad_ori_multiple_%s'%self.aff.optim)
      self.normalizer_grad_ori_multiple=self.obsidian_func(
        'normalizer_grad_ori_multiple_%s'%self.aff.optim)
    else:
      self.normalizer_grad_ori=lambda x: None

    self.descriptor=self.obsidian_func(
      'descriptor_%s_%s'%(self.desc.type,self.det.optim))

  def obsidian_func(self,funcname):
    return obsidian.__dict__[funcname]
      
  def parse_args(self,args):
    """ parse arguments, passes as a string in the form
    "pyr-min_size=45 ori-skip=0"
    NB that booleans must be passed as 0 or 1
    """
    argv=args.split()
    while argv:
      opt=argv.pop(0)
      [k,v]=opt.split('=')
      [stname,parname]=k.split('-')
      stage=getattr(self,stname)
      param=getattr(stage,parname)
      setattr(stage,parname,parse_as_type(type(param),v))

  def compute_descriptors(self,img,peek_fun=None):
    """ returns initial nb of descriptors & descriptor list """
    pstep=1/4.0
    if peek_fun: peek_fun(0)
    
    pyr=self.image_pyramid_build(
      img,self.pyr.scale_fact,self.pyr.min_size,self.pyr.smooth)
    if peek_fun: peek_fun(pstep)
    
    if self.det.type == 'file':
      descs=self.detector(self.det.fname)
    elif self.det.type == 'mser':
      descs=self.detector(pyr[0], self.det.mser_min_size, self.det.mser_max_size, 
                          self.det.mser_delta, self.det.mser_qmax, self.det.mser_exfac)
    else:
      descs=self.detector(pyr,self.det.threshold)
    n_interest=descs.size

    obsidian.local_desc_selector(descs,self.det.max_desc)
    if peek_fun: peek_fun(2*pstep)

    if not self.aff.skip:
      self.normalizer_second_moment(descs,pyr)
    if not self.ori.skip:
      self.normalizer_grad_ori_multiple(descs,pyr)
    if peek_fun: peek_fun(3*pstep)
    
    if self.desc.type=='cs_lbp':
      self.descriptor(descs,pyr,self.desc.window_size,self.desc.scale_one,
                 self.desc.index_size,self.desc.T,self.desc.norm,
                 self.desc.dim)
    elif self.desc.type=='sift':
      self.descriptor(descs,pyr,self.desc.window_size,self.desc.scale_one,
                 self.desc.index_size,self.desc.ori_size)
    if peek_fun: peek_fun(4*pstep)

    # convert from obsidian local_desc_list_t to bigimbaz pointset_t
    stringdescs = obsidian.local_desc_list_as_string(descs)
    pointset=self.string_to_pointset(stringdescs)
    
    assert descs.size==pointset.n
    
    return (n_interest,descs.size,pointset)
    
  def load_image(self,fname):
    if self.config.det_pre_filter:
      tmpname=os.tempnam()
      cmd=self.config.det_pre_filter%(fname,tmpname)
#      print "execing",cmd
      os.system(cmd)
      fname=tmpname
    else:
      tmpname=None
    try:
      if self.ins.max_pixel>0:
        im=obsidian.image_load_max(fname,self.ins.max_pixel)[0]
      else:
        im=obsidian.image_load(fname)
    finally:
      if tmpname: os.unlink(tmpname)
    if not im: raise IOError("bad image %s"%fname)
    return im

  def save_image(self,img,fname):
    " useful ?"
    obsidian.image_save(img,fname)

    
   

class MikoFeaturesComputer(FeaturesComputer):
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
      #pdb.set_trace()
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


#############################################################################
# Visual Word assignement 
# important function: compute_vw


class VWAssignement:
  """ assigns local descriptors to clusters (visual words) """

  def __init__(self,config,n_thread=1):
    print "reading clusters from ",config.clusterfile
    self.config=config
    self.ann_vw=None
    self.n_thread=n_thread
    
    if type(config.clusterfile)==type(''):
      (clustertab,n,d)=yael.fvecs_new_read(config.clusterfile)
    else:
      # ANN assignement
      (cf,coarseclust,edges)=config.clusterfile
      (clustertab,n,d)=yael.fvecs_new_read(cf)
      self.ann_vw=utils.load_ann_vw(coarseclust,edges,n,clustertab)        
      assert self.ann_vw

    if config.binsign_computer:
      self.binsign_computer=siftgeo_binarize_read(config.binsign_computer)      
      if not self.binsign_computer:
        raise IOError("could not read "+config.binsign_computer)
      assert d==self.binsign_computer.d and n==self.binsign_computer.nbvw     
    else:
      self.binsign_computer=None

    assert clustertab
    self.clusters=(clustertab,n,d)
    
    self.k = config.ma_k
    self.ma_disratio = config.ma_disratio
    

  def get_clusters(self):
    return self.clusters

  def load_descs(self, filename, ldc_filename = None):
    descs=pointset_read_cornerness(filename, self.config.vwsgeo_min_cornerness, 
                                   self.config.desc_format)
    if descs.n<0: raise IOError("bad desc file %s"%filename)

    # Update the cornerness
    if ldc_filename:
      (sal, tmp, n) = yael.fvecs_new_read (ldc_filename)

      assert n == descs.n
      sal = yael.FloatArray.acquirepointer (sal)
  
      for i in xrange (descs.n):
        descs[i].geom.cornerness = (descs[i].geom.cornerness * descs[i].geom.cornerness  / 5.0
                                    + sqrt (sal[i]))

      print [int(descs[i].geom.cornerness) for i in xrange (min (8, descs.n))]

    return descs


  def compute_vw(self, descs, peek_f = None):
    """ computes the visual words corresponding to descs """ 
    
    (centroids,nbvw,d)=self.get_clusters()
    k,ma_disratio=self.k,self.ma_disratio
    
    if descs and descs.n>0: # coords may be empty (None)
      npt=descs.n
      (fpoints,d2)=siftgeo_to_fvecs(descs.pts,descs.n)

      assert d==d2

      # output 
      vw=yael.IntArray(npt*k)

      if k==1:
        # only this part is threaded
        if not self.ann_vw:
          yael.quantize_codebook_thread(npt,nbvw,d,centroids,fpoints,vw,self.n_thread,peek_f)
        else:
          utils.quantize_codebook_annvw_thread(self.ann_vw,npt,fpoints,vw,self.n_thread,peek_f)
      else:
        # multiple assignement
        if not self.ann_vw:
          vw_centroid_dis = yael.quantize_codebook_multiple(npt,nbvw,d,k,centroids,fpoints,vw,peek_f)
        else:
          vw_centroid_dis = utils.quantize_codebook_annvw_multiple(self.ann_vw,npt,fpoints,k,vw)
        vw_centroid_dis_fa=yael.FloatArray.frompointer(vw_centroid_dis)
        vw_centroid_dis_fa.this.acquire()

      # build result
      vwgeo=pointset_t()
      vwgeo.pts = siftgeo_to_vwgeo(descs.pts,npt,vw,k)
      vwgeo.n=npt*k

      # add binary signatures for HE computation
      if self.binsign_computer:
        siftgeo_binarize_ffq(self.binsign_computer,fpoints,vwgeo,k)

      # ma_disratio=1 means that we keep all descriptors
      if ma_disratio > 1 and k > 1:
        vwgeoset_filter_ma (vwgeo, k, vw_centroid_dis, ma_disratio)        

      yael.FloatArray.frompointer(fpoints).this.acquire()

    else:
      vwgeo=pointset_t()
      
    return vwgeo
   
  def store_vw(self,pointset,outfname):
    pointset_write(outfname,pointset,self.binsign_computer and 2 or 1)

class DummyVWAssignement:

  def __init__(self,config):
    self.config=config
  
  def compute_vw(self, descs, peek_f = None):
    return descs

  
  def load_descs(self, filename, ldc_filename = None):

    descs=pointset_read_cornerness(filename, self.config.vwsgeo_min_cornerness, 
                                   self.config.desc_format)
    if descs.n<0: raise IOError("bad desc file %s"%filename)
    return descs
  

    
####################################################################
# Inverted files

class BigBase:  
  """ does searches of sparse vectors. Encapsulates a geometric
  inverted file"""
  def __init__(self,config,make_new=False):#,nbvw=0):    
    self.config=config    
    nbvw = self.config.nb_centroids
    if not make_new:
      print "reading inverted file %s"%self.config.invfilename
      self.ivf=ivfgeo.ivfgeo_read(self.config.invfilename,config.ivfgeo_format)
      if not self.ivf: raise RuntimeError("bad invfile")
    else:  
      self.ivf=ivfgeo.ivfgeo_new(nbvw,5)    
    self.ivfgeo_stat=None
    # set decent default parameters
    self.set_params(shortlist_len=20,wgc_type=0x50545,he_thres=22,
                    norm_type=None,scale_w=0.06,he_dist_w=1,
                    burstiness=ivfgeo.BURSTNORM_INTER | ivfgeo.BURSTNORM_INTRA,
                    tfidf=1)
    self.need_vw=True
    
  def __del__(self):
    ivfgeo.ivfgeo_delete(self.ivf)
    
  def set_params(self,shortlist_len=None,wgc_type=None,he_thres=None,
                 norm_type=None,scale_w=None,he_dist_w=None,
                 burstiness=None,tfidf=None):
    if shortlist_len!=None:
      self.shortlist_len=shortlist_len
    if wgc_type!=None:      
      ivfgeo.ivfgeo_set_wgc_type(self.ivf,wgc_type)
    if he_thres!=None:
      self.ivf.he_thres=he_thres
    if norm_type!=None and norm_type!=self.ivf.norm_type:
      self.ivf.norm_type=norm_type
      ivfgeo.ivfgeo_compute_norms(self.ivf)
    if scale_w!=None:
      self.ivf.scale_w=scale_w
      ivfgeo.ivfgeo_compute_scale_weights (self.ivf)
    if he_dist_w!=None:
      self.ivf.he_dist_w=he_dist_w
      ivfgeo.ivfgeo_compute_he_weights (self.ivf)
    if burstiness!=None:
      self.ivf.burstiness=burstiness
    if tfidf!=None:
      self.ivf.tfidf=tfidf
      ivfgeo.ivfgeo_compute_tfidf(self.ivf)
    
  def get_params(self):
    return {'shortlist_len':self.shortlist_len,
            'wgc_type': self.ivf.wgc_type,
            'he_thres':self.ivf.he_thres,
            'norm_type': self.ivf.norm_type,
            'scale_w':self.ivf.scale_w,
            'he_dist_w':self.ivf.he_dist_w,
            'burstiness':self.ivf.burstiness,
            'tfidf':self.ivf.tfidf,
            }
  
  def load_vw(self,vwf):
    pointset=pointset_read(vwf,self.config.binsign_computer and 2 or 1)      
    if pointset.n<0: raise IOError("bad vw file %s"%vwf)
    return pointset

  def add_label_to(self,label,ivf):
    vw=self.load_vw(self.config.vw_filename(label))
    self.add_vw_to(vw,label,ivf)
  
  def get_ivf(self,for_write=False):
    return self.ivf

  def max_label(self):
    return ivfgeo.ivfgeo_max_label(self.ivf)
    
  def release_ivf(self,for_write=False):
    """ useless in single thread context"""    
    pass
   
  def query(self,vw,peek_f=None,n_short=None):
    """ queries a visual word file and returns the corresponding
    image id's + distances """
    if n_short==None:
      n_short=self.shortlist_len
    vwgeoset_sort(vw)  # changes set!        
    ivf=self.get_ivf()    
    (res_i,dis_f)=ivfgeo.ivfgeo_query_peek(ivf,vw,n_short,peek_f)
    ressz=min(n_short,ivf.nbvec)              
    self.release_ivf()
    if not res_i: # when vw is empty
      return []
    res_i=yael.IntArray.frompointer(res_i)
    res_i.this.acquire() # transfer ownership to Python object
    dis_f=yael.FloatArray.frompointer(dis_f)
    dis_f.this.acquire()
    return [(res_i[i],dis_f[i]) for i in range(ressz)]

  def query_all(self,vw,peek_f=None):
    """ returns a C table with scores for each image """    
    vwgeoset_sort(vw)  # changes set!        
    ivf=self.get_ivf()
    (res_i,dis_f)=ivfgeo.ivfgeo_query_peek(ivf,vw,-1,peek_f)    
    self.release_ivf()
    if not dis_f:
      return None    
    assert not res_i
    dis_f=yael.FloatArray.frompointer(dis_f)
    dis_f.this.acquire()    
    return dis_f

  def begin_add(self):
    self.add_d=self.get_ivf(True)
    # ivf now blocked

  def add_vw_to(self,vw,n,ivf):
    """ changes input pointset """
    if self.config.ivfgeo_filter_n_cornerness_max >= 0:
      vwgeoset_filter_n_cornerness_max(vw, self.config.ivfgeo_filter_n_cornerness_max)
    if self.config.ivfgeo_min_cornerness >= 0:
      pointset_filter_cornerness(vw, self.config.ivfgeo_min_cornerness)
    vwgeoset_sort(vw)  
    ivfgeo.ivfgeo_add(ivf,vw,n)

  def add_vw(self,vw,n):
    self.add_vw_to(vw,n,self.add_d)
    
  def end_add(self):
    del self.add_d
    ivfgeo.ivfgeo_compute_tfidf(self.ivf)    
    self.release_ivf(True) 
  
  def write_invfile(self,fname):
    ivfgeo.ivfgeo_write(fname,self.ivf,self.config.ivfgeo_format)

  def display_short_stats(self):
    print "Unbalance factor:",ivfgeo.ivfgeo_unbalanced_factor(self.ivf)
    print "Total entries:",ivfgeo.ivfgeo_count_nbelems(self.ivf)

  def ivf_nim(self,ivf):
    return ivf.nbvec

  def ivf_npt(self,ivf):
    return ivfgeo.ivfgeo_count_nbelems(ivf)    
    
class PQBase:
  """ Version for PQ based queries """

  def __init__(self,config,make_new=False,nbvw=0):
    self.config=config
    (self.clustertab,n,d)=yael.fvecs_new_read(config.clusterfile)
    assert not make_new or n==nbvw,"want %d read %d"%(nbvw,n)
    self.clustertab=yael.FloatArray.acquirepointer(self.clustertab)    
    if make_new:
      print "reading ivf model from ",config.binsign_computer
      ivf=utils.ivfpq_fread(open(config.binsign_computer,"r"),self.clustertab)      
      self.ivf=ivfgeo.ivfgeo_pq_new(ivf)
      self.ivf.nbvw=nbvw
    else: 
      print "reading inverted file %s"%self.config.invfilename
      self.ivf=ivfgeo.ivfgeo_pq_fread(open(self.config.invfilename,"r"),self.clustertab)      
    self.set_params(wgc_type=0x50545,norm_type=None,scale_w=0.06,
                    burstiness=ivfgeo.BURSTNORM_INTER | ivfgeo.BURSTNORM_INTRA,
                    tfidf=1)
    self.need_vw=False
    
  def load_vw(self,filename):
    """ loads a siftgeo file in fact """     
    descs=pointset_read_cornerness(filename, self.config.vwsgeo_min_cornerness, 
                                   self.config.desc_format)
    if descs.n<0: raise IOError("bad desc file %s"%filename)
    return descs
    
  def set_params(self,wgc_type=None,
                 norm_type=None,scale_w=None,dist_w_type=None,
                 burstiness=None,tfidf=None,
                 ma=None,ma_disratio=None,nnn=None,sigma=None):
    if wgc_type!=None:      
      ivfgeo.ivfgeo_pq_set_wgc_type(self.ivf,wgc_type)
    if norm_type!=None and norm_type!=self.ivf.norm_type:
      self.ivf.norm_type=norm_type
      ivfgeo.ivfgeo_pq_compute_norms(self.ivf)
    if scale_w!=None:
      self.ivf.scale_w=scale_w
      ivfgeo.ivfgeo_pq_compute_scale_weights (self.ivf)
    if dist_w_type!=None:
      self.ivf.dist_w_type=dist_w_type
    #if burstiness!=None:
    #  self.ivf.burstiness=burstiness
    if tfidf!=None:
      self.ivf.tfidf=tfidf
      ivfgeo.ivfgeo_pq_compute_tfidf(self.ivf)
    if ma!=None and ma!=0: 
      self.ivf.ma=ma
    if ma_disratio!=None:
      self.ma_disratio=ma_disratio
    if nnn!=None:
      self.ivf.k=nnn
    if sigma!=None:
      self.ivf.sigma=sigma
        
  def get_params(self):
    return {'wgc_type': self.ivf.wgc_type,
            'norm_type': self.ivf.norm_type,
            'scale_w':self.ivf.scale_w,
            'dist_w_type':self.ivf.dist_w_type,
            'burstiness':self.ivf.burstiness,
            'tfidf':self.ivf.tfidf,
            'ma': self.ivf.ma,
            'ma_disratio': self.ivf.ma_disratio,
            'nnn': self.ivf.k,
            'sigma': self.ivf.sigma
            }
 
  def begin_add(self):
    pass

  def add_vw(self,vw,label):
    label2=ivfgeo.ivfgeo_pq_add(self.ivf,vw)
    assert label==label2

  def add_label_to(self,label,ivf):
    siftgeo=self.load_vw(self.config.desc_filename(label))
    vwf=self.config.vw_filename(label)
    if vwf==None:
      label2=ivfgeo.ivfgeo_pq_add(ivf,siftgeo)
      # assume label2 is ok
    else:
      # load wv file
      ps=pointset_read(vwf,self.config.binsign_computer and 2 or 1)      
      if ps.n<0: raise IOError("bad vw file %s"%vwf)
      # extract vw's
      vw=yael.IntArray.acquirepointer(vw_from_pointset(ps))
      label2=ivfgeo.ivfgeo_pq_add_with_vw(ivf,siftgeo,vw)          

  def end_add(self):
    ivfgeo.ivfgeo_pq_compute_tfidf(self.ivf)   

  def write_invfile(self,fname):
    ivfgeo.ivfgeo_pq_fwrite(self.ivf,open(fname,"w"))

  def get_ivf(self):
    return self.ivf
  
  def ivf_nim(self,ivf):
    return ivf.nim

  def ivf_npt(self,ivf):
    return utils.ivfpq_count_nbelems(ivf.ivfpq) 

  def release_ivf(self):
    pass

  def query_all(self,vw,peek_f=None):
    """ returns a C table with scores for each image """            
    ivf=self.get_ivf()
    scores=yael.FloatArray(ivf.nim)   
    ivfgeo.ivfgeo_pq_query(ivf,vw,scores)    
    self.release_ivf()
    return scores

  def __del__(self):
    # ivfgeo.ivfgeo_pq_delete(self.ivf)
    pass

  def display_short_stats(self):
    print "Unbalance factor:",ivfgeo.ivfgeo_pq_unbalanced_factor(bb.ivf)
    print "Total entries:",ivfgeo.ivfgeo_pq_count_nbelems(bb.ivf)

####################################################################
# Geometric verification
# important function: filter() 

class ImageMatch:
  """ mirrors a C imageMatch """
  def __init__(self,imm,includePtMatches=False):
    for a in ('stage1_nmatch','stage2_nmatch','stage2_votes','stage3_nmatch',
              'stage3_votes','final_votes'):
      setattr(self,a,getattr(imm,a))
    self.aff=self.aff_to_py(imm.affine)
    if includePtMatches!=False:      
      x = geom_filter.PointMatchArray.frompointer(imm.ptmatches)      
      # interest points matched list 
      self.matches = [(x[i].dbpt,x[i].qpt,x[i].score) for i in range(geom_filter.count_ptmatches(imm.ptmatches)) if x[i].score>0]
      assert len(self.matches)==self.stage3_nmatch #we have to include only won points after stage3 -> imm have to be updated
            
  def aff_to_py(self,sd):
    da=yael.DoubleArray.frompointer(sd)
    return [da[i] for i in range(6)]

  def __cmp__(self,other):
    " sort by increasing stage3_votes "
    return -cmp((self.stage3_votes,self.stage3_nmatch),
                (other.stage3_votes,other.stage3_nmatch))

  def __str__(self):
    # return "matches=%d,%d scores=%g,%g matrix=[%g %g %g; %g %g %g]"%(
    return "matches=%d,%d scores=%g,%g matrix=[%.2f %.2f %.0f; %.2f %.2f %.0f]"%(
      self.stage2_nmatch,self.stage3_nmatch,
      self.stage2_votes,self.stage3_votes,
      self.aff[0],self.aff[1],self.aff[2],
      self.aff[3],self.aff[4],self.aff[5])
    
class GeometricVerification:
  """ does point-based image matching for a short list of images """
  def __init__(self,config,n_thread=1):
    self.config=config
    self.lhp=geom_filter.lowehough_parameters_t()   
    #decent defaults
    self.set_params(geom_flags=0x408,thresh=22.0,fthresh=0.0)
    # fraction of computing time for loading (used for progress bar)
    self.fracs=Object()
    self.fracs.load=0.8
    self.fracs.match=0.1
    self.fracs.filter=0.1
    self.n_thread=n_thread
  
  def sync_with_geom_flags(self):
    """
    geom_flags is a bitfield:

    Point matching flags:
    2 :    crop pointsets to 1000 points
    4 :    do exact point matching
    8 :    match from vw's (query must be vwgeo), thresh is the Hamming distance threshold
    0x10:  match with thresh nearest neighbours
    0x20:  match with thresh nearest neighbours from the whole db

    Hough transform flags:
    0x40:  don't use Hough transform: test all hypotheses
    0xy00: Hough parameters
      y=0: default Hough parameters (good for 320*240 pixel images)
      y=4: for 1024*768
      y=5: more strict verification
      
    0x1000: use weighting 
    """
    lhp=self.lhp
    geom_filter.lowehough_parameters_default(lhp)    
    detail=(self.geom_flags>>8)&0xff
    if detail==4:
      lhp.pos_error_in_affine=25
      lhp.distinct_tolerance=15
      lhp.min_match_before_affine=3
      lhp.min_match_after_affine=3
      bin_sizes=yael.DoubleArray.frompointer(lhp.bin_sizes)
      bin_sizes[0]=math.log(6)
      bin_sizes[1]=math.pi/18.0        
      bin_sizes[2]=80
      bin_sizes[3]=80
    elif detail==5:
      lhp.pos_error_in_affine=10
      lhp.distinct_tolerance=15
      lhp.min_match_before_affine=3
      lhp.min_match_after_affine=4
      bin_sizes=yael.DoubleArray.frompointer(lhp.bin_sizes)
      bin_sizes[0]=math.log(6)
      bin_sizes[1]=math.pi/18.0        
      bin_sizes[2]=80
      bin_sizes[3]=80
    lhp.weight_deformation=self.geom_flags & 0x1000
    lhp.verbose=0
       
  def set_params(self,geom_flags=None,thresh=None,fthresh=None):
    """ set parameters that may be modified after __init__ """
    if thresh!=None:
      self.thresh=thresh
    if fthresh!=None:
      self.fthresh=fthresh
    if geom_flags!=None:
      self.geom_flags=geom_flags
      self.sync_with_geom_flags()
            
  def get_params(self):
    return {'geom_flags':self.geom_flags,
            'thresh':self.thresh,
            'fthresh':self.fthresh}    
  
  def crop_pointset(self,pointset,mem):
    nmax=1000	#TODO Amine, changes if working with logos (we use lots of interest points)
    if pointset.n>nmax:
      mem[pointset]=pointset.n
      pointset_sort_by_cornerness(pointset)
      pointset.n=nmax

  def undo_crop_pointset(self,mem):
    for k,v in mem.iteritems():
      k.n=v

  def make_sl(self,pa,n,mem):    
    if self.geom_flags & 2:
      for i in range(n):
        self.crop_pointset(pa[i],mem)
    return geom_filter.shortlist_new(pa,n)
  
  """def do_query(self,sl,n,query,peek_f=None):
    #if matching VWs, changes input query 
    if peek_f:
      match_peek_fun=lambda frac: peek_f(self.fracs.load+self.fracs.match*frac)
      filter_peek_fun=lambda frac: peek_f(self.fracs.load+self.fracs.match+self.fracs.filter*frac)
    else:
      match_peek_fun=filter_peek_fun=None
    if self.matches_vw():
      vwgeoset_sort(query)
    # match points    
    if self.geom_flags & 0x34:
      if self.geom_flags & 4:    method=0
      elif self.geom_flags&0x10: method=1
      elif self.geom_flags&0x20: method=2
      else: assert False
      imms=geom_filter.shortlist_match_points_exact(sl,query,method,self.thresh,match_peek_fun)
    elif self.geom_flags & 8:
      imms=geom_filter.shortlist_match_points_vw(sl,query,int(self.thresh),match_peek_fun)
    else:
      assert False,"KDTree deprecated"
    imms=geom_filter.ImageMatchArray.frompointer(imms)      
    # filter with Lowe's Hough transform
    the_filter=(self.geom_flags & 0x40 and geom_filter.shortlist_filter_allhyp or geom_filter.shortlist_filter_lowehough)
    the_filter(sl,query,imms,self.lhp,filter_peek_fun)
    # convert C structs to Python
    res2=[(ImageMatch(imms[i]),i) for i in range(n)]
    geom_filter.imagematches_delete(imms,n)
    # keep only successful images and sort by relevance
    res2=[(imm,i) for imm,i in res2 if imm.stage3_votes>0]
    res2.sort()
    return res2"""

  def load_imdesc(self,filename):
    if self.matches_vw():
      ps=pointset_read(filename,self.config.binsign_computer and 2 or 1)
      vwgeoset_sort(ps)
    else:
      ps=pointset_read(filename,self.config.desc_format)      
    if ps.n<0: raise IOError("bad descriptor file %s"%filename)
    return ps

  def ps_filename(self,imno):
    return self.matches_vw() and self.config.vw_filename(imno) or self.config.desc_filename(imno)

  def load_imno(self,imno):
    return self.load_imdesc(self.ps_filename(imno))
  
  def matches_vw(self):
    """ do we need VW's or descriptors at input?"""
    return bool(self.geom_flags & 8)

  def add_to_pa(self,i,j,pa):
    ps=self.load_imno(j)
    vwgeoset_sort(ps)
    pa[i]=ps    

  def load_dbimages(self,short_list,peek_f=None):
    n=len(short_list)
    pa=PointsetArray(n)
    if self.n_thread==1:
      for i,j in enumerate(short_list):
        self.add_to_pa(i,j,pa)
        if peek_f: peek_f(i*self.fracs.load/n)
    else:
      # multithreading load seem useful even when everything is on the same disk
      from utils.thread_utils import RunOnSetWithPeek
      sub_peek_f=peek_f and (lambda x: peek_f(x*self.fracs.load))
      RunOnSetWithPeek(self.n_thread,list(enumerate(short_list)), lambda (i,j): self.add_to_pa(i,j,pa),
                       peek_fun=sub_peek_f)        
    return pa


  # --------------------------------- Added ---------------------------------------
  # -------------------------------------------------------------------------------  
  
  def match_points(self,sl,query,peek_f=None):
    """ if matching VWs, changes input query """    
    if peek_f:
      match_peek_fun=lambda frac: peek_f(self.fracs.load+self.fracs.match*frac)
    else:
      match_peek_fun=None
    # match points    
    if self.geom_flags & 0x34:
      if self.geom_flags & 4:    method=0
      elif self.geom_flags&0x10: method=1
      elif self.geom_flags&0x20: method=2
      else: assert False
      imms=geom_filter.shortlist_match_points_exact(sl,query,method,self.thresh) #imms is list of image matches, len(imms)=sl->nim
    elif self.geom_flags & 8:
      vwgeoset_sort(query)
      imms=geom_filter.shortlist_match_points_vw(sl,query,int(self.thresh),match_peek_fun)
    else:
      assert False
    imms=geom_filter.ImageMatchArray.frompointer(imms)
    # has to be deallocated with geom_filter.imagematches_delete(imms,n)
    return imms

  # ----
  
  def match_images(self,sl,n,query,imms,peek_f=None,includePtMatches=False):
    if peek_f:
      filter_peek_fun=lambda frac: peek_f(self.fracs.load+self.fracs.match+self.fracs.filter*frac)
    else:
      filter_peek_fun=None          
    # filter with Lowe's Hough transform
    the_filter=(self.geom_flags & 0x40 and geom_filter.shortlist_filter_allhyp or geom_filter.shortlist_filter_lowehough)
    #TODO Amine, add function to geom_filter   
    
    the_filter(sl,query,imms,self.lhp,filter_peek_fun)
    
    # TODO imms should be updated after filter , each imm have to contain only winner points     
    
    # convert C structs to Python
    res2=[(ImageMatch(imms[i],includePtMatches),i) for i in range(n)]    
    
    # keep only successful images and sort by relevance
    res2=[(imm,i) for imm,i in res2 if imm.stage3_votes>0]
    res2.sort()
    return res2
        
  def filter(self, short_list, query, peek_f=None, includePtMatches=False):
    """ filters a result from query() with a geometric filter
    returns an ordered list of tuples (imno,nmatch,score,aff), where:
      imno   = image number
      nmatch = # of matching points
      score  = computed from matching points and distances
      aff    = affine matrix
    if matches_vw the order of the query changes!
    input: short_list: list of image indexes returned by inverted file (after quering)
    """
    n=len(short_list)
    if n==0: return []

    pa=self.load_dbimages(short_list,peek_f) #pa is the descs of short_list indexes
    mem={}
    if self.geom_flags & 2:
      self.crop_pointset(query,mem)     
    sl=self.make_sl(pa,n,mem) #sl is short list structure it contain descs of each image, npt per image see filter_shortlist.c 
    
    #  def match_points(self,sl,query,peek_f=None):
    imms=self.match_points(sl,query,peek_f)
    
    res2=self.match_images(sl,n,query,imms,peek_f,includePtMatches)

    geom_filter.imagematches_delete(imms,n)
    geom_filter.shortlist_delete(sl)

    ptmatches = []
    if includePtMatches:
      # include matched point coordinates (xdb,ydb) --> (xq,yq)
      res = []
      qpts = query # SWIG will wrap qpts[0] to query.pts[0]
      for (imm,imno) in res2:
        dbpts = pa[imno] # SWIG will wrap dbpts[0] to pa[imno].pts[0]
        ptmatches_xycoords = []
        for i,m in enumerate(imm.matches):      	
          ptmatches_xycoords.append( ((dbpts[m[0]].geom.x, dbpts[m[0]].geom.y), (qpts[m[1]].geom.x, qpts[m[1]].geom.y), m[2]) ) # last value = match score                    
        ptmatches.append(ptmatches_xycoords)        
    return [(short_list[imno],imm.stage3_nmatch,imm.stage3_votes,imm.aff) for (imm,imno) in res2], ptmatches
  # -------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------
        
############################################################################
# CachedData

class CachedData:
  """ Manages the processing stage objects. One instance of each is created,
  lazily, so that the ImBase object can be as shallow as possible.
  """
  def __init__(self,config,n_thread=1):
    self.config=config
    self.n_thread=n_thread    
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
      if self.config.he_type=='pq':
        self.the_VWAssignement=DummyVWAssignement(self.config)
      else: 
        self.the_VWAssignement=VWAssignement(self.config,self.n_thread)
    return self.the_VWAssignement

  def get_BigBase(self,**kwargs):
    if self.the_BigBase==None:
      if self.config.he_type=='pq':
        self.the_BigBase=PQBase(self.config,**kwargs)
        self.the_BigBase.vwa=self.get_VWAssignement()
      else: 
        self.the_BigBase=BigBase(self.config,**kwargs)        
    return self.the_BigBase

  def get_GeometricVerification(self):
    if self.the_GeometricVerification==None:
      self.the_GeometricVerification=GeometricVerification(self.config,self.n_thread)
    return self.the_GeometricVerification

  def alloc_nos(self,nf):
    """ allocates nf new image numbers """
    if self.config.nalloc<0: self.config.nalloc=self.config.nimg
    n0=self.config.nalloc
    self.config.nalloc+=nf
    return range(n0,self.config.nalloc)  
      
  def cache_all(self):
    "preload everything"
    self.fc = self.get_FeaturesComputer()
    self.vwa = self.get_VWAssignement()
    if not os.path.isfile(self.config.invfilename):
        self.bb = self.get_BigBase(make_new=True)
    else:
        self.bb = self.get_BigBase()
    self.gv = self.get_GeometricVerification()
  
  
class ImBase:
  """
  Implements all the image database functions. 
  """ 
  def __init__(self,cachedData, user_id):
    self.user_id = user_id
    self.cachedData=cachedData
    self.config=cachedData.config
    self.datadir="data/"

  def log(self,s):
    print "%s: %s\n"%(time.strftime("%H:%M:%S"),s)
  
  ########################################################
  # Step-per-step image search

  def compute_desc_scale(self,(xd,yd)):
    """ infer the scale applied to descriptors. Image pixels should be
    multiplied by this to obtain descriptor coordinates"""
    self.log("pre_scale=%s"%(self.config.det_pre_scale,))
    if self.config.det_pre_scale:
      assert self.config.det_pre_scale[0]=='maxpix'
      maxpix=self.config.det_pre_scale[1]
      self.log("maxpix=%d xd=%d yd=%d"%(maxpix,xd,yd))
      if xd*yd>maxpix:
        return math.sqrt(maxpix/float(xd*yd))
    return 1.0

  def loaded_image_size(self,fname):
    "image size after rescaling and loading"
    (xd,yd,ok)=obsidian.image_get_size(fname)
    if not ok: raise RuntimeError("%s bad image"%fname)
    scale=self.compute_desc_scale((xd,yd))    
    return (int(xd*scale),int(yd*scale))
  
  def draw_pts(self,image_name,bsiftgeo,out_name,max_size=-1,desc_format=0): 
    """ return an image with interest points drawn on it """
    im=obsidian.color_image_load(image_name)
    sz=max(im.width,im.height)
    descscale=self.compute_desc_scale((im.width,im.height))
    if max_size>0 and sz>max_size:
      scale=max_size/float(sz)
      im=obsidian.color_image_resize(im,1/scale);
    else:
      scale=1.0
    self.log("scale=%g descscale=%g"%(scale,descscale))
    desc_list = obsidian.local_desc_list_read_type(bsiftgeo,desc_format);
    if not desc_list:
      raise RuntimeError("descriptors not stored")
    obsidian.local_desc_list_draw_features(
      im,desc_list,scale/descscale, obsidian.color_new(255,0,0),2)
    obsidian.color_image_save(im, out_name)
  
  def draw_query_descs(self,image_name,bsiftgeo,max_size=-1):
    out_name=self.datadir+"query_descs.jpg"
    self.draw_pts(image_name,bsiftgeo,out_name,max_size)
    return out_name

  def draw_db_descs(self,dbno,max_size=-1):
    image_name=self.config.img_filename(dbno)
    bsiftgeo=self.config.desc_filename(dbno)
    out_name=self.datadir+"db_descs_%d.jpg"%dbno
    self.draw_pts(image_name,bsiftgeo,out_name,max_size)
    return out_name
      
  def draw_superpose(self,query_image_name,mname,aff,reverted=False):
    """ draw images superposed with the found affine transform """
    im0 = obsidian.color_image_load(query_image_name)
    im1 = obsidian.color_image_load(mname)    
    desc0scale=self.compute_desc_scale((im0.width,im0.height))    
    desc1scale=self.compute_desc_scale((im1.width,im1.height))
    aff=aff_pre_post_scale(aff,desc0scale,1.0/desc1scale)
    if reverted:
      (im1,im0)=(im0,im1)
      aff=invert_aff(aff)    
    da=yael.FloatArray(6)
    for i in range(6): da[i]=aff[i]  
    obsidian.color_image_draw_contrasted(im0,im1,da)
    outname=self.datadir+"superposed.jpg"
    obsidian.color_image_save(im0, outname)
    return outname


  ########################################################
  # Search in a single step

  # -- search 
  def search(self,im,state=None,includePtMatches=False):
    """ calls all the functions to do a query
    if im is a number, it is assumed to be a db image number
    """
    if state==None: state=DummyState()
    config=self.config
    # several ways of computing the visual words
    if type(im)==type(''):
      state.next('computing descriptors',image_file=im)
      #fc=self.cachedData.get_FeaturesComputer()
      fc = self.cachedData.fc
      image=fc.load_image(im)
      (n_int,n_desc,descs)=fc.compute_descriptors(image,state.set_frac)
      if self.datadir:
        # we store, to be able to display descriptors on query images
        siftgeo_name=self.datadir+"descs.siftgeo"
        fc.store_descriptors(descs,siftgeo_name)
      else:
        siftgeo_name="not_stored"      
      with basics.TicToc('assign vw'):
		  state.next('assign VW',n_int=n_int,n_desc=n_desc,siftgeo=siftgeo_name)	  
		  #vwa=self.cachedData.get_VWAssignement()   
		  vwa = self.cachedData.vwa
		  vw=vwa.compute_vw(descs,None)#state.set_frac)
    elif type(im)==type(1): # im is a db image number
      state.next('computing descriptors',image_file=config.img_filename(im))      
      #vwa=self.cachedData.get_VWAssignement()
      vwa = self.cachedData.vwa
      #gv=self.cachedData.get_GeometricVerification()      
      gv = self.cachedData.gv
      # load descs only if needed
      descs=not gv.matches_vw() and vwa.load_descs(config.desc_filename(im)) or None      
      state.next('assign VW',n_int=0,n_desc=descs and descs.n or 0,siftgeo=config.desc_filename(im))
      #bb=self.cachedData.get_BigBase()
      bb = self.cachedData.bb[self.user_id]
      vw=bb.load_vw(config.vw_filename(im))
    else: # im is an obsidian image struct
      state.next('computing descriptors',image_file="not_stored")
      #fc=self.cachedData.get_FeaturesComputer()
      fc = self.cachedData.fc
      (n_int,n_desc,descs)=fc.compute_descriptors(im,state.set_frac)
      state.next('assign VW',n_int=n_int,n_desc=n_desc,siftgeo="not_stored")
      #vwa=self.cachedData.get_VWAssignement()
      vwa = self.cachedData.vwa
      vw=vwa.compute_vw(descs,state.set_frac)      
    with basics.TicToc('query vw'):      
	    (res, res2, ptmatches) = self.query_vw(vw, descs, state, includePtMatches)    		#key function
    state.next('end',res2=res2)
    return (res,res2,ptmatches)
      

  # core function, from query's vw and descs compute matching
  def query_vw(self, vw, descs, state_obj=None,includePtMatches=False):
    if state_obj==None:
      state_obj=DummyState()  
    state_obj.next('query')
    #bb=self.cachedData.get_BigBase()
    bb = self.cachedData.bb[self.user_id]
    res=bb.query(vw,None)#state_obj.set_frac)
    res=[(i,dis) for i,dis in res]
    state_obj.next('filter',res=res)
    #gv=self.cachedData.get_GeometricVerification()
    gv = self.cachedData.gv
    gv_input=gv.matches_vw() and vw or descs    
    # filter using gem verification
    config=self.config
    res2,ptmatches = gv.filter([r[0] for r in res],gv_input,None,includePtMatches)		
    return (res, res2, ptmatches)
    

  ########################################################
  # TODO Remove images from db
  
  # Adding images to db
  def add_images_files(self, imnames, n_map=None, keep_files="idvt", state_obj=None):
    """ adds images from a list of files """
    #fc=self.cachedData.get_FeaturesComputer()
    fc = self.cachedData.fc
    imlist = []
    for imname in imnames:
      imlist.append(fc.load_image(imname))
    if state_obj:
      state_obj.next("loaded images",imlist=imnames)    
    n_map = self.add_images(imlist, n_map=n_map, keep_files=keep_files, state_obj=state_obj)
    #make_invfile(config,begin=None,end=None,n_thread=1)
    return n_map

  def add_images(self, imlist, n_map=None, keep_files="idvt", state_obj=None):
    """ adds images from a list of obsidian objects
    imlist must contain images that can be handled by the FeaturesComputer
    'i' in keep_files : save image file
    'd' in keep_files : save descriptor
    'v' in keep_files : save VW file
    't' in keep_files : save thumbnail
    """
    if state_obj==None:
      state_obj=DummyState()    
    config=self.config
    if not n_map:      
      n_map=self.cachedData.alloc_nos(len(imlist))
      pdb.set_trace()
      self.log("allocated n_map=%s"%n_map)
    #fc=self.cachedData.get_FeaturesComputer()
    fc = self.cachedData.fc
    #vwa=self.cachedData.get_VWAssignement()
    vwa = self.cachedData.vwa
    vws=[]
    state_obj.next("begin",n_map=n_map,thumbs_ok=0)
    self.log("begin enumerate")
    # compute descriptors and visual words
    for no,im in enumerate(imlist):
      n=n_map[no]
      self.log("handle %s %s"%(no,im))
      if 'i' in keep_files:
        imname=config.img_filename(n)
        fc.save_image(im,imname)      
      state_obj.next("compute descs",imno=no)
      (n_interest,n_desc,descs)=fc.compute_descriptors(im,state_obj.set_frac)
      if 'd' in keep_files:
        descname=config.desc_filename(n)
        fc.store_descriptors(descs,descname)
      state_obj.next("assign VW") 
      vw=vwa.compute_vw(descs,state_obj.set_frac) #sparse vw vector
      vws.append(vw)
      if 'v' in keep_files:
        vwname=config.vw_filename(n)
        vwa.store_vw(vw,vwname)
      if 't' in keep_files:
        state_obj.next("compute thumbnail")
        thumbname=config.thumb_filename(n)
        gen_thumbnail(imname,thumbname)
    state_obj.next("insert VWs in invfile")
    self.add_vws(n_map, vws, state_obj)    
    state_obj.next("end")
    return n_map

  def add_vws(self, n_map, vws, state_obj=None):    
    if state_obj==None:
      state_obj=DummyState()
    state_obj.next("lock invfile",thumbs_ok=1)
    #bb=self.cachedData.get_BigBase()            
    bb = self.cachedData.bb
    bb.begin_add()    
    # insert the sparse visual word vectors into the inverted file
    for no,vw in enumerate(vws):        
      n=n_map[no]
      state_obj.next("add image",imno=no)
      bb.add_vw(vw,n)
    self.config.nimg+=len(vws)
    # ivfm must be updated
    state_obj.next("update invfile metadata")
    bb.end_add()
    #TODO check it
    print "-> writing invfile ", self.config.invfilename
    prepare_dir(self.config.invfilename)
    bb.write_invfile(self.config.invfilename)
    	
  """def make_invfile(self,begin=None,end=None,n_thread=1):
    nbvw=self.config.nb_centroids
    bb=CachedData(self.config).get_BigBase(make_new=True,nbvw=nbvw) 
    print "nbvw=%d, adding images %d to %d (to be stored in %s)"%(
      bb.ivf.nbvw,begin,end,self.config.invfilename)
    bb.begin_add()
    # if n_thread>1:
    if True: 
      if n_thread>end-begin: n_thread=end-begin
      if self.config.he_type=='pq':
        res=[ivfgeo.ivfgeo_pq_dup(bb.ivf) for i in range(n_thread)]                    
      else: 
        res=[ivfgeo.ivfgeo_new(nbvw,10) for i in range(n_thread)]
      RunOnSet(n_thread,range(n_thread),lambda i: make_inv_slice(bb,begin,end,res,i))
      for b,e,ivf in res:
        print "merge %d-%d"%(b,e)
        if self.config.he_type=='pq':
          ivfgeo.ivfgeo_pq_merge(bb.ivf,ivf)
        else: 
          ivfgeo.ivfgeo_merge(bb.add_d,ivf,0)
    else:
      for label in range(begin,end):
        vwf=self.config.vw_filename(label) 
        vw=bb.load_vw(vwf)
        bb.add_vw(vw,label)
        print "\r",label,vwf[-40:],
        sys.stdout.flush()
    # ivfgeo.ivfgeo_pq_display(bb.ivf) 
    bb.end_add()
    print "writing invfile ",self.config.invfilename
    prepare_dir(self.config.invfilename)
    bb.write_invfile(self.config.invfilename)"""
  
  def store_imbase(self,basedir):
    """ store the database in a directory.
    Generates python ocde that can be added to config.py to use the
    stored data """
    config=self.config    
    bb=self.cachedData.get_BigBase()
    nim=config.nimg
    print "storing %d-%d"%(config.nfixed,nim)
    try_mkdir(basedir)
    datadir=basedir+"/data"
    try_mkdir(datadir)
    for i in range(config.nfixed,nim):
      if os.path.exists(config.img_filename(i)):
        copy_file(config.img_filename(i),  datadir+"/%09d.jpg"%i)
      if os.path.exists(config.desc_filename(i)):
        copy_file(config.desc_filename(i), datadir+"/%09d.bsifgeo"%i)
      if os.path.exists(config.vw_filename(i)):
        copy_file(config.vw_filename(i),   datadir+"/%09d.vw"%i)
      if os.path.exists(config.thumb_filename(i)):
        copy_file(config.thumb_filename(i),datadir+"/%09d.thumb.jpg"%i)
    parentclassname=config.__class__.__name__
    bb.write_invfile(basedir+"/extended.ivf")   
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
    



