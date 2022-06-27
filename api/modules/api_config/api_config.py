import app_config
import os
import pdb

class OverallConfig:
  """
  references to all datafiles used by the database.
  For each image there is (standard extension in brackets):
  - the image file (img)
  - the descriptor in Miko's -o4 format (bsifgeo)
  - the visual word vector (vwgeo)
  - the thumbnail (thumb.jpg)
  Images fall in 3 classes:
  - fixed initial: they were in the initial database
  - fixed extended: they were added in previous sessions
  (with an ExtendedXXXConfig class)
  - allocated: added in the current session. 
  The fixed images have numbers 0..n_fixed-1
  Allocated ones hase n_fixed..n_fixed+n_alloc
  """
  def __init__(self):
    self.dbname=self.__class__.__name__[:-6]
    # number of images to index
    self.nimg=0    
    # detector/descriptor parameters
    # det_prog=='obsidian': use obisidan detector/descriptor
    # det_prog=='compute_descriptors': use Miko's 
    self.det_prog='obsidian'
    # for compute_descriptors:
    #   copied verbatim to command line
    #   typical: "-hesaff -thres 100 -sift -max 1000"
    # for obsidian:
    #   cf. ObsidianFeaturesComputer
    #   typical: "pyr-min_size=45 ori-skip=0"
    self.det_args=''
    # pre-filter image with this command before loading
    # should take 2 arguments (as %s %s): in and out file
    self.det_pre_filter=''
    # used for display: if image is prescaled (in det_pre_filter or with
    # Obisidian's ins-max-pixel), set this to:
    # ('maxpix',12345): image is scaled uniformly so that # pixels < 12345
    self.det_pre_scale=None
    # where to find Krystian's compute_descriptors    
    self.miko_exec=app_config.APP_BASE_DIR+'/api/modules/compute_descriptors/compute_descriptors/compute_descriptors'
    # format of the descriptor files (may be 3 for .siftbin files)
    self.desc_format=0
    
    # clusters
    # cluster file formats:
    # 'fvecs': loaded with load_clusters_fvecs    
    self.cluster_format='fvecs'
    self.nb_centroids=20000    
    self.ma_k=1
    self.ma_disratio=1.2
    self.vwsgeo_min_cornerness = -1
    self.ldc = False

    # a siftgeo_binarize file to compute signature. 
    self.binsign_computer=''
    # Method used to generate the binsign file: 'pca'/'nopca'/'pq'
    self.he_type='nopca'
    
    # invfile parameters 
    self.ivfgeo_format=0
    self.invfilename=''
    # add at most this many points per image to invfile (selected by cornerness)
    self.ivfgeo_filter_n_cornerness_max=-1
    # idem, but select by cornerness threshold
    self.ivfgeo_min_cornerness=-1 

    # extensions: bookkeeping for added files
    self.extensions=[]
    self.add_dir=None
    self.nalloc=-1
    
    self.set_default_compute_descriptors()  #hessian affine parameters
    self.basedir = app_config.APP_BASE_DIR + app_config.DATA_FOLDER+"/"     
    self.clusterfile=self.basedir+"clust/clust_holidays_k20000.fvecs"    
    self.binsign_computer=self.basedir+"clust/params_holidays_k20000.binsign"

  def img_filename(self,n):
    return None
  
  def thumb_filename(self,n):
    return None

  def desc_filename(self,n):
    return None

  def vw_filename(self,n):
    return None 

  def img_url(self,n):
    return None

  def display_name(self,n):
    return self.img_filename(n).split('/')[-1]

  def set_default_compute_descriptors(self):
    self.det_prog='compute_descriptors'
    self.det_args="-hesaff -thres 500 -sift"
    """self.det_pre_filter=(
      'djpeg "%s" | ppmtopgm | '+
      'pnmnorm -bpercent=0.01 -wpercent=0.01 -maxexpand=400 |'+
      'pamscale -pixels $[1024*768] > "%s"')
	"""
    self.det_pre_filter=(
      'djpeg "%s" | ppmtopgm | '+
      'pnmnorm -bpercent=0.01 -wpercent=0.01 |'+
      'pnmscale -pixels '+str(1024*768)+' > "%s"')
    self.det_pre_scale=('maxpix',1024*768)    

  def set_obsidian_compute_descriptors(self):
    self.det_prog='obsidian'    
    self.det_args="det-type=hessian aff-skip=1 desc-type=cs_lbp desc-scale_one=15.0 desc-T=0.0001 pyr-min_size=12 det-threshold=500.0"
    self.det_pre_filter=(
      'djpeg "%s" | ppmtopgm | '+
      'pnmnorm -bpercent=0.01 -wpercent=0.01 -maxexpand=400 |'+
      'pnmscale -pixels $[1024*768] > "%s"')
    self.det_pre_scale=('maxpix',1024*768)

  def examples(self):
    # interesting queries
    return []    

# ======================================================================================

def select_config(dbname):
  """
  returns a config object from a configuration name. The config class
  name must end in Config (which needs not be specified in the
  dbname). It can be specified as:
  - XXXX or XXXXConfig which must be defined in this file
  - YYYY:XXXX loads XXXXConfig from imported module YYYY. To
  subclass the Config class in module YYYY, import config
  - YYYY.py:XXXX is expanded as XXXXConfig class defined in file
  YYYY.py. Nothing special is needed to access the Config class: it is
  passed with the module name
  XXXX can contain arguments to the class constructor:
  XXXX(12,13) will build class XXXXConfig(12,13)
  """
  
  if not dbname:
    sys.stderr.write("weird dbname %s\n"%dbname)
    sys.exit(1)

  if "(" in dbname and dbname[-1]==")":
    ba=dbname.index("(")
    dbsuff=dbname[ba:]
    dbargs=[eval(a) for a in dbname[ba+1:-1].split(',')]
    dbname=dbname[:ba]
  else:
    dbargs=()
    dbsuff=""
    
  if not dbname.endswith('Config'):
    dbname+="Config"
      
  if ':' not in dbname:
    try:
      # build config class and convert to object (=constructor)
      configclass=globals()[dbname]
    except KeyError:
      sys.stderr.write("unknown config %s\n"%dbname)
      sys.exit(1)
  else:
    colon=dbname.index(':')
    filename=dbname[:colon]
    classname=dbname[colon+1:]
    if '.' not in filename:
      # load from this directory
      locs={}
      mod=__import__(filename,globals(),locs)
      configclass=getattr(mod,classname)
    else:
      # load from arbitrary file
      globs=globals()
      execfile(filename,globs,globs)
      configclass=globs[classname]
  c=configclass(*dbargs)
  c.dbname=c.__class__.__name__[:-6]+dbsuff

  return c
  
# ======================================================================================

class Yari20kConfig(OverallConfig):  
  def __init__(self, user_id):
    self.user_id = user_id
    OverallConfig.__init__(self)      
    self.invfilename=self.basedir+'%s/ivf/invfile_k20k_binsign.ivfgeo'%self.user_id
    self.fetch_base()

  def fetch_base(self):
    project_dirs = os.listdir(self.basedir+"%s/"%self.user_id)
    project_dirs.remove('ivf')
    project_dirs.remove('mobile')    
    self.fnames=[]
    for project_dir in project_dirs:
        self.fnames.extend([fn[:-4] for fn in os.listdir(self.basedir+"%s/%s/images/"%(self.user_id, project_dir)) if (fn.endswith(".jpg") or fn.endswith(".jpeg") or fn.endswith(".png"))])
    self.fnames.sort()
    self.nimg=len(self.fnames)
    
  """def img_filename(self,n):
    return "%s/images/%s.jpg"%(self.basedir,self.fnames[n])"""
    
  """          
  def thumb_filename(self,n):
    return "%s/thumbnails/%s_thumb.jpg"%(self.basedir,self.fnames[n])

  def desc_filename(self,n):
    return "%s/siftgeo/hesaff_norm/%s.siftgeo"%(self.basedir,self.fnames[n])    
    
  def vw_filename(self,n):
    return "%s/vwsgeo/k20k/%s.vwsgeo"%(self.basedir,self.fnames[n])
  """      
