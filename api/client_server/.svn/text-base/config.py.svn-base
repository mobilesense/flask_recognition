"""
all the datafiles used by the image database.

To make a new database:

- add a new XXXXConfig from an example (SubNisterConfig). In the beginning
you only need the number of images (nfixed), their filenames (img_filename),
and a cluster file. All other files do not need to exist.

- run

> python prepare_imbase.py XXXX


"""

import sys

class Config:
  """
  references to all datafiles used by the database.
  For each image there is (standard extension in brackets):
  - the image file (img)
  - the descriptor in Miko's -o4 format (bsifgeo)
  - the visual word vector (vw)
  - the thumbnail (thumb.jpg)
  Images fall in 3 classes:
  - fixed initial: they were in the initial database
  - fixed extended: they were added in previous sessions
  (with an ExtendedXXXConfig class)
  - allocated: added in the current session. They may not be indexed yet.
  The fixed images have numbers 0..n_fixed-1
  Allocated ones hase n_fixed..n_fixed+n_alloc
  """

  def __init__(self):
    self.dbname=self.__class__.__name__[:-6]
    # store vw in sparse format
    self.sparse_vw=False
    # cluster file formats:
    # 'float': loaded with load_clusters
    # 'imat':  loaded with imat_new_fread
    # 'ivecs': loaded with load_clusters_ivecs
    self.cluster_format='float'
    # visual word allocation methods:
    # 'simple': simple exhaustive search
    # 'ann':    approximate nearest neighbour
    self.nn_method='simple'
    self.extensions=[]
    self.add_dir=None
    self.nalloc=0
    self.cdminvfilename=''
    self.cdmfactorsname=''
    self.ncdm=0
    # good for CS LBP
    self.geom_thr=0.3
    self.geom_fthr=0.15

  def added_file(self,n,ext):
    """ check if n refers to an allocated or fixed extended file"""
    if n>=self.nfixed+self.nalloc:
      raise ValueError("file %s %d does not exist"%(ext,n))
    if n>=self.nfixed:
      return '%s/%09d.%s'%(self.add_dir,n,ext)
    for (dirname,n0,n1) in self.extensions:
      if n0<=n<n1:
        return "%s/data/%09d.%s"%(dirname,n,ext)
    return None
    
  def extend(self,dirname,n0,n1):
    # add at beginning
    self.extensions=[(dirname,n0,n1)]+self.extensions
    self.nfixed=self.nimg=n1
    self.invfilename=dirname+"/extended.ivf"
    self.cdminvfilename=dirname+"/extended.cdmsup.ivf"
    self.cdmfactorsname=dirname+"/extended.cdmfac.hash"
    
  def vw_filename(self,n):
    fn=self.added_file(n,'vw')
    if fn: return fn
    return self.vwpattern%n

  def desc_filename(self,n):
    fn=self.added_file(n,'bsifgeo')
    if fn: return fn
    return self.descpattern%n

  def img_filename(self,n):
    fn=self.added_file(n,'img')
    if fn: return fn
    return self.imgpattern%n
  
  def thumb_filename(self,n):
    return "/home/cornwall/douze/img/generic_thumb.png"

    

class FotoliaConfig(Config):

  def __init__(self):
    Config.__init__(self)
    self.nimg=100000    
    self.nfixed=self.nimg
    self.imgpattern="/scratch2/curan/bourez/fotolia/originals/3%09d"
    self.thumbpattern="/scratch2/curan/bourez/fotolia/thumbnails/3%09d.jpg"
    self.det_prog='compute_descriptors'
    # self.clusterfile='/scratch2/cepheus/hjegou/clust/10000_chrystian.clust'
    # self.clusterfile='/scratch/cepheus/hjegou/nister_test/clust/10000_chrystian.clust'
    self.cluster_format='imat'
    self.clusterfile='/home/cornwall/douze/scratch/saved_clusters.imat'
    self.descpattern="/scratch2/curan/bourez/fotolia/bsifgeo/3%09d.bsifgeo"
    self.vwpattern='/scratch2/curan/bourez/fotolia/vw/3%09d.vw'
    self.q_nbcells=200000
    self.pcafile="/scratch2/curan/douze/fotolia.pca"
    self.invfilename="/scratch2/curan/douze/fotolia.ivf"
    self.det_thresh=100
    self.det_max_sift=1000
    self.add_dir='/scratch2/cepheus/douze/added_dir'

    self.geom_thr=0.2
    self.geom_fthr=0.12

  def thumb_filename(self,n):
    fn=self.added_file(n,'thumb.jpg')
    if fn: return fn
    return self.thumbpattern%n


class Nister2Config(Config):
  
  def __init__(self):
    Config.__init__(self)
    self.nimg=10200
    self.nfixed=self.nimg
    self.imgpattern='/home/clear/hjegou/dataset/nister/jpg/ukbench%05d.jpg'
    self.det_prog='compute_descriptors'
    self.cluster_format='ivecs'
    self.clusterfile='/home/clear/hjegou/dataset/nister/clust/nister_sample_thres200.clust'
    
    self.descpattern="/home/clear/douze/db/nister/bsiftgeo/ukbench%05d.bsiftgeo"
    self.vwpattern='/home/clear/douze/db/nister/vw/ukbench%05d.vw'
    self.sparse_vw=True
    self.nn_method='simple'
    self.q_nbcells=200000
    self.pcafile="/scratch2/curan/douze/nister_10200.pca"
    self.invfilename="/home/clear/douze/db/nister/nister_10200.ivf"
    # parameters to Miko's compute_descriptors
    self.det_thresh=100
    self.det_max_sift=1000
    self.thumbpattern="/scratch2/curan/douze/nister_thumbs/ukbench%05d.jpg"
    self.add_dir='/scratch2/cepheus/douze/added_dir'

  def thumb_filename(self,n):
    fn=self.added_file(n,'thumb.jpg')
    if fn: return fn
    return self.thumbpattern%n


class SubNisterConfig(Config):
  
  def __init__(self):
    Config.__init__(self)
    self.nimg=200
    self.nfixed=self.nimg
    self.imgpattern='/scratch2/cepheus/hjegou/dataset/nister/jpg/ukbench%05d.jpg'    
    self.det_prog='obsidian'
    self.det_obsidian_args=""
    
    self.cluster_format='ivecs'
    self.clusterfile='/scratch2/curan/douze/subnister/clusters.ivecs'
    
    self.descpattern="/scratch2/curan/douze/subnister/bsiftgeo/%05d.bsifgeo"
    self.thumbpattern="/scratch2/curan/douze/subnister/thumbs/%05d.jpg"

    self.vwpattern='/scratch2/curan/douze/subnister/vw/%05d.vw'
    self.q_nbcells=200000
    self.pcafile="/scratch2/curan/douze/subnister/subnister.pca"
    self.invfilename="/scratch2/curan/douze/subnister/subnister.ivf"
    
    self.add_dir='/scratch2/curan/douze/subnister/added_dir'

  def thumb_filename(self,n):
    fn=self.added_file(n,'thumb.jpg')
    if fn: return fn
    return self.thumbpattern%n

class NisterCDMConfig(Config):
  
  def __init__(self):
    Config.__init__(self)
    self.nimg=1000
    self.nfixed=self.nimg
    self.topdir='/scratch/adonis/mordelet/lear/benoit/cdm/run/db/nistercdm'
    self.imgpattern=self.topdir+'/jpg/cdm%05d.jpg'
    self.det_prog='compute_descriptors'
    self.cluster_format='ivecs'
    self.clusterfile=self.topdir+'/clust/nistercdm_sample_thres200.clust'
    
    self.descpattern=self.topdir+"/bsiftgeo/cdm%05d.bsiftgeo"
    self.vwpattern=self.topdir+'/vw/cdm%05d.vw'
    self.sparse_vw=True
    self.nn_method='simple'
    self.q_nbcells=20000
    self.pcafile=self.topdir+("/nistercdm_%d.pca"%self.nimg)
    self.invfilename=self.topdir+("/nistercdm_%d.ivf"%self.nimg)
    self.det_thresh=100
    self.det_max_sift=1000
    self.thumbpattern=self.topdir+"/nister_thumbs/cdm%05d.jpg"
    self.add_dir='/scratch/adonis/mordelet/lear/benoit/cdm/run/data/added_dir'

  def thumb_filename(self,n):
    fn=self.added_file(n,'thumb.jpg')
    if fn: return fn
    return self.thumbpattern%n


class DivPhotosConfig(Config):
  
  def __init__(self):
    Config.__init__(self)
    self.nimg=5216
    self.nfixed=self.nimg
    self.ncdm=0
    self.imgpattern="/home/clear/douze/db/div_photos/jpg/%05d.jpg"
    self.det_prog='compute_descriptors'
    self.cluster_format='ivecs'
    self.clusterfile='/home/clear/hjegou/dataset/nister/clust/nister_sample_thres200.clust'    
    self.descpattern="/home/clear/douze/db/div_photos/siftgeo/%05d.siftgeo"
    self.vwpattern='/home/clear/douze/db/div_photos/vw/%05d.vw'
    self.sparse_vw=True
    self.nn_method='simple'
    self.q_nbcells=200000
    self.pcafile="/scratch2/curan/douze/nister_10200.pca"
    self.invfilename="/home/clear/douze/db/div_photos/invfile.ivf"
    self.cdminvfilename=''
    self.cdmfactorsname=''
    # parameters to Miko's compute_descriptors
    self.det_thresh=100
    self.det_max_sift=1000
    self.thumbpattern="/home/clear/douze/db/div_photos/thumb/%05d.jpg"

  def thumb_filename(self,n):
    fn=self.added_file(n,'thumb.jpg')
    if fn: return fn
    return self.thumbpattern%n


def select_config(dbname):  
  try:
    # build config class and convert to object (=constructor)
    configclass=globals()[dbname+'Config']
  except KeyError:
    sys.stderr.write("unknown config %s\n"%dbname)
    sys.exit(1)
  return configclass()

