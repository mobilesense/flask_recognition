
from invfile import *

from config import select_config

from utils import *

from imbase import *

import array,pdb,os,sys


def load_ivec_to_sparse(vwf):
  v=ivec_new_fread(open(vwf,"r"))  
  (svi,sv)=ivec_to_sparse(v); ivec_delete(v)
  return (svi,sv)

def all_gen_thumb(config,im_offset=0,im_step=1):
  print "generating thumbnails"
  
  for label in range(im_offset,config.nimg,im_step):
    imname=config.img_filename(label)
    thumbname=config.thumb_filename(label)
    print "%d %s -> %s"%(label,imname,thumbname)
    gen_thumbnail(imname,thumbname)
  
def all_compute_descs(config,im_offset=0,im_step=1):
  print "computing descriptors"
  fc=CachedData(config).get_FeaturesComputer()
  for label in range(im_offset,config.nimg,im_step):
    imname=config.img_filename(label)
    descname=config.desc_filename(label)
    print "%d %s -> %s"%(label,imname,descname)
    img=fc.load_image(imname)
    (n_interest,n_desc,bsiftgeo)=fc.compute_descriptors(img)
    print "# pts %d %d"%(n_interest,n_desc)
    fc.store_descriptors(bsiftgeo)

  
def all_assign_vw(config,im_offset=0,im_step=1):
  print "assigning visual words"
  vwa=CachedData(config).get_VWAssignement()
  for label in range(im_offset,config.nimg,im_step):
    descname=config.desc_filename(label)
    vwname=config.vw_filename(label)
    print "%d %s -> %s"%(label,descname,vwname)
    descs=vwa.load_descs(descname)
    vw=vwa.compute_vw(descs)
    vwa.store_vw(vw,vwname)
 


def make_invfile(config):

  nbvw=10000
  hash_nbcells=100000

  ivf=invfile_new(nbvw)

  print "nbvw=%d nbvec=%d"%(ivf.nbvw,ivf.nbvec)

  missing_labels=[]

  for label in range(config.nimg):
    
    vwf=config.vw_filename(label)
    
    if label%1000==0:
      print "adding",vwf

    try:
      if not config.sparse_vw:
        v=ivec_new_fread(open(vwf,"r"))
        (svi,sv)=ivec_to_sparse(v); ivec_delete(v)
      else:
        f=open(vwf,"r")
        svi=ivec_new_fread(f)
        sv=ivec_new_fread(f)
    except IOError,e:
      if e.errno==2: # No such file or directory
        missing_labels.append(label)
        continue
      else:
        raise e

    invfile_add(ivf,label,svi,sv)

    ivec_delete(svi)
    ivec_delete(sv)

  print "missing vw files:",missing_labels

  print "writing invfile "

  invfile_write(open(config.invfilename,"w"),ivf)


def make_pcafile(config):
  pcafile=config.pcafile
  imlist="/tmp/prepare_imbase_imlist.dat"
  
  f=open(imlist,"w")

  missing_descs=[]

  if config.nimg<200:
    step=1
  else:
    step=config.nimg/200
  
  # use only 1 out of step descriptor files (too slow else)
  for i in range(0,config.nimg,step):
    fn=config.desc_filename(i)
    try:
      open(fn,'r')
    except IOError,e:     
      if e.errno==2: # No such file or directory
        missing_descs.append(i)
        continue
      else:
        raise e
      
    f.write(fn+"\n")

  del f

  print "missing desc files: ",missing_descs
    
  ret=os.system("../geom_filter/geom_filter_exec -verbose 4 -pcafile %s -desctype bsiftgeo -computepca @%s"%(
    pcafile,imlist))

  if ret!=0:
    raise RuntimeError("geom_filter_exec crashed ret=%d"%ret)

def make_clusters(config):
  clusterfile=config.clusterfile
  datapts="/tmp/prepare_imbase.siftgeo"
  
  f=open(datapts,"w")
  missing_vw=[]

  for i in range(0,config.nimg):
    fn=config.desc_filename(i)
    try:
      f2=open(fn,"r")
    except IOError,e:     
      if e.errno==2: # No such file or directory
        missing_vw.append(i)
        continue
      else:
        raise e
      
    f.write(f2.read())

  f.close()
    
  print "missing vw files: ",missing_vw

  nb_centroids=10000

  cmd="../siftgeo/cluster k=%d -i=%s -o=%s infmt=2 outfmt=0  max_iter=20"%(
    nb_centroids,datapts,clusterfile)

  print "execing",cmd
  sys.stdout.flush()
  
  ret=os.system(cmd)

  if ret!=0:
    raise RuntimeError("cluster crashed ret=%d"%ret)




def usage():
  sys.stderr.write("usage: %s [-db dbname] [-ofs offset] "%sys.argv[0]+
                   "[-step step] [desc] [vw] [thumb] [ivf] [pca] [clusters]\n")
  sys.exit(1)


if __name__=='__main__':

  args=sys.argv[1:]

  todo=[]
  dbname="none"
  im_offset=0
  im_step=1
  while args:
    if args[0] in ['-h','--help']:
      usage()
    elif args[0]=='-db':
      del args[0]
      dbname=args[0]
    elif args[0]=='-ofs':
      del args[0]
      im_offset=int(args[0])
    elif args[0]=='-step':
      del args[0]
      im_step=int(args[0])
    elif args[0] in ['thumb','desc','vw','ivf','pca','clusters']:
      todo.append(args[0])
    else:
      sys.stderr.write("unknown arg %s\n"%args[0])
      usage()
    del args[0]


  print "doing %s on db %s"%(todo,dbname)
  
  config=select_config(dbname)
  
  if 'thumb' in todo:
    all_gen_thumb(config,im_offset,im_step)
    
  if 'desc' in todo:
    all_compute_descs(config,im_offset,im_step)
    
  if 'vw' in todo:
    all_assign_vw(config,im_offset,im_step)
    
  if 'ivf' in todo:
    print "making inverted file"
    make_invfile(config)

  if 'pca' in todo:
    print "computing PCA"
    make_pcafile(config)
           
  if 'clusters' in todo:
    print "making clusters"
    make_clusters(config)
    
