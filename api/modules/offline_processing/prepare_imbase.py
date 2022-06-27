
import array,pdb,os,sys,random,thread,errno,traceback

from config.config import select_config

from basics import *

from imbase import *

from siftgeo.siftgeo import *

from utils.thread_utils import RunOnSet,parse_nt,clr_eol,parallel_map



def gen_thumb(config,label):
  imname=config.img_filename(label)
  thumbname=config.thumb_filename(label)
  prepare_dir(thumbname)
  print "%d %s -> %s"%(label,imname,thumbname)
  gen_thumbnail(imname,thumbname)

  
def compute_descs(config,fc,label):  
  imname=config.img_filename(label)
  descname=config.desc_filename(label)
  print "%d %s -> %s"%(label,imname,descname)
  prepare_dir(descname)
  img=fc.load_image(imname)
  (n_interest,n_desc,bsiftgeo)=fc.compute_descriptors(img)
  print "# pts %d %d"%(n_interest,n_desc)
  fc.store_descriptors(bsiftgeo,descname)
  # no need to dealloc
  
def assign_vw(config,vwa,label):
  descname=config.desc_filename(label)
  vwname=config.vw_filename(label)

  # *** should be added in compute_descs_and_vw if it works
  if config.ldc:
    ldc_filename = config.ldc_filename (label)
  else:
    ldc_filename = None

  prepare_dir(vwname)
  print "%d %s -> %s"%(label,descname,vwname)
  descs=vwa.load_descs(descname, ldc_filename)
  vw=vwa.compute_vw(descs)
  vwa.store_vw(vw,vwname)

def compute_descs_and_vw(config,fc,vwa,label,store_desc=False):
  imname=config.img_filename(label)
  vwname=config.vw_filename(label)
  descname=store_desc and config.desc_filename(label)  
  prepare_dir(vwname)
  print "%d %s -> %s %s"%(label,imname,vwname,descname)
  img=fc.load_image(imname)
  (n_interest,n_desc,descs)=fc.compute_descriptors(img)
  print "# pts %d %d"%(n_interest,n_desc)
  if descname:
    prepare_dir(descname)
    fc.store_descriptors(descs,descname)
  vw=vwa.compute_vw(descs)
  vwa.store_vw(vw,vwname)
  

def make_inv_slice(bb,begin,end,res,sliceno):
  n=len(res)
  b=begin+(end-begin)*sliceno/n
  e=begin+(end-begin)*(sliceno+1)/n
  ivf=res[sliceno]
  for label in range(b,e):
    bb.add_label_to(label,ivf)
    print "\r",label,config.display_name(label),clr_eol,
    sys.stdout.flush()

    if label%1000==900:
      tot_im=sum([bb.ivf_nim(i) for i in res if type(i)!=type(())])
      tot_pts=sum([bb.ivf_npt(i) for i in res if type(i)!=type(())])

      print "\ntotal %d images, %d pts"%(tot_im,tot_pts)
    
  res[sliceno]=(b,e,ivf)
  

def make_invfile(config,begin=None,end=None,n_thread=1):
  nbvw=config.nb_centroids

  bb=CachedData(config).get_BigBase(make_new=True,nbvw=nbvw)
  
  print "nbvw=%d, adding images %d to %d (to be stored in %s)"%(
    bb.ivf.nbvw,begin,end,config.invfilename)

  bb.begin_add()

  # if n_thread>1:
  if True: 
    if n_thread>end-begin: n_thread=end-begin

    if config.he_type=='pq':
      res=[ivfgeo.ivfgeo_pq_dup(bb.ivf) for i in range(n_thread)]                    
    else: 
      res=[ivfgeo.ivfgeo_new(nbvw,10) for i in range(n_thread)]
    RunOnSet(n_thread,range(n_thread),lambda i: make_inv_slice(bb,begin,end,res,i))
    for b,e,ivf in res:

      print "merge %d-%d"%(b,e)
      if config.he_type=='pq':
        ivfgeo.ivfgeo_pq_merge(bb.ivf,ivf)
      else: 
        ivfgeo.ivfgeo_merge(bb.add_d,ivf,0)
  else:
    for label in range(begin,end):
      vwf=config.vw_filename(label) 

      vw=bb.load_vw(vwf)

      bb.add_vw(vw,label)

      print "\r",label,vwf[-40:],
      sys.stdout.flush()


  # ivfgeo.ivfgeo_pq_display(bb.ivf) 

  bb.end_add()
  print "writing invfile ",config.invfilename
  prepare_dir(config.invfilename)
  bb.write_invfile(config.invfilename)

  
def do_count_points(config,i):
  if i%1000==0:
    print clr_eol,"%d\r"%i,
    sys.stdout.flush()
  return count_points(config.desc_filename(i),config.desc_format)

def make_desc_subset(config,npt_max,n_thread):
  """ get a subset of descriptors that fit in memory """
  npt_tot=0

  print "scanning %d files"%config.nimg
  nptperfile=parallel_map(n_thread,xrange(config.nimg),
                          lambda i: do_count_points(config,i))  
  npt_tot=sum(nptperfile)
  print "total %d pts"%npt_tot

  points,n=None,0
  
  if npt_tot<=npt_max:
    for i in range(0,config.nimg):
      (points,n)=read_points_add(config.desc_filename(i),points,n,config.desc_format)
    assert n==npt_tot
  else:    
    subset=random.sample(xrange(npt_tot), npt_max)
    subset.sort()
    print "selected subset of %d points"%len(subset)
    nr=0; j=0; prev_n=-1
    for i in range(config.nimg):
      fname=config.desc_filename(i)
      npt=count_points(fname,config.desc_format)
      nr_next=nr+npt;
           
      mask=yael.IntArray(npt)
      mask.clear(npt)
      j_prev=j;
      while j<npt_max and subset[j]<nr_next:
        mask[subset[j]-nr]=1
        j+=1
      nr=nr_next
      # print "taking %d pts from %s"%(j-j_prev,fname)
      if j==j_prev:
        continue
      (points,n)=read_points_add_with_mask(fname,points,n,config.desc_format,mask)
      if n/10000 != prev_n/10000:
        print clr_eol,"gethered %d pts\r"%n,
        sys.stdout.flush()        
      # print j,n,len([1 for k in range(npt) if mask[k]])
    assert n==npt_max
    
  return (points,n)



def make_clusters(config,only_he=False,n_thread=1):
  clusterfile=config.clusterfile

  print "Loading pts"
  
  nb_centroids=config.nb_centroids
  npt_max=1*1000*1000 #before it was 6*1000*1000
  
  (points,n)=make_desc_subset(config,npt_max,n_thread)
    
  (fpoints,d)=siftgeo_to_fvecs(points,n)
  delete_points(points,n);

  max_iter=40
  if not only_he:
    print "begin clustering/HE"
  else:
    print "begin HE"
  # Normalise vectors between iterations ? 
  nrlz=0
  nproj=64
  
  if not config.binsign_computer:
    centroids=yael.clustering_kmeans(n,d,fpoints,nb_centroids,max_iter,nrlz)
    prepare_dir(clusterfile)
    yael.fvecs_write(clusterfile,d,nb_centroids,centroids)    
  else:
    if not only_he:
      (centroids,clust_assign)=yael.clustering_kmeans_assign (n,d,fpoints,nb_centroids,max_iter, nrlz)      
      clust_assign=yael.IntArray.frompointer(clust_assign)
      clust_assign.this.acquire()
      prepare_dir(clusterfile)
      yael.fvecs_write(clusterfile,d,nb_centroids,centroids)      
    else:
      print "loading "+clusterfile
      (centroids,n2,d2)=yael.fvecs_new_read(clusterfile)
      assert n2==nb_centroids and d2==d
      print "compute assignement"
      clust_assign=yael.IntArray(n)
      yael.quantize_codebook_thread(n,nb_centroids,d,centroids,fpoints,clust_assign,n_thread,None)
    
    print "compute HE parameters"

    if config.he_type == 'nopca':
      sb = siftgeo_binarize_new (n, d, nb_centroids, nproj, fpoints,
                                 clust_assign)
    elif config.he_type in ('pca','pcanw'):
      sb = siftgeo_binarize_pca_new (n, d, nb_centroids, nproj,config.he_type=='pca' and 1 or 0,
                                     fpoints, clust_assign, centroids)
    else:
      assert False

    prepare_dir(config.binsign_computer)
    siftgeo_binarize_write (config.binsign_computer, sb)
    siftgeo_binarize_delete(sb)
  
  # dealloc
  yael.FloatArray.frompointer(centroids).this.acquire()
  yael.FloatArray.frompointer(fpoints).this.acquire()
  

   
def report_missing_1(name,fname,verb=False):
  print "%s:"%name,
  if verb: print fname,
  if not fname:
    print "undefined"
    return 2
  elif os.access(fname,os.R_OK):
    print "ok"
    return 0
  else:
    print "missing"
    return 1


def report_missing(name,(rmin,rmax),f,verb=False):
  if not verb: 
    print "Ranges of %s in [%d,%d]: "%(name,rmin,rmax),
  else:
    print "Checking %s in [%d,%d]: "%(name,rmin,rmax)
  sys.stdout.flush()
  last_missing=-1
  tot_missing=0  
  for r in range(rmin,rmax):
    fname=f(r)
    if verb: print r,fname,
    if not fname:
      print "file name for %d undefined"%r
      return 2
    if os.access(fname,os.R_OK):
      if verb: print "ok"
      if last_missing!=-1:
        if not verb: print "%d]"%r, ; sys.stdout.flush()
        tot_missing+=r-last_missing
        last_missing=-1
    else:
      if verb: print "missing"
      if last_missing==-1:
        last_missing=r
        if not verb: print "[%d,"%r, ;   sys.stdout.flush()
  if last_missing!=-1:
    if not verb: print "%d]"%rmax, ; sys.stdout.flush()
    tot_missing+=rmax-last_missing

  if tot_missing==0:
    print "ok"
    return 0
  elif tot_missing==rmax-rmin:
    print "all missing"
    return 1
  else:
    print "missing %d/%d files"%(tot_missing,rmax-rmin)
    return 4


def usage():
  sys.stderr.write("""usage: %s [-db dbname]
[-begin begin_index] [-end end_index] [-imnos list.dat] [-nt nthread]
[-check]
[desc] [vw] [thumb] [ivf] [cluster] [he] [desc+vw] [desc++vw]
[-ivfname]
"""%sys.argv[0])
  sys.stderr.write("""
given an image list, this computes:
   desc(*) the image descriptors 
   cluster: descriptor clusters + HE params file if needed (desc)
   he: HE param file (cluster,desc)
   vw(*): the visual words (desc,clusters)
   thumb(*): thumbnails for display
   ivf(*): the inverted file (vw)
   desc+vw(*): combined versions without storing intermediate results   
   desc++vw(*): combined versions also store descriptor   
(*): can be done in slices like [begin_index,end_index) or on selected imnos and decomposed in nthread threads
(in brackets): the dependencies
-ivfname: override inverted file name (to make in slices)
-check: don't compute, just check if outfiles exist
""")
  sys.exit(1)

def siftgeo_generator(config):
  for i in xrange(config.nimg):
    yield config.desc_filename(i)

def vwgeo_generator(config):
  for i in xrange(config.nimg):
    yield config.vw_filename(i)

if __name__=='__main__':

  args=sys.argv[1:]

  todo=[]
  dbname="none"
  im_begin=0
  im_end=-1
  n_thread=1
  ivfname=None
  docheck=False
  verbose=False
  imnos=None
  while args:
    a=args.pop(0)
    if a in ['-h','--help']:      usage()
    elif a=='-db':                dbname=args.pop(0)
    elif a=='-begin':             im_begin=int(args.pop(0))
    elif a=='-end':               im_end=int(args.pop(0))      
    elif a=='-imnos':             imnos=args.pop(0)
    elif a=='-nt':                n_thread=parse_nt(args.pop(0))
    elif a=='-ivfname':           ivfname=args.pop(0)
    elif a=='-check':             docheck=True
    elif a=='-v':                 verbose=True
    elif a in ['thumb','desc','vw','ivf','cluster','he',
               'desc+vw','desc++vw','img']:
      todo.append(a)
    else:
      sys.stderr.write("unknown arg %s\n"%a)
      usage()

  print "%s %s on db %s"%(docheck and "checking" or "doing",todo,dbname)
  
  config=select_config(dbname)

  if im_end==-1: im_end=config.nimg

  if im_end>config.nimg:
    print >> sys.stderr,"warn: forcing -end to",config.nimg
    im_end=config.nimg

  if docheck:
    err=0
    if not todo: todo=['img','thumb','desc','clusters','vw','ivf']
    if 'img' in todo:
      err|=report_missing('images',(im_begin,im_end),lambda im: config.img_filename(im),verbose)
    if 'thumb' in todo:
      err|=report_missing('thumb',(im_begin,im_end),lambda im: config.thumb_filename(im),verbose)
    if 'desc' in todo or 'desc++vw' in todo:
      err|=report_missing('desc',(im_begin,im_end),lambda im: config.desc_filename(im),verbose)
    if 'cluster' in todo:
      err|=report_missing_1('clusters',config.clusterfile,verbose)
      err|=report_missing_1('HE params',config.binsign_computer,verbose)
    if 'he' in todo:
      err|=report_missing_1('HE params',config.binsign_computer,verbose)
    if 'vw' in todo or 'desc+vw' in todo or 'desc++vw' in todo:
      err|=report_missing('vw',(im_begin,im_end),lambda im: config.vw_filename(im),verbose)
    if 'ivf' in todo:
      err|=report_missing_1('ivf',config.invfilename,verbose)
    sys.exit(err)

  if imnos:
    imnos=[int(l) for l in open(imnos,"r")]
  else:
    imnos=range(im_begin,im_end)
    
  
  if 'thumb' in todo:
    print "generating thumbnails"
    RunOnSet(n_thread,imnos,lambda im: gen_thumb(config,im))
    
  if 'desc' in todo:
    print "computing descriptors"
    fc=CachedData(config).get_FeaturesComputer()
    RunOnSet(n_thread,imnos,lambda im: compute_descs(config,fc,im))

  if 'cluster' in todo:
    print "making clusters"
    make_clusters(config,only_he=False,n_thread=n_thread)

  if 'he' in todo:
    print "making HE params"
    make_clusters(config,only_he=True,n_thread=n_thread)  

  if 'vw' in todo:
    print "assigning VWs"
    vwa=CachedData(config).get_VWAssignement()
    RunOnSet(n_thread,imnos,lambda im: assign_vw(config,vwa,im))
    
  if 'ivf' in todo:
    print "making inverted file"
    if ivfname: config.invfilename=ivfname
    make_invfile(config,im_begin,im_end,n_thread)
  
  if 'desc+vw' in todo or 'desc++vw' in todo:
    print "computing descriptors and assigning vw's"
    fc=CachedData(config).get_FeaturesComputer()
    vwa=CachedData(config).get_VWAssignement()    
    RunOnSet(n_thread,imnos,
             lambda im: compute_descs_and_vw(config,fc,vwa,im,'desc++vw' in todo))

    
