import sys,time,inspect,os,errno,types

class TicToc(object):
    """>>> with TicToc('test'):
        ...     time.sleep(2)
        Test Elapsed time is 2.000073 seconds.
    """
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %ims' % int(1000*(time.time() - self.tstart))


class State:
  """ state of a sequence of operations. Default implementation just displays
  stage on stdout"""

  def next(self,name,**kwargs):
    " step to a stage called name "
    print " "*70+"\rstage "+name,
    for k,v in kwargs.iteritems():
      print "%s=%s"%(k,v),
    print

  def set_frac(self,frac):
    " set fraction completed for current stage "
    n=int(frac*70)
    sys.stdout.write("+"*n+"-"*(70-n)+"\r")
    sys.stdout.flush()

class DummyState:
  def next(self,name,**kwargs):
    pass
  set_frac=None


def try_mkdir(dirname):
  """ returns False when the directory already existed """
  try:
    os.mkdir(dirname)
  except OSError,e:
    if e.errno==errno.EEXIST:
      return False
    else:
      raise e
  return True



def prepare_dir(fname):
  "make sure the enclosing directory of this file exists"
  dirname=fname[:fname.rfind('/')]
  if os.access(dirname,os.W_OK):
    return
  try:
    os.makedirs(dirname)
  except OSError,e:
    if e.errno!=errno.EEXIST:
      raise



def copy_file(src,dst):
  # print "copy %s -> %s"%(src,dst)
  open(dst,"w").write(open(src,"r").read())



def gen_thumbnail(infile,outfile):
  """ uses ImageMagick """
  os.system('convert -strip -quality 50  -geometry "150x150>" "%s" "%s"'%
            (infile,outfile))

def gen_thumbnail_string(infile,maxsz):
  """ uses ImageMagick """
  f=os.popen('convert -strip -quality 50  -geometry "%dx%d>" "%s" "-"'%
             (maxsz,maxsz,infile),"r")
  return f.read()



def invert_aff(a):
  """ invert the affine transform [a[0] a[1] a[2] ; a[3] a[4] a[5] ] """
  (a0,a1,a2,a3,a4,a5)=a
  det=a0 * a4 - a3 * a1
  return [
    a4 / det,
    -a1 / det,
   -(-a1 * a5 + a2 * a4) / det,
    -a3 / det,
    a0 / det,
    -(a0 * a5 - a2 * a3) / det]

  
def aff_pre_post_scale(a,pres,posts):
  """ pre- and post- scale an affine transfrom """
  p=pres*posts
  return [
    a[0]*p, a[1]*p, a[2]*posts,
    a[3]*p, a[4]*p, a[5]*posts]
    
  
def parse_as_type(ty,sval):
  """ interpret string sval as the same type as vty """
  if ty==types.BooleanType:
    if sval.lower() in ("0","false"): return False
    if sval.lower() in ("1","true"): return True
    raise ValueError("cannot interpret %s as boolean"%sval)
  else:
    return ty(sval)


class Object:
  """ only to add fields """
  pass

