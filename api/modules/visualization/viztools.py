import os, sys, pdb, time, cPickle, math, numpy as np, thread

import scripts.draw_postscript as drawps
from siftgeo import siftgeo

try:
  import Image, ImageDraw, ImageFont
except ImportError:
  pass

from yael import yael
try:
  import obsidian
except ImportError:
  pass

try:
  from imbase.config import select_config, ImageSizes
except ImportError:
  from config import select_config, ImageSizes

import warnings

colors = [(128,255,128),(128,128,255),(128,255,255),(0,128,0),(0,0,128),(0,128,128),(128,128,128),(128,0,128),(128,0,255),(64,64,64)]

def name_to_config_no( config, name ):
  no = [no for no in range(config.nimg) if config.display_name(no)==name]
  assert no, pdb.set_trace()#"error: name '%s' is not in config '%s'" % (name, config.dbname)
  assert len(no)==1
  no = no[0]
  print "translating name '%s' into config no %d for config %s" % (name, no, config.dbname)
  return no

def config_no_to_query_no( qconfig, cno ):
  no = [i for i,qno in enumerate(make_queryset(qconfig)) if qno==cno]
  assert no, pdb.set_trace()#"error: name '%s' is not in config '%s'" % (name, c.dbname)
  assert len(no)==1
  no = no[0]
  print "translating config no '%s' into query no %d for config %s" % (cno, no, qconfig.dbname)
  return no


def get_image_path( config, imno ):
  """retrieve an image from the base
     return (file_name, is_temporary)
  """
  path = config.prescaled_img_filename(imno)
  if path and os.path.isfile(path):
    return path, False
  else:
    if config.det_pre_scale and config.det_pre_scale[0]!='upscale':
      with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=RuntimeWarning)
        path = os.tmpnam()+".pgm" # assume unicity
      os.system(config.det_pre_filter%(config.img_filename(imno),path))
      return path, True
    else:
      path = config.img_filename(imno)
      return path, False

def get_image( config, imno ):
  fname, istemp = get_image_path(config,imno)
  img = Image.open(fname)
  if istemp:  os.remove(fname)
  return img

# print some specified text in an image
def print_text_in_img(img,txt,size=20,color=255,bbox=None,bbox_color=128):
  draw = ImageDraw.Draw(img)
  f= ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans.ttf", size)
  draw.text( (1,img.size[1]-size-4), txt, fill=0,     font=f)
  draw.text( (0,img.size[1]-size-5), txt, fill=color, font=f)
  if bbox: draw.rectangle( bbox, outline=bbox_color, fill=None )
  del draw

def get_homography( coefs ):
  homography = np.matrix(coefs+[1])
  homography.shape = (3,3)
  return homography

def apply_homo_pt( pt, mat ):
  res = mat * np.matrix(list(pt)+[1]).T
  res = np.array(res)
  res /= res[2,0]
  return [res[0,0],res[1,0]]

def apply_homo_pts( pts, mat ):
  res = []
  for pt in pts:
    res.append( apply_homo_pt(pt,mat) )
  return res

def bbox_to_corners( bbox ):
  return [(bbox[0],bbox[1]),(bbox[2],bbox[1]),(bbox[2],bbox[3]),(bbox[0],bbox[3])]

def bbox_center( bbox ):
  return (sum(bbox[::2])/2.0,sum(bbox[1::2])/2.0)

def bbox_width( bbox ):
  return abs(bbox[2]-bbox[0])

def bbox_height( bbox ):
  return abs(bbox[3]-bbox[1])

def concatenate_two_images_raw( qim, dbim, qtxt='', dbtxt='' ):
  print_text_in_img(qim,qtxt)
  print_text_in_img(dbim,dbtxt)
  qw = qim.size[0]
  width = dbim.size[0] + qw
  height = max( qim.size[1], dbim.size[1] )
  img = Image.new( "RGBA", (width, height) )
  img.paste(qim, (0, 0))
  img.paste(dbim, (qw, 0))
  return img

def concatenate_two_images( qconfig, qimno, bconfig, bimno, imms=None, ptmatches=None, selected=None, qbbox=None, dbbox=None, aff=None, 
                            show_null_matches=0, show_pt_scale=False, qpatch=None, bpatch=None, eps=None ):
  qim = get_image(qconfig,qimno)
  qw = qim.size[0]
  print_text_in_img(qim,"query '%s'" % qconfig.display_name(qimno))
  dbim = get_image( bconfig, bimno )
  print_text_in_img(dbim,"base '%s'" % bconfig.display_name(bimno))
  
  assert not(imms and (ptmatches or qbbox or dbbox or aff)), "error: imms can not be defined simultaneously with ptmatches"
  if imms==None:
    class Dummy: pass
    imm = Dummy()
    imm.ptmatches = ptmatches or []
    imm.db_bbox = dbbox
    imm.q_bbox = qbbox
    imm.aff = aff
    imm.selected=selected or []
    imms=[imm]
  
  width = dbim.size[0] + qw
  ph = max([0,qpatch and qpatch.height, bpatch and bpatch.height])
  height = max( qim.size[1], dbim.size[1] ) + ph
  img = Image.new( "RGBA", (width, height) )
  img.paste(qim, (0, 0))
  img.paste(dbim, (qw, 0))
  
  # EPS out
  if eps:
    (xs,ys) = img.size
    ps.ps_header(xs+20+100,ys+20,eps)
    eps.write("10 10 translate 1 0 0 setrgbcolor\n")
    ps.embed_jpeg(qconfig.img_filename(qimno),(0,0)+qim.size,eps)
    ps.embed_jpeg(bconfig.img_filename(bimno),(qw+100,0)+dbim.size,eps)
  
  # paste patches
  def paste_patch( patch, pos ):
    if patch:
      if type(patch)!=Image: patch=obsimg_to_Image(patch)
      img.paste(patch, pos )
  paste_patch(qpatch,(0,qim.size[1]))
  paste_patch(bpatch,(qw,dbim.size[1]))
  
  two_images_print_bboxes( img, qw, imms, eps )
  
  return two_images_print_ptmatches( img, qw, imms, show_null_matches, show_pt_scale, eps )

def print_bbox(draw,bbox,color):
  draw.rectangle( bbox, fill=None, outline=color )
  draw.line( bbox_to_corners(bbox)[2:4], fill=0 )

# complemetary with concatenate_two_images
def two_images_print_bboxes( img, qw, imms, eps ):
  draw = ImageDraw.Draw(img)
  for n,imm in enumerate(imms):
    aff = imm.aff
    qbbox = imm.q_bbox
    dbbox = imm.db_bbox
    
    if aff and qbbox:
      # draw query bounding box
      print_bbox(draw,qbbox,colors[n])
      
      # draw retrieved image bounding box
      trf = get_homography(aff)
      trf[0,2] += qw
      db_corners = apply_homo_pts( bbox_to_corners(qbbox), trf )
      for i in range(4):
        draw.line( db_corners[i-1]+db_corners[i], fill=colors[n] )
      draw.line( db_corners[2]+db_corners[3], fill=0 )
    
    if eps: 
      if qbbox: ps.ps_draw_rectangle(qbbox,eps,linewidth=10,dash=(50,50),rgb=(0,1,0))
      if dbbox: ps.ps_draw_rectangle((qw+100+dbbox[0],dbbox[1],qw+100+dbbox[2],dbbox[3]),eps,linewidth=10,dash=(50,50),rgb=(0,1,0))

# complemetary with concatenate_two_images
def two_images_print_ptmatches( img, qw, imms, show_null_matches, show_pt_scale, eps ):
  if not imms:  return img
  if not hasattr(imms[0],'ptmatches'): return img
  draw = ImageDraw.Draw(img)
  
  # non_selected over all groups
  non_selected = set(range(len(imms and imms[0].ptmatches or [])))-set([s for imm in imms for s in imm.selected])
  
  for n,imm in enumerate(imms):
    color = colors[n]
    
    if eps: ps.ps_prepare_draw_line(eps)
    
    if show_null_matches:
      old_img = img.copy()
      for m in non_selected:
        m=imm.ptmatches[m]
        draw.line( m[0][:2]+(qw+m[1][0],m[1][1]), fill=(192,0,0) )
        if show_pt_scale:
          s=show_pt_scale
          offx=0
          for i in range(2):
            draw.ellipse( (offx+m[i][0]-s*m[i][2],m[i][1]-s*m[i][2],offx+m[i][0]+s*m[i][2],m[i][1]+s*m[i][2]), outline=(128,0,0))
            if i==0:
              ang=m[i][3]+ (m[3] if i==0 else 0)
              draw.line((offx+m[i][0],m[i][1],offx+m[i][0]+s*m[i][2]*np.cos(ang),m[i][1]+s*m[i][2]*np.sin(ang)), fill=(128,0,0))
            offx=qw
      img = Image.blend( old_img, img, show_null_matches )
      draw = ImageDraw.Draw(img)
      show_null_matches = False # show null matches only once
    
    for m in imm.selected:
      m=imm.ptmatches[m]
      p1=m[0][:2]
      p2=(qw+m[1][0],m[1][1])
      draw.line( p1+p2, fill=color )
      if show_pt_scale:
        s=show_pt_scale
        offx=0
        for i in range(2):
          draw.ellipse( (offx+m[i][0]-s*m[i][2],m[i][1]-s*m[i][2],offx+m[i][0]+s*m[i][2],m[i][1]+s*m[i][2]), outline=color)
          ang=m[i][3]+ (m[3] if i==0 else 0)
          draw.line((offx+m[i][0],m[i][1],offx+m[i][0]+s*m[i][2]*np.cos(ang),m[i][1]+s*m[i][2]*np.sin(ang)), fill=color)
          offx=qw
        if eps: ps.ps_draw_line(p1,(p2[0]+100,p2[1]),eps,rgb=(1,1,0),linewidth=2)
  
  return img

tcos = [np.cos(ang/100.0) for ang in range(628)]
tsin = [np.sin(ang/100.0) for ang in range(628)]

def get_keypoint_ellipse(k):   # k is a geom struct
  (a,b,v1,v2)=drawps.mat22_svd(k.mi11,k.mi12,
                               k.mi21,k.mi22)
  det = np.sqrt(k.mi11*k.mi22-k.mi12*k.mi21)    # should be equal to one, otherwise we have to multiply the scale by it
  assert 0.999<det<1.001, 'error: weird determinant in ellipse definition'
  # (a,b) are the major and minor axis length
  # (v1,v2) is the major axis direction of the ellipse
  return (det*k.scale,a,b,v1,v2)

def draw_kpt_in_img(draw,d,color,verbose=0):
  assert (draw.mode=='RGB' or draw.mode=='RGBA') == (type(color)==tuple), "error in draw_kpt_in_img: image mode '%s' incompatible with color '%s'" % (draw.mode,str(color))
  k = d.geom
  s,a,b,v1,v2 = get_keypoint_ellipse(k)
  if verbose:
    print "(%g, %g) sc=%g ori=%g (a,b)=(%g %g) crn=%g dim/vw=%d" % (k.x,k.y,s,k.angle,a,b,k.cornerness,d.vw)
  if verbose>1:
    desc = [vr.get_desc_bin(i,j) for j in range(d.dim)]
    print "desc=%s (min=%d, max=%d)" % (str(desc), min(desc[2:]), max(desc[2:]))
  pts = [(k.x + s*(v1*a*tcos[ang]-v2*b*tsin[ang]), 
          k.y + s*(v1*b*tsin[ang]+v2*a*tcos[ang])) for ang in xrange(0,628,min(100,int(1+500/k.scale)))]   # one point every 5 pixel => npoint = 2*pi*r/5 => every 628/npoint
  draw.polygon( pts, outline=color )
  
  #draw.ellipse( (k.x-s,k.y-s,k.x+s,k.y+s), outline=color )
  #a *= k.scale
  #b *= k.scale
  #e,f = a*v1,a*v2;  draw.line( (k.x,k.y,k.x+e,k.y+f), fill=color )
  #e,f =-b*v2,b*v1;  draw.line( (k.x,k.y,k.x+e,k.y+f), fill=color )

# print keypoint position and scale in the image
def print_pointset_in_img(img,vr,group,bbox=None,recon=None,verbose=0,color=(0,255,0)):
  
  if img.mode=='L':
    img = img.convert('RGB')
  
  draw = ImageDraw.Draw(img)
  
  patch = obsidian.image_new(128,128)
  def show(self): obsimg_to_Image(self).show()
  setattr(patch.__class__, 'show', show)
  
  for i in group:
    try:
      draw_kpt_in_img( draw, vr[i], color, verbose )
    except IndexError:
      raise IndexError("pointset[0:%d] does not contain point #%d"%(vr.n,i))
    
    if recon:
      obsidian.reconstruct_zernike_cpu( patch, vr[i].desc, vr[i].dim, recon.det_args.find('desc-zer_polar=0')<0, vr.fmt )
      patch.show()
      if verbose: raw_input()
  
  if bbox:  draw.rectangle(bbox)
  return img

def show_groups( qconfig, qimno, qmatch, qconfig2=None, colors=None ):
  
  vw_to_group, ngroups, weights = qconfig.get_vw_to_group(qimno, qconfig2)
  
  img = get_image(qconfig,qimno)
  if img.mode=='L': img = img.convert('RGB')
  draw = ImageDraw.Draw(img)
  
  qdescs = siftgeo.pointset_read(qconfig.desc_filename(qimno),qconfig.desc_format)
  qdescs2 = qconfig2 and siftgeo.pointset_read(qconfig2.desc_filename(qimno),qconfig2.desc_format)
  
  for a in np.matrix(qmatch).A1.nonzero()[0]:
    best = int(weights[:,a].argmax())
    assert best<qdescs.n or qdescs2, pdb.set_trace()
    kpt = qdescs[best] if best<qdescs.n else qdescs2[best-qdescs.n]
    if colors:
      color = (128,128,128)
    else:
      color = (255,0,255) if best<qdescs.n else (0,255,0)
    draw_kpt_in_img(draw,kpt,color)
  
  return img

# return an image showing keypoint matches between two frames
def viz_keypoint_matches( qvidno, qvidfile, qframe_t, dbvidno, dbvidfile, dbframe_t, ptmatches=None, q_bbox=None, db_bbox=None, ttrans='', trf=None):
  
  qimg = Image.open(open(get_video_frame( qvidfile, qframe_t) or null_image_path,'r'))
  if q_bbox: ImageDraw.Draw(qimg).rectangle( apply_ttrans_query(q_bbox,ttrans,qimg.size[0]), outline=(255,128,0), fill=None)
  print_text_in_img(qimg,"query %d%s: %.2fs"%(qvidno,'f' in ttrans and 'f' or '',qframe_t))
  
  dbimg = Image.open(open(get_video_frame( dbvidfile, dbframe_t) or null_image_path,'r'))
  if db_bbox: ImageDraw.Draw(dbimg).rectangle( apply_ttrans_base(db_bbox,ttrans), outline=(255,128,0), fill=None)
  print_text_in_img(dbimg,"retrieved %d%s: %.2fs"%(dbvidno,'h' in ttrans and 'h' or '',dbframe_t))
  
  qw = qimg.size[0]
  qh = qimg.size[1]
  pair = concatenate_two_imgs(qimg,dbimg,'v')
  draw = ImageDraw.Draw(pair)
  
  if trf and trf[0]!=0:
    trf = get_homography(trf)
    if ptmatches:
      p_q = np.matrix( list(ptmatches[0][0:2])+[1] ).T
      p_b = np.matrix( list(ptmatches[0][2:4])+[1] ).T
      trf[0,2] += p_b[0,0]*(trf[2]*p_q) - trf[0]*p_q
      trf[1,2] += p_b[1,0]*(trf[2]*p_q) - trf[1]*p_q
    inv_trf = trf.I    
    def q_to_db( pt ):  return apply_homo_pt( pt, homography )
    def db_to_q( pt ):  return apply_homo_pt( pt, inv_homo )
    for t in range(1):
      t = 25*t
      corners = [(t,t),(dbimg.size[0]-t,t),(dbimg.size[0]-t,dbimg.size[1]-t),(t,dbimg.size[1]-t)]
      for i in range(4):
        draw.line( (corners[i-1][0],corners[i-1][1]+qh,corners[i][0],corners[i][1]+qh), fill=(255,2*t,255) )
        pq0 = apply_ttrans_query( db_to_q(apply_ttransinv_base( corners[i-1], ttrans )), ttrans, qw)
        pq1 = apply_ttrans_query( db_to_q(apply_ttransinv_base( corners[i  ], ttrans )), ttrans, qw)
        draw.line( pq0+pq1, fill=(255,2*t,255) )
  
  if ptmatches: # draw keypoint matches
    for (xq,yq,xb,yb,score,qs,bs) in ptmatches:
      c = int(255-255*score)
      pb = apply_ttrans_base( (xb,yb), ttrans )
      pq = apply_ttrans_query( (xq,yq), ttrans, qw) #db_to_q((xb,yb))
      draw.line( (pq,pb[0],pb[1]+qh), fill=(0,c,0) )
  
  del draw
  return pair


def extract_patch_affine( config, imno, q_bbox, homo, for_db, size, min_patch_size=0 ):
  center = bbox_center(q_bbox)
  aff22 = [max(min_patch_size,bbox_width(q_bbox))/float(size),0.0,0.0,max(min_patch_size,bbox_height(q_bbox))/float(size)]
  aff22 = np.matrix(aff22)
  aff22.shape=(2,2)
  if for_db:
    trf = get_homography( homo )
    center = apply_homo_pt( center, trf )
    tmp22 = np.matrix([homo[0],homo[1],homo[3],homo[4]])
    tmp22.shape=(2,2)
    aff22 = tmp22 * aff22
  img_name, istemp = get_image_path( config, imno )
  img = obsidian.image_load(img_name)
  if not img:
    print >>sys.stderr, "warning: obsidian could not load image '%s'" % img_name
    return None
  if istemp:  os.remove(img_name)
  pyr = obsidian.image_pyramid_build_cpu(img, 1.2, 11, True)
  patch = obsidian.image_new(size,size)
  patch.save = patch.save = lambda name: obsidian.image_save(patch,name) and name
  def show(self): obsimg_to_Image(self).show()
  setattr(patch.__class__, 'show', show)
  assert obsidian.image_pyramid_get_aff_patch( patch, pyr, yael.fvec(list(aff22.A1)), center[0], center[1])
  return patch

def concatenate_two_patch_affine( qconfig, qimno, q_bbox, bconfig, bimno, homo, size, min_patch_size=0 ):
  qp = extract_patch_affine( qconfig, qimno, q_bbox, homo, 0, size, min_patch_size )
  bp = extract_patch_affine( bconfig, bimno, q_bbox, homo, 1, size, min_patch_size )
  return concatenate_two_images_raw(obsimg_to_Image(qp),obsimg_to_Image(bp),"query patch","db patch")

def extract_patch_similarity( config, imno, q_bbox, db_bbox, homo, for_db, min_patch_size=128, out_size=256 ):
  center = bbox_center(q_bbox)
  size = np.array([max(min_patch_size,bbox_width(q_bbox)), max(min_patch_size,bbox_height(q_bbox))])
  aff22 = [1.0,0.0,0.0,1.0]
  aff22 = np.matrix(aff22)
  aff22.shape=(2,2)
  trf = get_homography( homo )
  if for_db:
    #center = bbox_center(db_bbox)
    trf = get_homography( homo )
    center = apply_homo_pt( center, trf )
    # correct scale
    size = np.array([max(min_patch_size,bbox_width(db_bbox)), max(min_patch_size,bbox_height(db_bbox))])
    # correct angle
    angle = np.arctan2(trf[1,0]+trf[1,1], trf[0,0]+trf[0,1]) 
    angle = ((angle+2*np.pi)//(np.pi/2))
    if np.fmod(angle+4,2)==1.0: size = np.array([size[1],size[0]])
    angle *= (np.pi/2)
    tmp22 = np.matrix([np.cos(angle),-np.sin(angle),np.sin(angle),np.cos(angle)])
    tmp22.shape=(2,2)
    aff22 = tmp22 * aff22
  
  img_name, istemp = get_image_path( config, imno )
  img = obsidian.image_load(img_name)
  if not img:
    print >>sys.stderr, "warning: obsidian could not load image '%s'" % img_name
    return None
  if istemp:  os.remove(img_name)
  pyr = obsidian.image_pyramid_build_cpu(img, 1.2, 11, True)
  
  maxdim = max(size)
  size *= out_size/maxdim
  aff22 /= out_size/maxdim
  
  patch = obsidian.image_new(int(size[0]),int(size[1]))
  patch.image_save = patch.save = lambda name: obsidian.image_save(patch,name) and name
  def show(self): obsimg_to_Image(self).show()
  setattr(patch.__class__, 'show', show)
  assert obsidian.image_pyramid_get_aff_patch( patch, pyr, yael.fvec(list(aff22.A1)), center[0], center[1])
  return patch

def concatenate_two_patch_similarity( qconfig, qimno, q_bbox, bconfig, bimno, db_bbox, homo, size, min_patch_size=0 ):
  qp = extract_patch_similarity( qconfig, qimno, q_bbox, db_bbox, homo, 0, min_patch_size, size )
  bp = extract_patch_similarity( bconfig, bimno, q_bbox, db_bbox, homo, 1, min_patch_size, size )
  return concatenate_two_images_raw(obsimg_to_Image(qp),obsimg_to_Image(bp),"query patch","db patch")


def obsimg_to_Image( obsimg ):
  with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=RuntimeWarning)
    name = os.tmpnam()+'.pgm'
  obsidian.image_save(obsimg,name)
  return Image.open(name)



def usage():
  print """viz_tools.py usage:
  
  visualize intermediary results of a search in bigimbaz
  
Options:

[-db name]               db name
[-querydb name]          query db name
[-no]                    image config no
[-name]                  image display name
[-i/-info]               information about this video according to config
[-vw]                    show keypoint visual words
[-v]                     increase verbosity  
"""
  sys.exit(1)
  
  
# Main program
if __name__=='__main__':
  args = sys.argv[1:]
  
  db=querydb=None
  config = qconfig = None
  info=False
  qno=qname=None
  no=name=None
  showimg=showvw=showdesc=showrecon=showgroups=False
  showframe=None
  ohyps=None
  verbose=0
  gt=False
  pos,rad=None,20
  
  while args:
    a = args.pop(0)
    if a == '-h':      usage()
    elif a=='-db':     db=args.pop(0)
    elif a=='-querydb':querydb=args.pop(0)
    elif a=='-no':     no=int(args.pop(0))
    elif a=='-qno':    qno=int(args.pop(0))
    elif a=='-name':   name=args.pop(0)
    elif a=='-qname':  qname=args.pop(0)
    elif a=='-ohyps':  ohyps=args.pop(0)
    elif a=='-gt':     gt=True
    elif a in {'-i','-info'}: info=True
    elif a in {'-img','-showimg'}: showimg=True
    elif a in {'-desc','-showdesc'}: showdesc=True
    elif a in {'-vw','-showvw'}: showvw=True
    elif a in {'-recon','-showrecon'}: showrecon=True
    elif a in {'-groups','-showgroups'}: showgroups=True
    elif a=='-pos':    pos=eval(args.pop(0))  # example: "45,89"
    elif a=='-rad':    rad=int(args.pop(0))
    elif a=='-v':      verbose+=1
    else:
      print >>sys.stderr,"unknown arg %s"%a
      usage()
  
  from scripts.groundtruth import make_queryset, make_groundtruth
  import random
  
  # select default config
  if db: config = select_config(db)
  if querydb: qconfig= select_config(querydb)
  
  if name: no = name_to_config_no( config, name )
  if qname: qno = name_to_config_no( qconfig, qname )
  if no and not config: qno = config_no_to_query_no( qconfig, no )
  
  if not config and not no: config, no = qconfig, qno
  assert config, "error: config is %s" % str(config)
  assert no!=None, "error: please set an image number with -no/-qno"
  
  # display some info about this video
  if info:
    c = config or qconfig
    if no==None:
      print "displaying info for config '%s'" % c.dbname
      print "  containing %d image:" % c.nimg
      nos = range(c.nimg)
    else:
      nos = [no]
    for i in nos:
      w,h = hasattr(c,'image_size') and c.image_size(i) or ('unkonwn','unkonwn')
      print "  #%6d, size=(%3d,%3d), name=%s, path=%s" % (i,w,h,c.display_name(i),c.img_filename(i))
  
  # show video summary 
  if showimg or showdesc or showvw or showrecon or showgroups:
    img=get_image(config,no)
    
    if hasattr(config,'fnames') and len(config.fnames[no])>=3 and config.fnames[no][2]:
      bbox = config.fnames[no][2]
      bbox = config.coordinates_img2desc(no,bbox[:2]) + config.coordinates_img2desc(no,bbox[2:])
      print_text_in_img(img,'',bbox=bbox,bbox_color=224)
    
    from siftgeo.siftgeo import *
    if showdesc or showrecon or showgroups:
      source = config.desc_filename(no)
      vr = pointset_read(source,config.desc_format+0)
    elif showvw:
      source = config.vw_filename(no)
      vr = pointset_read(source,config.binsign_computer and 2 or 1)
    else:
      source = config.img_filename(no)
      vr = pointset_t(0)
    
    if pos: pointset_crop(vr,pos[0]-rad,pos[1]-rad,pos[0]+rad,pos[1]+rad,0)
    
    if showgroups:
      groups = cPickle.load(open(config.groups_filename(no),'r'))[0]
    else:
      groups = [range(vr.n)]
      if not groups:  groups=[[]]
    
    print "    '%s'" % source
    color = (0,255,0)
    random.seed(0)
    for g,group in enumerate(groups):
      print "  # displaying %d keypoints from group %d/%d :" % (len(group), g+1, len(groups))
      img=print_pointset_in_img(img,vr,group,bbox=hasattr(config,'queries') and config.queries[no][2] or None,color=color,verbose=(1-showrecon)*verbose)
      if verbose: raw_input()
      color = (random.randint(0,255),random.randint(0,255),random.randint(0,255))
      if showrecon:
        print_pointset_in_img(img,vr,group,bbox=hasattr(config,'queries') and config.queries[no][2] or None,recon=config,verbose=verbose)
    img.show()
  
  if gt:
    assert qno!=None, "error: qno is none"
    gt,null = make_groundtruth(qconfig,qno,config,corrected=False)
    print "displaying ground truth for query #%d '%s'" % (qno, qconfig.display_name(qno))
    print "null =",null
    print "gt =",gt
  
  if ohyps:
    assert qno and no, "error: please define a query number/name and a db number/name"
    bnos_imms = cPickle.load(open(ohyps+'/query_%d.hyps.pickle'%qno,'r'))
    gt,null = make_groundtruth(qconfig,qno,config,corrected=True)
    for bno,imm in bnos_imms:
      if bno==no: 
        status = no in gt and "TP" or no in null and "NA" or "FP"
        print "found qno=%d (%s) --> bno=%d (%s); annotated as %s" %(qno,qconfig.display_name(qno),no,config.display_name(no),status)
        for name in sorted(dir(imm)):
          print name,"=",getattr(imm,name)
        print "imm.aff=",imm.aff
        pdb.set_trace()
        sys.exit(0)
    assert False, "not found !"


























