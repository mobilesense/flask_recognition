"""
Simple query test. Does not use the network.

"""

import sys

from imbase import *

from config.config import select_config

import time,sys

from basics import *

import pdb


def usage():
  print >> sys.stderr, "usage: test_imbase.py [-db dbname] [-add search]* [-store db] [-search img]* [-stat]"
  sys.exit(1)


qims=[]
toadd=[]
store_dirname=None
decompose=False
dbname=None
stat=False
draw=False

args=sys.argv[1:]

while args:
  arg=args.pop(0)
  if arg in ("-h","--help"):     usage()
  elif arg=='-add':              toadd.append(args.pop(0))
  elif arg=='-search':           qims.append(args.pop(0))
  elif arg=='-store':            store_dirname=args.pop(0)
  elif arg=='-decompose':        decompose=True
  elif arg=='-db':               dbname=args.pop(0)
  elif arg=='-stat':             stat=True
  elif arg=='-draw':             draw=True
  else:
    print >> sys.stderr, "unknown option %s "%arg
    
  
config=select_config(dbname)

imBase=ImBase(CachedData(config))
imBase.verbose=1

imBase.cachedData.cache_all()

#bb=imBase.cachedData.get_BigBase()
#gv=imBase.cachedData.get_GeometricVerification()

#wgc_type=0x50545
#bb.set_params(scale_w=0,he_dist_w=0,
#              burstiness=0x12,
#              tfidf=1,
#              wgc_type=wgc_type,
#              he_thres=22) 






print "datadir=",imBase.datadir


def sb(x):
  n=int(x*70)
  sys.stdout.write("+"*n+"-"*(70-n)+"\r")
  sys.stdout.flush()


if toadd!=[]:
  print "************** Add a images to db"

  n_map=imBase.add_images_files(toadd,keep_files="",state_obj=State())

  print "images got number ",n_map

if store_dirname:
  print "************** storing db to",store_dirname
  imBase.store_imbase(store_dirname)


if qims!=[]:

  for qim in qims:
    print "************** query %s"%qim
    (res,res2,_)=imBase.search(qim,state=State())    
    mnames = map(lambda x:config.img_filename(x[0]), res2)
    scores = map(lambda x:x[1], res2)
    print zip(mnames, scores)

if stat:
  print imBase.config.nimg
  for attr in dir(ivfgeo.cvar.query_stats):
    if attr!="this" and not attr.startswith("_") and not attr=="as_hist":
      print attr,"=",getattr(ivfgeo.cvar.query_stats,attr)


if draw:
  imBase.draw_pts("~/workspace/images/maison.jpg","~/workspace/obsidian/src/a.desc","test.jpg")
  
    
