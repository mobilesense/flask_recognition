"""
Simple query test. Does not use the network.

"""

import sys

from imbase import *

from config import select_config

import time

from utils import *


dbname=len(sys.argv)>=2 and sys.argv[1] or 'Nister2' 

# query
qimfile=len(sys.argv)>=3 and sys.argv[2] or None

config=select_config(dbname)

imBase=ImBase(CachedData(config))
imBase.verbose=1

print "datadir=",imBase.datadir

if not qimfile:
  if dbname in ['Nister', 'Nister2','SubNister']:
    qimfile="/home/cornwall/douze/img/ukbench00000.jpg"
    
  elif dbname=='Fotolia':
    qimfile=config.img_filename(0)



def sb(x):
  n=int(x*70)
  sys.stdout.write("+"*n+"-"*(70-n)+"\r")
  sys.stdout.flush()



todo=['add','desc','vw','query','filter','allin1']
# todo=['allin1']


if 'add' in todo and dbname=='Fotolia':
  print "************** Add a few images to db"
  
  iml=['/scratch2/cepheus/hjegou/dataset/nister/jpg/ukbench%05d.jpg'%i for
       i in range(2)]
  
  imBase.add_images(iml,sb,state_obj=State())

  qimfile="/scratch2/cepheus/hjegou/dataset/nister/jpg/ukbench00000.jpg"


if 'allin1' in todo:
  imBase.search(qimfile,state=State())


if 'desc' in todo:
  
  print "************** Compute descriptors"

  (n_interest,n_desc,bsiftgeo)=imBase.compute_descriptors(qimfile,sb)
  print

  print "got %d int pts -> %d descs in file %s"%(n_interest,n_desc,bsiftgeo)

  # print "drawn on: ",imBase.draw_query_descs(qimfile,bsiftgeo,200)
 
else:
  bsiftgeo="data/06297/descs.bsiftgeo"
  

if 'vw' in todo:

  print "************** Assign VW"

  t0=time.time()
  vwf=imBase.compute_vw(bsiftgeo,sb)
  t1=time.time()
  print
  print "assign: %.3f s"%(t1-t0)

else:

  vwf="data/00289/descs.vw"


if 'query' in todo:

  print "************** search"

  print "vw file %s"%vwf

  imBase.set_n_short(20)

  res=imBase.query(vwf,sb)

  open("test_imbase_res.py","w").write("res=%s\n"%res)

  print

else:
  res=[(0, 0.47531073169761506), (1, 1.1452782724057031), (3, 1.171089736894577), (2, 1.1961956488989949), (41, 1.474244194622806), (8, 1.4876234421758441), (43, 1.4979429686860444), (11, 1.5067992359097706), (42, 1.5069680666112766), (40, 1.5180471123751254), (79, 1.6364951077743255), (33, 1.6420074027063747), (131, 1.6491285125847848), (77, 1.6492662067681547), (9, 1.6539597766355163), (78, 1.6627002203087846), (76, 1.6654817598912046), (128, 1.6734370971701422), (129, 1.6787883078229346), (130, 1.6849774677027831)]
  

if 'filter' in todo:

  imBase.geom_thr=0.3
  imBase.geom_fthr=0.15

  res2=imBase.filter(bsiftgeo,[imno for (imno,dist) in res],sb)

  print

  print "filtered result:"

  for imno,npt,nmatch,score,aff in res2:
    print "imno %d (%s) %d %g aff=%s"%(
      imno,config.img_filename(imno),nmatch,score,aff)
    print "drawing on ",imBase.draw_superpose(qimfile,imno,aff)
