
from siftgeo import siftgeo

from imbase.config import select_config

from yael import yael




def desc_dists_stats((d1,n1,desc1),(d,n2,desc2)):
  assert d1==d
  k=1
  nn_indexes=yael.ivec(n1*k)
  nn_dists=yael.fvec(n1*k)
  yael.knn_full_thread(2,
                       n1,n2,d,k,
                       desc2,desc1,
                       None,
                       nn_indexes, nn_dists, 16,None)
  # stats on distances  
  # print "min,max=",yael.fvec_min(nn_dists,n1*k),yael.fvec_max(nn_dists,n1*k)
  print "avg=",yael.fvec_sum(nn_dists,n1*k)/(n1*k)

def load_desc(fname):
  ps=siftgeo.pointset_read(fname,0)
  (v, d) = siftgeo.siftgeo_to_fvecs (ps[0], ps.n)
  v=yael.fvec.acquirepointer(v)
  return (d,ps.n,v)

# config=select_config("holidays:HolidaysDenseArnau")
config=select_config("classemes:HolidaysDenseTX")


ims=range(1)

descs=[load_desc(config.desc_filename(i)) for i in ims]

print "nb descs: ",[n for (d,n,v) in descs]


for (d,n,desc) in descs:
  n0=0
  for i in range(n):
    print yael.fvec_norm(desc.plus(i*d),d,2)
    # if yael.fvec_all_0(desc.plus(i*d),d):
    #   n0+=1
  print "desc has %d 0s"%n0


if False:
  for i,desci in enumerate(descs):
    for j,descj in enumerate(descs[:i]):
      print i,j
      desc_dists_stats(desci,descj)


  

