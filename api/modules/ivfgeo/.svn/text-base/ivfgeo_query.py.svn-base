

import ivfgeo, siftgeo, sys, time

from yael import yael
#from utils import utils

ivfname=sys.argv[1]
qfilename=sys.argv[2:]
nq = len (qfilename);


#anglescale=len(sys.argv)>3
anglescale = 0

multiple_query = 0
vwfmt = 2
maxret=10

print "* Read invfile",ivfname

ivf=ivfgeo.ivfgeo_read(ivfname, 1)

print "* Compute the tf-idf weights"

ivfgeo.ivfgeo_compute_tfidf(ivf)

ivf.he_thres = 24

print "* Load query " #, qfilename

s={}
scompound = siftgeo.pointset_new_multiple (nq)
for i in range (nq):
  s[i] = siftgeo.pointset_read(qfilename[i], vwfmt)
  siftgeo.vwgeoset_sort(s[i])
  siftgeo.pointset_addto_multiple (scompound, s[i], i)

print "* querying, anglescale=",anglescale 

t0 = time.time()

if multiple_query:
  (res_i,dis_f)=ivfgeo.ivfgeo_query_multiple_peek(ivf,scompound,maxret,nq,None)

elif not anglescale:
  for i in range (nq):
    (res_i,dis_f)=ivfgeo.ivfgeo_query(ivf,s[i],maxret)

else:
  for i in range (nq):
    (res_i,dis_f)=ivfgeo.ivfgeo_query_scale_angle_peek (ivf,s[i],maxret,None)

ressz=min(maxret * nq,ivf.nbvec)

print time.time() - t0

# convert output to Python lists

res_i=yael.IntArray.frompointer(res_i)
res_i.this.acquire()

dis_f=yael.FloatArray.frompointer(dis_f)
dis_f.this.acquire()


res=[(res_i[i],dis_f[i]) for i in range(ressz)]

print "* results"

#print res

ivfgeo.ivfgeo_delete(ivf)

import _ivfgeo
_ivfgeo.delete_FloatArray(dis_f)
_ivfgeo.delete_IntArray(res_i)
