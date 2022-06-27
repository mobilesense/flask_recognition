import ivfgeo, siftgeo, sys

flistname = sys.argv[1]
nbvw = int(sys.argv[2])
fivfname = sys.argv[3]

binsign = len (sys.argv)>4 and sys.argv[4]=="-binsign"



print "* make new ivfgeo"

ivf = ivfgeo.ivfgeo_new (nbvw+1,5)

print "* add vwgeos"


i = 0
for fname in open (flistname,"r"):
  fname = fname.strip()
  s = siftgeo.pointset_read (fname, binsign and 2 or 1)
  print "  %d pts from %s"%(s.n,fname)

  siftgeo.vwgeoset_sort (s)
  ivfgeo.ivfgeo_add (ivf, s, i)

  i+=1

print "* storing", fivfname

ivfgeo.ivfgeo_write (fivfname, ivf)

ivfgeo.ivfgeo_delete (ivf)

