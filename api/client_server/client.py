import sys,socket,cPickle,pdb

from protocol import *

def usage():
  sys.stderr.write("usage: %s [-port port] [-search imname] [-ns n_short] [-add imname]* [-store dirname] [-cdm] [-cdmreset alpha] [-cdmloadivf file] [-cdmloadfac file] [-status] [-decompose] server_name\n"%sys.argv[0])
  sys.stderr.write("-decompose: decompose search steps")
  sys.exit(1)

HOST=None
qim=None
toadd=[]
port=PORT
n_short=20
store_dirname=None
get_status=False
cdm=False
cdm_alpha=0.0
cdm_loadivf=''
cdm_loadfac=''
args=sys.argv[1:]
decompose=False

while args:
  arg=args.pop(0)
  if arg in ("-h","--help"):     usage()
  elif arg=='-status':           get_status=True
  elif arg=='-add':              toadd.append(args.pop(0))
  elif arg=='-port':             port=int(args.pop(0))
  elif arg=='-ns':               n_short=int(args.pop(0))
  elif arg=='-search':           qim=args.pop(0)
  elif arg=='-store':            store_dirname=args.pop(0)
  elif arg=='-decompose':        decompose=True
  elif arg=='-cdm':              cdm=True
  elif arg=='-cdmreset':         cdm_alpha=float(args.pop(0))
  elif arg=='-cdmloadivf':       cdm_loadivf=args.pop(0)
  elif arg=='-cdmloadfac':       cdm_loadfac=args.pop(0)
  else:
    if not HOST:
      HOST=arg
    else:
      usage()
  

imBase=Client(HOST,port)

imBase.log('yo!')

print "rid=",imBase.rid

#print "get n=",imBase.get_n()

if get_status:
  print "ps:\n",imBase.get_ps_stats()
  print "rid:",imBase.rid
  print "running_rids:",imBase.get_server_pool()
  print "important_rids:",imBase.get_importants()


if toadd:
  radd=[]
  for i,fname in enumerate(toadd):
    radd.append(imBase.put_datafile(fname,'added_im_%d'%i))
  print "adding images"
  imBase.add_images(radd)
  

if cdm_loadivf or cdm_loadfac:
  imBase.cdm_load(cdm_loadivf, cdm_loadfac)

if cdm:
  imBase.cdm_op("iterate",(cdm_alpha,))
else:
  if cdm_alpha:
    imBase.cdm_reset(cdm_alpha)


if qim:

  imBase.set_n_short(n_short)

  image_file=imBase.put_datafile(qim,'image')

  print "image file is ",imBase.get_mimetype(image_file)

  if decompose:

    # bsiftgeo : descriptor file name
    # descimage : image (+ellipste) file name
    (n_int,n_desc,bsifgeo)=imBase.compute_descriptors(image_file)

    print "got %d int pts -> %d descs"%(n_int,n_desc)

    # Retrieve the image file containing the ellipse regions and store
    # it into a given target file name
    #imBase.get_datafile(descimage,"/tmp/descs.png")

    # Compute the visual words associated with these descriptors
    vwf = imBase.compute_vw(bsifgeo)

    print "vwf done, doing query"

    # Query the inverted file using this visual frequency vector
    res=imBase.query(vwf)


    # Display the result
    print res

    # Geometric filtering using on a subset of images
    print "filter:"

    res2=imBase.filter(bsifgeo,[imno for (imno,dist) in res])

  else:
    res2=imBase.process_search(image_file)

  print "res2=",res2


if store_dirname:
  imBase.store_imbase(store_dirname)
