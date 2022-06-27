
	def stat(self):
		print self.imBase.config.nimg
		for attr in dir(ivfgeo.cvar.query_stats):
			if attr!="this" and not attr.startswith("_") and not attr=="as_hist":
				print attr,"=",getattr(ivfgeo.cvar.query_stats,attr)


	def draw(self, qname, dbname, pts):
		""" 
		Draw matched points between the 2 images 
		"""
		img = self.concat(Image.open(qname), Image.open(dbname))
		qimg = cv2.imread(qname)
		dbimg = cv2.imread(dbname)
		for pt in pts:
			dbpt = pt[0]
			qpt = pt[1]			
			m = pt[2]
			cv2.circle(dbimg, (int(dbpt[0]),int(dbpt[1])), 4,(0,0,255),-1)			
			cv2.circle(qimg, (int(qpt[0]),int(qpt[1])), 4,(0,0,255),-1)			
			#cv2.line()
		cv2.imshow('qimg', qimg)
		cv2.imshow('dbimg', dbimg)		
		cv2.waitKey()

	def concat(self, qim, dbim):
	    """ 
	    Return a concatenated 2 images
	    """
		qw = qim.size[0]
		width = dbim.size[0] + qw
		height = max( qim.size[1], dbim.size[1] )
		img = Image.new( "RGBA", (width, height) )
		img.paste(qim, (0, 0))
		img.paste(dbim, (qw, 0))
		return img
		  
  def draw_pts(self,image_name,bsiftgeo,out_name,max_size=-1):
    """ return an image with interest points drawn on it """
    im=geom_filter.CRGBImage(image_name)
    sz=max(im.getxsize(),im.getysize())
    if max_size>0 and sz>max_size:
      scale=max_size/float(sz)
      im=geom_filter.downscale(im,scale); im.thisown=True
    else:
      scale=1.0
    imdesc=geom_filter.loadImdesc(bsiftgeo,True); imdesc.thisown=True
    geom_filter.drawInterestPoints(
      im,imdesc,scale,geom_filter.CRGBPixel(255,0,0),True)
    im.store(out_name)
  
  def draw_query_descs(self,image_name,bsiftgeo,max_size=-1):
    out_name=self.datadir+"query_descs.png"
    self.draw_pts(image_name,bsiftgeo,out_name,max_size)
    return out_name

  def draw_db_descs(self,dbno,max_size=-1):
    image_name=self.db_im_name(dbno)
    bsiftgeo=self.config.desc_filename(dbno)
    out_name=self.datadir+"db_descs_%d.png"%dbno
    self.draw_pts(image_name,bsiftgeo,out_name,max_size)
    return out_name
      
  def draw_superpose(self,query_image_name,dbno,aff):
    """ draw images superposed with the found affine transform """
    im0=geom_filter.CRGBImage(query_image_name)
    im1=geom_filter.CRGBImage(self.db_im_name(dbno))
    da=geom_filter.DoubleArray(6)
    for i in range(6): da[i]=aff[i]
    geom_filter.drawContrasted(im0,im1,da.cast())
    outname=self.datadir+"superposed_%d.png"%dbno
    im0.store(outname)
    return outname
