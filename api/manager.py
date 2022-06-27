from modules.imbase import *
from modules.api_config.api_config import select_config
from modules.basics import *
import Image
import cv2
import time
import sys
import pdb

class Manager:
    def __init__(self):
        #TODO shall we  create a copy for each user ?? 
        #or prevent from concurrent access, see: ThreadsafeCachedData in server.py
        self.overall_config = select_config("api_config:Overall")
        self.cached = CachedData(self.overall_config)
        self.cached.fc = self.cached.get_FeaturesComputer() #TODO see if there is shared things, exp: img_filename,..
        self.cached.vwa = self.cached.get_VWAssignement()
        self.cached.gv = self.cached.get_GeometricVerification()
        self.cached.bb = {}
        
class UserManager():
    def __init__(self, user_id, cached):  
        config = select_config('api_config:Yari20k(%d)'%user_id)
        cached.config = config #we adapt config to fc, vwa, gv
        if not os.path.isfile(config.invfilename):
            cached.bb[user_id] = cached.get_BigBase(make_new=True)
        else:
            cached.bb[user_id] = cached.get_BigBase()
        self.imBase = ImBase(cached, user_id)
        self.imBase.verbose = 1        
	
    def match(self, imname, includePtMatches=False):
        """ 
        Match image
        res = (imno,nmatch,score,aff)
        imno   = image number
        nmatch = # of matching points
        score  = computed from matching points and distances (vote)
        aff    = affine matrix	    
        """
        with TicToc("matching time"):
            (_, res, ptmatches) = self.imBase.search(imname, state=State(), includePtMatches=includePtMatches)     
            indexes = map(lambda x:x[0], res)					 
            scores = map(lambda x:x[1], res)					
            #mnames = map(lambda x:m.config.img_filename(x[0]), res2)		
        return indexes, scores
		
    def add(self, imnames):
        """
        Add images to ivfgeo index
        """  
        with TicToc("adding time"):      
            n_map = self.imBase.add_images_files(imnames, keep_files="", state_obj=State()) 
            print "[add] images got number ", n_map

    def remove(self, imnames):
        #TODO
        pass
						
if __name__=='__main__':
    action = sys.argv[1]
    user_id = 1
    m = manager('api_config:Yari20k(%d)'%user_id) 
    imname = '/home/amine/Bureau/yari2/api/server/processing/datasets/holidays/images/67.jpg' 
    if action == 'search'	:		
        indexes, scores = m.match(imname)				
        print indexes
        print scores
    elif action == 'add':
        m.add([imname])	



  
		    
