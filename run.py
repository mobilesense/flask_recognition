from app import app
from app import db
import gflags
import logging
import os
import sys

# tornado
from tornado.wsgi import WSGIContainer
from tornado.httpserver import HTTPServer
from tornado.ioloop import IOLoop

if __name__ == '__main__':
    # TODO tornado it threaded ??????????? for multiple clients
    #gflags.FLAGS(sys.argv)
    if not os.path.exists('db.sqlite'):
        db.create_all()    
    try:
        os.makedirs(UPLOAD_FOLDER)
    except Exception as err:
        pass
    logging.getLogger().setLevel(logging.INFO)
    #app.model = core.yarilabcore(ivf_file=FLAGS.ivf_file,
    #                          clust_file=FLAGS.clust_file)
    port = 5000
    print ('listening to port %d ..'%port)
    http_server = HTTPServer(WSGIContainer(app))    
    http_server.listen(port)
    IOLoop.instance().start()

