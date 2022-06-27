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
    if not os.path.exists('db.sqlite.subscription'):
        db.create_all()    
    port = 5000
    print 'listening to port %d ..'%port
    http_server = HTTPServer(WSGIContainer(app))    
    http_server.listen(port)
    IOLoop.instance().start()

