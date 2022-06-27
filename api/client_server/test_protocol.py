
from protocol import *

import sys,time,pdb

class MyRH(Server):

  id_counter=0

  def __init__(self,s):
    Server.__init__(self,s,MyRH.id_counter)
    MyRH.id_counter+=1
    
  def long_func(self,t):    
    time.sleep(t)
    return "finished long func"

  def add(self,x,y):
    return x+y

  def echo(self,s):
    return "echo "+s

  def alloc(self,sz):
    self.block=range(sz)

port=12025

if sys.argv[1]=='client':
  o=Client('localhost',port)
  
  o.log('coucou')

  rid=o.rid
  print "rid=",rid

  print "test add ",o.add(3,4)
  
  print "test several clients"

  for i in range(5):
    print Client('localhost',port).echo("coucou %d"%i)
  
  print "test put file"

  dist_name=o.put_datafile("/etc/fstab","fstab_%d"%rid)

  print "test get unknown file"
  
  try:
    o.get_datafile(dist_name+'totot',"/tmp/fstab_%d"%rid)
  except IOError,e:
    print "error ",e

  print "test get file"
  o.get_datafile(dist_name,"/tmp/fstab_%d"%rid)

  print "test async call"

  o.async_at_next_call()
  o.long_func(1)

  print "waiting..."

  print "result:",o.get_result()

  print "test async call with resume"
  
  o.detach_at_next_call()

  o.long_func(5)

  while True:
    time.sleep(1)
    o=Client('localhost',port)
    o.log("coucou")
    print "rid %d trying to attach"%o.rid
    try:
      output=o.attach(rid)
      break
    except CannotChangeRID,e:
      print "Still waiting..."
    o.sock.close()

  print "final result:",output
  
  print "re-test normal call"
  print o.echo("yo!")


if sys.argv[1]=='server':
  # pdb.set_trace()
  run_server(MyRH,port)

