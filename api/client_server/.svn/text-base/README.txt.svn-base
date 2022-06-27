
What is it?
===========

This is a client/server interface that combines Bigimbaz with a
PhotoMole geometric constraint check.

Prerequisites
=============

In addition to what comes in the CVS, you need a data/ directory that
will receive a lot of logs, temp files, etc (please link it to a
scratch).

If the database descriptors were computed with those, you also need
the older extract_features and the associated siftgeotobin.

To generate thumbnails (and to use compute_descriptors) you need
ImageMagick.

The images, descriptors and associated visual word files are defined
in config.py. The script prepare_imbase.py can prepare the data from a
set of images.


Scripts
=======

The scripts (and their dependencies) are:

* config.py: 
defines a select_config function that returns a Config objet from a
database name. The Config instance contains all the information about a
database (filenames, computation options, etc.)

* imbase.py (config): 
defines a ImBase class that does queries. It manages a data directory
where temporary information is stored. To use it, you should add
geom_filter to the PYTHONPATH:

export PYTHONPATH=$PWD/../geom_filter

* prepare_imbase.py (imbase):
computes the descriptors of an image base, prepares the visual words,
the PCA transform file and the inverted file.

* test_imbase.py (imbase):
test script that shows how to use the image database.

* protocol.py 
the communication protocol by which a Client object calls methods of a
Server object.

* test_protocol.py (protocol)
test script that shows how to use the client/server architecture

* server.py (imbase,protocol)
defines an ImbaseServer object and, used as a script, launches a
bigimbaz server

* client.py (protocol)
script that connects to an ImbaseServer and does a query.

The server is ready for use from a HTTP proxy (see the ../web/
subdirectory)

