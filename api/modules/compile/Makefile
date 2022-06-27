

# order matters

ALLDIRS=yael utils siftgeo ivfgeo geom_filter compute_descriptor

# ivfann

ALLINITPY=$(ALLDIRS)

all: $(addprefix compiledir.,$(ALLDIRS)) init_py

clean: $(addprefix cleandir.,$(ALLDIRS)) rm_init_py



compiledir.%: 
	cd ../; cd $(subst compiledir.,,$@); make

cleandir.%: 
	cd ../; cd $(subst cleandir.,,$@); make clean


# make empty __init__.py files to access the modules with dotted notation


init_py: __init__.py $(addsuffix /__init__.py,$(ALLINITPY)) 
	touch __init__.py

%__init__.py:
	touch ../$@

__init__.py:
	touch ../$@

rm_init_py:
	rm -f __init__.py $(addsuffix /__init__.py,$(ALLINITPY)) 

