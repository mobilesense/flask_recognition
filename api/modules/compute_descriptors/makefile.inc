# -*- makefile -*- 


IMAGEMAGICK=yes


ifeq "$(IMAGEMAGICK)" "yes"

  MAGICKCONFIG=Magick++-config
  MAGICKCPPFLAGS=`$(MAGICKCONFIG) --cppflags`
  MAGICKLDFLAGS= -L/usr/X11R6/lib `$(MAGICKCONFIG) --ldflags --libs` -ldl 

else 

  MAGICKCPPFLAGS=-DNOIMAGEMAGICK 
  MAGICKLDFLAGS=
endif





# CFLAGS = -g -O3 -c -Wall -fPIC  -Wno-deprecated 
CFLAGS = -g -O3 -c -Wall -fPIC  -Wno-deprecated 

# CFLAGS =  -g -c -Wall -fPIC  -Wno-deprecated  $$( Magick++-config --cflags )

CC = g++

.SUFFIXES = .c .cc .cpp 

#C_FILES = $(FILES:=.c)
CPP_FILES = $(FILES:=.cpp)
O_FILES = $(FILES:%=%.o)
SHAREDLIB = $(LIB:=.so)

%.o: %.cc   
	$(CC) $(CFLAGS) $< -o $@

%.o: %.cpp   
	$(CC) $(CFLAGS) $< -o $@ $(INCLUDES) $(EXTRACFLAGS)

$(LIB):  $(O_FILES)
	 ar rv $@ $?
	 ranlib $@

targets: $(LIB) 


