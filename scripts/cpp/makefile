UNAME := $(shell uname)
current_dir = $(shell pwd)

ifeq ($(UNAME), Linux)
TARGET =  Rfunc.so
endif
ifeq ($(UNAME), Darwin)
TARGET = Rfunc_mac.so
endif

$(TARGET): sim.h sim.cpp tree.h tree.cpp Rfunc.cpp
	g++ -std=c++11 -fpic  -O3 -pipe -I/$(current_dir) -Wall -Wp,-D_FORTIFY_SOURCE=2 -lstdc++ -fexceptions -fstack-protector --param=ssp-buffer-size=4 -shared $^ -o $@
