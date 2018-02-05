############################################################
# This makefile only works for the simulator target so far #
############################################################

all: hw4_opt.elf

VBXSIM=./mxp-simulator_only/repository/lib/vbxsim
VBXAPI=./mxp-simulator_only/repository/lib/vbxapi

hw4_opt.elf: hw4_opt.c $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a
	gcc -c -Wall -g -std=c99 -DVBX_SIMULATOR -I $(VBXAPI) -DDIM=4096 -o hw4_opt.o hw4_opt.c -O2 
	g++ hw4_opt.o  $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a  -o $@ -O2

$(VBXAPI)/libvbxapi.a  $(VBXSIM)/libvbxsim.a:
	make -C $(dir $@) SIMULATOR=true

clean:
	rm -rf hw4_opt.elf hw4_opt.o
	make -C $(VBXAPI) clean SIMULATOR=true
	make -C $(VBXSIM) clean SIMULATOR=true

.phony: all clean
