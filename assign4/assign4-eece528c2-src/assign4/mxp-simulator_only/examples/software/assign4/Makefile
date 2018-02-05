
############################################################
# This makefile only works for the simulator target so far #
############################################################

all: assignment4.elf

VBXSIM=../../../repository/lib/vbxsim
VBXAPI=../../../repository/lib/vbxapi


assignment4.elf: assignment4.c $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a
	gcc -O3 -c -Wall -g -std=c99 -DVBX_SIMULATOR -I$(VBXAPI) -o assignment4.o assignment4.c
	g++ -O3 assignment4.o  $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a  -o $@


$(VBXAPI)/libvbxapi.a  $(VBXSIM)/libvbxsim.a:
	make -C $(dir $@) SIMULATOR=true

clean:
	rm -rf assignment4.o
	rm -rf assignment4.elf
	make -C $(VBXAPI) clean SIMULATOR=true
	make -C $(VBXSIM) clean SIMULATOR=true

.phony: all clean
