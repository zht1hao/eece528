
############################################################
# This makefile only works for the simulator target so far #
############################################################

all: vectorAssign4.elf assignment4.elf

VBXSIM=../../../../repository/lib/vbxsim
VBXAPI=../../../../repository/lib/vbxapi


vectorAssign4.elf: vectorAssign4.c $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a
	gcc -c -Wall -g -std=c99 -DVBX_SIMULATOR -I$(VBXAPI) -o vectorAssign4.o vectorAssign4.c
	g++ vectorAssign4.o  $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a  -o $@

assignment4.elf: assignment4.c $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a
	gcc -c -Wall -g -std=c99 -DVBX_SIMULATOR -I$(VBXAPI) -o assignment4.o assignment4.c
	g++ assignment4.o  $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a  -o $@

$(VBXAPI)/libvbxapi.a  $(VBXSIM)/libvbxsim.a:
	make -C $(dir $@) SIMULATOR=true

clean:
	rm -rf vectorAssign4.o
	rm -rf vectorAssign4.elf
	rm -rf assignment4.o
	rm -rf assignment4.elf
	make -C $(VBXAPI) clean SIMULATOR=true
	make -C $(VBXSIM) clean SIMULATOR=true

.phony: all clean
