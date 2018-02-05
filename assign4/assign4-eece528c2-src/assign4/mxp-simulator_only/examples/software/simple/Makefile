
############################################################
# This makefile only works for the simulator target so far #
############################################################

all: test-cpp.elf test-c.elf

VBXSIM=../../../repository/lib/vbxsim
VBXAPI=../../../repository/lib/vbxapi


test-cpp.elf: test-cpp.cpp $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a
	g++ -Wall -g -DVBX_SIMULATOR -I$(VBXAPI) -o $@ $^

test-c.elf: test-c.c $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a
	gcc -c -Wall -g -std=c99 -DVBX_SIMULATOR -I$(VBXAPI) -o test-c.o test-c.c
	g++ test-c.o  $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a  -o $@


$(VBXAPI)/libvbxapi.a  $(VBXSIM)/libvbxsim.a:
	make -C $(dir $@) SIMULATOR=true

clean:
	rm -rf test-c.elf test-c.o
	rm -rf test-cpp.elf
	make -C $(VBXAPI) clean SIMULATOR=true
	make -C $(VBXSIM) clean SIMULATOR=true

.phony: all clean
