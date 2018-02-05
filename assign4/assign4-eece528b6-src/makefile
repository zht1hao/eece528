all: vector.elf

VBXSIM=vbxsim/
VBXAPI=vbxapi/

vector.elf: vector.c $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a
	gcc -c -Wall -g -std=c99 -DVBX_SIMULATOR -I$(VBXAPI) -o vector.o vector.c
	g++ vector.o  $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a  -o $@

$(VBXAPI)/libvbxapi.a  $(VBXSIM)/libvbxsim.a:
	make -C $(dir $@) SIMULATOR=true

clean:
	rm -rf vector.elf vector.o
	make -C $(VBXAPI) clean SIMULATOR=true
	make -C $(VBXSIM) clean SIMULATOR=true

.phony: all clean
