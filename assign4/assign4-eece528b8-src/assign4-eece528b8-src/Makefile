all: assign4.elf

VBXSIM=mxp-simulator_only/repository/lib/vbxsim/
VBXAPI=mxp-simulator_only/repository/lib/vbxapi/

assign4.elf: assign4.c $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a
	gcc -c -Wall -g -std=c99 -DVBX_SIMULATOR -I$(VBXAPI) -o assign4.o assign4.c
	g++ assign4.o  $(VBXAPI)/libvbxapi.a $(VBXSIM)/libvbxsim.a  -o $@

$(VBXAPI)/libvbxapi.a  $(VBXSIM)/libvbxsim.a:
	make -C $(dir $@) SIMULATOR=true

clean:
	rm -rf assign4.elf assign4.o
	make -C $(VBXAPI) clean SIMULATOR=true
	make -C $(VBXSIM) clean SIMULATOR=true

.phony: all clean
