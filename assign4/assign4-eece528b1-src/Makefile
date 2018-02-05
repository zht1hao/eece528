all: makelib GEV.elf
	@echo ==========================================
	@echo this is the VXM version of GE
	@echo run it with the matrix size
	@echo for example:
	@echo ./GEV.elf 4096
	@echo Notice: use make clean to clean the temps
	@echo ==========================================

GEV.elf: 
	gcc GEV.c -c -DVBX_SIMULATOR -Ivbxapi -O2 -o GEV.o 
	g++ GEV.o  vbxsim/libvbxsim.a vbxapi/libvbxapi.a -O2  -o GEV.elf


makelib:
	make -C ./vbxsim SIMULATOR=true
	make -C ./vbxapi SIMULATOR=true

clean:
	rm -f GEV.o
	rm -f GEV.elf
	make -C ./vbxsim clean SIMULATOR=true
	make -C ./vbxapi clean SIMULATOR=true

.phony: all clean
