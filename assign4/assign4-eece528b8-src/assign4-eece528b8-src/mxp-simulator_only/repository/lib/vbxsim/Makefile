
UNAME := $(shell uname)

SRC_LIB:=libvbxsim_linux.a

ifeq ($(UNAME),Linux)
SRC_LIB:=libvbxsim_linux.a
endif

ifeq ($(UNAME),Darwin)
SRC_LIB:=libvbxsim_mac.a
endif

ifdef ComSpec
SRC_LIB=libvbxsim_mingw32.a
endif
ifdef COMSPEC
SRC_LIB=libvbxsim_mingw32.a
endif

DST_LIB = libvbxsim.a

.PHONY:$(DST_LIB) clean
$(DST_LIB):$(SRC_LIB)
	cp $^ $@
clean:
	rm -f $(DST_LIB)
