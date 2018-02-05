GCC := gcc
GENERAL_FLAGS := -std=gnu99 -O3

SERIAL_SOURCE := gausElim.c

# this passes the DIM define
64_FLAG := -D DIM=64
256_FLAG := -D DIM=256
1024_FLAG := -D DIM=1024
2048_FLAG := -D DIM=2048
4096_FLAG := -D DIM=4096

all: s

s: s_64 s_256 s_1024 s_2048 s_4096

s_64:
	$(GCC) $(GENERAL_FLAGS) $(64_FLAG) $(SERIAL_SOURCE) -o serial_64
s_256:
	$(GCC) $(GENERAL_FLAGS) $(256_FLAG) $(SERIAL_SOURCE) -o serial_256
s_1024:
	$(GCC) $(GENERAL_FLAGS) $(1024_FLAG) $(SERIAL_SOURCE) -o serial_1024
s_2048:
	$(GCC) $(GENERAL_FLAGS) $(2048_FLAG) $(SERIAL_SOURCE) -o serial_2048
s_4096:
	$(GCC) $(GENERAL_FLAGS) $(4096_FLAG) $(SERIAL_SOURCE) -o serial_4096

s_debugging:
	$(GCC) $(GENERAL_FLAGS) $(SERIAL_SOURCE) -o serial_debugging
