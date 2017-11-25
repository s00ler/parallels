CC := g++-7
# replace this with your correct compiler as identified above

ARCH := core2 # Replace this with your CPU architecture.
# core2 is pretty safe for most modern machines.

CFLAGS := -march=$(ARCH) -O3 -fopenmp -m64 -std=c++11

COMPILE_COMMAND := $(CC) $(CFLAGS)

OUTPUT := test.out

all: main.cpp
	$(COMPILE_COMMAND) -o $(OUTPUT) main.cpp

clean:
	rm -f *.o $(OUTPUT).*
