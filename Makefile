EXECUTABLE_NAME=KPS
CC=g++-4.9
NVCC=nvcc
ARCH=sm_30
INC=-I/usr/local/cuda/include/
NVCCFLAGS=-Wall -Wextra -Werror -Wshadow -Ofast -mavx
CFLAGS=-Wall -Wextra -Werror -pedantic -Wshadow -Ofast -std=gnu++14 -mavx
CCLIBS=
LDFLAGS=-Wall -Wextra -Werror -pedantic -Wshadow -Ofast -std=gnu++14 -mavx
STATIC_LIBS=-L/usr/local/cuda/lib64 -ldl -lrt -lpthread -lcudart_static -lGeographic
LIBS=-L/usr/local/cuda/lib64 -lcudart -lGeographic
CPPSOURCES=$(wildcard *.cpp)
CUSOURCES=$(wildcard *.cu)

OBJECTS=$(CPPSOURCES:.cpp=.o) $(CUSOURCES:.cu=.o)

all: $(CPPSOURCES) $(CUSOURCES) $(EXECUTABLE_NAME)

$(EXECUTABLE_NAME) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

static: $(CPPSOURCES) $(CUSOURCES) $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE_NAME) $(STATIC_LIBS)
	
%.o:%.cpp
	$(CC) -c $(INC) $(CFLAGS) $< -o $@

%.o:%.cu
	$(NVCC) -arch=$(ARCH) -O3 -ccbin $(CC) -std=c++11 -c $(INC) -Xcompiler "$(NVCCFLAGS)" $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE_NAME)
