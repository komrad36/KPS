EXECUTABLE_NAME=KPS
CPP=g++
INC=
CPPFLAGS=-Wall -Wextra -Werror -pedantic -Wshadow -Ofast -std=gnu++14 -mavx
LDFLAGS=-Wall -Wextra -Werror -pedantic -Wshadow -Ofast -std=gnu++14 -mavx
LIBS=-lGeographic
CPPSOURCES=$(wildcard *.cpp)

OBJECTS=$(CPPSOURCES:.cpp=.o)

all: $(CPPSOURCES) $(EXECUTABLE_NAME)

$(EXECUTABLE_NAME) : $(OBJECTS)
	$(CPP) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)
	
%.o:%.cpp
	$(CPP) -c $(INC) $(CPPFLAGS) $< -o $@
clean:
	rm -rf *.o $(EXECUTABLE_NAME)
