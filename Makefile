CC=g++
#CPLUSPLUS=g++
SRC= util.cpp point.cpp node.cpp circuit.cpp net.cpp parser.cpp vec.cpp \
    main.cpp triplet.cpp algebra.cpp block.cpp 
#hash_mat.cpp map_mat.cpp 
HDR=$(SRC:.cpp=.h)
OBJ=$(SRC:.cpp=.o) 
BIN=pg
RELEASE=IPGS
CPPFLAGS=
CFLAGS=-Wall -Wextra -pipe -O2 -msse4.2 -mssse3 -mfpmath=sse -march=native
#CFLAGS=-Wall -g #-Wextra -pipe -O2 -msse4.2 -mssse3 -mfpmath=sse -march=core2
#LDFLAGS=-s -Wl,-O1,-hash-style=gnu
LDFLAGS=
CSCOPEFILES=cscope.files cscope.out cscope.po.out
CSCOPE=/usr/local/bin/cscope

PACKAGE= ./package_ck

UMFPACK=./umfpack
UMFPACK_LIB_DIR=$(UMFPACK)/lib
UMFPACK_INC_DIR=$(UMFPACK)/include
UMFPACK_LIB=$(UMFPACK_LIB_DIR)/libumfpack.a \
	    $(UMFPACK_LIB_DIR)/libamd.a \
	    $(UMFPACK_LIB_DIR)/libcholmod.a \
	    $(UMFPACK_LIB_DIR)/libcolamd.a \
            $(UMFPACK_LIB_DIR)/libccolamd.a \
            $(UMFPACK_LIB_DIR)/libcamd.a \
            $(UMFPACK_LIB_DIR)/libmetis.a \
            $(UMFPACK_LIB_DIR)/libgoto2.a

CHOLMOD= $(PACKAGE)/CHOLMOD
CHOLMOD_LIB_DIR=$(CHOLMOD)/Lib
CHOLMOD_INC_DIR=$(CHOLMOD)/Include
CHOLMOD_LIB=$(CHOLMOD_LIB_DIR)/libcholmod.a \
	    $(PACKAGE)/AMD/Lib/libamd.a

main: $(OBJ)
	@echo "Making project..."
	$(CC) $(LDFLAGS) -o $(BIN) $(OBJ) $(UMFPACK_LIB) $(CHOLMOD_LIB)

# test if cscope is installed
	@if test -f $(CSCOPE); then \
		make tags; \
	fi

release: $(OBJ)
	$(CC) $(LDFLAGS) -static -o $(BIN) $(OBJ) $(UMFPACK_LIB) $(CHOLMOD_LIB)

test: 
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -I$(UMFPACK_INC_DIR)\
	-I$(CHOLMOD_INC_DIR) -o test test.cpp $(UMFPACK_LIB) $(CHOLMOD_LIB)

all: main
	@echo "Making all..."

%.o: %.cpp  %.h global.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(UMFPACK_INC_DIR) -I$(CHOLMOD_INC_DIR) -c $<  -o $@

tags: $(SRC) $(HDR) main.cpp 
	@echo "Generating tags..."
	@find . -maxdepth 1 -name "*.h" -o -name "*.c" \
		-o -name "*.cpp" > cscope.files
	@cscope -bkq -i cscope.files
	@ctags -L cscope.files --sort=yes --c++-kinds=+p --fields=+iaS --extra=+q


.PHONY : clean
clean:
	@echo "Cleaning all..."
	rm -rf *.o $(OBJ) $(DBG) $(BIN) tags $(CSCOPEFILES)
