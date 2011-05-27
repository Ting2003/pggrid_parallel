##########################################################
# CUDA Makefile                                       #
#                                                            #
# Author      : Ting Yu: tingyu1@illinois.edu           #
# # Version     : 1.0                                          #
# # Date        : 05/03/2011                                  #
# # Discription : generic Makefile for making CUDA programs    #
# ##############################################################
#
BIN               := pg 
CXXFLAGS          :=-Wall -Wextra -pipe -O2 -msse4.2 -msse3 -mfpmath=sse -march=native -gstabs
#-O3 -g
CXX 		  :=g++

CUDA_INSTALL_PATH = /usr/local/cuda
CUDA_SDK_PATH = /home/ac/tingyu1/NVIDIA_GPU_Computing_SDK

LIBSUFFIX  :=_x86_64

NVCC ?= $(CUDA_INSTALL_PATH)/bin/nvcc
INCD = -I$(CUDA_SDK_PATH)/C/common/inc -I$(CUDA_INSTALL_PATH)/include -I$(CUDA_SDK_PATH)/shared/inc 

LIBS = -L/usr/local/libcuda -lcuda -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcublas -lcufft -L$(CUDA_SDK_PATH)/C/common/lib/linux -lcutil_x86_64 -L$(CUDA_SDK_PATH)/shared/lib/linux -L$(CUDA_SDK_PATH)/C/lib -lstdc++ -lpthread

CUDA_SDK?=3
COMMONFLAGS = -DCUDA_SDK=$(CUDA_SDK)
NVCCFLAGS := -G -g #--ptxas-options=-v -O3 -G -g 


# files
CPP_SOURCES       := util.cpp point.cpp node.cpp net.cpp parser.cpp\
		     vec.cpp main.cpp triplet.cpp algebra.cpp \
		     block.cpp circuit.cpp

CU_SOURCES        := circuit_host.cu circuit_kernel.cu
HEADERS           := $(wildcard *.h)
CPP_OBJS          := $(patsubst %.cpp, %.o, $(CPP_SOURCES))
CU_OBJS           := $(patsubst %.cu, %.cu_o, $(CU_SOURCES))

# packages
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

CHOLMOD=$(PACKAGE)/CHOLMOD
CHOLMOD_LIB_DIR=$(CHOLMOD)/Lib
CHOLMOD_INC_DIR=$(CHOLMOD)/Include
CHOLMOD_LIB=$(CHOLMOD_LIB_DIR)/libcholmod.a \
	    $(PACKAGE)/AMD/Lib/libamd.a

PACKAGE_LIB =$(UMFPACK_LIB) $(CHOLMOD_LIB)
PACKAGE_INC = -I$(UMFPACK_INC_DIR) -I$(CHOLMOD_INC_DIR)\

%.cu_o : %.cu
#circuit_host.o: circuit_host.cu
	$(NVCC) $(NVCCFLAGS) -c $(INCD) -I$(UMFPACK_INC_DIR) -I$(CHOLMOD_INC_DIR) -o $@ $<

%.o: %.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCD) -I$(UMFPACK_INC_DIR) -I$(CHOLMOD_INC_DIR) -o $@ $<

$(BIN): $(CPP_OBJS) $(CU_OBJS)
	$(CXX) -o $(BIN) $(CU_OBJS) $(CPP_OBJS) $(LDFLAGS) $(INCD)$(UMFPACK_LIB) $(CHOLMOD_LIB) $(LIBS)

clean:
	rm -f $(BIN) *.o *.cu_o
