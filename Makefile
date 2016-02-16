PROG=DVMcpp
SRC= pugixml.cpp BaseTypes.cpp DVMBase.cpp main.cpp

#Set compiler and version based on whether a HPC, normal PC with or without MATLAB
HOST=$(shell domainname | sed 's/\.//g')
HOSTHPC=cx1hpcicacuk
MATLAB=0

ifeq ($(HOST),$(HOSTHPC))
	OBJ=$(SRC:.cpp=.o)
	OBJ+=$(SRCF:.f90=.o)
	CXX=icpc
	#F90=ifort
	CXXFLAGS=-g3 -openmp -O3
	LDFLAGS=-openmp -O3
	CXXVER=$(shell icpc --version | head -1 | sed 's/(.*)//;s/  */_v/')
else
	ifeq ($(MATLAB),0)
		OBJ=$(SRC:.cpp=.o)
		OBJ+=$(SRCF:.f90=.o)
        CXX=g++
		CLIBS=""
		CXXFLAGS=-g3 -O3 -fopenmp 
		CFLAGS=""
        #Note: the framework accelerate is required under MAC and replaces lpack
		LDFLAGS=-O3 -I/opt/local/lib -fopenmp -larmadillo -framework Accelerate
		CXXVER=$(shell g++ --version | head -1 | sed 's/(.*)//;s/  */_v/')
	else
		#SRC += MatlabEngine.cpp MatlabKinematics.cpp
        SRC += MatlabEngine.cpp MatlabKinematics.cpp FloatingBodyHeaveMoored.cpp
		OBJ=$(SRC:.cpp=.o)

        #Path under Ubuntu
		#MATLAB_BASE=/usr/local/MATLAB/R2012a
        #LD_LIBRARY_PATH=$(MATLAB_BASE)/bin/glnx64

        #Path under MAC
        MATLAB_BASE =/Applications/MATLAB_R2013b.app
		LD_LIBRARY_PATH=$(MATLAB_BASE)/bin/maci64
        INCLUDE_PATH=$(MATLAB_BASE)/

        # Under Ubuntu
		#MATLABLIB= -leng -lmx -lut -lmat -lm -lgcc_s\
			   -lpthread -lc -lrt -lmwi18n -lmwfl -ldl -lexpat\
			   -lmwresource_core -lmwMATLAB_res -lz  -lcrypt
 
        # Under MAC
        MATLABLIB= -leng -lmx -lut -lmat -lm\
			   -lpthread -lc -lmwi18n -lmwfl -ldl -lexpat\
			   -lmwresource_core -lmwMATLAB_res -lz

		#-lstdc++ #do not link to lstdc++ - will blow up Matlab! 

		# CXX=g++-4.4 #requires g++-4.4 under linux
    
        CXX=g++
		CLIBS="$CLIBS -lstdc++"
		CXXFLAGS=-O2 -fopenmp -I$(MATLAB_BASE)/extern/include -ansi -D_GNU_SOURCE -D__BEMMATLAB__
		CFLAGS="$CFLAGS -fopenmp -fexceptions -fPIC -fno-omit-frame-pointer"
		
        # Under Ubuntu ...
        #LDFLAGS=-fopenmp -O2 -L$(MATLAB_BASE)/bin/glnxa64 $(MATLABLIB)\
			-Wl,-rpath,$(MATLAB_BASE)/bin/glnxa64
		
        # Under mac ...
        LDFLAGS=-fopenmp -O2 -L$(LD_LIBRARY_PATH) $(MATLABLIB)\
			-Wl,-rpath,$(LD_LIBRARY_PATH)

        
        CXXVER=$(shell g++ --version | head -1 | sed 's/(.*)//;s/  */_v/')
	endif
endif


#Set version number if revision controlled
GITSTAT=$(shell git status 2>&1 grep 'fatal' | grep -o 'On branch master')
HASGIT=On branch master
ifeq ($(GITSTAT),$(HASGIT))
	REVNO=$(shell git log | grep -c commit)
else
	REVNO="N/A"
endif

#Write various information to timestamp file
tstamp=$(shell date "+%A the %e %B %Y, %H:%M:%S")
$(shell echo "Revision $(REVNO) last compiled using $(CXXVER) on $(tstamp)" >timestamp)

# Compile any external Fortran code

# Begin C++ compilations

$(PROG):$(OBJ)
	@echo ""
	@echo "Building Binaries"
	@echo ""
	$(CXX) $(CXXFLAGS) -o $(PROG) $(OBJ) $(LDFLAGS)


BaseTypes.o: BaseTypes.cpp BaseTypes.hpp 
	$(CXX) $(CXXFLAGS) -c $<

DVMBase.o: DVMBase.cpp DVMBase.hpp BaseTypes.hpp
	$(CXX) $(CXXFLAGS) -c $<

main.o: main.cpp DVMBase.hpp DVMBase.cpp BaseTypes.hpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: cleanobj, clean

clean:
	touch $(OBJ) $(PROG)
	rm $(OBJ) $(PROG)
	@echo "Removed objects and exec"

cleanobj:
	touch $(OBJ)
	rm $(OBJ)
	@echo "Removed all objects"
