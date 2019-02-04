# C++ Compiler
CXX=g++

# Name of the main file and the target file
MAINFILE=analysis.cpp
TARGETFILE=analysis
OTHERFILES=ROOTFunctions.h DelphesFunctions.h #mt2_bisect.h  # Other files
#SOURCEFILES=mt2_bisect.cpp #mt2bl_bisect.cpp mt2w_bisect.cpp 

# Directory where HepMC is installed
HEPMCDIR=/usr/local/HEP/MG5_aMC_v2_6_3_2/HEPTools/hepmc/
HEPMCLIB=$(HEPMCDIR)/lib
HEPMCINC=$(HEPMCDIR)/include/

# Root flags and libs
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)  
ROOTGLIBS=$(shell root-config --glibs)  

# Directory of Exrootanalysis
EXROOTDIR=/home/aflorez/SIMULATIONS/MG5_aMC_v2_6_4/ExRootAnalysis
EXROOTLIB=$(EXROOTDIR)/lib
EXROOTINC=$(EXROOTDIR)/ExRootAnalysis


# Directory of Delphes
DELPHESDIR=/home/aflorez/SIMULATIONS/MG5_aMC_v2_6_4/Delphes-3.4.1
DELPHESLIB=${DELPHESDIR}
DELPHESINC=${DELPHESDIR}

all: 
	@echo " Usage: write one of the next possibilities:"
	@echo "	make compile    			To compile without Root nor HepMC"
	@echo "	make compile_ROOT			To compile with Root and without HepMC"
	@echo "	make compile_HepMC    			To compile without Root and with HepMC"
	@echo "	make compile_ROOT_HepMC    		To compile with Root and HepMC"
	@echo "	make compile_ROOT_Delphes		To compile with Root and Delphes"
	@echo " To edit make parameters type vi Makefile"


# To compile without Root nor HepMC
compile: $(MAINFILE) $(OTHERFILES)
	$(CXX) -o $(TARGETFILE) $(MAINFILE) $(SOURCEFILES)

# To compile with Root and without HepMC
compile_ROOT: $(MAINFILE) $(OTHERFILES)
	$(CXX) -o $(TARGETFILE) $(MAINFILE) $(SOURCEFILES)`root-config --cflags` \
          -I$(EXROOTINC) \
          -L$(EXROOTLIB) -lExRootAnalysis -L$(ROOTPHYLIB) -lPhysics

# To compile without Root and with HepMC
compile_HepMC: $(MAINFILE) $(OTHERFILES)
	$(CXX) -o $(TARGETFILE) $(MAINFILE) $(SOURCEFILES)\
	 -I$(HEPMCINC) \
	 -L$(HEPMCLIB) -lHepMC

# To compile with Root and HepMC
compile_ROOT_HepMC: $(MAINFILE) $(OTHERFILES)
	$(CXX) -o $(TARGETFILE) $(MAINFILE) $(SOURCEFILES) `root-config --cflags --glibs` \
         -I$(HEPMCINC) -I$(EXROOTINC) \
         -L$(HEPMCLIB) -lHepMC -L$(EXROOTLIB) -lExRootAnalysis -L$(ROOTLIB) -lPhysics

# To compile with Root and Delphes
compile_ROOT_Delphes: $(MAINFILE) $(OTHERFILES)
	$(CXX) -o $(TARGETFILE) $(MAINFILE) $(SOURCEFILES) -Wall $(ROOTCFLAGS) 	\
	-I$(DELPHESINC) -I$(EXROOTINC)	\
	-L$(DELPHESLIB) -lDelphes -L$(EXROOTLIB) -lExRootAnalysis $(ROOTLIBS)

