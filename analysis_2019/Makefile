# -------------------------------------------------
# -------     Universidad de los Andes      -------
# -------      Departamento de Física       -------
# -------    Proyecto Joven Investigador    -------
# -------  Andrés Felipe García Albarracín  -------
# -------    Juan Carlos Sanabria Arenas    -------
# -------          February - 2015          -------
# -------------------------------------------------
# This makefile serves to compile c++ programs using HepMC and Root libraries
#
# Usage of the program
# Write one of the next possibilities:
# make compile                    To compile without Root nor HepMC
# make compile_ROOT               To compile with Root and without HepMC
# make compile_HepMC              To compile without Root and with HepMC
# make compile_ROOT_HepMC         To compile with Root and HepMC

# Modify the following variables according to your c++ compiler and your 
# local Root and HepMC installation folders. If you do not use any of these
# programs, leave blank.

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
#ROOTFLAGLIB=root-config --cflags --glibs
ROOTLIB=$(ROOTSYS)/lib
ROOTLIBS='-lPhysics'
#ROOTPHYLIB=/usr/local/HEP/ROOT6.14.02/lib # Directory where is located libPhysics.so
#ROOTPHYLIB=/home/aflorez/Documents/root/lib # Directory where is located libPhysics.so

# Directory of Exrootanalysis
#EXROOTLIB=/usr/local/HEP/MG5_aMC_v2_6_3_2/ExRootAnalysis
#EXROOTINC=/usr/local/HEP/MG5_aMC_v2_6_3_2/ExRootAnalysis/ExRootAnalysis
#EXROOTINC=/home/aflorez/SIMULATIONS/MG5_aMC_v2_6_4/Delphes/external/ExRootAnalysis/
#EXROOTLIB=/home/aflorez/SIMULATIONS/MG5_aMC_v2_6_4/ExRootAnalysis/
EXROOTINC=/home/aflorez/SIMULATIONS/MG5_aMC_v2_6_4/Delphes-3.4.1/external/ExRootAnalysis

#EXROOTINC=/home/aflorez/SIMULATIONS/MG5_aMC_v2_6_4/ExRootAnalysis/ExRootAnalysis

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
	$(CXX) -o $(TARGETFILE) $(MAINFILE) $(SOURCEFILES)`root-config --cflags --glibs` \
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
	$(CXX) -o $(TARGETFILE) $(MAINFILE) $(SOURCEFILES) -Wall `root-config --ldflags --glibs --cflags` \
         -I$(DELPHESINC) -I$(HEPMCINC) -I$(EXROOTINC) \
         -L$(DELPHESLIB) -L$(HEPMCLIB) -lDelphes -L$(EXROOTLIB) -lExRootAnalysis -L$(ROOTLIB) $(ROOTLIBS)

