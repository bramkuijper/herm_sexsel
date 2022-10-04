EXE_HERM=herm_mutual_direct_benefits.exe
EXE_GONO=gonochorist_mutual_direct_benefits.exe
CPP_HERM=hermaphrodite_mutual_direct_benefit.cpp
CPP_GONO=gonochorist_mutual_direct_benefit.cpp

CXX=g++
CXXFLAGS=-Wall -O3

all : $(EXE_HERM) $(EXE_GONO)

$(EXE_HERM) : $(CPP_HERM)
	$(CXX) $(CXXFLAGS) -o $(EXE_HERM) $(CPP_HERM)

$(EXE_GONO) : $(CPP_GONO)
	$(CXX) $(CXXFLAGS) -o $(EXE_GONO) $(CPP_GONO)
