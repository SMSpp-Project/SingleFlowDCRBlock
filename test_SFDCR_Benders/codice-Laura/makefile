SHELL = /bin/bash

SYSTEM_32 = arm64_osx

LIBFORMAT  = static_pic

CPLEX_DIR = /Applications/CPLEX_STUDIO2211/cplex
GUROBI_DIR= /Library/gurobi1001

CC  = g++

CPLEX_LIB     = -L$(CPLEX_DIR)/lib/$(SYSTEM_32)/$(LIBFORMAT) -lcplex
CPLEX_INCLUDE = -I$(CPLEX_DIR)/include/ilcplex/

#GUROBI directories
GUROBIINCDIR = $(GUROBI_DIR)/macos_universal2/include/

#Include and libraries
GUROBI_LIB     = -L$(GUROBI_DIR)/macos_universal2/lib/ -lgurobi100
GUROBI_CPPLIB  = -L$(GUROBI_DIR)/macos_universal2/lib/ -lgurobi_c++ $(GUROBI_LIB)
GUROBI_INCLUDE = -I$(GUROBIINCDIR)

OBJECTS = main_ottimo.o Graph.o DCRGraph.o DCR_MFSP_SOCP_CPX.o DCR_MFSP_SOCP_GRB.o DCR_SPT.o DCR_INDI.o FlowSim.o DCR_3PRONG.o

main : $(OBJECTS)
	$(CC) $(COPT) -o main_ottimo $(OBJECTS) $(CPLEX_LIB) $(GUROBI_CPPLIB) -lpthread -lm 
	
#main.o: main.cpp 
main_ottimo.o: main_ottimo.cpp 
	$(CC) $(COPT) -c $(CPLEX_INCLUDE) $(GUROBI_INCLUDE) main_ottimo.cpp 
DCR_MFSP_SOCP_CPX.o: DCR_MFSP_SOCP_CPX.cpp DCR_MFSP_SOCP_CPX.h
	$(CC) $(COPT) -c $(CPLEX_INCLUDE) DCR_MFSP_SOCP_CPX.cpp 
DCR_MFSP_SOCP_GRB.o: DCR_MFSP_SOCP_GRB.cpp DCR_MFSP_SOCP_GRB.h
	$(CC) $(COPT) -c $(GUROBI_INCLUDE) DCR_MFSP_SOCP_GRB.cpp 
DCR_SPT.o: DCR_SPT.cpp DCR_SPT.h
	$(CC) $(COPT) -c DCR_SPT.cpp
DCR_INDI.o: DCR_INDI.cpp DCR_INDI.h
	$(CC) $(COPT) -c $(CPLEX_INCLUDE) DCR_INDI.cpp
FlowSim.o: FlowSim.cpp FlowSim.h
	$(CC) $(COPT) -c $(CPLEX_INCLUDE) $(GUROBI_INCLUDE) FlowSim.cpp
DCR_3PRONG.o: DCR_3PRONG.cpp DCR_3PRONG.h
	$(CC) $(COPT) -c $(CPLEX_INCLUDE) DCR_3PRONG.cpp
	
clean:
	rm *.o main
	
