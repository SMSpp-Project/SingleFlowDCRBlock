SHELL = /bin/bash

LIBFORMAT  = static_pic

SYSTEM_32 = x86-64_linux

CPLEX_DIR = /opt/ibm/ILOG/CPLEX_Studio201/cplex

CC  = g++

CPLEX_LIB     = -L$(CPLEX_DIR)/lib/$(SYSTEM_32)/$(LIBFORMAT) -lcplex
CPLEX_INCLUDE = -I$(CPLEX_DIR)/include/ilcplex

OBJECTS = Graph.o DCRGraph.o DCRLagrangianSolver.o SPT.o CPX_LAG.o BenBound.o CPL_BEND.o main.o #main_heur.o #raw_data_collector.o #test.o 

test: $(OBJECTS)
	$(CC) -O3 -o test $(OBJECTS) $(CPLEX_LIB) -lpthread -lm -ldl
#main: $(OBJECTS)
#	$(CC) -O3 -o main $(OBJECTS) $(CPLEX_LIB) -lpthread -lm -ldl

main_heur.o: main_heur.cpp 
	$(CC)  -c -O3 $(CPLEX_INCLUDE) main_heur.cpp 	
raw_data_collector.o: raw_data_collector.cpp 
	$(CC)  -c -O3 $(CPLEX_INCLUDE) raw_data_collector.cpp 	
test.o: test.cpp 
	$(CC)  -c -O3 $(CPLEX_INCLUDE) test.cpp 	
main.o: main.cpp 
	$(CC)  -c -O3 $(CPLEX_INCLUDE) main.cpp 
Graph.o: Graph.C Graph.h
	$(CC)  -c -O3 Graph.C 
DCRGraph.o: DCRGraph.cpp DCRGraph.h
	$(CC)  -c -O3 DCRGraph.cpp
DCRLagrangianSolver.o: DCRLagrangianSolver.cpp DCRLagrangianSolver.h
	$(CC)  -c -O3 DCRLagrangianSolver.cpp 
SPT.o: SPT.cpp SPT.h
	$(CC)  -c -O3 SPT.cpp 	
CPX_LAG.o: CPX_LAG.cpp CPX_LAG.h
	$(CC)  -c -O3 $(CPLEX_INCLUDE) CPX_LAG.cpp   
BenBound.o: BenBound.cpp BenBound.h
	$(CC)  -c -O3 BenBound.cpp 
CPL_BEND.o: CPL_BEND.cpp CPL_BEND.h
	$(CC)  -c -O3 $(CPLEX_INCLUDE) CPL_BEND.cpp 

clean:
	rm *.o main
	
