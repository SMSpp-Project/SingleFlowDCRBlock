##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SingleFlowDCRBlock                                           #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by SingleFlowDCRBlock, i.e., core SMS++. These are      #
#   *not* copied into $(SFDCRINC): adding those -I directives to the compile #
#   commands will have to done by whatever "main" makefile is using this.    #
#   Analogously, any external library and the corresponding -L< libdirs >    #
#   will have to be added to the final linking command by  whatever "main"   #
#   makefile is using this.                                                  #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by SingleFlowDCRBlock.                                              #
#                                                                            #
#   The makefile defines internally (cf. MCFClssSlvr below) which            #
#   MCFSolver< :MCFClass > will be available.                                #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(SFDCRSDR)    = the directory where the source is               #
#                                                                            #
#   Output: $(SFDCROBJ)    = the final object(s) / library                   #
#           $(SFDCRH)      = the .h files to include                         #
#           $(SFDCRINC)    = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################


# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

SFDCROBJ = $(SFDCRSDR)/obj/SingleFlowDCRBlock.o \
	$(SFDCRSDR)/obj/SingleFlowDCRBendersSolver.o \
	$(SFDCRSDR)/obj/MultiFlowDCRBlock.o \
	$(SFDCRSDR)/obj/BenBound.o \
	$(SFDCRSDR)/obj/DCRLagrangianSolver.o \
	$(SFDCRSDR)/obj/SPT.o \
	$(SFDCRSDR)/obj/DCR_SPT.o 

SFDCRINC = -I$(SFDCRSDR)/include -I/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include/ilcplex/
#-I$(SFDCRSDR)/../BundleSolver/NdoFiOracle/OPTUtils

SFDCRH   = $(SFDCRSDR)/include/SingleFlowDCRBlock.h \
	$(SFDCRSDR)/include/SingleFlowDCRBendersSolver.h \
	$(SFDCRSDR)/include/MultiFlowDCRBlock.h \
	$(SFDCRSDR)/include/BenBound.h \
	$(SFDCRSDR)/include/DCRLagrangianSolver.h \
	$(SFDCRSDR)/include/SPT.h \
	$(SFDCRSDR)/include/DCR.h \
	$(SFDCRSDR)/include/DCR_SPT.h 

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SFDCROBJ) $(SFDCRSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

all: $(SFDCROBJ)
#all : $(SFDCRSDR)/obj/SingleFlowDCRBlock.o

$(SFDCRSDR)/obj/SingleFlowDCRBlock.o: $(SFDCRSDR)/src/SingleFlowDCRBlock.cpp \
	$(SFDCRSDR)/include/SingleFlowDCRBlock.h $(SMS++OBJ)
	$(CC) -c $(SFDCRSDR)/src/SingleFlowDCRBlock.cpp -o $@ \
	$(SFDCRINC) $(SMS++INC) $(SW)

$(SFDCRSDR)/obj/MultiFlowDCRBlock.o: \
	$(SFDCRSDR)/src/MultiFlowDCRBlock.cpp $(SFDCRH) $(SMS++OBJ) $(MILPOBJ)
	$(CC) -c $(SFDCRSDR)/src/MultiFlowDCRBlock.cpp -o $@ \
	$(SFDCRINC) $(SMS++INC) $(libMCFClINC) -I$(MILPINC) $(MCFClssSlvr) $(SW)

$(SFDCRSDR)/obj/BenBound.o: \
	$(SFDCRSDR)/src/BenBound.cpp $(SFDCRH) $(SMS++OBJ) $(MILPOBJ)
	$(CC) -c $(SFDCRSDR)/src/BenBound.cpp -o $@ \
	$(SFDCRINC) $(SMS++INC) $(libMCFClINC) -I$(MILPINC) $(MCFClssSlvr) $(SW)

$(SFDCRSDR)/obj/DCRLagrangianSolver.o: \
	$(SFDCRSDR)/src/DCRLagrangianSolver.cpp $(SFDCRH) $(SMS++OBJ) $(MILPOBJ)
	$(CC) -c $(SFDCRSDR)/src/DCRLagrangianSolver.cpp -o $@ \
	$(SFDCRINC) $(SMS++INC) $(libMCFClINC) -I$(MILPINC) $(MCFClssSlvr) $(SW)

$(SFDCRSDR)/obj/SPT.o: \
	$(SFDCRSDR)/src/SPT.cpp $(SFDCRH) $(SMS++OBJ) $(MILPOBJ)
	$(CC) -c $(SFDCRSDR)/src/SPT.cpp -o $@ \
	$(SFDCRINC) $(SMS++INC) $(libMCFClINC) -I$(MILPINC) $(MCFClssSlvr) $(SW)

$(SFDCRSDR)/obj/DCR_SPT.o: \
	$(SFDCRSDR)/src/DCR_SPT.cpp $(SFDCRH) $(SMS++OBJ) $(MILPOBJ)
	$(CC) -c $(SFDCRSDR)/src/DCR_SPT.cpp -o $@ \
	$(SFDCRINC) $(SMS++INC) $(libMCFClINC) -I$(MILPINC) $(MCFClssSlvr) $(SW)

$(SFDCRSDR)/obj/SingleFlowDCRBendersSolver.o: \
	$(SFDCRSDR)/src/SingleFlowDCRBendersSolver.cpp $(SFDCRH) $(SMS++OBJ) $(MILPOBJ)
	$(CC) -c $(SFDCRSDR)/src/SingleFlowDCRBendersSolver.cpp -o $@ \
	$(SFDCRINC) $(SMS++INC) $(libMCFClINC) -I$(MILPINC) $(MCFClssSlvr) $(SW)

########################## End of makefile ###################################


