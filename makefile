##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SingleFlowDCRBlock                                           #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the     #
#   extent in which they are needed by the parts of SMS++ used by            #
#   SingleFlowDCRBlock.                                                      #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the libSMS++ library itself                     #
#           $(SFDCRBkSDR)  = the directory where the source is               #
#                                                                            #
#   Output: $(SFDCRBkOBJ)  = the final object(s) / library                   #
#           $(SFDCRBkH)    = the .h files to include                         #
#           $(SFDCRBkINC)  = the -I$( source directory )                     #
#                                                                            #
#                             Antonio Frangioni                              #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

SFDCRBkOBJ = $(SFDCRBkSDR)/obj/SingleFlowDCRBlock.o \
	$(SFDCRBkSDR)/obj/MultiFlowDCRBlock.o \
	$(SFDCRBkSDR)/obj/SingleFlowDCRBendersSolver.o \
	$(SFDCRBkSDR)/obj/BenBound.o \
	$(SFDCRBkSDR)/obj/DCRLagrangianSolver.o \
	$(SFDCRBkSDR)/obj/DCR_SPT.o \
	$(SFDCRBkSDR)/obj/SPT.o

SFDCRBkINC = -I$(SFDCRBkSDR)/include

SFDCRBkH   = $(SFDCRBkSDR)/include/SingleFlowDCRBlock.h \
	$(SFDCRBkSDR)/include/MultiFlowDCRBlock.h \
	$(SFDCRBkSDR)/include/SingleFlowDCRBendersSolver.h \
	$(SFDCRBkSDR)/include/BenBound.h \
	$(SFDCRBkSDR)/include/DCRLagrangianSolver.h \
	$(SFDCRBkSDR)/include/DCR_SPT.h \
	$(SFDCRBkSDR)/include/DCR.h \
	$(SFDCRBkSDR)/include/SPT.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SFDCRBkOBJ) $(SFDCRBkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(SFDCRBkSDR)/obj/SingleFlowDCRBlock.o: \
	$(SFDCRBkSDR)/src/SingleFlowDCRBlock.cpp $(SFDCRBkH) $(SMS++H) \
	$(SMS++OBJ)
	$(CC) -c $(SFDCRBkSDR)/src/SingleFlowDCRBlock.cpp -o $@ \
	$(SFDCRBkINC) $(SMS++INC) $(SW)

$(SFDCRBkSDR)/obj/MultiFlowDCRBlock.o: \
	$(SFDCRBkSDR)/src/MultiFlowDCRBlock.cpp $(SFDCRBkH) $(SMS++H) \
	$(SMS++OBJ)
	$(CC) -c $(SFDCRBkSDR)/src/MultiFlowDCRBlock.cpp -o $@ \
	$(SFDCRBkINC) $(SMS++INC) $(SW)

$(SFDCRBkSDR)/obj/SingleFlowDCRBendersSolver.o: \
	$(SFDCRBkSDR)/src/SingleFlowDCRBendersSolver.cpp $(SFDCRBkH) \
	$(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SFDCRBkSDR)/src/SingleFlowDCRBendersSolver.cpp -o $@ \
	$(SFDCRBkINC) $(SMS++INC) $(SW)

$(SFDCRBkSDR)/obj/BenBound.o: \
	$(SFDCRBkSDR)/src/BenBound.cpp $(SFDCRBkH)
	$(CC) -c $(SFDCRBkSDR)/src/BenBound.cpp -o $@ \
	$(SFDCRBkINC) $(SMS++INC) $(SW)

$(SFDCRBkSDR)/obj/DCRLagrangianSolver.o: \
	$(SFDCRBkSDR)/src/DCRLagrangianSolver.cpp $(SFDCRBkH)
	$(CC) -c $(SFDCRBkSDR)/src/DCRLagrangianSolver.cpp -o $@ \
	$(SFDCRBkINC) $(SMS++INC) $(SW)

$(SFDCRBkSDR)/obj/DCR_SPT.o: \
	$(SFDCRBkSDR)/src/DCR_SPT.cpp $(SFDCRBkH)
	$(CC) -c $(SFDCRBkSDR)/src/DCR_SPT.cpp -o $@ \
	$(SFDCRBkINC) $(SMS++INC) $(SW)

$(SFDCRBkSDR)/obj/SPT.o: \
	$(SFDCRBkSDR)/src/SPT.cpp $(SFDCRBkSDR)/include/SPT.h
	$(CC) -c $(SFDCRBkSDR)/src/SPT.cpp -o $@ \
	$(SFDCRBkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
