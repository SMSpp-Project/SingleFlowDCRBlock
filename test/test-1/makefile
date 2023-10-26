##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for DCRBlock / MCFSolver)                              #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
NAME = SFDCRBlock_test

# basic directory
DIR = .

# debug switches
#SW = -g -glldb -fno-inline -std=c++17 -ferror-limit=1
# debug switches with address sanitizer and extra pedantic warning
#SW = -g3 -glldb -fno-inline -std=c++17 -ferror-limit=1 -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter
# production switches
SW = -O3 -std=c++17 -DNDEBUG -DCLANG_1200_0_32_27_PATCH -ferror-limit=1

# compiler
CC = clang++ -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(DIR)/$(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/$(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# SFDCRBlock
SFDCRSDR = ..
include $(SFDCRSDR)/makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(SFDCROBJ)

# libraries
MLIB =  $(SFDCRLIB)

$(DIR)/$(NAME): $(MOBJ) $(DIR)/test.o
	$(CC) -o $(DIR)/$(NAME) $^ $(MLIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(SFDCRINC)

# includes
MH =    $(SFDCRH)

# compile command

$(DIR)/test.o: $(DIR)/test.cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

############################ End of makefile #################################
