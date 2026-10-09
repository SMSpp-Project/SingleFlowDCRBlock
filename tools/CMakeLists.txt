# --------------------------------------------------------------------------- #
#    CMake file for SingleFlowDCRBlock/tools                                  #
# --------------------------------------------------------------------------- #

# This is needed for setting the runtime path when installing.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ----- dcr2nc4 ------------------------------------------------------------- #
# the multi-file textual instances of MultiFlowDCRBlock, possibly reduced to
# their first flows, written as netCDF files, which is how the instances of
# the module are distributed [see data/README.md]
add_executable(dcr2nc4 dcr2nc4.cpp)
target_link_libraries(dcr2nc4 PRIVATE SMS++::SingleFlowDCRBlock)

# ----- Install instructions ------------------------------------------------ #
include(GNUInstallDirs)
install(TARGETS dcr2nc4
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# --------------------------------------------------------------------------- #
