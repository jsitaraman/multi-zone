PROJECT(mz)
ENABLE_LANGUAGE(Fortran C)
CMAKE_MINIMUM_REQUIRED(VERSION 3.8)
SET(CGNS_LOCATION /usr/local/lib)
SET(CMAKE_C_FLAGS "-fPIC -O2 -g")
SET(CMAKE_Fortran_FLAGS "-fbounds-check -fbacktrace")

set(DIST_SRCS
io.c
gridgen.c
utils.c
unique_nodes.c
ugrid_sub.c
write_cgns.c
boundaryface_interface.f90
main.c
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} /usr/local/include/)
LINK_DIRECTORIES(${CGNS_LOCATION}/lib)
ADD_EXECUTABLE(createmz ${DIST_SRCS})
TARGET_LINK_LIBRARIES(createmz cgns m)
