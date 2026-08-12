cmake_minimum_required(VERSION 3.25)
project(nbody_iv_conditioner VERSION 0.0.3 LANGUAGES C Fortran)

# ------------------------------------------------------------
# Python
# ------------------------------------------------------------

find_package(
    Python
    REQUIRED COMPONENTS
        Interpreter
        Development.Module
        NumPy
)
include(UseF2Py)

# ------------------------------------------------------------
# Fortran sources
# ------------------------------------------------------------

set(FORTRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../fortran)
set(FORTRAN_SOURCES
    ${FORTRAN_DIR}/utils.f90
    ${FORTRAN_DIR}/conditioners.f90
)

# ------------------------------------------------------------
# Fortran library
# ------------------------------------------------------------

add_library(nbody_iv_conditioner_fortran STATIC
    ${FORTRAN_SOURCES}
)
set_target_properties(
    nbody_iv_conditioner_fortran
    PROPERTIES
        Fortran_MODULE_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/mod
)
target_include_directories(
    nbody_iv_conditioner_fortran
    PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}/mod
)

# ------------------------------------------------------------
# F2PY module
# ------------------------------------------------------------
f2py_add_module(
    _core
    ${FORTRAN_DIR}/conditioners.f90
    ${FORTRAN_DIR}/utils.f90
)

add_dependencies(_core nbody_iv_conditioner_fortran)
target_link_libraries(_core PRIVATE nbody_iv_conditioner_fortran)
target_include_directories(_core PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/mod)

# ------------------------------------------------------------
# Install Python extension
# ------------------------------------------------------------
install(TARGETS _core DESTINATION nbody_iv_conditioner)

# ------------------------------------------------------------
# Install Python package
# ------------------------------------------------------------
install(
    DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/nbody_iv_conditioner/
    DESTINATION
        nbody_iv_conditioner
    FILES_MATCHING
        PATTERN "*.py"
)