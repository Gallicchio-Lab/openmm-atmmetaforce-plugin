set(WRAP_FILE ATMMetaForcePluginWrapper.cpp)
set(MODULE_NAME atmmetaforce)

# Execute SWIG to generate source code for the Python module.

add_custom_command(
    OUTPUT "${WRAP_FILE}"
    COMMAND "${SWIG_EXECUTABLE}"
        -python -c++
        -o "${WRAP_FILE}"
        "-I${OPENMM_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/atmmetaforceplugin.i"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/atmmetaforceplugin.i" "${CMAKE_CURRENT_SOURCE_DIR}/ATMMetaForceUtils.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Compile the Python module.

add_custom_target(PythonInstall DEPENDS "${WRAP_FILE}")
set(ATMMETAFORCEPLUGIN_HEADER_DIR "${CMAKE_SOURCE_DIR}/openmmapi/include")
set(ATMMETAFORCEPLUGIN_LIBRARY_DIR "${CMAKE_BINARY_DIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
add_custom_command(TARGET PythonInstall
    COMMAND "${PYTHON_EXECUTABLE}" -m pip install -e .
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

#tests
file(COPY ${CMAKE_SOURCE_DIR}/python/tests DESTINATION ${CMAKE_BINARY_DIR}/python )
add_custom_target(CopyTestFiles ALL DEPENDS ${CMAKE_BINARY_DIR}/python/tests )
add_test(NAME abfe
COMMAND pytest WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python/tests/
)

