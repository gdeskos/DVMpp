find_package(Doxygen)
if (DOXYGEN_FOUND)
    message(STATUS "Doxygen found: type \"make doc\" to build the documentation")
	add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOC_DIR}/Doxyfile
		WORKING_DIRECTORY ${DOC_DIR}
        COMMENT "Generating Documentation with Doxygen" VERBATIM)
else()
    message(STATUS "Doxygen not found: think about getting it to look at the documentation")
endif(DOXYGEN_FOUND)
