

ot_add_current_dir_to_include_dirs ()

ot_add_source_file ( Morris.cxx )
ot_add_source_file ( MorrisExperiment.cxx )
ot_add_source_file ( MorrisExperimentGrid.cxx )
ot_add_source_file ( MorrisExperimentLHS.cxx )

ot_install_header_file ( Morris.hxx )
ot_install_header_file ( MorrisExperiment.hxx )
ot_install_header_file ( MorrisExperimentGrid.hxx )
ot_install_header_file ( MorrisExperimentLHS.hxx )

include_directories (${INTERNAL_INCLUDE_DIRS})

add_library (otmorris ${SOURCEFILES})
set_target_properties (otmorris PROPERTIES POSITION_INDEPENDENT_CODE ON)
if (BUILD_SHARED_LIBS)
  set_target_properties (otmorris PROPERTIES COMPILE_DEFINITIONS "OTMORRIS_DLL_EXPORTS")
endif ()

if ( NOT DEFINED LIB_VERSION ) 
  set ( LIB_VERSION 0.0.0 )
endif ()
if ( NOT DEFINED LIB_SOVERSION ) 
  set ( LIB_SOVERSION 0 )
endif ()
set_target_properties ( otmorris PROPERTIES VERSION ${LIB_VERSION} )
set_target_properties ( otmorris PROPERTIES SOVERSION ${LIB_SOVERSION} )
target_link_libraries (otmorris ${OPENTURNS_LIBRARY})

# Add targets to the build-tree export set
export (TARGETS otmorris FILE ${PROJECT_BINARY_DIR}/OTMORRIS-Targets.cmake)

# Install the export set for use with the install-tree
install(EXPORT OTMORRIS-Targets
        DESTINATION "${OTMORRIS_CONFIG_CMAKE_PATH}"
        COMPONENT Development)

install(TARGETS otmorris
        EXPORT  OTMORRIS-Targets 
        RUNTIME DESTINATION ${OTMORRIS_BIN_PATH}
        LIBRARY DESTINATION ${OTMORRIS_LIBRARY_PATH}
        ARCHIVE DESTINATION ${OTMORRIS_LIBRARY_PATH}
)

