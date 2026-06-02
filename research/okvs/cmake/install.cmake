

#############################################
#            Install                        #
#############################################


configure_file("${CMAKE_CURRENT_LIST_DIR}/findDependancies.cmake" "findDependancies.cmake" COPYONLY)
configure_file("${CMAKE_CURRENT_LIST_DIR}/preamble.cmake" "preamble.cmake" COPYONLY)

# make cache variables for install destinations
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)


# generate the config file that is includes the exports
configure_package_config_file(
        "${CMAKE_CURRENT_LIST_DIR}/Config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/okvsConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/okvs
        NO_SET_AND_CHECK_MACRO
        NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

if (NOT DEFINED okvs_VERSION_MAJOR)
    message("\n\n\n\n warning, okvs_VERSION_MAJOR not defined ${okvs_VERSION_MAJOR}")
endif ()

set_property(TARGET bandokvs PROPERTY VERSION ${okvs_VERSION})

# generate the version file for the config file
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/okvsConfigVersion.cmake"
        VERSION "${okvs_VERSION_MAJOR}.${okvs_VERSION_MINOR}.${okvs_VERSION_PATCH}"
        COMPATIBILITY AnyNewerVersion
)

# install the configuration file
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/okvsConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/okvsConfigVersion.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/findDependancies.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/preamble.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/okvs
        )

# install library
install(
        TARGETS bandokvs
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        EXPORT okvsTargets)

# install headers
install(
        DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../frontend"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/"
        FILES_MATCHING PATTERN "*.h")


# install config
install(EXPORT okvsTargets
        FILE okvsTargets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/okvs
        NAMESPACE research::
        )
export(EXPORT okvsTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/okvsTargets.cmake"
        NAMESPACE research::
        )