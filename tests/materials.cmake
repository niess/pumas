include(ExternalProject)

set (MATERIALS_DIR "${CMAKE_BINARY_DIR}/materials")

ExternalProject_Add (Materials
    GIT_REPOSITORY "https://github.com/niess/pumas-materials"
    SOURCE_DIR "${CMAKE_BINARY_DIR}/materials"
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
