# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.
# Optionally glob, but only for CMake 3.12 or later:
file(GLOB HEADER_LIST CONFIGURE_DEPENDS
        "${NLPositivity_SOURCE_DIR}/include/nl-positivity/*.h"
        "${NLPositivity_SOURCE_DIR}/include/nl-positivity/*.hpp")

file(GLOB SOURCE_LIST CONFIGURE_DEPENDS
        "${NLPositivity_SOURCE_DIR}/src/*.h"
        "${NLPositivity_SOURCE_DIR}/src/*.hpp"
        "${NLPositivity_SOURCE_DIR}/src/*.cc"
        "${NLPositivity_SOURCE_DIR}/src/*.cpp")

add_library(nl-positivity ${SOURCE_LIST} ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(nl-positivity PUBLIC "${NLPositivity_SOURCE_DIR}/include/")

target_link_libraries(nl-positivity PUBLIC nlnum-lib PRIVATE ${Python3_LIBRARIES} cppitertools prettyprint OpenMP::OpenMP_CXX)

# IDEs should put the headers in a nice place
source_group(TREE "${PROJECT_SOURCE_DIR}/include" PREFIX "Header Files" FILES ${HEADER_LIST})
