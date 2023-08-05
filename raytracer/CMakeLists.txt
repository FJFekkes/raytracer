include(CMakePrintHelpers)

# fetches all .*cpp files
# puts them in a variable "sources"
file(GLOB_RECURSE sources CMAKE_CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(raytracer ${sources})
set_target_properties(raytracer PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_property(TARGET raytracer PROPERTY CXX_STANDARD 20)
