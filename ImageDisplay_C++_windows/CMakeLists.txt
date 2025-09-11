cmake_minimum_required(VERSION 3.10)
project(ImageDisplay_C_windows)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source files
set(SOURCES
    Main.cpp
    Image.cpp
)

# Header files (optional, for IDEs)
set(HEADERS
    Image.h
    resource.h
)

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Windows resource file
if(MSVC)
    target_sources(${PROJECT_NAME} PRIVATE Image.rc)
endif()

# Link libraries if needed (e.g., GDI32 for Windows graphics)
target_link_libraries(${PROJECT_NAME} PRIVATE gdi32 user32)

# Set output directory
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
