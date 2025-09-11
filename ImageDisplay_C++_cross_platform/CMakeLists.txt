cmake_minimum_required(VERSION 4.0)

project(ImageDisplay)

set(CMAKE_CXX_STANDARD 17)
set(WX_DOMAIN ./dependency/wxWidgets/include/wx)

include_directories(${WX_DOMAIN})
# Add source files
set(SOURCES
    src/Main.cpp
    manifest/MyImageApplication.rc
    ${WX_DOMAIN}/wx.h
    ${WX_DOMAIN}/dcbuffer.h

)


if(WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -static")
    set(wxBUILD_SHARED OFF)
    # or you can set it to ON, but you may face some problems about can't find dlls, in this case, you can copy dlls to your exe directory
else()
    set(wxBUILD_SHARED ON)
endif()

add_subdirectory(dependency/wxWidgets)

# Create executable
add_executable(MyImageApplication ${SOURCES})

# X11 required on Linux
if(UNIX AND NOT APPLE)
    find_package(X11 REQUIRED)
endif()

# Link libraries
target_link_libraries(MyImageApplication wx::core wx::base ${X11_LIBRARIES})