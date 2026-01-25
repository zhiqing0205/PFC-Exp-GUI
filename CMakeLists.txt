cmake_minimum_required(VERSION 3.20)

project(PFCExpGUI VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PFC_GRID_L "256" CACHE STRING "Grid size L (compile-time)")
set(PFC_GRID_M "256" CACHE STRING "Grid size M (compile-time)")
set(PFC_GRID_N "1" CACHE STRING "Grid size N (compile-time)")

set(PFC_BIN_DIR "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PFC_BIN_DIR}")
if(CMAKE_CONFIGURATION_TYPES)
  foreach(cfg ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER "${cfg}" cfg_upper)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${cfg_upper} "${PFC_BIN_DIR}/${cfg}")
  endforeach()
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)
find_package(FFTW3 REQUIRED)

add_executable(pfc_exp_cli
  src/experiment/misfit_uniform_time_seed_fftw_early.cpp
)
target_compile_definitions(pfc_exp_cli PRIVATE
  PFC_GRID_L=${PFC_GRID_L}
  PFC_GRID_M=${PFC_GRID_M}
  PFC_GRID_N=${PFC_GRID_N}
)
target_link_libraries(pfc_exp_cli PRIVATE FFTW3::fftw3)
if(UNIX AND NOT APPLE)
  target_link_libraries(pfc_exp_cli PRIVATE m)
endif()
set_target_properties(pfc_exp_cli PROPERTIES OUTPUT_NAME "pfc-exp-cli")

add_executable(pfc_exp_gui WIN32 MACOSX_BUNDLE
  src/gui/main.cpp
  src/gui/MainWindow.cpp
  src/gui/MainWindow.h
  src/gui/PlotWidget.cpp
  src/gui/PlotWidget.h
  src/gui/resources.qrc
)
target_link_libraries(pfc_exp_gui PRIVATE Qt6::Widgets)
set_target_properties(pfc_exp_gui PROPERTIES OUTPUT_NAME "PFC-Exp-GUI")
