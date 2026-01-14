cmake_minimum_required(VERSION 3.21)

project(QtExperimentRunner VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_GUI "Build Qt GUI client" ON)
option(BUILD_SIM "Build MPI+FFTW simulation" ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
include(GNUInstallDirs)

if(BUILD_SIM)
  if(NOT DEFINED MPI_SKIP_COMPILER_WRAPPER)
    set(MPI_SKIP_COMPILER_WRAPPER TRUE)
  endif()
  find_package(MPI REQUIRED COMPONENTS CXX)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

  pkg_check_modules(FFTW3_MPI QUIET IMPORTED_TARGET fftw3_mpi)
  if(NOT FFTW3_MPI_FOUND)
    find_path(FFTW3_MPI_INCLUDE_DIR NAMES fftw3-mpi.h)
    find_library(FFTW3_MPI_LIBRARY NAMES fftw3_mpi)
    if(NOT FFTW3_MPI_INCLUDE_DIR OR NOT FFTW3_MPI_LIBRARY)
      message(FATAL_ERROR "FFTW3 MPI not found. Install FFTW MPI development package (e.g. libfftw3-mpi-dev) or provide a fftw3_mpi pkg-config file.")
    endif()
  endif()

  add_executable(misfit_sim misfit_uniform_time_seed_fftw_early.cpp)
  target_compile_features(misfit_sim PRIVATE cxx_std_17)
  target_link_libraries(misfit_sim PRIVATE MPI::MPI_CXX PkgConfig::FFTW3 m)
  if(FFTW3_MPI_FOUND)
    target_link_libraries(misfit_sim PRIVATE PkgConfig::FFTW3_MPI)
  else()
    target_include_directories(misfit_sim PRIVATE "${FFTW3_MPI_INCLUDE_DIR}")
    target_link_libraries(misfit_sim PRIVATE "${FFTW3_MPI_LIBRARY}")
  endif()

  install(TARGETS misfit_sim RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()

if(BUILD_GUI)
  find_package(Qt6 6.4 REQUIRED COMPONENTS Quick Qml QuickControls2)

  qt_standard_project_setup()

  qt_add_executable(qt_gui_client
    src/gui/main.cpp
    src/gui/ExperimentController.cpp
    src/gui/RunQueueModel.cpp
  )

  target_include_directories(qt_gui_client PRIVATE src/gui)

  qt_add_qml_module(qt_gui_client
    URI ExperimentGui
    VERSION 1.0
    QML_FILES
      src/gui/qml/Main.qml
  )

  target_link_libraries(qt_gui_client PRIVATE Qt6::Quick Qt6::Qml Qt6::QuickControls2)

  if(WIN32)
    set_target_properties(qt_gui_client PROPERTIES WIN32_EXECUTABLE TRUE)
  endif()
  if(APPLE)
    set_target_properties(qt_gui_client PROPERTIES MACOSX_BUNDLE TRUE)
  endif()

  install(TARGETS qt_gui_client RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
  if(UNIX AND NOT APPLE)
    install(FILES packaging/linux/qt-experiment-runner.desktop DESTINATION "${CMAKE_INSTALL_DATADIR}/applications")
  endif()
endif()

set(CPACK_PACKAGE_NAME "qt-experiment-runner")
set(CPACK_PACKAGE_VENDOR "QtExperimentRunner")
set(CPACK_PACKAGE_CONTACT "maintainer@example.com")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Qt GUI to configure and batch-run misfit_sim experiments")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")

set(CPACK_GENERATOR "TGZ;DEB")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "${CPACK_PACKAGE_VENDOR}")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_PACKAGE_DEPENDS
  "openmpi-bin, libfftw3-mpi3, libfftw3-3, qml6-module-qtquick, qml6-module-qtquick-controls, qml6-module-qtquick-layouts, qml6-module-qtquick-dialogs"
)

include(CPack)
