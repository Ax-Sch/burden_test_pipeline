# Install script for directory: /media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/metal

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/bin/metal" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/bin/metal")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/bin/metal"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/bin/metal")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/bin" TYPE EXECUTABLE FILES "/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/metal/metal")
  if(EXISTS "$ENV{DESTDIR}/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/bin/metal" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/bin/metal")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tools/METAL-2020-05-05/build/bin/metal")
    endif()
  endif()
endif()

