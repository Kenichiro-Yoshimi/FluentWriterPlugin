cmake_minimum_required(VERSION 2.8)

# create a plugin that adds writers
# for Fluent files to the ParaView GUI
# it is added in the file dialog when doing saves.

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCUDE_DIRS}
  )
  INCLUDE_DIRECTORIES(SYSTEM
    ${MPI_INCLUDE_PATH}
  )

vtk_module_load(vtkzlib)

INCLUDE_DIRECTORIES(${vtkzlib_INCLUDE_DIRS})

ELSE()
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF ()

ADD_PARAVIEW_PLUGIN(FluentWriter "0.1"
  SERVER_MANAGER_SOURCES vtkFluentWriter.cxx
  SERVER_MANAGER_XML     FluentWriter.xml
  REQUIRED_ON_SERVER)

target_link_libraries(FluentWriter LINK_PRIVATE ${vtkzlib_LIBRARIES})

# Add testing if necessary
if (BUILD_TESTING)
  add_subdirectory(Testing)
endif()
