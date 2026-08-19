# ---------------------------------------------------------------------------
# packages/restG4 - optional Geant4 application package
# ---------------------------------------------------------------------------

# Si se compila este paquete de forma 100% aislada e independiente
if(CMAKE_PROJECT_NAME STREQUAL "restG4")
  cmake_minimum_required(VERSION 3.18)
  project(restG4 LANGUAGES CXX)

  set(ROOT_REQUIRED_LIBRARIES Core RIO Tree Graf Hist Physics Geom)
  find_package(ROOT 6.24 REQUIRED COMPONENTS ${ROOT_REQUIRED_LIBRARIES})
  include(${ROOT_USE_FILE})

  set(ROOT_TARGET_LIBRARIES "")
  foreach(_root_lib IN LISTS ROOT_REQUIRED_LIBRARIES)
    list(APPEND ROOT_TARGET_LIBRARIES ROOT::${_root_lib})
  endforeach()

  find_package(yaml-cpp REQUIRED)
  
  # Si es aislado, definimos los nombres canónicos que faltan del padre
  set(GEANT4_LIB RestGeant4) 
endif()

# Requisito obligatorio del paquete: Motor Geant4
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

if(Geant4_VERSION VERSION_LESS 11.0.0)
  add_compile_definitions(GEANT4_VERSION_LESS_11_0_0)
endif()

if(NOT EXISTS "${Geant4_INCLUDE_DIRS}/G4RunManagerFactory.hh")
  add_compile_definitions(GEANT4_WITHOUT_G4RunManagerFactory)
endif()

# Rastreo de fuentes internas de la aplicación
file(GLOB_RECURSE RESTG4_SOURCES src/*.cxx)
file(GLOB_RECURSE RESTG4_HEADERS include/*.h)

add_library(RestG4 SHARED ${RESTG4_SOURCES})

target_include_directories(
  RestG4
  SYSTEM
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/libraries/geant4/inc>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/source/core/inc>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/source/tools/inc>
    $<INSTALL_INTERFACE:include>)

# Enlace de dependencias de la aplicación (Core + Geant4_LIB unificados)
target_link_libraries(
  RestG4
  PUBLIC
    RestCore
    RestTools
    ${GEANT4_LIB}              # Ahora siempre vale "RestGeant4", garantizando el enlace
    ${Geant4_LIBRARIES}        # Librerías nativas de Geant4
    ${ROOT_TARGET_LIBRARIES}   # Librerías modernas ROOT::Core, ROOT::Tree, etc.
    yaml-cpp
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/ DESTINATION examples/restG4)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mac DESTINATION macros)
install(TARGETS RestG4 DESTINATION lib)
install(FILES ${RESTG4_HEADERS} DESTINATION include)

# --- Generación del Binario Ejecutable Final ---
add_executable(restG4 main.cxx)
target_link_libraries(restG4 PRIVATE RestG4)
install(TARGETS restG4 RUNTIME DESTINATION bin)

