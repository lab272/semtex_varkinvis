set (semtex_src
  analysis.cpp
  assemblymap.cpp
  auxfield.cpp
  bcmgr.cpp
  boundary.cpp
  boundarysys.cpp
  condition.cpp
  data2df.cpp
  domain.cpp
  edge.cpp
  element.cpp
  family.cpp
  feml.cpp
  field.cpp
  geometry.cpp
  history.cpp
  integration.cpp
  matrix.cpp
  mesh.cpp
  message.cpp
#  misc.cpp
  numbersys.cpp
  particle.cpp
  statistics.cpp
  svv.cpp
)

add_library (src STATIC ${semtex_src})
target_include_directories (src PRIVATE ${CMAKE_SOURCE_DIR}/src)

