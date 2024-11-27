

find_program(TETGEN_EXE tetgen NO_CACHE)
find_program(TRIANGLE_EXE triangle NO_CACHE)

if((NOT TETGEN_EXE) OR (NOT TRIANGLE_EXE))
  set(TRIFS OFF)
  if(NOT TETGEN_EXE)
    message(WARNING "tetgen not found: triFS.x utility will not be compiled!")
  endif()
  if(NOT TRIANGLE_EXE)
    message(WARNING "triangle not found: triFS.x utility will not be compiled!")
  endif()
else()
  set(TRIFS ON)
endif()
