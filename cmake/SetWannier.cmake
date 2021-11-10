#############################################################################
#
#############################################################################

unset(W_INCLUDE_DIRS CACHE)
find_path( W_INCLUDE_DIRS w90_constants.mod PATHS ${W_HOME} PATH_SUFFIXES src src/obj src/objp NO_DEFAULT_PATH )


set(W_LIBS ${W_HOME}/libwannier.a)

add_library(W_lib INTERFACE IMPORTED)
target_include_directories( W_lib INTERFACE ${W_INCLUDE_DIRS} )
target_link_libraries( W_lib INTERFACE ${W_LIBS})
