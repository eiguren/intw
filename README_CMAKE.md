MPI bilatzeakon gauza ezberdinek bilatzen ditu:
        message("${MPI_FOUND}")
        message("${MPI_VERSION}")
        message("${MPI_Fortran_FOUND}")
        message("${MPI_Fortran_COMPILER}")
        message("${MPI_Fortran_COMPILE_OPTIONS}")
        message("${MPI_Fortran_COMPILE_DEFINITIONS}")
        message("${MPI_Fortran_INCLUDE_DIRS}")
        message("${MPI_Fortran_LINK_FLAGS}")
        message("${MPI_Fortran_LIBRARIES}")
        message("${MPI_Fortran_COMPILER}")
        message("${MPI_Fortran_COMPILER_FLAGS}")
        message("${MPI_COMPILER_FLAGS}")
        message("${MPI_Fortran_INCLUDE_PATH}")

Oain interesatzen zaizkitenak hauek die: MPI_Fortran_LINK_FLAGS eta MPI_Fortran_LIBRARIES

MPI_Fortran_LINK_FLAGS: Hemen zuzenen gehitzen dio -rpath /liburutegin_direkzioa
MPI_Fortran_LIBRARIES: Hemen liburutegik bakarrik dare, path absolutoa erabilita

Ordun nere ustez hurrengoa gertatzen da:
link_libraries(${MPI_Fortran_LIBRARIES}) erabiltzeakon, -Wl,-rpath,/liburutegin_direkzioa gehitzen do, baño bakarrik liburutegin path absolutoa gehitzen do.
Aldiz, link_libraries(${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LINK_FLAGS}) edo link_libraries(MPI::MPI_Fortran) erabiltzeakon, MPI_Fortran_LINK_FLAGS-en ia agertzen da -rpath espezifikauta, eta hori re gitzen dio.

Azkenen gertatzen dana da, -rpath bi aldiz daola adierazita, batetik MPI_Fortran_LIBRARIES-eko liburutegin norabidekin, eta bestetik MPI_Fortran_LINK_FLAGS-en agertzen dana.


Nere ustez, .so liburutegi bat erabultzeakon, eta link_libraries erabiltzeakon, beti gehitzen do liburutegin ibilbidea rpath bitartez.

Proba:
Hau konprobatzeko .so liburutegi bat sortu, nik ezautzen deten ibilbide baten daona, eta link_libraries(/liburutegin_ibilbide_absolutoa) erabili.
Nere ustez liburutegie daon ibilbidea -rpath bitartez gehituko do.
   KONPROBAUTA: liburutegi bat gehitzeakon bere ibilbide absolutoa gehitzen do -rpath ekin

--------------------------------------------------

Beste asunto bat:

RUNPATH eta RPATH

Bi aldagaiek antzeko funtzioa due, baño inportantzi ezberdine. Programa bat ejekutatzeakon, ld-k orden hontan bilauko ditu liburutegik:
 1.RPATH-ekin adierazitako tokitan
 2.LD_LIBRARY_PATH-ekin adierazitako tokitan
 3.RUNPATH-ekin adierazitako tokitan

Baino, RPATH eta RUNPATH bik balio bat badue, RPATH ignorau eiten da, eta ordun lehenengo LD_LIBRARY_PATH beitzen da, ta gero RUNPATH.

Hau berrie dala uste det, ta oaingo linkerrak beti RUNPATH adierazten dola. Alare, --enable-new-dtags erabilita beti RUNPATH idatziko da, eta --disable-new-dtags erabilita RPATH.
