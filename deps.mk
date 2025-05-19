obj/autocalc.o : src/autocalc.f90 obj/mod_geninfo.o obj/mod_intrainfo.o
obj/becke.o : src/becke.f90
obj/c1hole.o : src/c1hole.f90 obj/mod_wfxinfo.o obj/mod_geninfo.o obj/mod_fractions.o
obj/centercalc.o : src/centercalc.f90 obj/mod_geninfo.o obj/mod_intrainfo.o
obj/derive.o : src/derive.f90 obj/mod_geninfo.o obj/mod_wfxinfo.o
obj/functions.o : src/functions.f90 obj/mod_geninfo.o obj/mod_loginfo.o obj/mod_wfxinfo.o
obj/gauss_hermite.o : src/gauss_hermite.f90 obj/mod_quadratures.o
obj/gaussleg.o : src/gaussleg.f90
obj/gen_cube.o : src/gen_cube.f90 obj/mod_geninfo.o obj/mod_cubeinfo.o
obj/gridpoints.o : src/gridpoints.f90 obj/mod_quadratures.o
obj/gridpoints2.o : src/gridpoints2.f90 obj/mod_quadratures.o
obj/gridpoints3.o : src/gridpoints3.f90 obj/mod_quadratures.o
obj/intracule.o : src/intracule.f90 obj/mod_geninfo.o obj/mod_intrastuff.o obj/mod_intrainfo.o obj/mod_quadratures.o
obj/mod_cubeinfo.o : src/mod_cubeinfo.f90
obj/mod_fractions.o : src/mod_fractions.f90
obj/mod_geninfo.o : src/mod_geninfo.f90
obj/mod_inputdat.o : src/mod_inputdat.f90
obj/mod_intrainfo.o : src/mod_intrainfo.f90
obj/mod_intrastuff.o : src/mod_intrastuff.f90 obj/mod_geninfo.o obj/mod_quadratures.o
obj/mod_located.o : src/mod_located.f90
obj/mod_loginfo.o : src/mod_loginfo.f90
obj/mod_quadratures.o : src/mod_quadratures.f90
obj/mod_radis.o : src/mod_radis.f90
obj/mod_wfxinfo.o : src/mod_wfxinfo.f90
obj/pdint.o : src/pdint.f90 obj/mod_geninfo.o obj/mod_intrainfo.o
obj/read_files.o : src/read_files.f90 obj/mod_wfxinfo.o obj/mod_geninfo.o obj/mod_located.o obj/mod_loginfo.o
obj/roda.o : src/roda.f90 obj/mod_inputdat.o obj/mod_cubeinfo.o obj/mod_located.o obj/mod_intrainfo.o obj/mod_geninfo.o obj/mod_radis.o
