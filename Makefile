LDFLAGS=-lgsl -lm

oneloopcalc: oneloopcalc.o integrator.o hardfactors_position.o hardfactors_momentum.o dss_pinlo.o mstwpdf.o cubature.o gluondist.o integrationcontext.o integrationtype.o libinterp2d.a
	g++ $(CPPFLAGS) -o oneloopcalc oneloopcalc.o dss_pinlo.o mstwpdf.o cubature.o gluondist.o integrationcontext.o integrationtype.o hardfactors_position.o hardfactors_momentum.o integrator.o libinterp2d.a $(LDFLAGS)
gbwoneloopcalc: gbwoneloopcalc.cpp dss_pinlo.cpp mstwpdf.cc cubature.c libinterp2d.a
	g++ $(CPPFLAGS) -o gbwoneloopcalc gbwoneloopcalc.cpp dss_pinlo.cpp mstwpdf.cc cubature.c libinterp2d.a -lgsl -lm
mvgluondist: gluondist.cpp gluondist.h interp2d.h
	g++ $(CPPFLAGS) -DGLUON_DIST_DRIVER -o mvgluondist gluondist.cpp libinterp2d.a -lgsl -lm

.PHONY: clean

clean:
	rm *.o

dss_pinlo.o: dss_pinlo.cpp dss_pinlo.h interp2d.h
dss_pinlo_test.o: dss_pinlo_test.cpp dss_pinlo.h interp2d.h
gbwoneloopcalc.o: gbwoneloopcalc.cpp cubature.h mstwpdf.h dss_pinlo.h \
 interp2d.h
gluondist.o: gluondist.cpp gluondist.h
hardfactors_momentum.o: hardfactors_momentum.cpp hardfactor.h \
 integrationcontext.h context.h mstwpdf.h dss_pinlo.h interp2d.h \
 coupling.h gluondist.h integrationtype.h hardfactors_position.h
hardfactors_position.o: hardfactors_position.cpp hardfactor.h \
 integrationcontext.h context.h mstwpdf.h dss_pinlo.h interp2d.h \
 coupling.h gluondist.h integrationtype.h hardfactors_position.h
integrationcontext.o: integrationcontext.cpp mstwpdf.h dss_pinlo.h \
 interp2d.h integrationcontext.h context.h coupling.h gluondist.h
integrationtype.o: integrationtype.cpp integrationtype.h \
 integrationcontext.h context.h mstwpdf.h dss_pinlo.h interp2d.h \
 coupling.h gluondist.h
integrator.o: integrator.cpp integrator.h context.h mstwpdf.h dss_pinlo.h \
 interp2d.h coupling.h gluondist.h integrationcontext.h integrationtype.h \
 hardfactor.h
oneloopcalc.o: oneloopcalc.cpp cubature.h mstwpdf.h dss_pinlo.h \
 interp2d.h coupling.h gluondist.h context.h integrationcontext.h \
 hardfactors_position.h hardfactors_momentum.h hardfactor.h integrationtype.h \
 integrator.h
