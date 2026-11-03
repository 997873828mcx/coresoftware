AUTOMAKE_OPTIONS = foreign
# List of shared libraries to produce
lib_LTLIBRARIES = \
  libtpctrackreco_io.la \
  libtpctrackreco.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -fopenmp

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  TpcHitsetPedestalCorrection.h \
  Tpc_ModuleTrackReco.h \
  Tpc_ModuleTrack.h \
  Tpc_ModuleTrackv1.h \
  Tpc_ModuleTrackContainer.h \
  Tpc_ModuleTrackContainerv1.h \
  Tpc_AssembledTrack.h \
  Tpc_AssembledTrackv1.h \
  Tpc_AssembledTrackContainer.h \
  Tpc_AssembledTrackContainerv1.h \
  Tpc_AssembledTrackReco.h \
  TpcCrossingDecision.h \
  TpcCrossingDecisionv1.h \
  TpcCrossingDecisionContainer.h \
  TpcCrossingDecisionContainerv1.h \
  TpcCrossingFinder.h \
  Tpc_PolyTrack.h \
  Tpc_PolyTrackv1.h \
  Tpc_PolyTrackContainer.h \
  Tpc_PolyTrackContainerv1.h \
  Tpc_PolyTrackVertex.h \
  Tpc_PolyTrackVertexv1.h \
  Tpc_PolyTrackVertexContainer.h \
  Tpc_PolyTrackVertexContainerv1.h \
  Tpc_PolyTrackVertexer.h \
  Tpc_PolyCluster.h \
  Tpc_PolyClusterv1.h \
  Tpc_PolyClusterContainer.h \
  Tpc_PolyClusterContainerv1.h \
  Tpc_PolyClusterizer.h \
  Tpc_PolyTrackReco.h \
  TpcTrackFit.h \
  TpcTrackHelixFitter.h \
  TpcTrackKalmanFitter.h \
  TpcTrackDiagnosticFitter.h \
  TpcDriftTimePolicy.h \
  TpcPolyTrackSeedConverter.h \
  TpcPolyClusterTrkrClusterConverter.h \
  Tpc_FittingTools.h \
  IdealPadMap.h


ROOTDICTS = \
  Tpc_AssembledTrack_Dict.cc \
  Tpc_AssembledTrackContainer_Dict.cc \
  Tpc_AssembledTrackContainerv1_Dict.cc \
  Tpc_AssembledTrackv1_Dict.cc \
  TpcCrossingDecision_Dict.cc \
  TpcCrossingDecisionv1_Dict.cc \
  TpcCrossingDecisionContainer_Dict.cc \
  TpcCrossingDecisionContainerv1_Dict.cc \
  Tpc_ModuleTrack_Dict.cc \
  Tpc_ModuleTrackContainer_Dict.cc \
  Tpc_ModuleTrackContainerv1_Dict.cc \
  Tpc_ModuleTrackv1_Dict.cc \
  Tpc_PolyCluster_Dict.cc \
  Tpc_PolyClusterContainer_Dict.cc \
  Tpc_PolyClusterContainerv1_Dict.cc \
  Tpc_PolyClusterv1_Dict.cc \
  Tpc_PolyTrack_Dict.cc \
  Tpc_PolyTrackContainer_Dict.cc \
  Tpc_PolyTrackContainerv1_Dict.cc \
  Tpc_PolyTrackVertex_Dict.cc \
  Tpc_PolyTrackVertexContainer_Dict.cc \
  Tpc_PolyTrackVertexContainerv1_Dict.cc \
  Tpc_PolyTrackVertexv1_Dict.cc \
  Tpc_PolyTrackv1_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libtpctrackreco_io_la_SOURCES = \
  Tpc_ModuleTrackv1.cc \
  Tpc_ModuleTrackContainerv1.cc \
  Tpc_AssembledTrackv1.cc \
  Tpc_AssembledTrackContainerv1.cc \
  TpcCrossingDecisionv1.cc \
  TpcCrossingDecisionContainerv1.cc \
  Tpc_PolyTrackv1.cc \
  Tpc_PolyTrackContainerv1.cc \
  Tpc_PolyTrackVertexv1.cc \
  Tpc_PolyTrackVertexContainerv1.cc \
  Tpc_PolyClusterv1.cc \
  Tpc_PolyClusterContainerv1.cc \
  $(ROOTDICTS)

libtpctrackreco_la_SOURCES = \
  TpcHitsetPedestalCorrection.cc \
  Tpc_ModuleTrackReco.cc \
  Tpc_AssembledTrackReco.cc \
  TpcCrossingFinder.cc \
  Tpc_PolyTrackVertexer.cc \
  Tpc_PolyClusterizer.cc \
  Tpc_PolyTrackReco.cc \
  TpcPolyTrackSeedConverter.cc \
  TpcPolyClusterTrkrClusterConverter.cc \
  TpcTrackHelixFitter.cc \
  TpcTrackKalmanFitter.cc \
  TpcTrackDiagnosticFitter.cc \
  Tpc_FittingTools.cc \
  IdealPadMap.cc

libtpctrackreco_io_la_LIBADD = \
  -lphool

libtpctrackreco_la_LIBADD = \
  libtpctrackreco_io.la \
  -lphool \
  -lpdbcalBase \
  -lfun4all \
  -lffamodules \
  -lcdbobjects \
  -lg4eval \
  -lg4eval_io \
  -lphg4hit \
  -lmbd_io \
  -lg4detectors_io \
  -lmvtx \
  -lmvtx_io \
  -lintt \
  -lintt_io \
  -lmicromegas_io \
  -lSubsysReco \
  -ltpc_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lphfield \
  -lPHGarfield \
  -lTpcConditions

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_tpctrackreco_io \
  testexternals_tpctrackreco


testexternals_tpctrackreco_SOURCES = testexternals.cc
testexternals_tpctrackreco_LDADD = libtpctrackreco.la

testexternals_tpctrackreco_io_SOURCES = testexternals.cc
testexternals_tpctrackreco_io_LDADD = libtpctrackreco_io.la

check_PROGRAMS = test_TpcTrackDiagnosticFitter test_TpcDriftTimePolicy test_TpcHitsetPedestalCorrection
TESTS = $(check_PROGRAMS)
test_TpcTrackDiagnosticFitter_SOURCES = tests/test_TpcTrackDiagnosticFitter.cc
test_TpcTrackDiagnosticFitter_CPPFLAGS = $(AM_CPPFLAGS) -I$(srcdir)
test_TpcTrackDiagnosticFitter_LDADD = libtpctrackreco.la
test_TpcDriftTimePolicy_SOURCES = tests/test_TpcDriftTimePolicy.cc
test_TpcDriftTimePolicy_CPPFLAGS = -I$(srcdir)
test_TpcHitsetPedestalCorrection_SOURCES = tests/test_TpcHitsetPedestalCorrection.cc
test_TpcHitsetPedestalCorrection_CPPFLAGS = $(AM_CPPFLAGS) -I$(srcdir)
test_TpcHitsetPedestalCorrection_LDFLAGS = $(AM_LDFLAGS) -L$(ROOTSYS)/lib
test_TpcHitsetPedestalCorrection_LDADD = libtpctrackreco.la -ltrack_io -lpdbcalBase -lphool -lCore -lRIO

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
