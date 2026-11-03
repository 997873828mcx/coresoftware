AUTOMAKE_OPTIONS = foreign subdir-objects

EXTRA_DIST = \
  tpc_frame_geometry.csv

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

noinst_HEADERS = PHGarfieldRootOutput.h

pkginclude_HEADERS = \
  PHGarfield.h \
  PHGarfieldForwardMap.h \
  PHGarfieldAdaptiveMap.h \
  PHGarfieldCompactMap.h \
  PHGarfieldRossegger.h

lib_LTLIBRARIES = \
  libPHGarfield.la

libPHGarfield_la_LIBADD = \
  -lffamodules \
  -lffarawobjects \
  -lcdbobjects \
  -lg4detectors \
  -lEvent \
  -lphool \
  -lphfield \
  -lGarfield \
  -lSubsysReco
libPHGarfield_la_LIBADD += -lcrypto

libPHGarfield_la_SOURCES = \
  PHGarfield.cc \
  PHGarfieldForwardMap.cc \
  PHGarfieldAdaptiveMap.cc \
  PHGarfieldAdaptiveBuilder.cc \
  PHGarfieldAdaptiveGenerator.cc \
  PHGarfieldCompactMap.cc \
  PHGarfieldRossegger.cc

bin_PROGRAMS = \
  GasModel \
  MergeGasFiles \
  PHGarfieldForwardMapGenerator \
  PHGarfieldCompactMapConverter

PHGarfieldCompactMapConverter_SOURCES = PHGarfieldCompactMapConverter.cc
PHGarfieldCompactMapConverter_LDADD = libPHGarfield.la

PHGarfieldForwardMapGenerator_SOURCES = \
  PHGarfieldForwardMapGenerator.cc

PHGarfieldForwardMapGenerator_LDADD = \
  libPHGarfield.la

MergeGasFiles_SOURCES = MergeGasFiles.cc

MergeGasFiles_LDADD = \
  -lGarfield \
  -lphool

GasModel_SOURCES = GasModel.cc

GasModel_LDADD = \
  -lGarfield

################################################
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

check_PROGRAMS = \
  test_forward_map_formula \
  test_adaptive_map \
  test_compact_map

TESTS = \
  test_forward_map_formula \
  test_adaptive_map \
  test_compact_map

test_forward_map_formula_SOURCES = \
  tests/test_forward_map_formula.cc

test_forward_map_formula_LDADD = \
  libPHGarfield.la

test_adaptive_map_SOURCES = tests/test_adaptive_map.cc
test_adaptive_map_LDADD = libPHGarfield.la
test_compact_map_SOURCES = tests/test_compact_map.cc
test_compact_map_LDADD = libPHGarfield.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libPHGarfield.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
