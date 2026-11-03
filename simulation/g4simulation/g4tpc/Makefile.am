##############################################
# please add new classes/includes/everything else in alphabetical order

AUTOMAKE_OPTIONS = foreign subdir-objects

AM_CXXFLAGS = `geant4-config --cflags`

lib_LTLIBRARIES = \
  libg4tpc.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

libg4tpc_la_LIBADD = \
  -lphool \
  -lpdbcalBase \
  -lcdbobjects \
  -lffamodules \
  -lg4detectors \
  -lg4tracking_io \
  -lPHGarfield \
  -lphg4hit \
  -lphparameter \
  -ltrack \
  -ltrackbase_historic_io \
  -ltpc_io

pkginclude_HEADERS = \
  PHG4TpcCentralMembrane.h \
  TpcClusterBuilder.h \
  PHG4TpcDigitizer.h \
  PHG4TpcDirectLaser.h \
  PHG4TpcDistortion.h \
  PHG4TpcElectronDrift.h \
  PHG4TpcEndCapSubsystem.h \
  PHG4TpcPadBaselineShift.h \
  PHG4TpcPadPlane.h \
  PHG4TpcPadPlaneReadout.h \
  PHG4TpcSubsystem.h \
  TpcChannelMask.h \
  TpcSampaResponse.h

noinst_HEADERS = \
  TpcPrimaryIonizationModel.h \
  TpcRootQaTree.h

libg4tpc_la_SOURCES = \
  PHG4TpcCentralMembrane.cc \
  TpcClusterBuilder.cc \
  PHG4TpcDetector.cc \
  PHG4TpcDigitizer.cc \
  PHG4TpcDirectLaser.cc \
  PHG4TpcDisplayAction.cc \
  PHG4TpcDistortion.cc \
  PHG4TpcElectronDrift.cc \
  PHG4TpcEndCapDetector.cc \
  PHG4TpcEndCapDisplayAction.cc \
  PHG4TpcEndCapSteppingAction.cc \
  PHG4TpcEndCapSubsystem.cc \
  PHG4TpcPadBaselineShift.cc \
  PHG4TpcPadPlane.cc \
  PHG4TpcPadPlaneReadout.cc \
  PHG4TpcSteppingAction.cc \
  PHG4TpcSubsystem.cc \
  TpcChannelMask.cc \
  TpcPrimaryIonizationModel.cc

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

check_PROGRAMS = \
  test_PHG4TpcDigitizerMask \
  test_PHG4TpcDigitizerPedestal \
  test_PHG4TpcGarfieldLut \
  test_PHG4TpcPadPlaneReadoutSerf \
  test_PHG4TpcSampaTimeBounds \
  test_PHG4TpcTransverseDiffusion \
  test_TpcChannelMask \
  test_TpcRootQaTree \
  test_TpcSampaResponse

TESTS = \
  test_PHG4TpcDigitizerMask \
  test_PHG4TpcDigitizerPedestal \
  test_PHG4TpcGarfieldLut \
  test_PHG4TpcPadPlaneReadoutSerf \
  test_PHG4TpcSampaTimeBounds \
  test_PHG4TpcTransverseDiffusion \
  test_TpcChannelMask \
  test_TpcRootQaTree \
  test_TpcSampaResponse

test_PHG4TpcDigitizerMask_SOURCES = \
  tests/test_PHG4TpcDigitizerMask.cc

test_PHG4TpcDigitizerMask_LDADD = \
  libg4tpc.la \
  -lcdbobjects \
  -lg4detectors_io \
  -lphool \
  -ltrack_io

test_PHG4TpcDigitizerPedestal_SOURCES = \
  tests/test_PHG4TpcDigitizerPedestal.cc

test_PHG4TpcDigitizerPedestal_CPPFLAGS = $(AM_CPPFLAGS) -I$(srcdir)

test_PHG4TpcDigitizerPedestal_LDADD = \
  libg4tpc.la \
  -lg4detectors_io \
  -lphparameter_io \
  -lphool \
  -ltrack_io \
  -lpdbcalBase

test_PHG4TpcGarfieldLut_SOURCES = \
  tests/test_PHG4TpcGarfieldLut.cc

test_PHG4TpcGarfieldLut_LDADD = \
  libg4tpc.la \
  -lPHGarfield

test_PHG4TpcPadPlaneReadoutSerf_SOURCES = \
  tests/test_PHG4TpcPadPlaneReadoutSerf.cc

test_PHG4TpcPadPlaneReadoutSerf_LDADD = \
  libg4tpc.la \
  -lg4detectors_io \
  -lphg4hit \
  -lphool \
  -ltrack_io \
  -lgsl

test_TpcChannelMask_SOURCES = \
  tests/test_TpcChannelMask.cc

test_TpcChannelMask_LDADD = \
  libg4tpc.la \
  -ltrack_io

test_PHG4TpcSampaTimeBounds_SOURCES = \
  tests/test_PHG4TpcSampaTimeBounds.cc

test_PHG4TpcSampaTimeBounds_LDADD = \
  libg4tpc.la \
  -lg4detectors_io \
  -lphool

test_PHG4TpcTransverseDiffusion_SOURCES = \
  tests/test_PHG4TpcTransverseDiffusion.cc

test_PHG4TpcTransverseDiffusion_CPPFLAGS = $(AM_CPPFLAGS) -I$(srcdir)

test_PHG4TpcTransverseDiffusion_LDADD = \
  libg4tpc.la \
  -lg4detectors_io \
  -lphg4hit \
  -lphparameter \
  -lphparameter_io \
  -lphool \
  -ltrack_io \
  -lpdbcalBase

test_TpcRootQaTree_SOURCES = \
  tests/test_TpcRootQaTree.cc

test_TpcRootQaTree_LDADD =

test_TpcSampaResponse_SOURCES = \
  tests/test_TpcSampaResponse.cc

test_TpcSampaResponse_CPPFLAGS = $(AM_CPPFLAGS) -I$(srcdir)

test_TpcSampaResponse_LDADD =

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libg4tpc.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
