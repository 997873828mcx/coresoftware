// Exercise the production electron-drift sampler, not a duplicate of its
// formulas. Synthetic PAI seeds specify the electron count; a capture pad plane
// observes the endpoints before avalanche gain, SAMPA shaping or ADC cuts.
#include "PHG4TpcElectronDrift.h"
#include "PHG4TpcPadPlane.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <g4detectors/PHG4TpcGeomv1.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>
#include <trackbase/ActsGeometry.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace
{
  constexpr unsigned int ElectronCount = 100000;
  constexpr double MaxDriftLengthCm = 100.0;
  constexpr double MembraneHalfWidthCm = 0.25;
  constexpr double DriftVelocityCmPerNs = 0.01;
  constexpr double StartX = 50.0;
  constexpr double StartY = 8.0;
  constexpr double StartTimeNs = 100.0;
  constexpr double PhysicalDiffusion = 0.005313;

  void require(bool condition, const std::string& message)
  {
    if (!condition)
    {
      std::cerr << "TPC transverse diffusion regression: " << message << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  struct Moments
  {
    unsigned int count{0};
    double sumX{0.0};
    double sumY{0.0};
    double sumXX{0.0};
    double sumYY{0.0};
    double sumXY{0.0};
    double sumXXXX{0.0};
    double sumYYYY{0.0};
    double maxAbsDisplacement{0.0};

    void add(double x, double y)
    {
      ++count;
      sumX += x;
      sumY += y;
      sumXX += x * x;
      sumYY += y * y;
      sumXY += x * y;
      sumXXXX += x * x * x * x;
      sumYYYY += y * y * y * y;
      maxAbsDisplacement = std::max({maxAbsDisplacement, std::abs(x), std::abs(y)});
    }
  };

  class CapturePadPlane final : public PHG4TpcPadPlane
  {
   public:
    CapturePadPlane(unsigned int expectedSide, double expectedTime)
      : PHG4TpcPadPlane("DiffusionCapturePadPlane")
      , m_expectedSide(expectedSide)
      , m_expectedTime(expectedTime)
    {
    }

    int InitRun(PHCompositeNode*) override
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

    void SetDefaultParameters() override {}

    void MapToPadPlane(TpcClusterBuilder&, TrkrHitSetContainer*,
                       TrkrHitSetContainer*, TrkrHitTruthAssoc*,
                       double x, double y, double time, unsigned int side,
                       PHG4HitContainer::ConstIterator, TNtuple*, TNtuple*) override
    {
      require(std::isfinite(x) && std::isfinite(y), "nonfinite captured endpoint");
      require(side == m_expectedSide, "wrong TPC side reached the pad plane");
      require(std::abs(time - m_expectedTime) < 1.e-9,
              "transverse diffusion changed the nominal drift time");
      moments.add(x - StartX, y - StartY);
      // No detector hits are made: this test stops at the electron endpoint.
      // Consequently truth clustering requires no Acts detector surfaces.
    }

    Moments moments;

   private:
    unsigned int m_expectedSide;
    double m_expectedTime;
  };

  void checkCase(const std::string& name, double driftLength,
                 double physicalCoefficient, double addedSigma,
                 double minimumDriftLength, unsigned int side,
                 unsigned int seed)
  {
    auto top = std::make_unique<PHCompositeNode>("TOP");
    auto* dst = new PHCompositeNode("DST");
    auto* run = new PHCompositeNode("RUN");
    auto* par = new PHCompositeNode("PAR");
    top->addNode(dst);
    top->addNode(run);
    top->addNode(par);

    auto* geometry = new PHG4TpcGeomContainer;
    auto* layer = new PHG4TpcGeomv1;
    layer->set_layer(20);
    layer->set_radius(50.0);
    layer->set_thickness(1.0);
    layer->set_max_driftlength(MaxDriftLengthCm);
    layer->set_CM_halfwidth(MembraneHalfWidthCm);
    layer->set_drift_velocity_sim(DriftVelocityCmPerNs);
    layer->set_extended_readout_time(1000.0);
    require(geometry->AddLayerCellGeom(20, layer) == 0, name + ": geometry initialization failed");
    run->addNode(new PHIODataNode<PHObject>(geometry, "TPCGEOMCONTAINER", "PHObject"));

    // Only the envelope transform is used by drift. Initialize it explicitly:
    // an otherwise empty ActsGeometry has no valid default affine transform.
    auto* acts = new ActsGeometry;
    acts->set_tpc_world_envelope_transform(Acts::Transform3::Identity());
    dst->addNode(new PHDataNode<ActsGeometry>(acts, "ActsGeometry"));

    auto* tpcParNode = new PHCompositeNode("TPC");
    par->addNode(tpcParNode);
    auto* gasContainer = new PHParametersContainer("TPC");
    auto* gas = new PHParameters("TPC");
    gas->set_double_param("Ne_frac", 0.0);
    gas->set_double_param("Ar_frac", 0.75);
    gas->set_double_param("CF4_frac", 0.20);
    gas->set_double_param("N2_frac", 0.0);
    gas->set_double_param("isobutane_frac", 0.05);
    gas->set_int_param("use_primary_cluster_ionization", 0);
    gasContainer->AddPHParameters(0, gas);
    tpcParNode->addNode(new PHDataNode<PHParametersContainer>(gasContainer, "G4GEO_TPC"));

    // InitRun requires the ordinary G4 input node even when explicit PAI
    // ionization seeds supply the drift sources. Both are synthetic here.
    auto* ordinaryHits = new PHG4HitContainer("G4HIT_TPC");
    dst->addNode(new PHIODataNode<PHObject>(ordinaryHits, "G4HIT_TPC", "PHObject"));
    auto* seeds = new PHG4HitContainer("G4HIT_TPC_PAI_CLUSTER");
    dst->addNode(new PHIODataNode<PHObject>(seeds, "G4HIT_TPC_PAI_CLUSTER", "PHObject"));
    auto* hit = new PHG4Hitv1;
    const double absZ = MaxDriftLengthCm + MembraneHalfWidthCm - driftLength;
    for (unsigned int endpoint = 0; endpoint < 2; ++endpoint)
    {
      hit->set_x(endpoint, StartX);
      hit->set_y(endpoint, StartY);
      hit->set_z(endpoint, side == 0 ? -absZ : absZ);
      hit->set_t(endpoint, StartTimeNs);
    }
    hit->set_trkid(1);
    hit->set_layer(20);
    hit->set_edep(1.e-6);
    hit->set_eion(1.e-6);
    hit->set_index_i(ElectronCount);
    seeds->AddHit(20, hit);

    PHG4TpcElectronDrift drift("DiffusionTest_" + name);
    drift.Detector("TPC");
    drift.set_seed(seed);
    drift.set_do_ElectronDriftQAHistos(false);
    drift.set_use_pai_cluster_seeds(true);
    drift.set_double_param("diffusion_trans", physicalCoefficient);
    drift.set_double_param("added_smear_trans", addedSigma);
    drift.set_double_param("diffusion_long", 0.0);
    drift.set_double_param("added_smear_long", 0.0);
    drift.set_force_min_trans_drift_length(minimumDriftLength);
    auto* capture = new CapturePadPlane(side, StartTimeNs + driftLength / DriftVelocityCmPerNs);
    drift.registerPadPlane(capture);  // PHG4TpcElectronDrift owns the pad plane.

    require(drift.Init(top.get()) == Fun4AllReturnCodes::EVENT_OK, name + ": Init failed");
    require(drift.InitRun(top.get()) == Fun4AllReturnCodes::EVENT_OK, name + ": InitRun failed");
    require(drift.process_event(top.get()) == Fun4AllReturnCodes::EVENT_OK, name + ": process_event failed");

    const auto& result = capture->moments;
    require(result.count == ElectronCount, name + ": electron endpoints were dropped or duplicated");
    const double expectedVariance = physicalCoefficient * physicalCoefficient *
                                        std::max(driftLength, minimumDriftLength) +
                                    addedSigma * addedSigma;
    if (expectedVariance == 0.0)
    {
      require(result.maxAbsDisplacement == 0.0, name + ": zero diffusion displaced an electron");
      std::cout << name << ": " << result.count << " exactly unchanged endpoints" << std::endl;
      return;
    }

    const double meanX = result.sumX / result.count;
    const double meanY = result.sumY / result.count;
    const double varianceX = result.sumXX / result.count - meanX * meanX;
    const double varianceY = result.sumYY / result.count - meanY * meanY;
    const double covariance = result.sumXY / result.count - meanX * meanY;
    const double meanTolerance = 8.0 * std::sqrt(expectedVariance / result.count);
    require(std::abs(meanX) < meanTolerance && std::abs(meanY) < meanTolerance,
            name + ": nonzero transverse mean");
    // 2.5% of variance is about 5.6 statistical standard deviations for N=1e5.
    // The old sampler has half the required variance and fails decisively.
    require(std::abs(varianceX / expectedVariance - 1.0) < 0.025 &&
                std::abs(varianceY / expectedVariance - 1.0) < 0.025,
            name + ": per-axis width differs from the configured coefficient");
    require(std::abs(covariance / expectedVariance) < 0.025,
            name + ": x/y diffusion is correlated");
    const double expectedFourth = 3.0 * expectedVariance * expectedVariance;
    require(std::abs(result.sumXXXX / result.count / expectedFourth - 1.0) < 0.06 &&
                std::abs(result.sumYYYY / result.count / expectedFourth - 1.0) < 0.06,
            name + ": transverse fourth moments are inconsistent with Gaussian axes");

    std::cout << std::setprecision(8) << name << ": N=" << result.count
              << " sigma_expected=" << std::sqrt(expectedVariance)
              << " sigma_x=" << std::sqrt(varianceX)
              << " sigma_y=" << std::sqrt(varianceY)
              << " covariance/variance=" << covariance / expectedVariance << std::endl;
  }
}  // namespace

int main()
{
  checkCase("zero", 40.0, 0.0, 0.0, 0.0, 0, 11001);
  checkCase("physical_side0", 40.0, PhysicalDiffusion, 0.0, 0.0, 0, 11002);
  checkCase("physical_side1", 40.0, PhysicalDiffusion, 0.0, 0.0, 1, 11003);
  checkCase("short_drift", 10.0, PhysicalDiffusion, 0.0, 0.0, 0, 11004);
  checkCase("added_only", 40.0, 0.0, 0.025, 0.0, 1, 11005);
  checkCase("forced_minimum", 10.0, PhysicalDiffusion, 0.0, 80.0, 0, 11006);
  checkCase("physical_plus_added", 40.0, PhysicalDiffusion, 0.020, 0.0, 1, 11007);
  checkCase("forced_plus_added", 10.0, PhysicalDiffusion, 0.020, 80.0, 1, 11008);
  checkCase("minimum_below_actual", 80.0, PhysicalDiffusion, 0.0, 10.0, 0, 11009);
  std::cout << "TPC production transverse diffusion regression passed" << std::endl;
  return EXIT_SUCCESS;
}
