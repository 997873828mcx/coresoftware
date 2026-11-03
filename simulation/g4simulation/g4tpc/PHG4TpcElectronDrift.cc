// this is the new containers version
// it uses the same MapToPadPlane as the old containers version

#include "PHG4TpcElectronDrift.h"
#include "PHG4TpcDistortion.h"
#include "PHG4TpcPadPlane.h"  // for PHG4TpcPadPlane
#include "TpcClusterBuilder.h"
#include "TpcPrimaryIonizationModel.h"
#include "TpcRootQaTree.h"

#include <trackbase/ClusHitsVerbosev1.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>  // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssoc.h>  // for TrkrHitTruthA...
#include <trackbase/TrkrHitTruthAssocv1.h>
#include <trackbase/TrkrHitv2.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <g4tracking/TrkrTruthTrackContainerv1.h>
#include <g4tracking/TrkrTruthTrackv1.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Particlev3.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phparameter/PHParameterInterface.h>  // for PHParameterIn...
#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <phgarfield/PHGarfield.h>

#include <phfield/PHFieldConfig.h>
#include <phfield/PHFieldUtility.h>

#include <pdbcalbase/PdbParameterMapContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>  // for PHDataNode
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TMD5.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TTree.h>

// for event id tagging in debug prints
#include <ffaobjects/EventHeader.h>

#include <boost/format.hpp>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_alloc

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>    // for sqrt, abs, NAN
#include <cstdlib>  // for exit
#include <exception>
#include <format>
#include <iostream>
#include <map>  // for _Rb_tree_cons...
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>  // for pair
#include <vector>

namespace
{
  // Stores one track-layer truth intersection (treated as a "true cluster")
  // per event entry as a PHG4Hit.
  static constexpr const char *kTpcTruthIntersectionNodeName = "G4HIT_TPC_TRUECLUSTER";
  static constexpr std::size_t kMaxTruthIntersectionsPerTrackLayer = 3;

  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }

  bool fatal_garfield_failure(const PHGarfield::ForwardDriftStatus status)
  {
    using Status = PHGarfield::ForwardDriftStatus;
    return status == Status::NotRun ||
           status == Status::NotInitialized ||
           status == Status::InvalidInput ||
           status == Status::InvalidVelocity ||
           status == Status::WrongDirection ||
           status == Status::LeftMagneticFieldMap ||
           status == Status::MaxSteps;
  }
}  // namespace

PHG4TpcElectronDrift::PHG4TpcElectronDrift(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , temp_hitsetcontainer(new TrkrHitSetContainerv1)
  , single_hitsetcontainer(new TrkrHitSetContainerv1)
{
  InitializeParameters();
  RandomGenerator.reset(gsl_rng_alloc(gsl_rng_mt19937));
  set_seed(PHRandomSeed());
}

PHG4TpcElectronDrift::~PHG4TpcElectronDrift()
{
  delete driftXY;
  delete driftStepQA;
  delete ionizationEventQA;
  delete ionizationClusterQA;
  delete ionizationSegmentQA;
}

//_____________________________________________________________
int PHG4TpcElectronDrift::Init(PHCompositeNode *topNode)
{
  padplane->Init(topNode);
  event_num = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
int PHG4TpcElectronDrift::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  auto runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  auto parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  const std::string paramnodename = "G4CELLPARAM_" + detector;
  const std::string geonodename = "G4CELLPAR_" + detector;
  const std::string tpcgeonodename = "G4GEO_" + detector;
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    std::cout << Name() << " Could not locate G4HIT node " << hitnodename << std::endl;
    topNode->print();
    gSystem->Exit(1);
    exit(1);
  }
  // new containers
  hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    hitsetcontainer = new TrkrHitSetContainerv1;
    auto newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    DetNode->addNode(newNode);
  }

  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    hittruthassoc = new TrkrHitTruthAssocv1;
    auto newNode = new PHIODataNode<PHObject>(hittruthassoc, "TRKR_HITTRUTHASSOC", "PHObject");
    DetNode->addNode(newNode);
  }

  truthtracks = findNode::getClass<TrkrTruthTrackContainer>(topNode, "TRKR_TRUTHTRACKCONTAINER");
  if (!truthtracks)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    truthtracks = new TrkrTruthTrackContainerv1();
    auto newNode = new PHIODataNode<PHObject>(truthtracks, "TRKR_TRUTHTRACKCONTAINER", "PHObject");
    DetNode->addNode(newNode);
  }

  truthclustercontainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
  if (!truthclustercontainer)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    truthclustercontainer = new TrkrClusterContainerv4;
    auto newNode = new PHIODataNode<PHObject>(truthclustercontainer, "TRKR_TRUTHCLUSTERCONTAINER", "PHObject");
    DetNode->addNode(newNode);
  }

  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (Verbosity() > 0 && m_track_map)
  {
    std::cout << Name() << ": found SvtxTrackMap with " << m_track_map->size() << " entries" << std::endl;
  }

  seggeonodename = "TPCGEOMCONTAINER";  // + detector;
  seggeo = findNode::getClass<PHG4TpcGeomContainer>(topNode, seggeonodename);
  assert(seggeo);

  // The layer geometry carries the z geometry and timing used by the complete
  // simulation chain. Drift and readout must use the same values.
  auto *z_geometry = seggeo->GetLayerCellGeom(20);
  if (!z_geometry)
  {
    std::cerr << "PHG4TpcElectronDrift: missing TPC layer geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const double max_drift_length = z_geometry->get_max_driftlength();
  const double central_membrane_half_width = z_geometry->get_CM_halfwidth();
  drift_velocity = z_geometry->get_drift_velocity_sim();
  tpc_length = 2.0 * (max_drift_length + central_membrane_half_width);

  if (max_drift_length <= 0.0 ||
      central_membrane_half_width < 0.0 ||
      drift_velocity <= 0.0 ||
      !std::isfinite(tpc_length))
  {
    std::cerr << "PHG4TpcElectronDrift: invalid TPC drift geometry"
              << " (max drift length " << max_drift_length
              << ", central membrane half width " << central_membrane_half_width
              << ", drift velocity " << drift_velocity << ")" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_garfield && m_distortionMap)
  {
    std::cerr << Name()
              << ": PHGarfield forward drift and PHG4TpcDistortion cannot be "
                 "enabled together; that would apply the distortion twice."
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (m_distortionMap && !m_distortionMap->is_initialized())
  {
    std::cerr << Name()
              << ": PHG4TpcDistortion was not successfully initialized"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (!std::isfinite(m_garfield_max_nominal_time_mismatch) ||
      m_garfield_max_nominal_time_mismatch < 0.0)
  {
    std::cerr << Name()
              << ": invalid maximum PHGarfield nominal-time mismatch "
              << m_garfield_max_nominal_time_mismatch << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_garfield)
  {
    if (!std::isfinite(m_garfield_step_ns) || m_garfield_step_ns <= 0.0)
    {
      std::cerr << Name() << ": invalid PHGarfield forward-drift step "
                << m_garfield_step_ns << " ns" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    const double nominal_max_drift_time_ns =
        max_drift_length / drift_velocity;
    if (m_garfield_step_ns * PHGarfield::MaxForwardDriftSteps <
        nominal_max_drift_time_ns)
    {
      std::cerr << Name()
                << ": PHGarfield forward-drift step "
                << m_garfield_step_ns
                << " ns cannot cover the nominal maximum drift time "
                << nominal_max_drift_time_ns << " ns within "
                << PHGarfield::MaxForwardDriftSteps << " steps"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    // PHG4 drift positions use the local physical TPC-envelope frame. Keep
    // PHGarfield in that frame as well, and use the geometry placement only
    // when sampling the global magnetic-field map.
    m_garfield->SetCentralMembraneHalfWidth(central_membrane_half_width);
    m_garfield->MoveTpc(
        z_geometry->get_place_x(),
        z_geometry->get_place_y(),
        z_geometry->get_place_z());
    m_garfield->SetTpcRotation(
        z_geometry->get_rot_x(),
        z_geometry->get_rot_y(),
        z_geometry->get_rot_z());
    m_garfield->SetZeroMagneticField(zero_bfield);

    if (m_garfield->InitRun(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        !m_garfield->IsReady())
    {
      std::cerr << Name() << ": PHGarfield failed to initialize its field or gas"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    const bool has_nonzero_electric_correction_scale =
        m_garfield->GetSpaceChargeScaleSide0() != 0.0 ||
        m_garfield->GetSpaceChargeScaleSide1() != 0.0;
    if ((m_require_garfield_electric_field_map ||
         has_nonzero_electric_correction_scale) &&
        !m_garfield->HasElectricFieldCorrections())
    {
      std::cerr << Name()
                << ": PHGarfield forward drift requires a valid electric-field "
                   "correction map when it is required explicitly or either "
                   "correction scale is nonzero, but none was loaded"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  }

  m_tpc_layer_radii.clear();
  {
    const auto layer_range = seggeo->get_begin_end();
    m_tpc_layer_radii.reserve(std::distance(layer_range.first, layer_range.second));
    for (auto layeriter = layer_range.first; layeriter != layer_range.second; ++layeriter)
    {
      auto *layergeom = layeriter->second;
      if (!layergeom)
      {
        continue;
      }
      m_tpc_layer_radii.push_back({static_cast<unsigned int>(layergeom->get_layer()), layergeom->get_radius()});
    }
    std::sort(
        m_tpc_layer_radii.begin(),
        m_tpc_layer_radii.end(),
        [](const TpcLayerRadius &lhs, const TpcLayerRadius &rhs)
        { return lhs.layer < rhs.layer; });
  }

  m_truth_intersection_hits = findNode::getClass<PHG4HitContainer>(topNode, kTpcTruthIntersectionNodeName);
  if (!m_truth_intersection_hits)
  {
    m_truth_intersection_hits = new PHG4HitContainer(kTpcTruthIntersectionNodeName);
    auto *newNode = new PHIODataNode<PHObject>(m_truth_intersection_hits, kTpcTruthIntersectionNodeName, "PHObject");
    dstNode->addNode(newNode);
  }

  UpdateParametersWithMacro();

  m_density_enabled = false;
  m_density_layer = get_int_param("density_layer");
  m_density_bin_width_cm = get_double_param("density_bin_width_cm");
  m_density_window_cm = get_double_param("density_window_cm");
  m_uniform_density_test = (get_int_param("uniform_density_test") != 0);
  m_use_pai_cluster_seeds = (get_int_param("use_pai_cluster_seeds") != 0);
  fixed_primary_clusters_per_cm = get_double_param("fixed_primary_clusters_per_cm");
  if (fixed_primary_clusters_per_cm <= 0.0)
  {
    fixed_primary_clusters_per_cm = get_double_param("fixed_primary_electrons_per_cm");
  }
  if (m_density_layer >= 0 && m_density_bin_width_cm > 0.0 && m_density_window_cm > 0.0)
  {
    if (auto *layer_geom = seggeo->GetLayerCellGeom(m_density_layer))
    {
      m_density_layer_radius = layer_geom->get_radius();
      const double half_thickness = 0.5 * layer_geom->get_thickness();
      m_density_layer_rlow = m_density_layer_radius - half_thickness;
      m_density_layer_rhigh = m_density_layer_radius + half_thickness;
      m_density_radial_margin_cm = std::max(0.0, get_double_param("density_radial_margin_cm"));
      m_density_capture_rlow = m_density_layer_rlow - m_density_radial_margin_cm;
      m_density_capture_rhigh = m_density_layer_rhigh + m_density_radial_margin_cm;
      m_density_hist_half_range = (m_density_window_cm > 0.0) ? m_density_window_cm : (half_thickness + m_density_radial_margin_cm);
      const double low_edge = -m_density_hist_half_range;
      double high_edge = m_density_hist_half_range;
      int nbins = static_cast<int>(std::ceil((high_edge - low_edge) / m_density_bin_width_cm));
      if (nbins < 1)
      {
        nbins = 1;
      }
      high_edge = low_edge + nbins * m_density_bin_width_cm;
      if (electronDensityProfile)
      {
        delete electronDensityProfile;
        electronDensityProfile = nullptr;
      }
      const std::string hname = Name() + std::string("_ElectronDensityLayer") + std::to_string(m_density_layer);
      const std::string htitle = "Ionization start vs. layer " + std::to_string(m_density_layer) + ";r - R_{layer} (cm);Electrons";
      electronDensityProfile = new TH1F(hname.c_str(), htitle.c_str(), nbins, low_edge, high_edge);
      electronDensityProfile->SetDirectory(nullptr);
      m_density_enabled = true;
    }
    else if (Verbosity() > 0)
    {
      std::cout << Name() << ": density_layer " << m_density_layer << " not present in geometry; disabling density histogram" << std::endl;
    }
  }
  m_avg_x_enabled = false;
  m_avg_layer = get_int_param("average_x_layer");
  if (m_avg_layer >= 0)
  {
    if (auto *layer_geom = seggeo->GetLayerCellGeom(m_avg_layer))
    {
      m_avg_layer_radius = layer_geom->get_radius();
      const double half_thickness = 0.5 * layer_geom->get_thickness();
      m_avg_layer_rlow = m_avg_layer_radius - half_thickness;
      m_avg_layer_rhigh = m_avg_layer_radius + half_thickness;
      m_avg_x_enabled = true;
      if (Verbosity() > 0)
      {
        const double full_thickness = layer_geom->get_thickness();
        std::cout << Name() << ": averaging x for layer " << m_avg_layer
                  << " at R=" << m_avg_layer_radius << " cm"
                  << " thickness=" << full_thickness << " cm (rlow=" << m_avg_layer_rlow
                  << ", rhigh=" << m_avg_layer_rhigh << ")" << std::endl;
      }
    }
    else if (Verbosity() > 0)
    {
      std::cout << Name() << ": average_x_layer " << m_avg_layer << " not present in geometry; disabling average-x residual histogram" << std::endl;
    }
  }
  PHNodeIterator runIter(runNode);
  auto RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(detector);
    runNode->addNode(RunDetNode);
  }
  SaveToNodeTree(RunDetNode, paramnodename);

  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  auto ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", detector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(detector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, geonodename);

  // find Tpc Geo
  PHNodeIterator tpcpariter(ParDetNode);
  auto tpcparams = findNode::getClass<PHParametersContainer>(ParDetNode, tpcgeonodename);
  if (!tpcparams)
  {
    const std::string runparamname = "G4GEOPARAM_" + detector;
    auto tpcpdbparams = findNode::getClass<PdbParameterMapContainer>(RunDetNode, runparamname);
    if (tpcpdbparams)
    {
      tpcparams = new PHParametersContainer(detector);
      if (Verbosity())
      {
        tpcpdbparams->print();
      }
      tpcparams->CreateAndFillFrom(tpcpdbparams, detector);
      ParDetNode->addNode(new PHDataNode<PHParametersContainer>(tpcparams, tpcgeonodename));
    }
    else
    {
      std::cout << "PHG4TpcElectronDrift::InitRun - failed to find " << runparamname << " in order to initialize " << tpcgeonodename << ". Aborting run ..." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  assert(tpcparams);

  if (Verbosity())
  {
    tpcparams->Print();
  }
  const PHParameters *tpcparam = tpcparams->GetParameters(0);
  assert(tpcparam);

  if (m_garfield)
  {
    const double gas_inner_radius =
        tpcparam->get_double_param("gas_inner_radius");
    const double gas_outer_radius =
        tpcparam->get_double_param("gas_outer_radius");
    const double configured_tpc_length =
        tpcparam->get_double_param("tpc_length");
    const double central_membrane_width =
        tpcparam->get_double_param("window_thickness");
    if (!std::isfinite(gas_inner_radius) ||
        !std::isfinite(gas_outer_radius) ||
        !std::isfinite(configured_tpc_length) ||
        !std::isfinite(central_membrane_width) ||
        gas_inner_radius < 0.0 ||
        gas_outer_radius <= gas_inner_radius ||
        std::abs(configured_tpc_length - tpc_length) > 1.0e-6 ||
        std::abs(0.5 * central_membrane_width -
                 central_membrane_half_width) >
            1.0e-6)
    {
      std::cerr << Name()
                << ": inconsistent TPC geometry for PHGarfield: gas radii=("
                << gas_inner_radius << ", " << gas_outer_radius
                << ") cm, tpc_length=" << configured_tpc_length
                << " cm (drift geometry gives " << tpc_length
                << " cm), window_thickness=" << central_membrane_width
                << " cm (CM half-width=" << central_membrane_half_width
                << " cm)"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    m_garfield->SetDriftVolumeRadii(gas_inner_radius, gas_outer_radius);

    if (m_garfield->HasElectricFieldCorrections() &&
        !m_garfield->ElectricFieldMapCovers(
            gas_inner_radius,
            gas_outer_radius,
            max_drift_length))
    {
      std::cerr << Name()
                << ": PHGarfield electric-field map does not cover the full "
                   "TPC gas radius and drift length"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (m_garfield_full_transport)
    {
      const double reference_radius =
          0.5 * (gas_inner_radius + gas_outer_radius);
      const double nominal_time_ns = max_drift_length / drift_velocity;
      double maximum_relative_time_mismatch = 0.0;
      for (unsigned int side = 0; side < 2; ++side)
      {
        PHGarfield::ForwardDriftResult reference_result;
        const double reference_z =
            side == 0
                ? -central_membrane_half_width
                : central_membrane_half_width;
        if (!m_garfield->ForwardDrift(
                reference_radius,
                0.0,
                reference_z,
                tpc_length / 2.0,
                m_garfield_step_ns,
                reference_result,
                false) ||
            !reference_result.reached_readout)
        {
          std::cerr
              << Name()
              << ": failed PHGarfield scale-zero timing calibration drift "
              << "on side " << side << ": "
              << PHGarfield::ForwardDriftStatusName(
                     reference_result.status)
              << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }

        const double relative_time_mismatch =
            std::abs(reference_result.drift_time_ns - nominal_time_ns) /
            nominal_time_ns;
        maximum_relative_time_mismatch =
            std::max(
                maximum_relative_time_mismatch,
                relative_time_mismatch);
        std::cout << Name()
                  << ": PHGarfield full-transport scale-zero time side "
                  << side << " = " << reference_result.drift_time_ns
                  << " ns; PHG4 nominal time " << nominal_time_ns
                  << " ns; fractional mismatch "
                  << relative_time_mismatch << std::endl;
      }

      if (maximum_relative_time_mismatch >
          m_garfield_max_nominal_time_mismatch)
      {
        std::cerr << Name()
                  << ": PHGarfield full-transport timing mismatch exceeds "
                  << m_garfield_max_nominal_time_mismatch
                  << "; make the Garfield gas/E field and PHG4 drift velocity "
                     "consistent, or use distortion-only mode"
                  << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
  }

  if (m_distortionMap &&
      m_distortionMap->application_mode() ==
          PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport)
  {
    const auto& metadata = m_distortionMap->garfield_lut_metadata();
    const double gas_inner_radius =
        tpcparam->get_double_param("gas_inner_radius");
    const double gas_outer_radius =
        tpcparam->get_double_param("gas_outer_radius");
    const double configured_tpc_length =
        tpcparam->get_double_param("tpc_length");
    const double central_membrane_width =
        tpcparam->get_double_param("window_thickness");
    constexpr double geometry_tolerance_cm = 1.0e-6;
    const auto differs =
        [](const double first, const double second)
    {
      return !std::isfinite(first) || !std::isfinite(second) ||
             std::abs(first - second) > geometry_tolerance_cm;
    };

    if (!std::isfinite(gas_inner_radius) ||
        !std::isfinite(gas_outer_radius) ||
        !std::isfinite(configured_tpc_length) ||
        !std::isfinite(central_membrane_width) ||
        gas_inner_radius <= 0.0 ||
        gas_outer_radius <= gas_inner_radius ||
        differs(configured_tpc_length, tpc_length) ||
        differs(
            0.5 * central_membrane_width,
            central_membrane_half_width) ||
        differs(metadata.gas_inner_radius_cm, gas_inner_radius) ||
        differs(metadata.gas_outer_radius_cm, gas_outer_radius) ||
        differs(
            metadata.central_membrane_half_width_cm,
            central_membrane_half_width) ||
        differs(metadata.readout_abs_z_cm, tpc_length / 2.0) ||
        differs(metadata.tpc_place_x_cm, z_geometry->get_place_x()) ||
        differs(metadata.tpc_place_y_cm, z_geometry->get_place_y()) ||
        differs(metadata.tpc_place_z_cm, z_geometry->get_place_z()) ||
        differs(metadata.tpc_rotation_x_rad, z_geometry->get_rot_x()) ||
        differs(metadata.tpc_rotation_y_rad, z_geometry->get_rot_y()) ||
        differs(metadata.tpc_rotation_z_rad, z_geometry->get_rot_z()) ||
        metadata.zero_magnetic_field != zero_bfield)
    {
      std::cerr
          << Name()
          << ": Garfield LUT geometry does not match runtime TPC geometry"
          << " (map r=[" << metadata.gas_inner_radius_cm << ","
          << metadata.gas_outer_radius_cm << "], CM half-width="
          << metadata.central_membrane_half_width_cm
          << ", readout |z|=" << metadata.readout_abs_z_cm
          << "; runtime r=[" << gas_inner_radius << ","
          << gas_outer_radius << "], CM half-width="
          << central_membrane_half_width << ", readout |z|="
          << tpc_length / 2.0 << ")" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    const double velocity_mismatch =
        std::abs(
            metadata.nominal_phg4_drift_velocity_cm_per_ns -
            drift_velocity) /
        drift_velocity;
    if (!std::isfinite(velocity_mismatch) ||
        velocity_mismatch > m_garfield_max_nominal_time_mismatch)
    {
      std::cerr
          << Name()
          << ": Garfield LUT expects PHG4 drift velocity "
          << metadata.nominal_phg4_drift_velocity_cm_per_ns
          << " cm/ns, runtime uses " << drift_velocity
          << " cm/ns (fractional mismatch " << velocity_mismatch
          << ", allowed " << m_garfield_max_nominal_time_mismatch << ")"
          << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (metadata.scale_zero_maximum_time_mismatch >
        m_garfield_max_nominal_time_mismatch)
    {
      std::cerr
          << Name()
          << ": Garfield LUT measured scale-zero timing mismatch "
          << metadata.scale_zero_maximum_time_mismatch
          << ", exceeding the runtime allowance "
          << m_garfield_max_nominal_time_mismatch
          << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (!metadata.zero_magnetic_field)
    {
      const PHFieldConfig* field_config =
          PHFieldUtility::GetFieldConfigNode(nullptr, topNode);
      // Both readers consume the same Cartesian xyz payload. The interpolated
      // reader is the mapped Geant4/Kalman baseline; keep the payload and scale
      // checks below for either reader instead of rejecting that configuration.
      if (!field_config ||
          (field_config->get_field_config() != PHFieldConfig::Field3DCartesian &&
           field_config->get_field_config() != PHFieldConfig::FieldInterpolated) ||
          !std::isfinite(field_config->get_magfield_rescale()))
      {
        std::cerr
            << Name()
            << ": a nonzero-B Garfield LUT requires the active simulation "
               "field to use a Cartesian or interpolated Cartesian PHField map"
            << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      const double field_scale = field_config->get_magfield_rescale();
      const double scale_normalization =
          std::max(
              {1.0,
               std::abs(field_scale),
               std::abs(metadata.magnetic_field_scale)});
      if (std::abs(field_scale - metadata.magnetic_field_scale) >
          1.0e-6 * scale_normalization)
      {
        std::cerr
            << Name()
            << ": Garfield LUT magnetic-field scale "
            << metadata.magnetic_field_scale
            << " does not match the active simulation scale "
            << field_scale << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      std::unique_ptr<TMD5> active_field_md5(
          TMD5::FileChecksum(field_config->get_filename().c_str()));
      const std::string active_field_digest =
          active_field_md5 ? active_field_md5->AsString() : "";
      if (active_field_digest.empty() ||
          active_field_digest != metadata.magnetic_field_map_md5)
      {
        std::cerr
            << Name()
            << ": Garfield LUT magnetic-field payload does not match the "
               "active simulation field "
            << field_config->get_filename() << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
  }

  diffusion_long = get_double_param("diffusion_long");
  added_smear_sigma_long = get_double_param("added_smear_long");
  diffusion_trans = get_double_param("diffusion_trans");
  if (zero_bfield)
  {
    diffusion_trans *= zero_bfield_diffusion_factor;
  }
  added_smear_sigma_trans = get_double_param("added_smear_trans");
  force_min_trans_drift_length = get_double_param("force_min_trans_drift_length");
  if (force_min_trans_drift_length < 0.)
  {
    force_min_trans_drift_length = 0.;
  }

  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html

  double Ne_dEdx = 1.56;  // keV/cm
  double Ne_primary_clusters_per_cm = 13;  // Primary ionization clusters/cm
  double Ne_frac = tpcparam->get_double_param("Ne_frac");

  double Ar_dEdx = 2.44;  // keV/cm
  double Ar_primary_clusters_per_cm = 23;  // Primary ionization clusters/cm
  double Ar_ionization_electrons_per_cm = 97;  // Total ionized electrons/cm for a MIP
  double Ar_frac = tpcparam->get_double_param("Ar_frac");

  double CF4_dEdx = 6.38;     // keV/cm
  double CF4_primary_clusters_per_cm = 51;  // Primary ionization clusters/cm
  double CF4_ionization_electrons_per_cm = 120;  // Total ionized electrons/cm for a MIP
  double CF4_frac = tpcparam->get_double_param("CF4_frac");

  double N2_dEdx = 2.127;  // keV/cm https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/nitrogen_gas.html
  double N2_primary_clusters_per_cm = 25;   // Primary ionization clusters/cm (approximate, small impact here)
  double N2_frac = tpcparam->get_double_param("N2_frac");

  double isobutane_dEdx = 5.67;  // keV/cm
  double isobutane_primary_clusters_per_cm = 84;  // Primary ionization clusters/cm
  double isobutane_ionization_electrons_per_cm = 220;  // Total ionized electrons/cm for a MIP
  double isobutane_frac = tpcparam->get_double_param("isobutane_frac");

  m_use_primary_cluster_ionization =
      tpcparam->exist_int_param("use_primary_cluster_ionization") &&
      tpcparam->get_int_param("use_primary_cluster_ionization") != 0;
  m_primaryIonizationModel.reset();
  if (m_use_primary_cluster_ionization)
  {
    // For the electron/positron energy-based fallback of the opt-in model,
    // use the internally consistent stopping powers and total ionization
    // yields in Table 1 of arXiv:1008.3736. Keep the established values above
    // unchanged when the optional model is disabled.
    CF4_dEdx = 6.38;         // keV/cm
    CF4_ionization_electrons_per_cm = 120;        // Number/cm
    isobutane_dEdx = 5.67;   // keV/cm
    isobutane_ionization_electrons_per_cm = 220;  // Number/cm

    constexpr double fraction_tolerance = 1.e-12;
    if (!std::isfinite(Ne_frac) || !std::isfinite(N2_frac) ||
        std::abs(Ne_frac) > fraction_tolerance ||
        std::abs(N2_frac) > fraction_tolerance)
    {
      std::cerr
          << Name()
          << ": the optional primary-cluster ionization model supports only "
             "Ar/CF4/iC4H10 mixtures; Ne_frac="
          << Ne_frac << " and N2_frac=" << N2_frac << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    try
    {
      m_primaryIonizationModel =
          std::make_unique<TpcPrimaryIonizationModel>(
              TpcPrimaryIonizationModel::GasFractions{
                  Ar_frac, CF4_frac, isobutane_frac});
    }
    catch (const std::exception &error)
    {
      std::cerr << Name()
                << ": failed to configure the optional primary-cluster "
                   "ionization model: "
                << error.what() << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (Verbosity() > 0)
    {
      std::cout
          << Name()
          << ": using path-length primary-cluster ionization for "
             "non-electron charged particles; electrons and positrons use "
             "the energy-deposit response; "
          << m_primaryIonizationModel->primary_clusters_per_cm()
          << " primary clusters/cm, mean cluster size "
          << m_primaryIonizationModel->mean_cluster_size() << std::endl;
    }
  }

  if (m_use_PDG_gas_params)
  {
    Ne_dEdx = 1.446;
    Ne_primary_clusters_per_cm = 13;

    Ar_dEdx = 2.525;

    CF4_dEdx = 6.382;
  }

  double Tpc_primary_clusters_per_cm = (Ne_primary_clusters_per_cm * Ne_frac) + (Ar_primary_clusters_per_cm * Ar_frac) + (CF4_primary_clusters_per_cm * CF4_frac) + (N2_primary_clusters_per_cm * N2_frac) + (isobutane_primary_clusters_per_cm * isobutane_frac);

  double Tpc_gas_mixture_dEdx = (Ne_dEdx * Ne_frac) + (Ar_dEdx * Ar_frac) + (CF4_dEdx * CF4_frac) + (N2_dEdx * N2_frac) + (isobutane_dEdx * isobutane_frac);

  double Tpc_supported_ionization_gases_dEdx = (Ar_dEdx * Ar_frac) + (CF4_dEdx * CF4_frac) + (isobutane_dEdx * isobutane_frac);
  ionization_electrons_per_cm = (Ar_ionization_electrons_per_cm * Ar_frac) + (CF4_ionization_electrons_per_cm * CF4_frac) + (isobutane_ionization_electrons_per_cm * isobutane_frac);
  avg_ionization_energy_kev_per_electron = (ionization_electrons_per_cm > 0.0) ? (Tpc_supported_ionization_gases_dEdx / ionization_electrons_per_cm) : std::numeric_limits<double>::quiet_NaN();
  ionization_electrons_per_gev = (Tpc_supported_ionization_gases_dEdx > 0.0) ? (ionization_electrons_per_cm / Tpc_supported_ionization_gases_dEdx) * 1e6 : std::numeric_limits<double>::quiet_NaN();

  if (Ne_frac > 0.0 || N2_frac > 0.0)
  {
    std::cout << "PHG4TpcElectronDrift::InitRun - Warning: total-ionization-electron dE/dx QA uses only Ar/CF4/isobutane constants; "
              << "Ne_frac=" << Ne_frac << " and N2_frac=" << N2_frac
              << " are excluded from the ionization-electron energy scale." << std::endl;
  }

  primary_clusters_per_cm = Tpc_primary_clusters_per_cm;
  if (fixed_primary_clusters_per_cm > 0.0)
  {
    primary_clusters_per_cm = fixed_primary_clusters_per_cm;
  }
  primary_clusters_per_gev = (Tpc_gas_mixture_dEdx > 0.0) ? (primary_clusters_per_cm / Tpc_gas_mixture_dEdx) * 1e6 : std::numeric_limits<double>::quiet_NaN();

  std::cout << "PHG4TpcElectronDrift::InitRun - primary clusters/cm = " << primary_clusters_per_cm
            << ", primary clusters/GeV (from dE/dx) = " << primary_clusters_per_gev << std::endl;
  std::cout << "PHG4TpcElectronDrift::InitRun - ionization electrons/cm = " << ionization_electrons_per_cm
            << ", average ionization energy = " << avg_ionization_energy_kev_per_electron << " keV/electron"
            << ", ionization electrons/GeV (from dE/dx) = " << ionization_electrons_per_gev << std::endl;
  if (fixed_primary_clusters_per_cm > 0.0)
  {
    std::cout << "PHG4TpcElectronDrift::InitRun - overriding gas-derived primary cluster density with fixed_primary_clusters_per_cm = "
              << fixed_primary_clusters_per_cm << std::endl;
  }

  // Preserve the upstream energy-deposit response (also used for the
  // electron/positron fallback of the opt-in path-length model). The local
  // primary-cluster QA energy scale above is an independent diagnostic.
  const double energy_ne_total = m_use_PDG_gas_params ? 40.0 : 43.0;
  const double energy_ar_total = m_use_PDG_gas_params ? 97.0 : 94.0;
  const double energy_cf4_total = (m_use_primary_cluster_ionization || m_use_PDG_gas_params) ? 120.0 : 100.0;
  const double energy_isobutane_total = m_use_primary_cluster_ionization ? 220.0 : 195.0;
  const double energy_cf4_dedx = m_use_PDG_gas_params ? 6.382 : (m_use_primary_cluster_ionization ? 6.38 : 7.0);
  const double energy_isobutane_dedx = m_use_primary_cluster_ionization ? 5.67 : 5.93;
  const double energy_total = energy_ne_total * Ne_frac + energy_ar_total * Ar_frac +
                              energy_cf4_total * CF4_frac + 25.0 * N2_frac + energy_isobutane_total * isobutane_frac;
  const double energy_dedx = Ne_dEdx * Ne_frac + Ar_dEdx * Ar_frac + energy_cf4_dedx * CF4_frac +
                             N2_dEdx * N2_frac + energy_isobutane_dedx * isobutane_frac;
  if (!(energy_dedx > 0.0) || !std::isfinite(energy_total))
  {
    std::cerr << Name() << ": invalid ionization energy normalization" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  electrons_per_gev = energy_total / energy_dedx * 1.e6;

  // Initialize cluster size CDF for sPHENIX Ar/CF4/iC4H10 mixture.
  constexpr int kClusterSizeCutoff = 384;
  // Argon probabilities from the provided table (column b), converted to fractions.
  const std::map<int, double> ar_cluster_prob = {
      {1, 0.656},
      {2, 0.150},
      {3, 0.064},
      {4, 0.035},
      {5, 0.0225},
      {6, 0.0155},
      {7, 0.0105},
      {8, 0.0081},
      {9, 0.0061},
      {10, 0.0049},
      {11, 0.0039},
      {12, 0.0030},
      {13, 0.0025},
      {14, 0.0020},
      {15, 0.0016},
      {16, 0.0012},
      {17, 0.00095},
      {18, 0.00075},
      {19, 0.00063}};
  // CH4 probabilities from the provided table (column b), converted to fractions.
  const std::map<int, double> ch4_cluster_prob = {
      {1, 0.786},
      {2, 0.120},
      {3, 0.034},
      {4, 0.016},
      {5, 0.0095},
      {6, 0.0060},
      {7, 0.0044},
      {8, 0.0034},
      {9, 0.0027},
      {10, 0.0021},
      {11, 0.0017},
      {12, 0.0013},
      {13, 0.0010},
      {14, 0.0008},
      {15, 0.0006},
      {16, 0.00050},
      {17, 0.00042},
      {18, 0.00037},
      {19, 0.00033}};

  auto build_cluster_prob = [&](const std::map<int, double> &discrete, int tail_start)
  {
    std::vector<double> prob(kClusterSizeCutoff + 1, 0.0);
    double sum_discrete = 0.0;
    for (auto const &[n, p] : discrete)
    {
      if (n > kClusterSizeCutoff)
      {
        continue;
      }
      prob[n] = p;
      sum_discrete += p;
    }
    double remaining = 1.0 - sum_discrete;
    if (remaining > 0.0 && tail_start <= kClusterSizeCutoff)
    {
      double sum_inv_sq = 0.0;
      for (int n = tail_start; n <= kClusterSizeCutoff; ++n)
      {
        sum_inv_sq += 1.0 / (n * n);
      }
      double C = (sum_inv_sq > 0.0) ? remaining / sum_inv_sq : 0.0;
      for (int n = tail_start; n <= kClusterSizeCutoff; ++n)
      {
        if (prob[n] == 0.0)
        {
          prob[n] = C / (n * n);
        }
      }
    }
    double norm = 0.0;
    for (int n = 1; n <= kClusterSizeCutoff; ++n)
    {
      norm += prob[n];
    }
    if (norm > 0.0)
    {
      for (int n = 1; n <= kClusterSizeCutoff; ++n)
      {
        prob[n] /= norm;
      }
    }
    return prob;
  };

  const auto prob_ar = build_cluster_prob(ar_cluster_prob, 20);
  const auto prob_ch4 = build_cluster_prob(ch4_cluster_prob, 20);

  const double mix_primary = (Ar_frac * Ar_primary_clusters_per_cm) + (CF4_frac * CF4_primary_clusters_per_cm) + (isobutane_frac * isobutane_primary_clusters_per_cm);
  const double weight_ar = (mix_primary > 0.0) ? (Ar_frac * Ar_primary_clusters_per_cm / mix_primary) : 0.0;
  const double weight_cf4 = (mix_primary > 0.0) ? (CF4_frac * CF4_primary_clusters_per_cm / mix_primary) : 0.0;
  const double weight_iso = (mix_primary > 0.0) ? (isobutane_frac * isobutane_primary_clusters_per_cm / mix_primary) : 0.0;

  std::vector<double> mix_prob(kClusterSizeCutoff + 1, 0.0);
  if (mix_primary > 0.0)
  {
    for (int n = 1; n <= kClusterSizeCutoff; ++n)
    {
      mix_prob[n] = (weight_ar + weight_cf4) * prob_ar[n] + weight_iso * prob_ch4[n];
    }
  }
  else
  {
    mix_prob = prob_ar;
    std::cout << "PHG4TpcElectronDrift::InitRun - Warning: gas mixture has no Ar, CF4, or isobutane content, defaulting to pure Argon cluster size distribution." << std::endl;
  }

  double mix_norm = 0.0;
  for (int n = 1; n <= kClusterSizeCutoff; ++n)
  {
    mix_norm += mix_prob[n];
  }
  if (mix_norm > 0.0)
  {
    for (int n = 1; n <= kClusterSizeCutoff; ++n)
    {
      mix_prob[n] /= mix_norm;
    }
  }

  cluster_size_cdf.clear();
  cluster_size_cdf.resize(kClusterSizeCutoff + 1, 0.0);
  double current_cdf = 0.0;
  for (int n = 1; n <= kClusterSizeCutoff; ++n)
  {
    current_cdf += mix_prob[n];
    cluster_size_cdf[n] = current_cdf;
  }
  cluster_size_cdf[kClusterSizeCutoff] = 1.0;

  // min_time to max_time is the time window for accepting drifted electrons after the trigger
  min_time = 0.0;
  const bool use_garfield_full_lut =
      m_distortionMap &&
      m_distortionMap->application_mode() ==
          PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport;
  const bool use_garfield_full_direct =
      m_garfield && m_garfield_full_transport;
  const double maximum_deterministic_drift_time =
      use_garfield_full_lut
          ? m_distortionMap->garfield_lut_metadata().maximum_drift_time_ns
          : (use_garfield_full_direct
                 ? m_garfield_step_ns *
                       static_cast<double>(
                           PHGarfield::MaxForwardDriftSteps)
                 : max_drift_length / drift_velocity);
  max_time =
      maximum_deterministic_drift_time +
      z_geometry->get_extended_readout_time();
  min_active_radius = get_double_param("min_active_radius");
  max_active_radius = get_double_param("max_active_radius");

  if (Verbosity() > 0)
  {
    std::cout << Name() << " gas mixture fractions (Ne/Ar/CF4/N2/iC4H10): "
              << Ne_frac << "/" << Ar_frac << "/" << CF4_frac << "/" << N2_frac << "/" << isobutane_frac << std::endl;
    std::cout << Name() << " primary ionization summary: primary clusters per cm " << primary_clusters_per_cm
              << ", dE/dx (keV/cm) " << Tpc_gas_mixture_dEdx
              << ", primary clusters per GeV " << primary_clusters_per_gev
              << ", ionization electrons per cm " << ionization_electrons_per_cm
              << ", average ionization energy (keV/electron) " << avg_ionization_energy_kev_per_electron << std::endl;
    std::cout << Name() << " diffusion sigmas (long/trans) [cm^0.5]: " << diffusion_long << "/" << diffusion_trans
              << " with additional smearing (long/trans): " << added_smear_sigma_long << "/" << added_smear_sigma_trans << std::endl;
    std::cout << Name() << " drift window [min,max] (ns): " << min_time << ", " << max_time << std::endl;
    std::cout << Name() << " acceptance radii [min,max] (cm): " << min_active_radius << ", " << max_active_radius << std::endl;
    if (force_min_trans_drift_length > 0.)
    {
      std::cout << Name() << " forcing minimum transverse drift length (cm): " << force_min_trans_drift_length << std::endl;
    }
    std::cout << PHWHERE << " drift velocity " << drift_velocity
              << " extended_readout_time " << z_geometry->get_extended_readout_time()
              << " max time cutoff " << max_time << std::endl;
  }

  auto se = Fun4AllServer::instance();
  if (do_ElectronDriftQAHistos)
  {
    dlong = new TH1F("difflong", "longitudinal diffusion", 100, diffusion_long - diffusion_long / 2., diffusion_long + diffusion_long / 2.);
    se->registerHisto(dlong);
    dtrans = new TH1F("difftrans", "transversal diffusion", 100, diffusion_trans - diffusion_trans / 2., diffusion_trans + diffusion_trans / 2.);
    se->registerHisto(dtrans);
  }

  if (do_ElectronDriftQAHistos)
  {
    diffDistance = new TH1F("diffDistance", "Transverse diffusion displacement;#Delta r (cm);Counts", 300, 0.0, 3.0);
    diffDX = new TH1F("diffDX", "Transverse diffusion #Delta x;#Delta x (cm);Counts", 400, -3.0, 3.0);
    diffDY = new TH1F("diffDY", "Transverse diffusion #Delta y;#Delta y (cm);Counts", 400, -3.0, 3.0);
    nPrimaryClusters = new TH1F("nPrimaryClusters", "Sampled primary clusters per step;N_{primary clusters};Counts", 400, -0.5, 399.5);
    primaryClusterMean = new TH1F("primaryClusterMean", "Poisson mean primary clusters per step;#bar{N}_{primary clusters};Counts", 400, 0.0, 200.0);
    nPrimaryClustersVsMean = new TH2F("nPrimaryClustersVsMean", "Sampled primary clusters vs. mean;#bar{N}_{primary clusters};N_{primary clusters}", 200, 0.0, 200.0, 200, -0.5, 199.5);
    diffPerSqrtL = new TH1F("diffPerSqrtL", "Transverse diffusion normalized by #sqrt{L};#Delta r/#sqrt{L} (cm^{0.5});Counts", 400, 0.0, 0.2);
    diffDXPerSqrtL = new TH1F("diffDXPerSqrtL", "#Delta x normalized by #sqrt{L};#Delta x/#sqrt{L} (cm^{0.5});Counts", 400, -0.2, 0.2);
    diffDYPerSqrtL = new TH1F("diffDYPerSqrtL", "#Delta y normalized by #sqrt{L};#Delta y/#sqrt{L} (cm^{0.5});Counts", 400, -0.2, 0.2);
    nPrimaryClustersPerCm = new TH1F("nPrimaryClustersPerCm", "Sampled primary clusters per cm;N_{primary clusters}/cm;Counts", 400, 0.0, 400.0);
    nIonizationElectrons = new TH1F("nIonizationElectrons", "Sampled ionization electrons per step;N_{e};Counts", 500, -0.5, 999.5);
    nIonizationElectronsPerCm = new TH1F("nIonizationElectronsPerCm", "Sampled ionization electrons per cm;N_{e}/cm;Counts", 500, 0.0, 1000.0);
    nIonizationElectronsVsPrimaryClusters = new TH2F("nIonizationElectronsVsPrimaryClusters", "Ionization electrons vs. primary clusters;N_{primary clusters};N_{e}", 200, -0.5, 199.5, 500, -0.5, 999.5);
    g4StepDedx = new TH1F("g4StepDedx", "Geant4 step dE/dx;dE/dx (keV/cm);Counts", 500, 0.0, 50.0);
    ionizedElectronsDedx = new TH1F("ionizedElectronsDedx", "Ionized electrons #times average ionization energy;dE/dx (keV/cm);Counts", 500, 0.0, 50.0);
    ionizedElectronsDedxVsG4StepDedx = new TH2F("ionizedElectronsDedxVsG4StepDedx", "Ionized-electron dE/dx vs. Geant4 dE/dx;Geant4 dE/dx (keV/cm);Ionized-electron dE/dx (keV/cm)", 300, 0.0, 50.0, 300, 0.0, 50.0);
    diffVsDrift = new TH2F("diffVsDrift", "Transverse diffusion vs. drift length;Drift length L (cm);#Delta r (cm)", 200, 0.0, tpc_length / 2., 300, 0.0, 3.0);
    diffDXVsDrift = new TH2F("diffDXVsDrift", "#Delta x vs. drift length;Drift length L (cm);#Delta x (cm)", 200, 0.0, tpc_length / 2., 400, -3.0, 3.0);
    diffDYVsDrift = new TH2F("diffDYVsDrift", "#Delta y vs. drift length;Drift length L (cm);#Delta y (cm)", 200, 0.0, tpc_length / 2., 400, -3.0, 3.0);
    diffPerSqrtLVsDrift = new TH2F("diffPerSqrtLVsDrift", "#Delta r/#sqrt{L} vs. drift length;Drift length L (cm);#Delta r/#sqrt{L} (cm^{0.5})", 200, 0.0, tpc_length / 2., 300, 0.0, 0.2);
    diffDXPerSqrtLVsDrift = new TH2F("diffDXPerSqrtLVsDrift", "#Delta x/#sqrt{L} vs. drift length;Drift length L (cm);#Delta x/#sqrt{L} (cm^{0.5})", 200, 0.0, tpc_length / 2., 400, -0.2, 0.2);
    diffDYPerSqrtLVsDrift = new TH2F("diffDYPerSqrtLVsDrift", "#Delta y/#sqrt{L} vs. drift length;Drift length L (cm);#Delta y/#sqrt{L} (cm^{0.5})", 200, 0.0, tpc_length / 2., 400, -0.2, 0.2);
    hitmapstart = new TH2F("hitmapstart", "g4hit starting X-Y locations", 1560, -78, 78, 1560, -78, 78);
    hitmapend = new TH2F("hitmapend", "g4hit final X-Y locations", 1560, -78, 78, 1560, -78, 78);
    hitmapstart_z = new TH2F("hitmapstart_z", "g4hit starting Z-R locations", 2000, -tpc_length / 2.0, tpc_length / 2.0, 780, 0, 78);
    hitmapend_z = new TH2F("hitmapend_z", "g4hit final Z-R locations", 2000, -tpc_length / 2.0, tpc_length / 2.0, 780, 0, 78);
    z_startmap = new TH2F("z_startmap", "g4hit starting Z vs. R locations", 2000, -tpc_length / 2.0, tpc_length / 2.0, 780, 0, 78);
    deltaphi = new TH2F("deltaphi", "Total delta phi; phi (rad);#Delta phi (rad)", 600, -M_PI, M_PI, 1000, -.2, .2);
    deltaRphinodiff = new TH2F("deltaRphinodiff", "Total delta R*phi, no diffusion; r (cm);#Delta R*phi (cm)", 600, 20, 80, 1000, -3, 5);
    deltaphivsRnodiff = new TH2F("deltaphivsRnodiff", "Total delta phi vs. R; phi (rad);#Delta phi (rad)", 600, 20, 80, 1000, -.2, .2);
    deltaz = new TH2F("deltaz", "Total delta z; z (cm);#Delta z (cm)", 1000, -tpc_length / 2.0, tpc_length / 2.0, 1000, -10, 10);
    deltaphinodiff = new TH2F("deltaphinodiff", "Total delta phi (no diffusion, only SC distortion); phi (rad);#Delta phi (rad)", 600, -M_PI, M_PI, 1000, -.2, .2);
    deltaphinodist = new TH2F("deltaphinodist", "Total delta phi (no SC distortion, only diffusion); phi (rad);#Delta phi (rad)", 600, -M_PI, M_PI, 1000, -.2, .2);
    deltar = new TH2F("deltar", "Total Delta r; r (cm);#Delta r (cm)", 580, 20, 78, 1000, -3, 5);
    deltarnodiff = new TH2F("deltarnodiff", "Delta r (no diffusion, only SC distortion); r (cm);#Delta r (cm)", 580, 20, 78, 1000, -2, 5);
    deltarnodist = new TH2F("deltarnodist", "Delta r (no SC distortion, only diffusion); r (cm);#Delta r (cm)", 580, 20, 78, 1000, -2, 5);
    ratioElectronsRR = new TH1F("ratioElectronsRR", "Ratio of electrons reach readout vs all in acceptance", 1561, -0.0325, 1.0465);
    TpcRootQaTree::reset(driftXY, "driftXY", "Electron start/end positions");
    driftXY->Branch("start_x", &m_drift_start_x, "start_x/F");
    driftXY->Branch("start_y", &m_drift_start_y, "start_y/F");
    driftXY->Branch("end_x", &m_drift_end_x, "end_x/F");
    driftXY->Branch("end_y", &m_drift_end_y, "end_y/F");
    TpcRootQaTree::reset(driftStepQA, "driftStepQA", "Per-step energy loss and dE/dx");
    driftStepQA->Branch("event", &m_stepqa_event, "event/I");
    driftStepQA->Branch("track_id", &m_stepqa_track_id, "track_id/I");
    driftStepQA->Branch("step_length_cm", &m_stepqa_step_length_cm, "step_length_cm/F");
    driftStepQA->Branch("edep_kev", &m_stepqa_edep_kev, "edep_kev/F");
    driftStepQA->Branch("eion_kev", &m_stepqa_eion_kev, "eion_kev/F");
    driftStepQA->Branch("dedx_kev_per_cm", &m_stepqa_dedx_kev_per_cm, "dedx_kev_per_cm/F");
    driftStepQA->Branch("deiondx_kev_per_cm", &m_stepqa_deiondx_kev_per_cm, "deiondx_kev_per_cm/F");
    driftStepQA->Branch("mean_primary_clusters", &m_stepqa_mean_primary_clusters, "mean_primary_clusters/F");
    driftStepQA->Branch("n_primary_clusters", &m_stepqa_n_primary_clusters, "n_primary_clusters/I");
    driftStepQA->Branch("ionization_source_model", &m_stepqa_source_model, "ionization_source_model/I");
    driftStepQA->Branch("n_ionization_sources", &m_stepqa_n_sources, "n_ionization_sources/I");
    driftStepQA->Branch("avg_ionization_energy_kev_per_electron", &m_stepqa_avg_ionization_energy_kev_per_electron, "avg_ionization_energy_kev_per_electron/F");
    driftStepQA->Branch("mean_ionization_electrons", &m_stepqa_mean_ionization_electrons, "mean_ionization_electrons/F");
    driftStepQA->Branch("n_ionization_electrons", &m_stepqa_n_ionization_electrons, "n_ionization_electrons/I");
    driftStepQA->Branch("ionized_electrons_edep_kev", &m_stepqa_ionized_electrons_edep_kev, "ionized_electrons_edep_kev/F");
    driftStepQA->Branch("ionized_electrons_dedx_kev_per_cm", &m_stepqa_ionized_electrons_dedx_kev_per_cm, "ionized_electrons_dedx_kev_per_cm/F");
    TpcRootQaTree::reset(ionizationEventQA, "ionizationEventQA", "Event-level primary-ionization QA");
    ionizationEventQA->Branch("event", &m_eventqa_event, "event/I");
    ionizationEventQA->Branch("use_pai_cluster_seeds", &m_eventqa_use_pai_cluster_seeds, "use_pai_cluster_seeds/I");
    ionizationEventQA->Branch("source_hits", &m_eventqa_source_hits, "source_hits/I");
    ionizationEventQA->Branch("path_length_cm", &m_eventqa_path_length_cm, "path_length_cm/F");
    ionizationEventQA->Branch("n_primary_clusters", &m_eventqa_n_primary_clusters, "n_primary_clusters/I");
    ionizationEventQA->Branch("n_energy_fallback_sources", &m_eventqa_n_energy_sources, "n_energy_fallback_sources/I");
    ionizationEventQA->Branch("n_ionization_electrons", &m_eventqa_n_ionization_electrons, "n_ionization_electrons/I");
    ionizationEventQA->Branch("primary_clusters_per_cm", &m_eventqa_primary_clusters_per_cm, "primary_clusters_per_cm/F");
    ionizationEventQA->Branch("ionization_electrons_per_cm", &m_eventqa_ionization_electrons_per_cm, "ionization_electrons_per_cm/F");

    TpcRootQaTree::reset(ionizationClusterQA, "ionizationClusterQA", "Primary-ionization cluster QA");
    ionizationClusterQA->Branch("event", &m_clusterqa_event, "event/I");
    ionizationClusterQA->Branch("use_pai_cluster_seeds", &m_clusterqa_use_pai_cluster_seeds, "use_pai_cluster_seeds/I");
    ionizationClusterQA->Branch("ionization_source_model", &m_clusterqa_source_model, "ionization_source_model/I");
    ionizationClusterQA->Branch("track_id", &m_clusterqa_track_id, "track_id/I");
    ionizationClusterQA->Branch("g4hit_id", &m_clusterqa_g4hit_id, "g4hit_id/I");
    ionizationClusterQA->Branch("cluster_index", &m_clusterqa_cluster_index, "cluster_index/I");
    ionizationClusterQA->Branch("segment_index", &m_clusterqa_segment_index, "segment_index/I");
    ionizationClusterQA->Branch("cluster_size", &m_clusterqa_cluster_size, "cluster_size/I");
    ionizationClusterQA->Branch("track_path_cm", &m_clusterqa_track_path_cm, "track_path_cm/F");
    ionizationClusterQA->Branch("x_cm", &m_clusterqa_x_cm, "x_cm/F");
    ionizationClusterQA->Branch("y_cm", &m_clusterqa_y_cm, "y_cm/F");
    ionizationClusterQA->Branch("z_cm", &m_clusterqa_z_cm, "z_cm/F");
    ionizationClusterQA->Branch("t_ns", &m_clusterqa_t_ns, "t_ns/F");

    TpcRootQaTree::reset(ionizationSegmentQA, "ionizationSegmentQA", "One-centimeter primary-ionization segment QA");
    ionizationSegmentQA->Branch("event", &m_segmentqa_event, "event/I");
    ionizationSegmentQA->Branch("use_pai_cluster_seeds", &m_segmentqa_use_pai_cluster_seeds, "use_pai_cluster_seeds/I");
    ionizationSegmentQA->Branch("track_id", &m_segmentqa_track_id, "track_id/I");
    ionizationSegmentQA->Branch("segment_index", &m_segmentqa_segment_index, "segment_index/I");
    ionizationSegmentQA->Branch("segment_start_cm", &m_segmentqa_segment_start_cm, "segment_start_cm/F");
    ionizationSegmentQA->Branch("segment_length_cm", &m_segmentqa_segment_length_cm, "segment_length_cm/F");
    ionizationSegmentQA->Branch("n_primary_clusters", &m_segmentqa_n_primary_clusters, "n_primary_clusters/I");
    ionizationSegmentQA->Branch("n_energy_fallback_sources", &m_segmentqa_n_energy_sources, "n_energy_fallback_sources/I");
    ionizationSegmentQA->Branch("n_ionization_electrons", &m_segmentqa_n_ionization_electrons, "n_ionization_electrons/I");
    ionizationSegmentQA->Branch("primary_clusters_per_cm", &m_segmentqa_primary_clusters_per_cm, "primary_clusters_per_cm/F");
    ionizationSegmentQA->Branch("ionization_electrons_per_cm", &m_segmentqa_ionization_electrons_per_cm, "ionization_electrons_per_cm/F");
  }

  if (m_avg_x_enabled)
  {
    m_avgOutf.reset(new TFile(m_avg_output_file.c_str(), "recreate"));
    if (!m_avgOutf || m_avgOutf->IsZombie())
    {
      std::cout << Name() << ": failed to open avg output file " << m_avg_output_file << std::endl;
    }
    if (m_avgOutf)
    {
      m_avgOutf->cd();
    }
    m_avgXResidualTree = new TTree("avgXResidual", "Track-level <x>-x_int residuals");
    m_avgXResidualTree->Branch("trackid", &m_avgTree_trackid, "trackid/F");
    m_avgXResidualTree->Branch("residual", &m_avgTree_residual, "residual/F");
    m_avgXResidualTree->Branch("avgx", &m_avgTree_avgx, "avgx/F");
    m_avgXResidualTree->Branch("xint", &m_avgTree_xint, "xint/F");
    m_avgXResidualTree->Branch("residual_primary", &m_avgTree_residual_primary, "residual_primary/F");
    m_avgXResidualTree->Branch("avgx_primary", &m_avgTree_avgx_primary, "avgx_primary/F");
    m_avgXResidualTree->Branch("count_primary", &m_avgTree_count_primary, "count_primary/F");
    m_avgXResidualTree->Branch("residual_end", &m_avgTree_residual_end, "residual_end/F");
    m_avgXResidualTree->Branch("avgx_end", &m_avgTree_avgx_end, "avgx_end/F");
    m_avgXResidualTree->Branch("count_end", &m_avgTree_count_end, "count_end/F");
    m_avgXResidualTree->Branch("residual_end_primary", &m_avgTree_residual_end_primary, "residual_end_primary/F");
    m_avgXResidualTree->Branch("avgx_end_primary", &m_avgTree_avgx_end_primary, "avgx_end_primary/F");
    m_avgXResidualTree->Branch("count_end_primary", &m_avgTree_count_end_primary, "count_end_primary/F");
    m_avgXResidualTree->Branch("residual_x", &m_avgTree_residual_x, "residual_x/F");
    m_avgXResidualTree->Branch("avgx_cart", &m_avgTree_avgx_cart, "avgx_cart/F");
    m_avgXResidualTree->Branch("residual_x_primary", &m_avgTree_residual_x_primary, "residual_x_primary/F");
    m_avgXResidualTree->Branch("avgx_cart_primary", &m_avgTree_avgx_cart_primary, "avgx_cart_primary/F");
    m_avgXResidualTree->Branch("residual_x_end", &m_avgTree_residual_x_end, "residual_x_end/F");
    m_avgXResidualTree->Branch("avgx_cart_end", &m_avgTree_avgx_cart_end, "avgx_cart_end/F");
    m_avgXResidualTree->Branch("residual_x_end_primary", &m_avgTree_residual_x_end_primary, "residual_x_end_primary/F");
    m_avgXResidualTree->Branch("avgx_cart_end_primary", &m_avgTree_avgx_cart_end_primary, "avgx_cart_end_primary/F");
  }

  if (Verbosity())
  {
    // eval tree only when verbosity is on
    m_outf.reset(new TFile("nt_out.root", "recreate"));
    nt = new TNtuple("nt", "electron drift stuff", "hit:ts:tb:tsig:rad:zstart:zfinal");
    nthit = new TNtuple("nthit", "TrkrHit collecting", "layer:phipad:zbin:neffelectrons");
    ntfinalhit = new TNtuple("ntfinalhit", "TrkrHit collecting", "layer:phipad:zbin:neffelectrons");
    ntpad = new TNtuple("ntpad", "electron by electron pad centroid", "layer:phigem:phiclus:zgem:zclus");
    se->registerHisto(nt);
    se->registerHisto(nthit);
    se->registerHisto(ntpad);
  }

  const int padplane_status = padplane->InitRun(topNode);
  if (padplane_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << PHWHERE
              << "TPC pad-plane initialization failed with status "
              << padplane_status << std::endl;
    return padplane_status;
  }

  // print all layers radii
  if (Verbosity())
  {
    const auto range = seggeo->get_begin_end();
    std::cout << "PHG4TpcElectronDrift::InitRun - layers: " << std::distance(range.first, range.second) << std::endl;
    int counter = 0;
    for (auto layeriter = range.first; layeriter != range.second; ++layeriter)
    {
      const auto radius = layeriter->second->get_radius();
      std::cout << boost::str(boost::format("%.3f ") % radius);
      if (++counter == 8)
      {
        counter = 0;
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
  }

  if (record_ClusHitsVerbose)
  {
    // get the node
    mClusHitsVerbose = findNode::getClass<ClusHitsVerbosev1>(topNode, "Trkr_TruthClusHitsVerbose");
    if (!mClusHitsVerbose)
    {
      PHNodeIterator dstiter(dstNode);
      auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode("TRKR");
        dstNode->addNode(DetNode);
      }
      mClusHitsVerbose = new ClusHitsVerbosev1();
      auto newNode = new PHIODataNode<PHObject>(mClusHitsVerbose, "Trkr_TruthClusHitsVerbose", "PHObject");
      DetNode->addNode(newNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcElectronDrift::process_event(PHCompositeNode *topNode)
{
  padplane->BeginEvent(event_num);

  truth_track = nullptr;  // track to which truth clusters are built

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "ActsGeometry not found on node tree. Exiting" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_track_layer_data.clear();

  if (truth_clusterer.needs_input_nodes())
  {
    truth_clusterer.set_input_nodes(truthclustercontainer, m_tGeometry,
                                    seggeo, mClusHitsVerbose);
  }

  static constexpr unsigned int print_layer = 49;

  // tells m_distortionMap which event to look at
  if (m_distortionMap)
  {
    m_distortionMap->load_event(event_num);
  }

  // g4hits
  auto g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << hitnodename << std::endl;
    gSystem->Exit(1);
  }
  PHG4HitContainer *pai_cluster_hits = nullptr;
  if (m_use_pai_cluster_seeds)
  {
    pai_cluster_hits = findNode::getClass<PHG4HitContainer>(topNode, m_pai_cluster_node_name);
    if (!pai_cluster_hits)
    {
      std::cout << Name() << " requested PAI cluster seeds, but could not locate node "
                << m_pai_cluster_node_name << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  PHG4TruthInfoContainer *truthinfo =
      findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  PHG4HitContainer *drift_source_hits = m_use_pai_cluster_seeds ? pai_cluster_hits : g4hit;
  PHG4HitContainer::ConstRange hit_begin_end = drift_source_hits->getHits();
  const std::size_t drift_source_size = drift_source_hits->size();
  bool event_has_tpc_secondaries = false;
  auto is_secondary_track = [&](const int track_id) -> bool
  {
    if (track_id < 0)
    {
      return true;
    }

    /* if (truthinfo)
    {
      if (const auto *particle = truthinfo->GetParticle(track_id))
      {
        return !truthinfo->is_primary(particle);
      }
    } */

    return false;
  };

  if (m_qa_write_only_with_secondaries)
  {
    if (m_use_pai_cluster_seeds && drift_source_size > 0)
    {
      event_has_tpc_secondaries = true;
      m_seen_event_with_secondaries = true;
    }
    for (auto hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      if (is_secondary_track(hiter->second->get_trkid()))
      {
        event_has_tpc_secondaries = true;
        m_seen_event_with_secondaries = true;
        break;
      }
    }
  }
  const bool fill_qa_hists_this_event =
      do_ElectronDriftQAHistos &&
      (!m_qa_write_only_with_secondaries || event_has_tpc_secondaries);

  struct OriginalTrackPathSegment
  {
    double track_path_start_cm{0.0};
    double length_cm{0.0};
    double x0{0.0};
    double y0{0.0};
    double z0{0.0};
    double dx{0.0};
    double dy{0.0};
    double dz{0.0};
  };

  struct IonizationSegmentAccumulator
  {
    double segment_start_cm{0.0};
    double segment_length_cm{0.0};
    int n_primary_clusters{0};
    int n_energy_sources{0};
    int n_ionization_electrons{0};
  };

  static constexpr double ionization_qa_segment_length_cm = 1.0;
  std::unordered_map<int, std::vector<OriginalTrackPathSegment>> original_track_path_segments;
  std::unordered_map<int, double> original_track_path_lengths;
  std::map<std::pair<int, int>, IonizationSegmentAccumulator> ionization_segment_accumulators;

  auto clamp_track_path_cm = [&](const int track_id, double track_path_cm) -> double
  {
    auto length_iter = original_track_path_lengths.find(track_id);
    if (length_iter != original_track_path_lengths.end() && length_iter->second > 0.0)
    {
      track_path_cm = std::clamp(track_path_cm, 0.0, length_iter->second);
      if (track_path_cm >= length_iter->second)
      {
        track_path_cm = std::nextafter(length_iter->second, 0.0);
      }
    }
    return track_path_cm;
  };

  auto add_segment_length = [&](const int track_id, const double track_path_start_cm, const double length_cm)
  {
    if (length_cm <= 0.0)
    {
      return;
    }

    double cursor = track_path_start_cm;
    double remaining = length_cm;
    while (remaining > 1e-9)
    {
      const int segment_index = static_cast<int>(std::floor(cursor / ionization_qa_segment_length_cm));
      const double segment_start_cm = segment_index * ionization_qa_segment_length_cm;
      const double segment_end_cm = segment_start_cm + ionization_qa_segment_length_cm;
      const double step = std::min(remaining, segment_end_cm - cursor);
      if (step <= 0.0)
      {
        break;
      }

      auto &segment = ionization_segment_accumulators[std::make_pair(track_id, segment_index)];
      segment.segment_start_cm = segment_start_cm;
      segment.segment_length_cm += step;

      cursor += step;
      remaining -= step;
    }
  };

  auto path_from_position = [&](const int track_id, const PHG4Hit *hit, const double fallback_track_path_cm) -> double
  {
    auto segment_iter = original_track_path_segments.find(track_id);
    if (segment_iter == original_track_path_segments.end() || segment_iter->second.empty())
    {
      return clamp_track_path_cm(track_id, fallback_track_path_cm);
    }

    const double x = hit->get_x(0);
    const double y = hit->get_y(0);
    const double z = hit->get_z(0);
    double best_distance2 = std::numeric_limits<double>::max();
    double best_track_path_cm = fallback_track_path_cm;

    for (const auto &segment : segment_iter->second)
    {
      if (segment.length_cm <= 0.0)
      {
        continue;
      }

      const double length2 = segment.length_cm * segment.length_cm;
      double f = ((x - segment.x0) * segment.dx + (y - segment.y0) * segment.dy + (z - segment.z0) * segment.dz) / length2;
      f = std::clamp(f, 0.0, 1.0);

      const double x_projected = segment.x0 + f * segment.dx;
      const double y_projected = segment.y0 + f * segment.dy;
      const double z_projected = segment.z0 + f * segment.dz;
      const double distance2 = square(x - x_projected) + square(y - y_projected) + square(z - z_projected);
      if (distance2 < best_distance2)
      {
        best_distance2 = distance2;
        best_track_path_cm = segment.track_path_start_cm + f * segment.length_cm;
      }
    }

    return clamp_track_path_cm(track_id, best_track_path_cm);
  };

  auto add_cluster_to_segment = [&](const int track_id, const double track_path_cm, const int cluster_size, const bool is_primary_source) -> int
  {
    if (!fill_qa_hists_this_event || !std::isfinite(track_path_cm) || track_path_cm < 0.0 || cluster_size <= 0)
    {
      return -1;
    }

    const double clamped_track_path_cm = clamp_track_path_cm(track_id, track_path_cm);
    const int segment_index = static_cast<int>(std::floor(clamped_track_path_cm / ionization_qa_segment_length_cm));
    if (segment_index < 0)
    {
      return -1;
    }

    auto &segment = ionization_segment_accumulators[std::make_pair(track_id, segment_index)];
    segment.segment_start_cm = segment_index * ionization_qa_segment_length_cm;
    if (is_primary_source) ++segment.n_primary_clusters;
    else ++segment.n_energy_sources;
    segment.n_ionization_electrons += cluster_size;
    return segment_index;
  };

  double event_path_length_cm = 0.0;
  if (fill_qa_hists_this_event)
  {
    std::unordered_map<int, double> original_track_path_offsets;
    PHG4HitContainer::ConstRange g4hit_begin_end = g4hit->getHits();
    for (auto path_iter = g4hit_begin_end.first; path_iter != g4hit_begin_end.second; ++path_iter)
    {
      const double t0 = std::fmax(path_iter->second->get_t(0), path_iter->second->get_t(1));
      if (t0 > max_time)
      {
        continue;
      }
      const double dx_hit = path_iter->second->get_x(1) - path_iter->second->get_x(0);
      const double dy_hit = path_iter->second->get_y(1) - path_iter->second->get_y(0);
      const double dz_hit = path_iter->second->get_z(1) - path_iter->second->get_z(0);
      const double step_length = std::sqrt(square(dx_hit) + square(dy_hit) + square(dz_hit));
      event_path_length_cm += step_length;

      const int track_id = path_iter->second->get_trkid();
      const double track_path_start_cm = original_track_path_offsets[track_id];
      if (step_length > 0.0)
      {
        OriginalTrackPathSegment segment;
        segment.track_path_start_cm = track_path_start_cm;
        segment.length_cm = step_length;
        segment.x0 = path_iter->second->get_x(0);
        segment.y0 = path_iter->second->get_y(0);
        segment.z0 = path_iter->second->get_z(0);
        segment.dx = dx_hit;
        segment.dy = dy_hit;
        segment.dz = dz_hit;
        original_track_path_segments[track_id].push_back(segment);
        add_segment_length(track_id, track_path_start_cm, step_length);
      }
      original_track_path_offsets[track_id] = track_path_start_cm + step_length;
      original_track_path_lengths[track_id] = track_path_start_cm + step_length;
    }

    m_eventqa_event = event_num;
    m_eventqa_use_pai_cluster_seeds = m_use_pai_cluster_seeds ? 1 : 0;
    m_eventqa_source_hits = static_cast<int>(drift_source_size);
    m_eventqa_path_length_cm = static_cast<float>(event_path_length_cm);
    m_eventqa_n_primary_clusters = 0;
    m_eventqa_n_energy_sources = 0;
    m_eventqa_n_ionization_electrons = 0;
    m_eventqa_primary_clusters_per_cm = std::numeric_limits<float>::quiet_NaN();
    m_eventqa_ionization_electrons_per_cm = std::numeric_limits<float>::quiet_NaN();
  }

  if (m_truth_intersection_hits)
  {
    m_truth_intersection_hits->Reset();
  }

  struct TruthLayerIntersection
  {
    int track_id{0};
    unsigned int layer{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double t{0.0};
    //double dirx{0.0};
    //double diry{0.0};
    //double dirz{0.0};
    double path{0.0};
    //double eion{0.0};
    //double edep{0.0};
  };
  std::map<std::pair<int, unsigned int>, std::vector<TruthLayerIntersection>> truth_layer_intersections;

  unsigned int count_g4hits = 0;
  //  int count_electrons = 0;

  m_track_path_offset.clear();
  m_track_anchor_set.clear();
  m_track_anchor_length.clear();

  //  double ecollectedhits = 0.0;
  //  int ncollectedhits = 0;
  double ihit = 0;
  unsigned int dump_interval = 5000;  // dump temp_hitsetcontainer to the node tree after this many g4hits
  unsigned int dump_counter = 0;

  const auto flush_temp_hitsets = [&]()
  {
    const auto temp_hitset_range = temp_hitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
    if (temp_hitset_range.first == temp_hitset_range.second)
    {
      return;
    }

    double eg4hit = 0.0;
    for (auto temp_hitset_iter = temp_hitset_range.first;
         temp_hitset_iter != temp_hitset_range.second;
         ++temp_hitset_iter)
    {
      const auto node_hitsetkey = temp_hitset_iter->first;
      const auto layer = TrkrDefs::getLayer(node_hitsetkey);
      const auto sector = TpcDefs::getSectorId(node_hitsetkey);
      const auto side = TpcDefs::getSide(node_hitsetkey);
      if (Verbosity() > 100)
      {
        std::cout << "PHG4TpcElectronDrift: temp_hitset with key: " << node_hitsetkey << " in layer " << layer
                  << " with sector " << sector << " side " << side << std::endl;
      }

      auto node_hitsetit = hitsetcontainer->findOrAddHitSet(node_hitsetkey);
      const auto temp_hit_range = temp_hitset_iter->second->getHits();
      for (auto temp_hit_iter = temp_hit_range.first;
           temp_hit_iter != temp_hit_range.second;
           ++temp_hit_iter)
      {
        const auto temp_hitkey = temp_hit_iter->first;
        auto *temp_tpchit = temp_hit_iter->second;
        if (Verbosity() > 10 && layer == print_layer)
        {
          std::cout << "      temp_hitkey " << temp_hitkey << " layer " << layer << " pad " << TpcDefs::getPad(temp_hitkey)
                    << " z bin " << TpcDefs::getTBin(temp_hitkey)
                    << "  energy " << temp_tpchit->getEnergy() << " eg4hit " << eg4hit << std::endl;
          eg4hit += temp_tpchit->getEnergy();
        }

        auto *node_hit = node_hitsetit->second->getHit(temp_hitkey);
        if (!node_hit)
        {
          node_hit = new TrkrHitv2();
          node_hitsetit->second->addHitSpecificKey(temp_hitkey, node_hit);
        }
        node_hit->addEnergy(temp_tpchit->getEnergy());
      }

      if (Verbosity() > 100 && layer == print_layer)
      {
        std::cout << "  ihit " << ihit << " collected energy = " << eg4hit << std::endl;
      }
    }

    temp_hitsetcontainer->Reset();
  };

  int trkid = -1;

  PHG4Hit *prior_g4hit = nullptr;  // used to check for jumps in g4hits;
  // if there is a big jump (such as crossing into the INTT area or out of the TPC)
  // then cluster the truth clusters before adding a new hit. This prevents
  // clustering loopers in the same HitSetKey surfaces in multiple passes
  for (auto hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    ++count_g4hits;
    dump_counter++;

    const double t0 = std::fmax(hiter->second->get_t(0), hiter->second->get_t(1));
    if (t0 > max_time)
    {
      continue;
    }

    int trkid_new = hiter->second->get_trkid();
    if (trkid != trkid_new)
    {  // starting a new track
      prior_g4hit = nullptr;
      if (truth_track)
      {
        truth_clusterer.cluster_hits(truth_track);
      }
      trkid = trkid_new;

      if (Verbosity() > 1000)
      {
        std::cout << " New track : " << trkid << " is embed? : ";
      }

      if (truthinfo && truthinfo->isEmbeded(trkid))
      {
        truth_track = truthtracks->getTruthTrack(trkid, truthinfo);
        truth_clusterer.b_collect_hits = true;
        if (Verbosity() > 1000)
        {
          std::cout << " YES embedded" << std::endl;
        }
      }
      else
      {
        truth_track = nullptr;
        truth_clusterer.b_collect_hits = false;
        if (Verbosity() > 1000)
        {
          std::cout << " NOT embedded" << std::endl;
        }
      }

      if (m_avg_x_enabled && m_track_map && trkid >= 0)
      {
        auto it = m_track_map->find(static_cast<unsigned int>(trkid));
        if (it != m_track_map->end())
        {
          const SvtxTrack *trk = it->second;
          const double vx = trk->get_px();
          const double vy = trk->get_py();
          const double vz = trk->get_pz();
          const double vmag = std::sqrt(square(vx) + square(vy) + square(vz));
          if (vmag > 0.)
          {
            const Acts::Vector3 base_world(trk->get_x(), trk->get_y(), trk->get_z());
            const Acts::Vector3 ahead_world(
                trk->get_x() + vx / vmag,
                trk->get_y() + vy / vmag,
                trk->get_z() + vz / vmag);
            const Acts::Vector3 base_envelope =
                m_tGeometry->transformTpcWorldToEnvelope(base_world);
            const Acts::Vector3 ahead_envelope =
                m_tGeometry->transformTpcWorldToEnvelope(ahead_world);
            const Acts::Vector3 direction_envelope = ahead_envelope - base_envelope;
            const double direction_norm = direction_envelope.norm();
            if (direction_norm > 0.0)
            {
              auto &layer_data = m_track_layer_data[trkid];
              layer_data.base_x = base_envelope.x();
              layer_data.base_y = base_envelope.y();
              layer_data.base_z = base_envelope.z();
              layer_data.dir_x = direction_envelope.x() / direction_norm;
              layer_data.dir_y = direction_envelope.y() / direction_norm;
              layer_data.dir_z = direction_envelope.z() / direction_norm;
              layer_data.have_line = true;
            }
          }
        }
      }
    }

    // see if there is a jump in x or y relative to previous PHG4Hit
    if (truth_clusterer.b_collect_hits)
    {
      if (prior_g4hit)
      {
        // if the g4hits jump in x or y by > max_g4hit_jump, cluster the truth tracks
        if (std::abs(prior_g4hit->get_x(0) - hiter->second->get_x(0)) > max_g4hitstep || std::abs(prior_g4hit->get_y(0) - hiter->second->get_y(0)) > max_g4hitstep)
        {
          if (truth_track)
          {
            truth_clusterer.cluster_hits(truth_track);
          }
        }
      }
      prior_g4hit = hiter->second;
    }

    // for very high occupancy events, accessing the TrkrHitsets on the node tree
    // for every drifted electron seems to be very slow
    // Instead, use a temporary map to accumulate the charge from all
    // drifted electrons, then copy to the node tree later

    const double dx_hit = hiter->second->get_x(1) - hiter->second->get_x(0);
    const double dy_hit = hiter->second->get_y(1) - hiter->second->get_y(0);
    const double dz_hit = hiter->second->get_z(1) - hiter->second->get_z(0);
    const double step_length = std::sqrt(square(dx_hit) + square(dy_hit) + square(dz_hit));

    const Acts::Vector3 hit_start_world(
        hiter->second->get_x(0),
        hiter->second->get_y(0),
        hiter->second->get_z(0));
    const Acts::Vector3 hit_end_world(
        hiter->second->get_x(1),
        hiter->second->get_y(1),
        hiter->second->get_z(1));
    const Acts::Vector3 hit_start_envelope =
        m_tGeometry->transformTpcWorldToEnvelope(hit_start_world);
    const Acts::Vector3 hit_end_envelope =
        m_tGeometry->transformTpcWorldToEnvelope(hit_end_world);
    const Acts::Vector3 hit_delta_envelope = hit_end_envelope - hit_start_envelope;

    const int track_id = hiter->second->get_trkid();
    double track_segment_start = m_track_path_offset[track_id];

    const double eion = hiter->second->get_eion();
    const double edep = hiter->second->get_edep();
    // Keep upstream's opt-in path-length model and established energy-deposit
    // fallback, together with the local explicit PAI/uniform test modes.
    bool use_primary_clusters_for_hit = false;
    double poisson_mean = 0.0;
    if (m_use_pai_cluster_seeds)
    {
      poisson_mean = hiter->second->get_index_i() > 0 ? 1.0 : 0.0;
    }
    else if (m_use_primary_cluster_ionization)
    {
      const double path_length_cm = hiter->second->get_path_length();
      if (!std::isfinite(path_length_cm))
      {
        std::cerr << Name() << ": G4 hit " << hiter->first
                  << " lacks path_length required by primary-cluster ionization; regenerate the G4 input"
                  << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
      use_primary_clusters_for_hit = path_length_cm > 0.0;
      poisson_mean = use_primary_clusters_for_hit
                         ? (fixed_primary_clusters_per_cm > 0.0
                                ? fixed_primary_clusters_per_cm * path_length_cm
                                : m_primaryIonizationModel->mean_primary_clusters(path_length_cm))
                         : std::max(0.0, eion * electrons_per_gev);
    }
    else
    {
      poisson_mean = std::max(0.0, eion * electrons_per_gev);
    }
    unsigned int n_primary_clusters = 0;
    if (m_use_pai_cluster_seeds || m_uniform_density_test)
    {
      n_primary_clusters = static_cast<unsigned int>(std::lround(poisson_mean));
    }
    else
    {
      n_primary_clusters = gsl_ran_poisson(RandomGenerator.get(), poisson_mean);
    }
    // Source model: 0 energy-deposit electron, 1 path-length primary
    // cluster, 2 explicit PAI seed. Keep physical-primary QA distinct.
    const int source_model = m_use_pai_cluster_seeds ? 2 : (use_primary_clusters_for_hit ? 1 : 0);
    const bool is_primary_source = source_model != 0;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double dedx_kev_per_cm = (step_length > 0.) ? (edep * 1e6 / step_length) : nan;
    const double deiondx_kev_per_cm = (step_length > 0.) ? (eion * 1e6 / step_length) : nan;

    auto fill_step_qa = [&](const int n_ionization_electrons)
    {
      if (!fill_qa_hists_this_event)
      {
        return;
      }

      const double mean_ionization_electrons = m_use_pai_cluster_seeds
          ? std::max(0, hiter->second->get_index_i())
          : poisson_mean * (use_primary_clusters_for_hit && m_enable_cluster_size_fluctuations
                                ? m_primaryIonizationModel->mean_cluster_size() : 1.0);
      const double ionized_electrons_edep_kev = std::isfinite(avg_ionization_energy_kev_per_electron) ? n_ionization_electrons * avg_ionization_energy_kev_per_electron : nan;
      const double ionized_electrons_dedx_kev_per_cm =
          (step_length > 0. && std::isfinite(ionized_electrons_edep_kev)) ? ionized_electrons_edep_kev / step_length : nan;

      if (primaryClusterMean && is_primary_source)
      {
        primaryClusterMean->Fill(poisson_mean);
      }
      if (nPrimaryClusters && is_primary_source)
      {
        nPrimaryClusters->Fill(static_cast<double>(n_primary_clusters));
      }
      if (nPrimaryClustersPerCm && is_primary_source && step_length > 0.)
      {
        nPrimaryClustersPerCm->Fill(static_cast<double>(n_primary_clusters) / step_length);
      }
      if (nPrimaryClustersVsMean && is_primary_source)
      {
        nPrimaryClustersVsMean->Fill(poisson_mean, static_cast<double>(n_primary_clusters));
      }
      if (nIonizationElectrons)
      {
        nIonizationElectrons->Fill(static_cast<double>(n_ionization_electrons));
      }
      if (nIonizationElectronsPerCm && step_length > 0.)
      {
        nIonizationElectronsPerCm->Fill(static_cast<double>(n_ionization_electrons) / step_length);
      }
      if (nIonizationElectronsVsPrimaryClusters && is_primary_source)
      {
        nIonizationElectronsVsPrimaryClusters->Fill(static_cast<double>(n_primary_clusters), static_cast<double>(n_ionization_electrons));
      }
      if (g4StepDedx && std::isfinite(dedx_kev_per_cm))
      {
        g4StepDedx->Fill(dedx_kev_per_cm);
      }
      if (ionizedElectronsDedx && std::isfinite(ionized_electrons_dedx_kev_per_cm))
      {
        ionizedElectronsDedx->Fill(ionized_electrons_dedx_kev_per_cm);
      }
      if (ionizedElectronsDedxVsG4StepDedx && std::isfinite(dedx_kev_per_cm) && std::isfinite(ionized_electrons_dedx_kev_per_cm))
      {
        ionizedElectronsDedxVsG4StepDedx->Fill(dedx_kev_per_cm, ionized_electrons_dedx_kev_per_cm);
      }

      if (driftStepQA)
      {
        m_stepqa_event = event_num;
        m_stepqa_track_id = track_id;
        m_stepqa_step_length_cm = static_cast<float>(step_length);
        m_stepqa_edep_kev = static_cast<float>(edep * 1e6);
        m_stepqa_eion_kev = static_cast<float>(eion * 1e6);
        m_stepqa_dedx_kev_per_cm = static_cast<float>(dedx_kev_per_cm);
        m_stepqa_deiondx_kev_per_cm = static_cast<float>(deiondx_kev_per_cm);
        m_stepqa_source_model = source_model;
        m_stepqa_n_sources = static_cast<int>(n_primary_clusters);
        m_stepqa_mean_primary_clusters = is_primary_source ? static_cast<float>(poisson_mean) : 0.F;
        m_stepqa_n_primary_clusters = is_primary_source ? static_cast<int>(n_primary_clusters) : 0;
        m_stepqa_avg_ionization_energy_kev_per_electron = static_cast<float>(avg_ionization_energy_kev_per_electron);
        m_stepqa_mean_ionization_electrons = static_cast<float>(mean_ionization_electrons);
        m_stepqa_n_ionization_electrons = n_ionization_electrons;
        m_stepqa_ionized_electrons_edep_kev = static_cast<float>(ionized_electrons_edep_kev);
        m_stepqa_ionized_electrons_dedx_kev_per_cm = static_cast<float>(ionized_electrons_dedx_kev_per_cm);
        driftStepQA->Fill();
      }
    };
    if (m_avg_x_enabled)
    {
      auto &layer_data = m_track_layer_data[trkid_new];
      if (!layer_data.have_line)
      {
        const double norm = hit_delta_envelope.norm();
        if (norm > 0.)
        {
          layer_data.base_x = hit_start_envelope.x();
          layer_data.base_y = hit_start_envelope.y();
          layer_data.base_z = hit_start_envelope.z();
          layer_data.dir_x = hit_delta_envelope.x() / norm;
          layer_data.dir_y = hit_delta_envelope.y() / norm;
          layer_data.dir_z = hit_delta_envelope.z() / norm;
          layer_data.have_line = true;
        }
      }
    }

    if (m_truth_intersection_hits && !m_tpc_layer_radii.empty() && step_length > 1e-9)
    {
      const double x0 = hit_start_envelope.x();
      const double y0 = hit_start_envelope.y();
      const double z0 = hit_start_envelope.z();
      const double dx_envelope = hit_delta_envelope.x();
      const double dy_envelope = hit_delta_envelope.y();
      const double dz_envelope = hit_delta_envelope.z();
      const double t0_hit = hiter->second->get_t(0);
      const double dt = hiter->second->get_t(1) - t0_hit;
      const double r0 = std::sqrt(square(x0) + square(y0));
      const double r1 = std::hypot(hit_end_envelope.x(), hit_end_envelope.y());
      const double rmin = std::min(r0, r1);
      const double rmax = std::max(r0, r1);
      const double a = square(dx_envelope) + square(dy_envelope);
      const double b = 2.0 * (dx_envelope * x0 + dy_envelope * y0);
      const double c0 = square(x0) + square(y0);
      static constexpr double radial_tol = 5e-3;
      static constexpr double t_tol = 1e-6;

      if (a > 1e-12)
      {
        for (const auto &layer_info : m_tpc_layer_radii)
        {
          const double radius = layer_info.radius;
          if (radius < rmin - radial_tol || radius > rmax + radial_tol)
          {
            continue;
          }

          const double c = c0 - square(radius);
          double disc = b * b - 4.0 * a * c;
          if (disc < -1e-10)
          {
            continue;
          }

          disc = std::max(0.0, disc);
          const double sqrt_disc = std::sqrt(disc);
          const std::array<double, 2> t_candidates = {
              (-b - sqrt_disc) / (2.0 * a),
              (-b + sqrt_disc) / (2.0 * a)};

          bool have_candidate = false;
          double best_t = 0.0;
          double best_residual = std::numeric_limits<double>::max();

          for (double t_candidate : t_candidates)
          {
            if (t_candidate < -t_tol || t_candidate > 1.0 + t_tol)
            {
              continue;
            }
            const double t = std::clamp(t_candidate, 0.0, 1.0);
            const double xi = x0 + t * dx_envelope;
            const double yi = y0 + t * dy_envelope;
            const double residual = std::abs(std::sqrt(square(xi) + square(yi)) - radius);
            if (!have_candidate || residual < best_residual)
            {
              best_t = t;
              best_residual = residual;
              have_candidate = true;
            }
          }

          if (!have_candidate)
          {
            continue;
          }

          const Acts::Vector3 intersection_envelope(
              x0 + best_t * dx_envelope,
              y0 + best_t * dy_envelope,
              z0 + best_t * dz_envelope);
          const Acts::Vector3 intersection_world =
              m_tGeometry->transformTpcEnvelopeToWorld(intersection_envelope);
          const double ti = t0_hit + best_t * dt;
          const double path = track_segment_start + best_t * step_length;

          TruthLayerIntersection intersection;
          intersection.track_id = track_id;
          intersection.layer = layer_info.layer;
          intersection.x = intersection_world.x();
          intersection.y = intersection_world.y();
          intersection.z = intersection_world.z();
          intersection.t = ti;
          // intersection.dirx = dx_hit;
          // intersection.diry = dy_hit;
          // intersection.dirz = dz_hit;
          intersection.path = path;
          // intersection.eion = hiter->second->get_eion();
          // intersection.edep = hiter->second->get_edep();
          const std::pair<int, unsigned int> key(track_id, layer_info.layer);
          truth_layer_intersections[key].push_back(intersection);
        }
      }
    }

    if (Verbosity() > 100)
    {
      std::cout << "  new hit with t0, " << t0 << " g4hitid " << hiter->first
                << " eion " << eion << " step length " << step_length
                << " poisson mean " << poisson_mean
                << " n_primary_clusters " << n_primary_clusters
                << " entry z " << hiter->second->get_z(0) << " exit z "
                << hiter->second->get_z(1) << " avg z"
                << (hiter->second->get_z(0) + hiter->second->get_z(1)) / 2.0
                << std::endl;
    }

    if (Verbosity() > 0 && count_g4hits <= 5)
    {
      std::cout << Name() << " hit " << count_g4hits
                << " eion " << eion
                << " step length " << step_length
                << " poisson mean " << poisson_mean
                << " sampled primary clusters " << n_primary_clusters
                << " track dz " << hiter->second->get_z(1) - hiter->second->get_z(0)
                << " entry radius " << std::sqrt(square(hiter->second->get_x(0)) + square(hiter->second->get_y(0)))
                << " exit radius " << std::sqrt(square(hiter->second->get_x(1)) + square(hiter->second->get_y(1)))
                << std::endl;
    }

    if (n_primary_clusters == 0)
    {
      fill_step_qa(0);
      m_track_path_offset[track_id] += step_length;
      continue;
    }

    if (Verbosity() > 100)
    {
      std::cout << std::endl
                << "electron drift: g4hit " << hiter->first << " created primary clusters: "
                << n_primary_clusters << " over step length " << step_length << " cm"
                << " (eion " << eion * 1000000 << " keV)" << std::endl;
      std::cout << " entry x,y,z = " << hiter->second->get_x(0) << "  "
                << hiter->second->get_y(0) << "  " << hiter->second->get_z(0)
                << " radius " << sqrt(pow(hiter->second->get_x(0), 2) + pow(hiter->second->get_y(0), 2)) << std::endl;
      std::cout << " exit x,y,z = " << hiter->second->get_x(1) << "  "
                << hiter->second->get_y(1) << "  " << hiter->second->get_z(1)
                << " radius " << sqrt(pow(hiter->second->get_x(1), 2) + pow(hiter->second->get_y(1), 2)) << std::endl;
    }

    int notReachingReadout = 0;
    int n_secondary_electrons = 0;
    //    int notInAcceptance = 0;

    // Loop over primary ionization clusters
    for (unsigned int i = 0; i < n_primary_clusters; i++)
    {
      // Sample cluster size
      int cluster_size = m_use_pai_cluster_seeds ? hiter->second->get_index_i() : 1;
      if (!m_use_pai_cluster_seeds && use_primary_clusters_for_hit &&
          m_enable_cluster_size_fluctuations)
      {
        cluster_size = static_cast<int>(m_primaryIonizationModel->sample_cluster_size(
            gsl_ran_flat(RandomGenerator.get(), 0.0, 1.0)));
      }
      if (cluster_size <= 0)
      {
        continue;
      }

      // We choose the electron starting position at random from a flat
      // distribution along the path length the parameter t is the fraction of
      // the distance along the path betwen entry and exit points, it has
      // values between 0 and 1
      double f = gsl_ran_flat(RandomGenerator.get(), 0.0, 1.0);
      if (m_use_pai_cluster_seeds)
      {
        f = 0.0;
      }
      else if (m_uniform_density_test && n_primary_clusters > 0)
      {
        f = (static_cast<double>(i) + 0.5) / static_cast<double>(n_primary_clusters);
        if (f >= 1.0)
        {
          f = std::nextafter(1.0, 0.0);
        }
      }

      const double cluster_x = m_use_pai_cluster_seeds
                                   ? hiter->second->get_x(0)
                                   : hiter->second->get_x(0) + f * (hiter->second->get_x(1) - hiter->second->get_x(0));
      const double cluster_y = m_use_pai_cluster_seeds
                                   ? hiter->second->get_y(0)
                                   : hiter->second->get_y(0) + f * (hiter->second->get_y(1) - hiter->second->get_y(0));
      const double cluster_z = m_use_pai_cluster_seeds
                                   ? hiter->second->get_z(0)
                                   : hiter->second->get_z(0) + f * (hiter->second->get_z(1) - hiter->second->get_z(0));
      const double cluster_t = m_use_pai_cluster_seeds
                                   ? hiter->second->get_t(0)
                                   : hiter->second->get_t(0) + f * (hiter->second->get_t(1) - hiter->second->get_t(0));
      const Acts::Vector3 cluster_world(cluster_x, cluster_y, cluster_z);
      const Acts::Vector3 cluster_envelope =
          m_tGeometry->transformTpcWorldToEnvelope(cluster_world);
      const double fallback_track_path_cm = track_segment_start + f * step_length;
      const double cluster_track_path_cm = m_use_pai_cluster_seeds
                                               ? path_from_position(track_id, hiter->second, fallback_track_path_cm)
                                               : clamp_track_path_cm(track_id, fallback_track_path_cm);
      const int cluster_segment_index = add_cluster_to_segment(track_id, cluster_track_path_cm, cluster_size, is_primary_source);

      if (fill_qa_hists_this_event)
      {
        if (is_primary_source) ++m_eventqa_n_primary_clusters;
        else ++m_eventqa_n_energy_sources;
        m_eventqa_n_ionization_electrons += cluster_size;
        if (ionizationClusterQA)
        {
          m_clusterqa_event = event_num;
          m_clusterqa_use_pai_cluster_seeds = m_use_pai_cluster_seeds ? 1 : 0;
          m_clusterqa_source_model = source_model;
          m_clusterqa_track_id = track_id;
          m_clusterqa_g4hit_id = static_cast<int>(hiter->first);
          m_clusterqa_cluster_index = static_cast<int>(i);
          m_clusterqa_segment_index = cluster_segment_index;
          m_clusterqa_cluster_size = cluster_size;
          m_clusterqa_track_path_cm = static_cast<float>(cluster_track_path_cm);
          m_clusterqa_x_cm = static_cast<float>(cluster_x);
          m_clusterqa_y_cm = static_cast<float>(cluster_y);
          m_clusterqa_z_cm = static_cast<float>(cluster_z);
          m_clusterqa_t_ns = static_cast<float>(cluster_t);
          ionizationClusterQA->Fill();
        }
      }

      const double x_start = cluster_envelope.x();
      const double y_start = cluster_envelope.y();
      const double z_start = cluster_envelope.z();
      const double t_start = cluster_t;
      const unsigned int side = z_start > 0.0 ? 1U : 0U;
      const double radstart = std::hypot(x_start, y_start);
      const double phistart = std::atan2(y_start, x_start);

      double drift_distance = tpc_length / 2.0 - std::abs(z_start);
      if (drift_distance < 0.0)
      {
        drift_distance = 0.0;
      }

      // The deterministic part of the transport is common to all secondary
      // electrons in this primary cluster. In Garfield mode, calculate it once
      // here and retain the existing stochastic diffusion inside the loop.
      double deterministic_x = x_start;
      double deterministic_y = y_start;
      double deterministic_t_path = drift_distance / drift_velocity;
      const bool use_garfield_full_lut =
          m_distortionMap &&
          m_distortionMap->application_mode() ==
              PHG4TpcDistortion::ApplicationMode::GarfieldFullTransport;
      if (use_garfield_full_lut)
      {
        PHG4TpcDistortion::GarfieldTransportSample sample;
        const auto query_status =
            m_distortionMap->query_garfield_transport(
                radstart, phistart, z_start, sample);
        if (query_status ==
            PHG4TpcDistortion::GarfieldTransportQueryStatus::OutsideDriftVolume)
        {
          // A Geant4 boundary point or the straight chord used to sample a
          // curved step can lie just outside the physical gas volume. Such a
          // source cannot drift to the readout, but it must not abort the rest
          // of the event.
          n_secondary_electrons += cluster_size;
          notReachingReadout += cluster_size;
          continue;
        }
        if (query_status !=
            PHG4TpcDistortion::GarfieldTransportQueryStatus::Success)
        {
          std::cerr
              << Name()
              << ": invalid Garfield full-transport LUT "
                 "query at (r,phi,z)=("
              << radstart << ", " << phistart << ", " << z_start << ")"
              << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }

        if (sample.reach !=
            PHG4TpcDistortion::ReachCellClass::AllReach)
        {
          n_secondary_electrons += cluster_size;
          notReachingReadout += cluster_size;
          continue;
        }

        if (radstart <= 0.0)
        {
          std::cerr << Name()
                    << ": invalid zero-radius Garfield LUT query"
                    << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
        const double deterministic_radius =
            radstart + sample.delta_r_cm;
        const double deterministic_phi =
            phistart + sample.delta_rphi_cm / radstart;
        deterministic_x =
            deterministic_radius * std::cos(deterministic_phi);
        deterministic_y =
            deterministic_radius * std::sin(deterministic_phi);
        deterministic_t_path = sample.drift_time_ns;
      }

      const double side_space_charge_scale =
          side == 0U ? (m_garfield ? m_garfield->GetSpaceChargeScaleSide0() : 0.0)
                     : (m_garfield ? m_garfield->GetSpaceChargeScaleSide1() : 0.0);
      if (m_garfield &&
          (m_garfield_full_transport ||
           (m_garfield->HasElectricFieldCorrections() &&
            side_space_charge_scale != 0.0)))
      {
        PHGarfield::ForwardDriftResult distorted_result;
        const bool distorted_ok = m_garfield->ForwardDrift(
            x_start,
            y_start,
            z_start,
            tpc_length / 2.0,
            m_garfield_step_ns,
            distorted_result,
            true);
        if (!distorted_ok || !distorted_result.reached_readout)
        {
          if (fatal_garfield_failure(distorted_result.status))
          {
            std::cerr << Name() << ": fatal PHGarfield forward-drift failure: "
                      << PHGarfield::ForwardDriftStatusName(
                             distorted_result.status)
                      << " at (" << x_start << ", " << y_start << ", "
                      << z_start << ") cm" << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
          }
          n_secondary_electrons += cluster_size;
          notReachingReadout += cluster_size;
          continue;
        }

        if (m_garfield_full_transport)
        {
          deterministic_x = distorted_result.x_cm;
          deterministic_y = distorted_result.y_cm;
          deterministic_t_path = distorted_result.drift_time_ns;
        }
        else
        {
          PHGarfield::ForwardDriftResult reference_result;
          const bool reference_ok = m_garfield->ForwardDrift(
              x_start,
              y_start,
              z_start,
              tpc_length / 2.0,
              m_garfield_step_ns,
              reference_result,
              false,
              false);
          if (!reference_ok || !reference_result.reached_readout)
          {
            if (fatal_garfield_failure(reference_result.status))
            {
              std::cerr
                  << Name()
                  << ": fatal PHGarfield scale-zero reference failure: "
                  << PHGarfield::ForwardDriftStatusName(
                         reference_result.status)
                  << " at (" << x_start << ", " << y_start << ", "
                  << z_start << ") cm" << std::endl;
              return Fun4AllReturnCodes::ABORTRUN;
            }
            n_secondary_electrons += cluster_size;
            notReachingReadout += cluster_size;
            continue;
          }

          // Inject only the field-distortion delta into PHG4's calibrated
          // nominal clock and endpoint. This makes scale=0 an exact no-op and
          // matches the contract of the existing displacement-map mode.
          deterministic_x += distorted_result.x_cm - reference_result.x_cm;
          deterministic_y += distorted_result.y_cm - reference_result.y_cm;
          deterministic_t_path +=
              distorted_result.drift_time_ns - reference_result.drift_time_ns;
        }

        if (!std::isfinite(deterministic_x) ||
            !std::isfinite(deterministic_y) ||
            !std::isfinite(deterministic_t_path) ||
            deterministic_t_path < 0.0)
        {
          n_secondary_electrons += cluster_size;
          notReachingReadout += cluster_size;
          continue;
        }
      }

      // Loop over secondary electrons in the cluster
      for (int j = 0; j < cluster_size; ++j)
      {
        ++n_secondary_electrons;
        // Diffusion coefficients specify the width in each transverse axis.
        // A single Gaussian radius at a random angle would make both the x
        // and y widths smaller by sqrt(2), and would not give a 2D Gaussian.
        const double transverse_sigma = diffusion_trans * std::sqrt(drift_distance);
        const double diffusion_x = gsl_ran_gaussian(RandomGenerator.get(), transverse_sigma);
        const double diffusion_y = gsl_ran_gaussian(RandomGenerator.get(), transverse_sigma);
        const double added_smear_x = gsl_ran_gaussian(RandomGenerator.get(), added_smear_sigma_trans);
        const double added_smear_y = gsl_ran_gaussian(RandomGenerator.get(), added_smear_sigma_trans);

        const double t_path = deterministic_t_path;
        const double t_sigma = diffusion_long * std::sqrt(drift_distance) / drift_velocity;
        const double rantime =
            gsl_ran_gaussian(RandomGenerator.get(), t_sigma) +
            gsl_ran_gaussian(RandomGenerator.get(), added_smear_sigma_long) / drift_velocity;
        double t_final = t_start + t_path + rantime;

        if (t_final < min_time || t_final > max_time)
        {
          continue;
        }

        double z_final;
        if (z_start < 0)
        {
          z_final = -tpc_length / 2. + t_final * drift_velocity;
        }
        else
        {
          z_final = tpc_length / 2. - t_final * drift_velocity;
        }

        if (m_avg_x_enabled && radstart >= m_avg_layer_rlow && radstart <= m_avg_layer_rhigh)
        {
          // collect undiffused electron start positions for the configured layer
          auto &layer_data = m_track_layer_data[track_id];
          layer_data.sum_rphi += m_avg_layer_radius * phistart;
          layer_data.sum_x += x_start;
          layer_data.count++;
          if (j == 0)
          {
            layer_data.sum_rphi_primary += m_avg_layer_radius * phistart;
            layer_data.sum_x_primary += x_start;
            layer_data.count_primary++;
          }
        }
        double x_final = deterministic_x + diffusion_x + added_smear_x;
        double y_final = deterministic_y + diffusion_y + added_smear_y;

        if (force_min_trans_drift_length > drift_distance)
        {
          const double additional_length = force_min_trans_drift_length - drift_distance;
          if (additional_length > 0.)
          {
            const double extra_sigma = diffusion_trans * std::sqrt(additional_length);
            x_final += gsl_ran_gaussian(RandomGenerator.get(), extra_sigma);
            y_final += gsl_ran_gaussian(RandomGenerator.get(), extra_sigma);
          }
        }

        // These are diffusion-only displacements even when the deterministic
        // starting point at the readout was supplied by Garfield.
        const double delta_x = x_final - deterministic_x;
        const double delta_y = y_final - deterministic_y;
        const double rantrans = std::sqrt(square(delta_x) + square(delta_y));

        double rad_final = sqrt(square(x_final) + square(y_final));
        double phi_final = atan2(y_final, x_final);
        if (m_avg_x_enabled && rad_final >= m_avg_layer_rlow && rad_final <= m_avg_layer_rhigh)
        {
          auto &layer_data = m_track_layer_data[track_id];
          // Use the nominal layer radius so residuals are purely angular and not
          // inflated by radius fluctuations from diffusion.
          layer_data.sum_rphi_end += m_avg_layer_radius * phi_final;
          layer_data.sum_x_end += x_final;
          layer_data.count_end++;
          if (j == 0)
          {
            layer_data.sum_rphi_end_primary += m_avg_layer_radius * phi_final;
            layer_data.sum_x_end_primary += x_final;
            layer_data.count_end_primary++;
          }
        }

        if (fill_qa_hists_this_event)
        {
          if (diffDistance)
          {
            diffDistance->Fill(rantrans);
          }
          if (diffDX)
          {
            diffDX->Fill(delta_x);
          }
          if (diffDY)
          {
            diffDY->Fill(delta_y);
          }
          if (diffVsDrift)
          {
            diffVsDrift->Fill(drift_distance, rantrans);
          }
          if (diffDXVsDrift)
          {
            diffDXVsDrift->Fill(drift_distance, delta_x);
          }
          if (diffDYVsDrift)
          {
            diffDYVsDrift->Fill(drift_distance, delta_y);
          }

          if (drift_distance > 0.)
          {
            const double sqrtL = std::sqrt(drift_distance);
            const double norm = rantrans / sqrtL;
            if (diffPerSqrtL)
            {
              diffPerSqrtL->Fill(norm);
            }
            if (diffPerSqrtLVsDrift)
            {
              diffPerSqrtLVsDrift->Fill(drift_distance, norm);
            }
            const double dxnorm = delta_x / sqrtL;
            const double dynorm = delta_y / sqrtL;
            if (diffDXPerSqrtL)
            {
              diffDXPerSqrtL->Fill(dxnorm);
            }
            if (diffDYPerSqrtL)
            {
              diffDYPerSqrtL->Fill(dynorm);
            }
            if (diffDXPerSqrtLVsDrift)
            {
              diffDXPerSqrtLVsDrift->Fill(drift_distance, dxnorm);
            }
            if (diffDYPerSqrtLVsDrift)
            {
              diffDYPerSqrtLVsDrift->Fill(drift_distance, dynorm);
            }
          }
          z_startmap->Fill(z_start, radstart);                   // map of starting location in Z vs. R
          deltaphinodist->Fill(phistart, rantrans / rad_final);  // delta phi no distortion, just diffusion+smear
          deltarnodist->Fill(radstart, rantrans);                // delta r no distortion, just diffusion+smear
        }

        if (m_distortionMap && !use_garfield_full_lut)
        {
          // zhangcanyu
          const double reaches = m_distortionMap->get_reaches_readout(radstart, phistart, z_start);
          if (reaches < thresholdforreachesreadout)
          {
            notReachingReadout++;
            continue;
          }

          const double r_distortion = m_distortionMap->get_r_distortion(radstart, phistart, z_start);
          const double phi_distortion = m_distortionMap->get_rphi_distortion(radstart, phistart, z_start) / radstart;
          const double z_distortion = m_distortionMap->get_z_distortion(radstart, phistart, z_start);

          rad_final += r_distortion;
          phi_final += phi_distortion;
          z_final += z_distortion;
          if (z_start < 0)
          {
            t_final = (z_final + tpc_length / 2.0) / drift_velocity;
          }
          else
          {
            t_final = (tpc_length / 2.0 - z_final) / drift_velocity;
          }

          x_final = rad_final * std::cos(phi_final);
          y_final = rad_final * std::sin(phi_final);

          //	if(i < 1)
          //{std::cout << " electron " << i << " r_distortion " << r_distortion << " phi_distortion " << phi_distortion << " rad_final " << rad_final << " phi_final " << phi_final << " r*dphi distortion " << rad_final * phi_distortion << " z_distortion " << z_distortion << std::endl;}

          if (fill_qa_hists_this_event)
          {
            const double phi_final_nodiff = phistart + phi_distortion;
            const double rad_final_nodiff = radstart + r_distortion;
            deltarnodiff->Fill(radstart, rad_final_nodiff - radstart);    // delta r no diffusion, just distortion
            deltaphinodiff->Fill(phistart, phi_final_nodiff - phistart);  // delta phi no diffusion, just distortion
            deltaphivsRnodiff->Fill(radstart, phi_final_nodiff - phistart);
            deltaRphinodiff->Fill(radstart, radstart * phi_distortion);
            deltar->Fill(radstart, rad_final - radstart);    // total delta r
            deltaphi->Fill(phistart, phi_final - phistart);  // total delta phi
            deltaz->Fill(z_start, z_distortion);             // map of distortion in Z (time)
          }
        }

        if ((m_garfield || use_garfield_full_lut) &&
            fill_qa_hists_this_event)
        {
          const double rad_final_nodiff =
              std::hypot(deterministic_x, deterministic_y);
          const double phi_final_nodiff =
              std::atan2(deterministic_y, deterministic_x);
          const double phi_distortion =
              std::remainder(phi_final_nodiff - phistart, 2.0 * M_PI);
          const double apparent_z_nodiff =
              z_start < 0.0
                  ? -tpc_length / 2.0 + t_path * drift_velocity
                  : tpc_length / 2.0 - t_path * drift_velocity;
          const double z_distortion = apparent_z_nodiff - z_start;

          deltarnodiff->Fill(radstart, rad_final_nodiff - radstart);
          deltaphinodiff->Fill(phistart, phi_distortion);
          deltaphivsRnodiff->Fill(radstart, phi_distortion);
          deltaRphinodiff->Fill(
              radstart,
              radstart * phi_distortion);
          deltar->Fill(radstart, rad_final - radstart);
          deltaphi->Fill(
              phistart,
              std::remainder(phi_final - phistart, 2.0 * M_PI));
          deltaz->Fill(z_start, z_distortion);
        }

        if (fill_qa_hists_this_event)
        {
          // Fill start/end maps regardless of whether distortion maps are enabled.
          hitmapstart->Fill(x_start, y_start);
          hitmapend->Fill(x_final, y_final);
          hitmapstart_z->Fill(z_start, radstart);
          hitmapend_z->Fill(z_final, rad_final);
        }

        if (fill_qa_hists_this_event && driftXY)
        {
          m_drift_start_x = static_cast<float>(x_start);
          m_drift_start_y = static_cast<float>(y_start);
          m_drift_end_x = static_cast<float>(x_final);
          m_drift_end_y = static_cast<float>(y_final);
          driftXY->Fill();
        }

        if (m_density_enabled && electronDensityProfile)
        {
          if (rad_final >= m_density_capture_rlow && rad_final <= m_density_capture_rhigh)
          {
            const double delta_r = radstart - m_density_layer_radius;
            if (delta_r >= -m_density_hist_half_range && delta_r <= m_density_hist_half_range)
            {
              const double s_global = cluster_track_path_cm;
              if (!m_track_anchor_set[track_id])
              {
                m_track_anchor_set[track_id] = true;
                m_track_anchor_length[track_id] = s_global;
              }
              electronDensityProfile->Fill(s_global - m_track_anchor_length[track_id]);
            }
          }
        }

        // remove electrons outside of our acceptance. Careful though, electrons from just inside 30 cm can contribute in the 1st active layer readout, so leave a little margin
        if (rad_final < min_active_radius - 2.0 || rad_final > max_active_radius + 1.0)
        {
          //        notInAcceptance++;
          continue;
        }

        if (Verbosity() > 1000)
        //      if(i < 1)
        {
          std::cout << "electron " << i << " g4hitid " << hiter->first << " f " << f << std::endl;
          std::cout << "radstart " << radstart << " x_start: " << x_start
                    << ", y_start: " << y_start
                    << ",z_start: " << z_start
                    << " t_start " << t_start
                    << " t_path " << t_path
                    << " t_sigma " << t_sigma
                    << " rantime " << rantime
                    << std::endl;

          std::cout << "       rad_final " << rad_final << " x_final " << x_final
                    << " y_final " << y_final
                    << " z_final " << z_final << " t_final " << t_final
                    << " zdiff " << z_final - z_start << std::endl;
        }

        if (Verbosity() > 0)
        {
          assert(nt);
          nt->Fill(ihit, t_start, t_final, t_sigma, rad_final, z_start, z_final);
        }
        padplane->MapToPadPlane(truth_clusterer, single_hitsetcontainer.get(),
                                temp_hitsetcontainer.get(), hittruthassoc, x_final, y_final, t_final,
                                side, hiter, ntpad, nthit);
      }  // end loop over secondary electrons
    }  // end loop over electrons for this g4hit

    m_track_path_offset[track_id] = track_segment_start + step_length;

    if (fill_qa_hists_this_event && n_secondary_electrons > 0)
    {
      ratioElectronsRR->Fill(static_cast<double>(n_secondary_electrons - notReachingReadout) / static_cast<double>(n_secondary_electrons));
    }
    fill_step_qa(n_secondary_electrons);

    TrkrHitSetContainer::ConstRange single_hitset_range = single_hitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
    for (TrkrHitSetContainer::ConstIterator single_hitset_iter = single_hitset_range.first;
         single_hitset_iter != single_hitset_range.second;
         ++single_hitset_iter)
    {
      // we have an itrator to one TrkrHitSet for the Tpc from the single_hitsetcontainer
      TrkrDefs::hitsetkey node_hitsetkey = single_hitset_iter->first;
      const unsigned int layer = TrkrDefs::getLayer(node_hitsetkey);
      const int sector = TpcDefs::getSectorId(node_hitsetkey);
      const int side = TpcDefs::getSide(node_hitsetkey);

      if (Verbosity() > 8)
      {
        std::cout << " hitsetkey " << node_hitsetkey << " layer " << layer << " sector " << sector << " side " << side << std::endl;
      }
      // get all of the hits from the single hitset
      TrkrHitSet::ConstRange single_hit_range = single_hitset_iter->second->getHits();
      for (TrkrHitSet::ConstIterator single_hit_iter = single_hit_range.first;
           single_hit_iter != single_hit_range.second;
           ++single_hit_iter)
      {
        TrkrDefs::hitkey single_hitkey = single_hit_iter->first;

        // Add the hit-g4hit association
        // no need to check for duplicates, since the hit is new
        hittruthassoc->addAssoc(node_hitsetkey, single_hitkey, hiter->first);
        if (Verbosity() > 100)
        {
          std::cout << "        adding assoc for node_hitsetkey " << node_hitsetkey << " single_hitkey " << single_hitkey << " g4hitkey " << hiter->first << std::endl;
        }
      }
    }

    if (dump_counter >= dump_interval)
    {
      flush_temp_hitsets();
      dump_counter = 0;
    }

    ++ihit;

    single_hitsetcontainer->Reset();

  }  // end loop over g4hits

  // A hit near the end of the event can be followed only by G4 hits rejected
  // by an early-continue condition above. Flush unconditionally so accepted
  // charge is neither dropped nor carried into the next event.
  flush_temp_hitsets();

  if (m_truth_intersection_hits)
  {
    std::size_t nwritten_truth_intersections = 0;
    for (const auto &entry : truth_layer_intersections)
    {
      auto points = entry.second;
      std::sort(points.begin(), points.end(),
                [](const TruthLayerIntersection &lhs, const TruthLayerIntersection &rhs)
                { return lhs.path < rhs.path; });

      const std::size_t n_to_write = std::min(points.size(), kMaxTruthIntersectionsPerTrackLayer);
      for (std::size_t i = 0; i < n_to_write; ++i)
      {
        const auto &point = points[i];
        auto *intersection_hit = new PHG4Hitv1();
        intersection_hit->set_trkid(point.track_id);
        intersection_hit->set_layer(point.layer);
        intersection_hit->set_x(0, point.x);
        intersection_hit->set_y(0, point.y);
        intersection_hit->set_z(0, point.z);
        intersection_hit->set_t(0, point.t);
        intersection_hit->set_x(1, point.x);
        intersection_hit->set_y(1, point.y);
        intersection_hit->set_z(1, point.z);
        intersection_hit->set_t(1, point.t);
        // Optional direction/energy metadata intentionally omitted.
        // intersection_hit->set_px(0, point.dirx);
        // intersection_hit->set_py(0, point.diry);
        // intersection_hit->set_pz(0, point.dirz);
        // intersection_hit->set_px(1, point.dirx);
        // intersection_hit->set_py(1, point.diry);
        // intersection_hit->set_pz(1, point.dirz);
        // intersection_hit->set_edep(point.edep);
        // intersection_hit->set_eion(point.eion);
        intersection_hit->set_path_length(point.path);
        m_truth_intersection_hits->AddHit(point.layer, intersection_hit);
        ++nwritten_truth_intersections;
      }
    }

    if (Verbosity() > 1)
    {
      std::cout << Name() << ": wrote " << nwritten_truth_intersections
                << " truth layer intersections to " << kTpcTruthIntersectionNodeName << std::endl;
    }
  }

  if (m_avg_x_enabled && m_avgXResidualTree)
  {
    for (const auto &entry : m_track_layer_data)
    {
      const auto &data = entry.second;
      if (!data.have_line || data.count == 0)
      {
        continue;
      }
      double t_out = 0.0;
      double xi = 0.0;
      double yi = 0.0;
      double zi = 0.0;
      if (!cylinder_intersection(data.base_x, data.base_y, data.base_z,
                                 data.dir_x, data.dir_y, data.dir_z,
                                 m_avg_layer_radius, t_out, xi, yi, zi))
      {
        continue;
      }
      const double avg_rphi = data.sum_rphi / static_cast<double>(data.count);
      const double rphi_int = m_avg_layer_radius * std::atan2(yi, xi);
      m_avgTree_trackid = static_cast<float>(entry.first);
      m_avgTree_residual = static_cast<float>(avg_rphi - rphi_int);
      m_avgTree_avgx = static_cast<float>(avg_rphi);
      m_avgTree_xint = static_cast<float>(rphi_int);
      const double avg_x_cart = data.sum_x / static_cast<double>(data.count);
      m_avgTree_residual_x = static_cast<float>(avg_x_cart - xi);
      m_avgTree_avgx_cart = static_cast<float>(avg_x_cart);
      if (data.count_primary > 0)
      {
        const double avg_rphi_primary = data.sum_rphi_primary / static_cast<double>(data.count_primary);
        m_avgTree_residual_primary = static_cast<float>(avg_rphi_primary - rphi_int);
        m_avgTree_avgx_primary = static_cast<float>(avg_rphi_primary);
        m_avgTree_count_primary = static_cast<float>(data.count_primary);
        const double avg_x_primary = data.sum_x_primary / static_cast<double>(data.count_primary);
        m_avgTree_residual_x_primary = static_cast<float>(avg_x_primary - xi);
        m_avgTree_avgx_cart_primary = static_cast<float>(avg_x_primary);
      }
      else
      {
        m_avgTree_residual_primary = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_avgx_primary = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_count_primary = 0.0f;
        m_avgTree_residual_x_primary = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_avgx_cart_primary = std::numeric_limits<float>::quiet_NaN();
      }
      if (data.count_end > 0)
      {
        const double avg_rphi_end = data.sum_rphi_end / static_cast<double>(data.count_end);
        m_avgTree_residual_end = static_cast<float>(avg_rphi_end - rphi_int);
        m_avgTree_avgx_end = static_cast<float>(avg_rphi_end);
        m_avgTree_count_end = static_cast<float>(data.count_end);
        const double avg_x_end = data.sum_x_end / static_cast<double>(data.count_end);
        m_avgTree_residual_x_end = static_cast<float>(avg_x_end - xi);
        m_avgTree_avgx_cart_end = static_cast<float>(avg_x_end);
      }
      else
      {
        m_avgTree_residual_end = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_avgx_end = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_count_end = 0.0f;
        m_avgTree_residual_x_end = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_avgx_cart_end = std::numeric_limits<float>::quiet_NaN();
      }
      if (data.count_end_primary > 0)
      {
        const double avg_rphi_end_primary = data.sum_rphi_end_primary / static_cast<double>(data.count_end_primary);
        m_avgTree_residual_end_primary = static_cast<float>(avg_rphi_end_primary - rphi_int);
        m_avgTree_avgx_end_primary = static_cast<float>(avg_rphi_end_primary);
        m_avgTree_count_end_primary = static_cast<float>(data.count_end_primary);
        const double avg_x_end_primary = data.sum_x_end_primary / static_cast<double>(data.count_end_primary);
        m_avgTree_residual_x_end_primary = static_cast<float>(avg_x_end_primary - xi);
        m_avgTree_avgx_cart_end_primary = static_cast<float>(avg_x_end_primary);
      }
      else
      {
        m_avgTree_residual_end_primary = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_avgx_end_primary = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_count_end_primary = 0.0f;
        m_avgTree_residual_x_end_primary = std::numeric_limits<float>::quiet_NaN();
        m_avgTree_avgx_cart_end_primary = std::numeric_limits<float>::quiet_NaN();
      }
      m_avgXResidualTree->Fill();
    }
  }

  if (truth_track)
  {
    truth_clusterer.cluster_hits(truth_track);
  }
  truth_clusterer.clear_hitsetkey_cnt();

  if (Verbosity() > 20)
  {
    int evtseq = -1;
    if (auto *eh = findNode::getClass<EventHeader>(topNode, "EventHeader"))
    {
      evtseq = eh->get_EvtSequence();
    }
    std::cout << "From PHG4TpcElectronDrift: hitsetcontainer printout at end: evt=" << evtseq << std::endl;
    // We want all hitsets for the Tpc
    TrkrHitSetContainer::ConstRange hitset_range = hitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
    for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
         hitset_iter != hitset_range.second;
         ++hitset_iter)
    {
      // we have an itrator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
      TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
      const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
      // if (layer != print_layer)
      //{
      //  continue;
      // }
      const int sector = TpcDefs::getSectorId(hitsetkey);
      const int side = TpcDefs::getSide(hitsetkey);

      std::cout << "PHG4TpcElectronDrift: hitset with key: " << hitsetkey << " in layer " << layer << " with sector " << sector << " side " << side << std::endl;

      // get all of the hits from this hitset
      TrkrHitSet *hitset = hitset_iter->second;
      TrkrHitSet::ConstRange hit_range = hitset->getHits();

      for (TrkrHitSet::ConstIterator hit_iter = hit_range.first;
           hit_iter != hit_range.second;
           ++hit_iter)
      {
        TrkrDefs::hitkey hitkey = hit_iter->first;
        TrkrHit *tpchit = hit_iter->second;
        std::cout << "      hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey) << " z bin " << TpcDefs::getTBin(hitkey)
                  << "  energy " << tpchit->getEnergy() << std::endl;
      }
    }
  }

  if (Verbosity() > 1000)
  {
    std::cout << "From PHG4TpcElectronDrift: hittruthassoc dump:" << std::endl;
    hittruthassoc->identify();

    hittruthassoc->identify();
  }

  if (fill_qa_hists_this_event && ionizationSegmentQA)
  {
    for (const auto &segment_iter : ionization_segment_accumulators)
    {
      const int track_id = segment_iter.first.first;
      const int segment_index = segment_iter.first.second;
      const auto &segment = segment_iter.second;
      const double length_cm = segment.segment_length_cm;

      m_segmentqa_event = event_num;
      m_segmentqa_use_pai_cluster_seeds = m_use_pai_cluster_seeds ? 1 : 0;
      m_segmentqa_track_id = track_id;
      m_segmentqa_segment_index = segment_index;
      m_segmentqa_segment_start_cm = static_cast<float>(segment.segment_start_cm);
      m_segmentqa_segment_length_cm = static_cast<float>(length_cm);
      m_segmentqa_n_primary_clusters = segment.n_primary_clusters;
      m_segmentqa_n_energy_sources = segment.n_energy_sources;
      m_segmentqa_n_ionization_electrons = segment.n_ionization_electrons;
      m_segmentqa_primary_clusters_per_cm = (length_cm > 0.0) ? static_cast<float>(segment.n_primary_clusters / length_cm) : std::numeric_limits<float>::quiet_NaN();
      m_segmentqa_ionization_electrons_per_cm = (length_cm > 0.0) ? static_cast<float>(segment.n_ionization_electrons / length_cm) : std::numeric_limits<float>::quiet_NaN();
      ionizationSegmentQA->Fill();
    }
  }

  if (fill_qa_hists_this_event && ionizationEventQA)
  {
    if (event_path_length_cm > 0.0)
    {
      m_eventqa_primary_clusters_per_cm = static_cast<float>(m_eventqa_n_primary_clusters / event_path_length_cm);
      m_eventqa_ionization_electrons_per_cm = static_cast<float>(m_eventqa_n_ionization_electrons / event_path_length_cm);
    }
    ionizationEventQA->Fill();
  }

  padplane->EndEvent(event_num);
  ++event_num;  // if doing more than one event, event_num will be incremented.

  if (Verbosity() > 500)
  {
    std::cout << " TruthTrackContainer results at end of event in PHG4TpcElectronDrift::process_event " << std::endl;
    truthtracks->identify();
  }

  if (Verbosity() > 800)
  {
    truth_clusterer.print(truthtracks);
    truth_clusterer.print_file(truthtracks, "drift_clusters.txt");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcElectronDrift::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    assert(m_outf);
    assert(nt);
    assert(ntpad);
    assert(nthit);
    assert(ntfinalhit);

    m_outf->cd();
    nt->Write();
    ntpad->Write();
    nthit->Write();
    ntfinalhit->Write();
    m_outf->Close();
  }
  auto write_ionization_qa_trees = [&]()
  {
    if (ionizationEventQA)
    {
      ionizationEventQA->Write();
    }
    if (ionizationClusterQA)
    {
      ionizationClusterQA->Write();
    }
    if (ionizationSegmentQA)
    {
      ionizationSegmentQA->Write();
    }
  };

  const bool skip_qa_output = m_qa_write_only_with_secondaries && !m_seen_event_with_secondaries;
  if (do_ElectronDriftQAHistos && m_write_legacy_qa_histograms)
  {
    if (skip_qa_output)
    {
      if (Verbosity() > 0)
      {
        std::cout << Name() << ": skip writing " << m_qa_output_file
                  << " (no secondary TPC g4hits found)." << std::endl;
      }
    }
    else
    {
      EDrift_outf.reset(new TFile(m_qa_output_file.c_str(), "recreate"));
      EDrift_outf->cd();
      deltar->Write();
      deltaphi->Write();
      deltaz->Write();
      deltarnodist->Write();
      deltaphinodist->Write();
      deltarnodiff->Write();
      deltaphinodiff->Write();
      deltaRphinodiff->Write();
      deltaphivsRnodiff->Write();
      hitmapstart->Write();
      hitmapend->Write();
      hitmapstart_z->Write();
      hitmapend_z->Write();
      z_startmap->Write();
      ratioElectronsRR->Write();
      if (diffDistance)
      {
        diffDistance->Write();
      }
      if (diffDX)
      {
        diffDX->Write();
      }
      if (diffDY)
      {
        diffDY->Write();
      }
      if (diffPerSqrtL)
      {
        diffPerSqrtL->Write();
      }
      if (nPrimaryClustersPerCm)
      {
        nPrimaryClustersPerCm->Write();
      }
      if (nIonizationElectrons)
      {
        nIonizationElectrons->Write();
      }
      if (nIonizationElectronsPerCm)
      {
        nIonizationElectronsPerCm->Write();
      }
      if (nIonizationElectronsVsPrimaryClusters)
      {
        nIonizationElectronsVsPrimaryClusters->Write();
      }
      if (g4StepDedx)
      {
        g4StepDedx->Write();
      }
      if (ionizedElectronsDedx)
      {
        ionizedElectronsDedx->Write();
      }
      if (ionizedElectronsDedxVsG4StepDedx)
      {
        ionizedElectronsDedxVsG4StepDedx->Write();
      }
      if (diffDXPerSqrtL)
      {
        diffDXPerSqrtL->Write();
      }
      if (diffDYPerSqrtL)
      {
        diffDYPerSqrtL->Write();
      }
      if (nPrimaryClusters)
      {
        nPrimaryClusters->Write();
      }
      if (primaryClusterMean)
      {
        primaryClusterMean->Write();
      }
      if (nPrimaryClustersVsMean)
      {
        nPrimaryClustersVsMean->Write();
      }
      if (diffVsDrift)
      {
        diffVsDrift->Write();
      }
      if (diffPerSqrtLVsDrift)
      {
        diffPerSqrtLVsDrift->Write();
      }
      if (diffDXVsDrift)
      {
        diffDXVsDrift->Write();
      }
      if (diffDYVsDrift)
      {
        diffDYVsDrift->Write();
      }
      if (diffDXPerSqrtLVsDrift)
      {
        diffDXPerSqrtLVsDrift->Write();
      }
      if (diffDYPerSqrtLVsDrift)
      {
        diffDYPerSqrtLVsDrift->Write();
      }
      if (driftXY)
      {
        driftXY->Write();
      }
      if (driftStepQA && m_step_qa_output_file.empty())
      {
        driftStepQA->Write();
      }
      if (m_ionization_qa_output_file.empty())
      {
        write_ionization_qa_trees();
      }
      if (electronDensityProfile)
      {
        electronDensityProfile->Write();
      }
      EDrift_outf->Close();
    }
  }
  if (do_ElectronDriftQAHistos && !m_ionization_qa_output_file.empty())
  {
    if (skip_qa_output)
    {
      if (Verbosity() > 0)
      {
        std::cout << Name() << ": skip writing " << m_ionization_qa_output_file
                  << " (no secondary TPC g4hits found)." << std::endl;
      }
    }
    else
    {
      m_ionizationQAOutf.reset(new TFile(m_ionization_qa_output_file.c_str(), "recreate"));
      m_ionizationQAOutf->cd();
      write_ionization_qa_trees();
      m_ionizationQAOutf->Close();
    }
  }
  if (do_ElectronDriftQAHistos && driftStepQA && !m_step_qa_output_file.empty())
  {
    if (skip_qa_output)
    {
      if (Verbosity() > 0)
      {
        std::cout << Name() << ": skip writing " << m_step_qa_output_file
                  << " (no secondary TPC g4hits found)." << std::endl;
      }
    }
    else
    {
      m_stepQAOutf.reset(new TFile(m_step_qa_output_file.c_str(), "recreate"));
      m_stepQAOutf->cd();
      driftStepQA->Write();
      m_stepQAOutf->Close();
    }
  }
  if (m_avgOutf && m_avgXResidualTree)
  {
    m_avgOutf->cd();
    m_avgXResidualTree->Write();
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": wrote avgXResidual tree to " << m_avg_output_file << std::endl;
    }
    m_avgOutf->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcElectronDrift::set_seed(const unsigned int seed)
{
  gsl_rng_set(RandomGenerator.get(), seed);
}

void PHG4TpcElectronDrift::SetDefaultParameters()
{
  // longitudinal diffusion for 50:50 Ne:CF4 is 0.012, transverse is 0.004, drift velocity is 0.008
  // longitudinal diffusion for 60:40 Ar:CF4 is 0.012, transverse is 0.004, drift velocity is 0.008 (chosen to be the same at 50/50 Ne:CF4)
  // longitudinal diffusion for 65:25:10 Ar:CF4:N2 is 0.013613, transverse is 0.005487, drift velocity is 0.006965
  // longitudinal diffusion for 75:20:05 Ar:CF4:i-C4H10 is 0.014596, transverse is 0.005313, drift velocity is 0.007550

  set_default_double_param("diffusion_long", 0.014596);   // cm/SQRT(cm)
  set_default_double_param("diffusion_trans", 0.005313);  // cm/SQRT(cm)
  set_default_double_param("Ne_frac", 0.00);
  set_default_double_param("Ar_frac", 0.75);
  set_default_double_param("CF4_frac", 0.20);
  set_default_double_param("N2_frac", 0.00);
  set_default_double_param("isobutane_frac", 0.05);
  set_default_double_param("min_active_radius", 30.);        // cm
  set_default_double_param("max_active_radius", 78.);        // cm

  // These are purely fudge factors, used to increase the resolution to 150 microns and 500 microns, respectively
  // override them from the macro to get a different resolution
  set_default_double_param("added_smear_trans", 0.0);             // cm (used to be 0.085 before sims got better)
  set_default_double_param("added_smear_long", 0.0);              // cm (used to be 0.105 before sims got better)
  set_default_double_param("force_min_trans_drift_length", 0.0);  // cm
  set_default_int_param("density_layer", -1);
  set_default_double_param("density_bin_width_cm", 0.05);
  set_default_double_param("density_window_cm", -1.0);
  set_default_double_param("density_radial_margin_cm", 0.2);
  set_default_double_param("fixed_primary_clusters_per_cm", -1.0);
  set_default_double_param("fixed_primary_electrons_per_cm", -1.0);  // backward-compatible alias
  set_default_int_param("uniform_density_test", 0);
  set_default_int_param("use_pai_cluster_seeds", 0);
  set_default_int_param("average_x_layer", -1);

  return;
}

bool PHG4TpcElectronDrift::cylinder_intersection(const double x0, const double y0, const double z0,
                                                 const double vx, const double vy, const double vz,
                                                 const double R, double &t_out,
                                                 double &xi, double &yi, double &zi) const
{
  const double a = square(vx) + square(vy);
  const double b = 2.0 * (vx * x0 + vy * y0);
  const double c = square(x0) + square(y0) - square(R);
  const double disc = b * b - 4.0 * a * c;
  if (a == 0.0 || disc < 0.0)
  {
    return false;
  }
  const double s = std::sqrt(disc);
  const double t1 = (-b + s) / (2.0 * a);
  const double t2 = (-b - s) / (2.0 * a);
  double t = std::numeric_limits<double>::max();
  if (t1 > 0.0)
  {
    t = std::min(t, t1);
  }
  if (t2 > 0.0)
  {
    t = std::min(t, t2);
  }
  if (!std::isfinite(t) || t <= 0.0)
  {
    return false;
  }
  t_out = t;
  xi = x0 + t * vx;
  yi = y0 + t * vy;
  zi = z0 + t * vz;
  return true;
}

void PHG4TpcElectronDrift::setTpcDistortion(PHG4TpcDistortion *distortionMap)
{
  m_distortionMap.reset(distortionMap);
}

void PHG4TpcElectronDrift::setPHGarfield(
    PHGarfield *garfield,
    const double step_ns)
{
  m_garfield.reset(garfield);
  m_garfield_step_ns = step_ns;
}

void PHG4TpcElectronDrift::registerPadPlane(PHG4TpcPadPlane *inpadplane)
{
  if (Verbosity())
  {
    std::cout << "Registering padplane " << std::endl;
  }
  padplane.reset(inpadplane);
  padplane->Detector(Detector());
  padplane->UpdateInternalParameters();
  if (Verbosity())
  {
    std::cout << "padplane registered and parameters updated" << std::endl;
  }

  return;
}

void PHG4TpcElectronDrift::set_flag_threshold_distortion(bool setflag, float setthreshold)
{
  std::cout << boost::str(boost::format("The logical status of threshold is now %d! and the value is set to %f") % setflag % setthreshold)
            << std::endl
            << std::endl
            << std::endl;
  do_getReachReadout = setflag;
  thresholdforreachesreadout = setthreshold;
}
