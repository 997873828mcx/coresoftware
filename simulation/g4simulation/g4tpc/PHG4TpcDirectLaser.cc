#include "PHG4TpcDirectLaser.h"

#include <phparameter/PHParameterInterface.h>  // for PHParameterInterface

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for get_volume_id
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Particlev3.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPointv1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrack_v2.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>  // for TVector3, operator*
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <phool/PHRandomSeed.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_const_mksa.h>  // for the speed of light

#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>  // for operator<<, basic_os...
#include <optional>
#include <cstdlib>  // for std::getenv

namespace
{
  using PHG4Particle_t = PHG4Particlev3;
  using PHG4VtxPoint_t = PHG4VtxPointv1;
  using PHG4Hit_t = PHG4Hitv1;

  // utility
  template <class T>
  constexpr T square(const T& x)
  {
    return x * x;
  }

  // unique detector id for all direct lasers
  const int detId = PHG4HitDefs::get_volume_id("PHG4TpcDirectLaser");

  ///@name units
  //@{
  constexpr double cm = 1.0;
  //@}

  /// speed of light, in cm per ns
  constexpr double speed_of_light = GSL_CONST_MKSA_SPEED_OF_LIGHT * 1e-7;

  /// length of generated G4Hits along laser track
  constexpr double maxHitLength = 1. * cm;

  /// TPC half length
  constexpr double halflength_tpc = 105.5 * cm;

  // inner and outer radii of field cages/TPC
  constexpr double begin_CM = 20. * cm;
  constexpr double end_CM = 78. * cm;

  // half the thickness of the CM;
  constexpr double halfwidth_CM = 0.5 * cm;

  //_____________________________________________________________
  std::optional<TVector3> central_membrane_intersection(const TVector3& start, const TVector3& direction)
  {
    const double end = start.z() > 0 ? halfwidth_CM : -halfwidth_CM;
    const double dist = end - start.z();

    // if line is vertical, it will never intercept the endcap
    if (direction.z() == 0)
    {
      return std::nullopt;
    }

    // check that distance and direction have the same sign
    if (dist * direction.z() < 0)
    {
      return std::nullopt;
    }

    const double direction_scale = dist / direction.z();
    return start + direction * direction_scale;
  }

  //_____________________________________________________________
  std::optional<TVector3> endcap_intersection(const TVector3& start, const TVector3& direction)
  {
    const double end = start.z() > 0 ? halflength_tpc : -halflength_tpc;
    const double dist = end - start.z();

    // if line is vertical, it will never intercept the endcap
    if (direction.z() == 0)
    {
      return std::nullopt;
    }

    // check that distance and direction have the same sign
    if (dist * direction.z() < 0)
    {
      return std::nullopt;
    }

    const double direction_scale = dist / direction.z();
    return start + direction * direction_scale;
  }

  //_____________________________________________________________
  std::optional<TVector3> cylinder_line_intersection(const TVector3& s, const TVector3& v, double radius)
  {
    const double R2 = square(radius);

    // Generalized Parameters for collision with cylinder of radius R:
    // from quadratic formula solutions of when a vector intersects a circle:
    const double a = square(v.x()) + square(v.y());
    const double b = 2 * (v.x() * s.x() + v.y() * s.y());
    const double c = square(s.x()) + square(s.y()) - R2;

    const double rootterm = square(b) - 4 * a * c;

    /*
     * if a==0 then we are parallel and will have no solutions.
     * if the rootterm is negative, we will have no real roots,
     * we are outside the cylinder and pointing skew to the cylinder such that we never cross.
     */
    if (rootterm < 0 || a == 0)
    {
      return std::nullopt;
    }

    // Find the (up to) two points where we collide with the cylinder:
    const double sqrtterm = std::sqrt(rootterm);
    const double t1 = (-b + sqrtterm) / (2 * a);
    const double t2 = (-b - sqrtterm) / (2 * a);

    /*
     * keep only intersections in front of the start point. For lasers that
     * originate outside the reference cylinder and point outward, both roots
     * are negative and we should report "no intersection" so that callers can
     * skip the tilt instead of tilting around a point behind the origin.
     */
    double min_t = 0.0;
    bool have_solution = false;
    if (t1 >= 0.0)
    {
      min_t = t1;
      have_solution = true;
    }
    if (t2 >= 0.0 && (!have_solution || t2 < min_t))
    {
      min_t = t2;
      have_solution = true;
    }
    if (!have_solution)
    {
      return std::nullopt;
    }
    return s + v * min_t;
  }

  //_____________________________________________________________
  std::optional<TVector3> field_cage_intersection(const TVector3& start, const TVector3& direction)
  {
    auto ofc_strike = cylinder_line_intersection(start, direction, end_CM);
    auto ifc_strike = cylinder_line_intersection(start, direction, begin_CM);

    // if either of the two intersection is invalid, return the other
    if (!ifc_strike)
    {
      return ofc_strike;
    }
    if (!ofc_strike)
    {
      return ifc_strike;
    }

    // both intersection are valid, calculate signed distance to start z
    const auto ifc_dist = (ifc_strike->Z() - start.Z()) / direction.Z();
    const auto ofc_dist = (ofc_strike->Z() - start.Z()) / direction.Z();

    if (ifc_dist < 0)
    {
      return (ofc_dist > 0) ? ofc_strike : std::nullopt;
    }
    if (ofc_dist < 0)
    {
      return ifc_strike;
    }
    else
    {
      return (ifc_dist < ofc_dist) ? ifc_strike : ofc_strike;
    }
  }

  /// TVector3 stream
  inline std::ostream& operator<<(std::ostream& out, const TVector3& vector)
  {
    out << "( " << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out;
  }

}  // namespace

//_____________________________________________________________
PHG4TpcDirectLaser::PHG4TpcDirectLaser(const std::string& name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

//_____________________________________________________________
int PHG4TpcDirectLaser::InitRun(PHCompositeNode* topNode)
{
  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_g4truthinfo)
  {
    std::cout << "Fun4AllDstPileupMerger::load_nodes - creating node G4TruthInfo" << std::endl;

    PHNodeIterator iter(topNode);
    auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_g4truthinfo = new PHG4TruthInfoContainer();
    dstNode->addNode(new PHIODataNode<PHObject>(m_g4truthinfo, "G4TruthInfo", "PHObject"));
  }

  // load and check G4Hit node
  hitnodename = "G4HIT_" + detector;
  auto* g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    std::cout << Name() << " Could not locate G4HIT node " << hitnodename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // find or create track map
  /* it is used to store laser parameters on a per event basis */
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name);
  if (!m_track_map)
  {
    // find DST node and check
    PHNodeIterator iter(topNode);
    auto* dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // find or create SVTX node
    iter = PHNodeIterator(dstNode);
    auto* node = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
    if (!node)
    {
      dstNode->addNode(node = new PHCompositeNode("SVTX"));
    }

    // add track node
    m_track_map = new SvtxTrackMap_v2;
    node->addNode(new PHIODataNode<PHObject>(m_track_map, m_track_map_name, "PHObject"));
  }

  // setup parameters
  UpdateParametersWithMacro();
  electrons_per_cm = get_double_param("electrons_per_cm");
  electrons_per_gev = get_double_param("electrons_per_gev");
  m_launch_offset_cm = std::max(0.0, get_double_param("launch_offset_cm"));

  m_tilt_layer = get_int_param("tilt_layer");
  m_tilt_steps = std::max(1, get_int_param("tilt_steps"));
  m_tilt_min_deg = get_double_param("tilt_min_deg");
  m_tilt_max_deg = get_double_param("tilt_max_deg");
  const double tilt_angle_deg = get_double_param("tilt_angle_deg");
  m_tilt_angle_rad = tilt_angle_deg * M_PI / 180.;
  m_active_tilt_angle_rad = m_tilt_angle_rad;
  m_current_tilt_step = 0;
  const bool use_tilt_range = (m_tilt_steps > 1) && (std::abs(m_tilt_max_deg - m_tilt_min_deg) > 1e-6);
  if (use_tilt_range)
  {
    m_active_tilt_angle_rad = m_tilt_min_deg * M_PI / 180.;
  }
  m_enable_tilt = (m_tilt_layer >= 0) && (use_tilt_range || std::abs(m_active_tilt_angle_rad) > 0.0);
  m_tilt_reference_radius = std::numeric_limits<double>::quiet_NaN();

  m_refine_layer = get_int_param("refine_layer");
  m_refine_halfwidth_cm = get_double_param("refine_halfwidth_cm");
  m_refine_step_cm = get_double_param("refine_step_cm");
  m_refine_active = false;

  m_tpc_geom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  if (m_enable_tilt)
  {
    if (!m_tpc_geom)
    {
      std::cout << Name() << ": CYLINDERCELLGEOM_SVTX node not found; disabling transverse tilt" << std::endl;
      m_enable_tilt = false;
    }
    else
    {
      if (auto* geom = m_tpc_geom->GetLayerCellGeom(m_tilt_layer))
      {
        m_tilt_reference_radius = geom->get_radius();
      }
      else
      {
        std::cout << Name() << ": geometry for tilt layer " << m_tilt_layer << " not found; disabling transverse tilt" << std::endl;
        m_enable_tilt = false;
      }
    }
  }

  if (m_refine_layer >= 0 && m_refine_halfwidth_cm > 0.0 && m_refine_step_cm > 0.0)
  {
    if (!m_tpc_geom)
    {
      std::cout << Name() << ": refine-layer request ignored (CYLINDERCELLGEOM_SVTX missing)" << std::endl;
    }
    else if (auto* geom = m_tpc_geom->GetLayerCellGeom(m_refine_layer))
    {
      const double radius = geom->get_radius();
      const double thickness = geom->get_thickness();
      const double rlow = std::max(0.0, radius - 0.5 * thickness);
      const double rhigh = radius + 0.5 * thickness;
      m_refine_rmin = std::max(0.0, rlow - m_refine_halfwidth_cm);
      m_refine_rmax = rhigh + m_refine_halfwidth_cm;
      m_refine_active = true;
    }
    else
    {
      std::cout << Name() << ": refine-layer " << m_refine_layer << " not found in geometry; disabling refinement" << std::endl;
    }
  }

  // setup lasers
  SetupLasers();

  // allocate and seed private GSL RNG for random-phi sampling
  m_rng.reset(gsl_rng_alloc(gsl_rng_mt19937));
  if (m_rng)
  {
    gsl_rng_set(m_rng.get(), PHRandomSeed());
  }

  // print configuration
  if (m_steppingpattern == true)
  {
    std::cout << "PHG4TpcDirectLaser::InitRun - m_steppingpattern: " << m_steppingpattern << std::endl;
    std::cout << "PHG4TpcDirectLaser::InitRun - nTotalSteps: " << nTotalSteps << std::endl;
  }
  else
  {
    std::cout << "PHG4TpcDirectLaser::InitRun - m_autoAdvanceDirectLaser: " << m_autoAdvanceDirectLaser << std::endl;
    std::cout << "PHG4TpcDirectLaser::InitRun - phi steps: " << nPhiSteps << " min: " << minPhi << " max: " << maxPhi << std::endl;
    std::cout << "PHG4TpcDirectLaser::InitRun - theta steps: " << nThetaSteps << " min: " << minTheta << " max: " << maxTheta << std::endl;
    std::cout << "PHG4TpcDirectLaser::InitRun - nTotalSteps: " << nTotalSteps << std::endl;
  }
  std::cout << "PHG4TpcDirectLaser::InitRun - electrons_per_cm: " << electrons_per_cm << std::endl;
  std::cout << "PHG4TpcDirectLaser::InitRun - electrons_per_gev " << electrons_per_gev << std::endl;
  if (m_enable_tilt)
  {
    std::cout << "PHG4TpcDirectLaser::InitRun - transverse tilt layer: " << m_tilt_layer
              << " angle(deg): ";
    if (use_tilt_range)
    {
      std::cout << m_tilt_min_deg << " -> " << m_tilt_max_deg
                << " in " << m_tilt_steps << " steps" << std::endl;
    }
    else
    {
      std::cout << tilt_angle_deg << std::endl;
    }
  }

  // If using pattern stepping from file, load angles from CALIBRATIONROOT only then
  if (m_steppingpattern)
  {
    const char* calibroot_env = std::getenv("CALIBRATIONROOT");
    if (!calibroot_env)
    {
      std::cout << Name() << ": CALIBRATIONROOT is not set; disabling file-based pattern stepping" << std::endl;
      m_steppingpattern = false;
    }
    else
    {
      const std::string LASER_ANGLES_ROOTFILE = std::string(calibroot_env) + "/TPC/DirectLaser/theta_phi_laser.root";
      TFile* infile1 = TFile::Open(LASER_ANGLES_ROOTFILE.c_str());
      if (!infile1 || infile1->IsZombie())
      {
        std::cout << Name() << ": cannot open " << LASER_ANGLES_ROOTFILE << "; disabling file-based pattern stepping" << std::endl;
        m_steppingpattern = false;
      }
      else
      {
        pattern = dynamic_cast<TNtuple*>(infile1->Get("angles"));
        if (!pattern)
        {
          std::cout << Name() << ": TNtuple 'angles' not found in " << LASER_ANGLES_ROOTFILE << "; disabling file-based pattern stepping" << std::endl;
          m_steppingpattern = false;
        }
        else
        {
          pattern->SetBranchAddress("#theta", &theta_p);
          pattern->SetBranchAddress("#phi", &phi_p);
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
int PHG4TpcDirectLaser::process_event(PHCompositeNode* topNode)
{
  // g4 input event
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  assert(m_g4truthinfo);

  // load g4hit container
  m_g4hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  assert(m_g4hitcontainer);

  // load track map
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name);
  assert(m_track_map);

  // if random phi mode is enabled, override stepping and pick a uniform phi in [minPhi,maxPhi]
  if (m_use_random_phi)
  {
    // use configured theta; if a range was set, pick the lower edge (common usage is a fixed theta)
    const double theta = (nThetaSteps > 0) ? minTheta : 0.0;
    // draw phi parameter uniformly in [minPhi, maxPhi]
    double phi = minPhi;
    if (maxPhi > minPhi)
    {
      // uniform draw from [minPhi,maxPhi) using module-private GSL RNG
      const double width = (maxPhi - minPhi);
      phi = minPhi + gsl_ran_flat(m_rng.get(), 0.0, width);
    }
    AimToThetaPhi(theta, phi);
  }
  else if (m_autoAdvanceDirectLaser || m_steppingpattern)
  {
    AimToNextPatternStep();
  }
  //_________________________________________________
  else
  {
    // use arbitrary direction
    AimToThetaPhi(arbitrary_theta, arbitrary_phi);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::SetDefaultParameters()
{
  // same gas parameters as in PHG4TpcElectronDrift::SetDefaultParameters

  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html
/*   static constexpr double Ne_dEdx = 1.56;    // keV/cm
  static constexpr double CF4_dEdx = 7.00;   // keV/cm
  static constexpr double Ne_NTotal = 43;    // Number/cm
  static constexpr double CF4_NTotal = 100;  // Number/cm
  static constexpr double Tpc_NTot = 0.5 * Ne_NTotal + 0.5 * CF4_NTotal;
  static constexpr double Tpc_dEdx = 0.5 * Ne_dEdx + 0.5 * CF4_dEdx;
  static constexpr double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx; */


  //static constexpr double Ne_dEdx = 1.56;    // keV/cm
  static constexpr double CF4_dEdx = 7.00;   // keV/cm
  static constexpr double Ar_dEdx = 2.44;// keV/cm
  static constexpr double isobutane_dEdx = 5.93;// keV/cm

  //static constexpr double Ne_NTotal = 43;    // Number/cm
  static constexpr double CF4_NTotal = 100;  // Number/cm
  static constexpr double Ar_NTotal = 94; // Number/cm
  static constexpr double isobutane_NTotal = 195; // Number/cm
  
  static constexpr double Tpc_NTot = 0.75 * Ar_NTotal + 0.20 * CF4_NTotal + 0.05 * isobutane_NTotal;
  static constexpr double Tpc_dEdx = 0.75 * Ar_dEdx + 0.20 * CF4_dEdx + 0.05 * isobutane_dEdx;
  static constexpr double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  // number of electrons per deposited GeV in TPC gas
  set_default_double_param("electrons_per_gev", Tpc_ElectronsPerKeV * 1e6);

  // number of electrons deposited by laser per cm
  set_default_double_param("electrons_per_cm", 100.25);

  // optional transverse tilt (disabled by default)
  set_default_int_param("tilt_layer", -1);
  set_default_double_param("tilt_angle_deg", 0.0);
  set_default_double_param("tilt_min_deg", 0.0);
  set_default_double_param("tilt_max_deg", 0.0);
  set_default_int_param("tilt_steps", 1);
  set_default_int_param("refine_layer", -1);
  set_default_double_param("refine_halfwidth_cm", 0.0);
  set_default_double_param("refine_step_cm", 0.0);
  set_default_double_param("launch_offset_cm", 50.0);
}

//_____________________________________________________________
void PHG4TpcDirectLaser::SetPhiStepping(int n, double min, double max)
{
  if (n < 0 || max < min)
  {
    std::cout << PHWHERE << " - invalid" << std::endl;
    return;
  }
  nPhiSteps = n;
  minPhi = min;
  maxPhi = max;
  nTotalSteps = nThetaSteps * nPhiSteps;
  return;
}
//_____________________________________________________________
void PHG4TpcDirectLaser::SetThetaStepping(int n, double min, double max)
{
  if (n < 0 || max < min)
  {
    std::cout << PHWHERE << " - invalid" << std::endl;
    return;
  }
  nThetaSteps = n;
  minTheta = min;
  maxTheta = max;
  nTotalSteps = nThetaSteps * nPhiSteps;

  return;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::SetFileStepping(int n)
{
  if (n < 0 || n > 13802)  // 13802 = hard coded number of tuple entries
  {
    std::cout << PHWHERE << " - invalid" << std::endl;
    return;
  }
  nTotalSteps = n;

  return;
}

//_____________________________________________________________

void PHG4TpcDirectLaser::SetupLasers()
{
  // clear previous lasers
  m_lasers.clear();

  // position of first laser at positive z
  const double launch_offset = m_launch_offset_cm * cm;
  const double launch_radius = 60.0 * cm;
  const TVector3 position_base(launch_radius, 0., halflength_tpc - launch_offset);

  // add lasers
  for (int i = 0; i < 8; ++i)
  {
    Laser laser;

    // set laser direction
    /*
     * first four lasers are on positive z readout plane, and shoot towards negative z
     * next four lasers are on negative z readout plane and shoot towards positive z
     */
    laser.m_position = position_base;
    if (i < 4)
    {
      laser.m_position.SetZ(position_base.z());
      laser.m_direction = -1;
      laser.m_phi = M_PI / 2 * i - (15 * M_PI / 180);  // additional offset of 15 deg.
    }
    else
    {
      laser.m_position.SetZ(-(halflength_tpc - launch_offset));
      laser.m_direction = 1;
      laser.m_phi = M_PI / 2 * i + (15 * M_PI / 180);  // additional offset of 15 deg.
    }

    // rotate around z
    laser.m_position.RotateZ(laser.m_phi);

    // append: either all lasers or a selected single index
    if (m_selected_laser_index < 0 || i == m_selected_laser_index)
    {
      m_lasers.push_back(laser);
    }
    //  if(i==0) m_lasers.push_back(laser);//Only laser 1
    //  if(i==3) m_lasers.push_back(laser);// Laser 4
    // if(i<4) m_lasers.push_back(laser);//Lasers 1, 2, 3, 4
  }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToNextPatternStep()
{
  if (nTotalSteps >= 1)
  {
    if (m_steppingpattern)
    {
      AimToPatternStep_File(currentPatternStep);
      ++currentPatternStep;
    }
    else
    {
      AimToPatternStep(currentPatternStep);
      ++currentPatternStep;
    }
  }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToThetaPhi(double theta, double phi)
{
  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AimToThetaPhi - theta: " << theta << " phi: " << phi << std::endl;
  }

  if (m_enable_tilt)
  {
    UpdateActiveTiltAngle();
    if (Verbosity())
    {
      std::cout << "PHG4TpcDirectLaser::AimToThetaPhi - active tilt angle (deg): "
                << m_active_tilt_angle_rad * 180. / M_PI << std::endl;
    }
  }

  // all lasers
  for (const auto& laser : m_lasers)
  {
    AppendLaserTrack(theta, phi, laser);
  }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::UpdateActiveTiltAngle()
{
  if (!m_enable_tilt)
  {
    m_active_tilt_angle_rad = m_tilt_angle_rad;
    return;
  }

  if (m_tilt_steps <= 1 || std::abs(m_tilt_max_deg - m_tilt_min_deg) < 1e-6)
  {
    m_active_tilt_angle_rad = m_tilt_angle_rad;
    return;
  }

  const double fraction = (m_tilt_steps > 1) ? static_cast<double>(m_current_tilt_step) / (m_tilt_steps - 1) : 0.0;
  const double angle_deg = m_tilt_min_deg + fraction * (m_tilt_max_deg - m_tilt_min_deg);
  m_active_tilt_angle_rad = angle_deg * M_PI / 180.;
  m_current_tilt_step = (m_current_tilt_step + 1) % m_tilt_steps;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToPatternStep(int n)
{
  // trim against overflows
  n = n % nTotalSteps;

  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AimToPatternStep - step: " << n << "/" << nTotalSteps << std::endl;
  }

  // store as current pattern
  currentPatternStep = n;

  // calculate theta
  const int thetaStep = n / nPhiSteps;
  const double theta = minTheta + thetaStep * (maxTheta - minTheta) / nThetaSteps;

  // calculate phi
  const int phiStep = n % nPhiSteps;
  const double phi = minPhi + phiStep * (maxPhi - minPhi) / nPhiSteps;

  // generate laser tracks
  AimToThetaPhi(theta, phi);

  return;
}

//_____________________________________________________________

void PHG4TpcDirectLaser::AimToPatternStep_File(int n)
{
  // trim against overflows
  n = n % nTotalSteps;

  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AimToPatternStep_File - step: " << n << "/" << nTotalSteps << std::endl;
  }

  // store as current pattern
  currentPatternStep = n;

  if (!pattern)
  {
    if (Verbosity())
    {
      std::cout << Name() << ": pattern TNtuple not available; skipping file-based pattern step" << std::endl;
    }
    return;
  }

  pattern->GetEntry(n);

  // calculate theta
  std::cout << "From file, current entry = " << n << " Theta: " << theta_p << " Phi: " << phi_p << std::endl;

  const double theta = theta_p * M_PI / 180.;

  // calculate phi
  const double phi = phi_p * M_PI / 180.;

  // generate laser tracks
  AimToThetaPhi(theta, phi);

  return;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::ApplyTilt(TVector3& pos, TVector3& dir) const
{
  if (!m_enable_tilt)
  {
    return;
  }

  const double original_dir_mag = dir.Mag();
  if (original_dir_mag <= 0)
  {
    return;
  }

  const TVector3 original_pos = pos;
  const TVector3 original_dir = dir;

  const double dir_xy_mag = std::hypot(original_dir.x(), original_dir.y());
  if (dir_xy_mag <= 0)
  {
    return;
  }

  std::optional<TVector3> pivot;
  if (std::isfinite(m_tilt_reference_radius))
  {
    pivot = cylinder_line_intersection(original_pos, original_dir, m_tilt_reference_radius);
  }

  if (!pivot)
  {
    return;
  }

  TVector3 radial = *pivot;
  radial.SetZ(0.0);

  const double radial_mag = radial.Perp();
  if (radial_mag <= 0)
  {
    return;
  }
  radial *= 1.0 / radial_mag;

  TVector3 phi_hat = TVector3(0.0, 0.0, 1.0).Cross(radial);
  const double phi_mag = phi_hat.Mag();
  if (phi_mag <= 0)
  {
    return;
  }
  phi_hat *= 1.0 / phi_mag;

  const double cos_tilt = std::cos(m_active_tilt_angle_rad);
  const double sin_tilt = std::sin(m_active_tilt_angle_rad);

  const double dir_z = original_dir.z();
  TVector3 desired_dir = dir_xy_mag * (cos_tilt * radial + sin_tilt * phi_hat) + TVector3(0.0, 0.0, dir_z);
  const double desired_mag = desired_dir.Mag();
  if (desired_mag <= 0)
  {
    return;
  }
  TVector3 dir_unit = desired_dir * (1.0 / desired_mag);
  TVector3 new_dir = dir_unit * original_dir_mag;

  if (Verbosity())
  {
    const double cos_alpha = radial.Dot(dir_unit);
    const double clamped = std::max(-1.0, std::min(1.0, cos_alpha));
    const double alpha_deg = std::acos(clamped) * 180. / M_PI;
    std::cout << Name() << " tilt layer " << m_tilt_layer
              << " cos(alpha)=" << cos_alpha
              << " alpha(deg)=" << alpha_deg << std::endl;
  }

  const double original_distance = (*pivot - original_pos).Mag();
  TVector3 new_pos = *pivot - dir_unit * original_distance;

  pos = new_pos;
  dir = new_dir;
}

//_____________________________________________________________

void PHG4TpcDirectLaser::AppendLaserTrack(double theta, double phi, const PHG4TpcDirectLaser::Laser& laser)
{
  if (!m_g4hitcontainer)
  {
    std::cout << PHWHERE << "invalid g4hit container. aborting" << std::endl;
    return;
  }

  // store laser position; optionally rotate origin with phi so that
  // the origin azimuth matches the direction azimuth (radial emission)
  TVector3 pos = laser.m_position;

  // define track direction
  const auto& direction = laser.m_direction;
  TVector3 dir(0, 0, direction);

  // adjust direction
  dir.RotateY(theta * direction);

  if (laser.m_direction == -1)
  {
    dir.RotateZ(phi);  // if +z facing -z
    if (m_lock_origin_to_phi) pos.RotateZ(phi);
  }
  else
  {
    dir.RotateZ(-phi);  // if -z facting +z
    if (m_lock_origin_to_phi) pos.RotateZ(-phi);
  }

  // also rotate by laser azimuth
  dir.RotateZ(laser.m_phi);

  // apply optional transverse tilt relative to reference layer
  ApplyTilt(pos, dir);

  // print
  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AppendLaserTrack - position: " << pos << " direction: " << dir << std::endl;
  }

  // dummy momentum
  static constexpr double total_momentum = 1;

  // mc track id
  int trackid = -1;

  // create truth vertex and particle
  if (m_g4truthinfo)
  {
    // add vertex
    const auto vtxid = m_g4truthinfo->maxvtxindex() + 1;
    auto* const vertex = new PHG4VtxPoint_t(pos.x(), pos.y(), pos.z(), 0, vtxid);
    m_g4truthinfo->AddVertex(vtxid, vertex);

    // increment track id
    trackid = m_g4truthinfo->maxtrkindex() + 1;

    // create new g4particle
    auto* particle = new PHG4Particle_t();
    particle->set_track_id(trackid);
    particle->set_vtx_id(vtxid);
    particle->set_parent_id(0);
    particle->set_primary_id(trackid);
    particle->set_px(total_momentum * dir.x());
    particle->set_py(total_momentum * dir.y());
    particle->set_pz(total_momentum * dir.z());

    m_g4truthinfo->AddParticle(trackid, particle);
  }

  // store in SvtxTrack map
  if (m_track_map)
  {
    // allocate on heap; SvtxTrackMap takes ownership
    auto* track = new SvtxTrack_v2();
    track->set_x(pos.x());
    track->set_y(pos.y());
    track->set_z(pos.z());

    // total momentum is irrelevant. What matters is the direction
    track->set_px(total_momentum * dir.x());
    track->set_py(total_momentum * dir.y());
    track->set_pz(total_momentum * dir.z());

    // insert in map
    m_track_map->insert(track);

    if (Verbosity())
    {
      std::cout << "PHG4TpcDirectLaser::AppendLaserTrack - position: " << pos << " direction: " << dir << std::endl;
    }
  }

  // find collision point
  /*
   * intersection to either central membrane or endcaps
   * if the position along beam and laser direction have the same sign, it will intercept the endcap
   * otherwise will intercept the central membrane
   */
  const auto plane_strike = (pos.z() * dir.z() > 0) ? endcap_intersection(pos, dir) : central_membrane_intersection(pos, dir);

  // field cage intersection
  const auto fc_strike = field_cage_intersection(pos, dir);

  // if none of the strikes is valid, there is no valid information found.
  if (!(plane_strike || fc_strike))
  {
    return;
  }

  // decide relevant end of laser
  /* chose field cage intersection if valid, and if either plane intersection is invalid or happens on a larger z along the laser direction) */
  const TVector3& strike = (fc_strike && (!plane_strike || fc_strike->z() / dir.z() < plane_strike->z() / dir.z())) ? *fc_strike : *plane_strike;

  // find length
  const double fullLength = (strike - pos).Mag();
  int nHitSteps = fullLength / maxHitLength + 1;

  TVector3 start = pos;
  TVector3 end = start;
  TVector3 step = dir * (maxHitLength / (dir.Mag()));

  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AppendLaserTrack -"
              << " fullLength: " << fullLength
              << " nHitSteps: " << nHitSteps
              << std::endl;
  }

  const double refine_step = m_refine_step_cm * cm;

  for (int i = 0; i < nHitSteps; i++)
  {
    start = end;  // new starting point is the previous ending point.
    if (i + 1 == nHitSteps)
    {
      // last step is the remainder size
      end = strike;
    }
    else
    {
      // all other steps are uniform length
      end = start + step;
    }

    const bool do_refine = m_refine_active && m_refine_step_cm > 0.0 && segmentNeedsRefinement(start, end);
    if (do_refine)
    {
      const TVector3 segment = end - start;
      const double length = segment.Mag();
      const int n_sub = std::max(1, static_cast<int>(std::ceil(length / refine_step)));
      for (int isub = 0; isub < n_sub; ++isub)
      {
        const double frac0 = static_cast<double>(isub) / n_sub;
        const double frac1 = static_cast<double>(isub + 1) / n_sub;
        const TVector3 sub_start = start + segment * frac0;
        const TVector3 sub_end = start + segment * frac1;
        emitLaserHit(trackid, sub_start, sub_end, dir, pos);
      }
    }
    else
    {
      emitLaserHit(trackid, start, end, dir, pos);
    }
  }

  return;
}

bool PHG4TpcDirectLaser::segmentNeedsRefinement(const TVector3& start, const TVector3& end) const
{
  if (!m_refine_active)
  {
    return false;
  }

  const double r_start = start.Perp();
  const double r_end = end.Perp();

  if (std::max(r_start, r_end) >= m_refine_rmin && std::min(r_start, r_end) <= m_refine_rmax)
  {
    return true;
  }

  // Radial laser shots are effectively straight towards or away from the origin.
  // The endpoint check above is therefore sufficient; skip the more expensive
  // closest-approach test unless we re-enable it for non-radial studies.
  // (Keep the original logic below for quick reactivation.)
  /*
  const TVector3 delta = end - start;
  const double a = delta.X() * delta.X() + delta.Y() * delta.Y();
  if (a <= 0.0)
  {
    return false;
  }

  const double b = start.X() * delta.X() + start.Y() * delta.Y();
  const double s_ext = -b / a;
  if (s_ext > 0.0 && s_ext < 1.0)
  {
    const TVector3 point = start + delta * s_ext;
    const double r_mid = point.Perp();
    if (r_mid >= m_refine_rmin && r_mid <= m_refine_rmax)
    {
      return true;
    }
  }
  */
  return false;
}

void PHG4TpcDirectLaser::emitLaserHit(int trackid, const TVector3& start, const TVector3& end, const TVector3& dir, const TVector3& origin)
{
  const double stepLength = (end - start).Mag();

  // from phg4tpcsteppingaction.cc
  auto* hit = new PHG4Hit_t;
  hit->set_trkid(trackid);
  hit->set_layer(99);

  // here we set the entrance values in cm
  hit->set_x(0, start.X() / cm);
  hit->set_y(0, start.Y() / cm);
  hit->set_z(0, start.Z() / cm);
  hit->set_t(0, (start - origin).Mag() / speed_of_light);

  hit->set_x(1, end.X() / cm);
  hit->set_y(1, end.Y() / cm);
  hit->set_z(1, end.Z() / cm);
  hit->set_t(1, (end - origin).Mag() / speed_of_light);

  // momentum
  hit->set_px(0, dir.X());  // GeV
  hit->set_py(0, dir.Y());
  hit->set_pz(0, dir.Z());

  hit->set_px(1, dir.X());
  hit->set_py(1, dir.Y());
  hit->set_pz(1, dir.Z());

  const double totalE = electrons_per_cm * stepLength / electrons_per_gev;
  if (Verbosity() > 0)
  {
    const double r_probe = 0.5 * (start.Perp() + end.Perp());
    const int layer_guess = findLayerForRadius(r_probe);
    std::cout << Name() << ": laser step length " << stepLength
              << " cm, electrons/cm " << electrons_per_cm
              << ", electrons/GeV " << electrons_per_gev
              << ", assigned eion " << totalE << " GeV";
    if (layer_guess >= 0)
    {
      std::cout << " (layer " << layer_guess << ", r ~ " << r_probe << " cm)";
    }
    std::cout << std::endl;
  }

  hit->set_eion(totalE);
  hit->set_edep(totalE);
  m_g4hitcontainer->AddHit(detId, hit);
}

int PHG4TpcDirectLaser::findLayerForRadius(double radius) const
{
  if (!m_tpc_geom)
  {
    return -1;
  }

  const int nLayers = m_tpc_geom->get_NLayers();
  for (int ilayer = 0; ilayer < nLayers; ++ilayer)
  {
    if (auto* geom = m_tpc_geom->GetLayerCellGeom(ilayer))
    {
      const double rcen = geom->get_radius();
      const double half = 0.5 * geom->get_thickness();
      if (radius >= rcen - half && radius < rcen + half)
      {
        return geom->get_layer();
      }
    }
  }
  return -1;
}
