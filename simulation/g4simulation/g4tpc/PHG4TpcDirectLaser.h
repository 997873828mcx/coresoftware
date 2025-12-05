#ifndef G4TPC_PHG4TPCDIRECTLASER_H
#define G4TPC_PHG4TPCDIRECTLASER_H

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameterInterface.h>

#include <TNtuple.h>
#include <TVector3.h>

#include <cmath>
#include <limits>
#include <string>  // for string, allocator
#include <vector>  // for vector
#include <memory>

#include <gsl/gsl_rng.h>

class PHG4HitContainer;
class SvtxTrackMap;
class PHG4TruthInfoContainer;
class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;

class PHG4TpcDirectLaser : public SubsysReco, public PHParameterInterface
{
 public:
  /// constructor
  PHG4TpcDirectLaser(const std::string &name = "PHG4TpcDirectLaser");

  /// destructor
  ~PHG4TpcDirectLaser() override = default;

  /// run initialization
  int InitRun(PHCompositeNode *) override;

  /// per event processing
  int process_event(PHCompositeNode *) override;

  /// default parameters
  void SetDefaultParameters() override;

  /// detector name
  void Detector(const std::string &d)
  {
    detector = d;
  }

  /// define steps along phi
  void SetPhiStepping(int n, double min, double max);

  /// define steps along theta
  void SetThetaStepping(int n, double min, double max);

  /// define steps for file
  void SetFileStepping(int n);

  /// get total number of steps
  int GetNpatternSteps() const
  {
    return nPhiSteps * nThetaSteps;
  };

  /// set current patter step
  void SetCurrentPatternStep(int value)
  {
    currentPatternStep = value;
  }

  /// advance automatically through patterns
  void SetDirectLaserAuto(bool value)
  {
    m_autoAdvanceDirectLaser = value;
  };

  /// advance automatically through pattern from file
  void SetDirectLaserPatternfromFile(bool value)
  {
    m_steppingpattern = value;
  };

  void SetArbitraryThetaPhi(double theta, double phi)
  {
    arbitrary_theta = theta;
    arbitrary_phi = phi;
  }

  /// select a single laser index (0..7). Set to -1 to use all.
  void SetSingleLaserIndex(int idx)
  {
    m_selected_laser_index = idx;
  }

  /// when true, rotate the laser origin around z by the same phi step
  /// used for the direction, so that the global azimuth of the origin
  /// matches the direction (radial emission).
  void SetLockOriginToPhi(bool value)
  {
    m_lock_origin_to_phi = value;
  }

  /// enable random phi per event within [minPhi,maxPhi]
  /// when enabled, process_event ignores stepping and draws a uniform phi each event
  void EnableRandomPhi(bool value)
  {
    m_use_random_phi = value;
  }

  /// convenience: set random phi range and enable
  void SetRandomPhiRange(double min, double max)
  {
    minPhi = min;
    maxPhi = max;
    m_use_random_phi = true;
  }

  /// set number of laser tracks to fire per event (default 1)
  void set_tracks_per_event(int n)
  {
    m_tracks_per_event = n;
    set_int_param("tracks_per_event", n);
  }

 private:
  /// define lasers
  /* by default there are 4 lasers on each side of the TPC */
  void SetupLasers();

  /// aim lasers to a given theta and phi angle
  void AimToThetaPhi(double theta, double phi);

  /// aim lasers to a give step
  void AimToPatternStep(int n);

  /// aim lasers to a give step from file
  void AimToPatternStep_File(int n);

  float theta_p{0};
  float phi_p{0};
  TNtuple *pattern{nullptr};

  /// aim to next step
  void AimToNextPatternStep();

  /// stores laser position and direction along z
  class Laser
  {
   public:
    /// laser position
    TVector3 m_position;

    /// laser phi position
    double m_phi{0};

    /// laser direction along z
    int m_direction{1};
  };

  /// append track in given angular direction and for a given laser
  void AppendLaserTrack(double theta, double phi, const Laser &);

  /// apply transverse tilt around the radial direction for the configured layer
  void ApplyTilt(TVector3 &pos, TVector3 &dir) const;
  void UpdateActiveTiltAngle();

  /// detector name
  std::string detector{"TPC"};

  /// g4hitnode name
  std::string hitnodename;

  /// lasers
  std::vector<Laser> m_lasers;

  /// number of electrons deposited per cm laser track
  double electrons_per_cm{300.0};

  // number of electrons per deposited GeV in TPC gas
  /**
   * it is used to convert a given number of electrons into an energy
   * as expected by G4Hit. The energy is then converted back to a number of electrons
   * inside PHG4TpcElectronDrift
   */
  double electrons_per_gev{std::numeric_limits<double>::signaling_NaN()};
  double m_launch_offset_cm{50.0};

  double arbitrary_theta{-30.0};  // degrees
  double arbitrary_phi{-30.0};    // degrees

  ///@name default phi and theta steps
  //@{
  int nPhiSteps{1};
  int nThetaSteps{1};
  int nTotalSteps{1};
  double minPhi{0};
  double maxPhi{0};
  double minTheta{0};
  double maxTheta{0};
  //@}

  // current patter step
  int currentPatternStep{0};

  /// set to true to change direct laser tracks from one event to the other
  bool m_autoAdvanceDirectLaser{false};

  /// set to true to get stepping patern from file
  bool m_steppingpattern{false};

  /// g4hit container
  PHG4HitContainer *m_g4hitcontainer{nullptr};

  //! truth information
  PHG4TruthInfoContainer *m_g4truthinfo{nullptr};

  /// track map, used to store track parameters
  std::string m_track_map_name{"SvtxTrackMap"};
  SvtxTrackMap *m_track_map{nullptr};

  /// single-laser selection; -1 means all lasers
  int m_selected_laser_index{-1};

  /// if true, rotate origin with phi so direction is radial wrt layers
  bool m_lock_origin_to_phi{false};

  /// if true, choose a random phi each event between [minPhi,maxPhi]
  bool m_use_random_phi{false};

  /// number of tracks to fire per event
  int m_tracks_per_event{1};

  // GSL RNG for random-phi sampling (private to this module)
  struct GslDeleter { void operator()(gsl_rng* p) const { if(p) gsl_rng_free(p); } };
  std::unique_ptr<gsl_rng, GslDeleter> m_rng;

  /// geometry container (needed for per-layer tilt reference)
  PHG4TpcCylinderGeomContainer *m_tpc_geom{nullptr};
  int m_refine_layer{-1};
  double m_refine_halfwidth_cm{0.0};
  double m_refine_step_cm{0.0};
  double m_refine_rmin{std::numeric_limits<double>::quiet_NaN()};
  double m_refine_rmax{std::numeric_limits<double>::quiet_NaN()};
  bool m_refine_active{false};
  bool segmentNeedsRefinement(const TVector3& start, const TVector3& end) const;
  void emitLaserHit(int trackid, const TVector3& start, const TVector3& end, const TVector3& dir, const TVector3& origin);
  int findLayerForRadius(double radius) const;

  ///@name transverse tilt configuration
  //@{
  bool m_enable_tilt{false};
  int m_tilt_layer{-1};
  double m_tilt_angle_rad{0};
  double m_active_tilt_angle_rad{0};
  double m_tilt_min_deg{0};
  double m_tilt_max_deg{0};
  int m_tilt_steps{1};
  int m_current_tilt_step{0};
  double m_tilt_reference_radius{std::numeric_limits<double>::quiet_NaN()};
  //@}
};

#endif
