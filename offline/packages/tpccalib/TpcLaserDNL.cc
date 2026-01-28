#include "TpcLaserDNL.h"
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/ActsGeometry.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <g4tracking/TrkrTruthTrackContainer.h>
#include <g4tracking/TrkrTruthTrack.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <Eigen/Dense>

// For EventHeader event sequence tagging
#include <ffaobjects/EventHeader.h>

TpcLaserDNL::TpcLaserDNL(const std::string& name)
: SubsysReco(name)
{}

int TpcLaserDNL::Init(PHCompositeNode*)
{
m_tf.reset(TFile::Open(m_outfile.c_str(),"RECREATE"));
m_tt = new TTree("dnl","laser dnl");
m_tt->Branch("event",&m_event,"event/I");
m_tt->Branch("trkid",&m_trkid,"trkid/I");
m_tt->Branch("layer",&m_layer,"layer/i");
m_tt->Branch("side",&m_side,"side/I");
m_tt->Branch("r",&m_r,"r/D");
m_tt->Branch("phi_true",&m_phi_true,"phi_true/D");
m_tt->Branch("phi_reco",&m_phi_reco,"phi_reco/D");
m_tt->Branch("dphi",&m_dphi,"dphi/D");
m_tt->Branch("dRphi",&m_dRphi,"dRphi/D");
m_tt->Branch("nused",&m_nused,"nused/I");
m_tt->Branch("nhit_scanned",&m_nhit_scanned,"nhit_scanned/I");
m_tt->Branch("adcsum",&m_adcsum,"adcsum/D");
m_tt->Branch("xtrue",&m_xtrue,"xtrue/D");
m_tt->Branch("ytrue",&m_ytrue,"ytrue/D");
m_tt->Branch("ztrue",&m_ztrue,"ztrue/D");
m_tt->Branch("xreco",&m_xreco,"xreco/D");
m_tt->Branch("yreco",&m_yreco,"yreco/D");
m_tt->Branch("zreco",&m_zreco,"zreco/D");
m_tt->Branch("npad_used",&m_npad_used,"npad_used/I");
m_tt->Branch("ntbin_used",&m_ntbin_used,"ntbin_used/I");
m_tt->Branch("nbins_used",&m_nbins_used,"nbins_used/I");
m_tt->Branch("phi_pad_max",&m_phi_pad_max,"phi_pad_max/D");
m_tt->Branch("phase",&m_phase,"phase/D");
m_tt->Branch("phase_reco",&m_phase_reco,"phase_reco/D");
m_tt->Branch("pad_phi_center",&m_pad_phi_centers);
// fitted straight-line prediction at layer 44 (NaN elsewhere)
m_tt->Branch("xfit_layer44",&m_xfit_layer44,"xfit_layer44/D");
m_tt->Branch("yfit_layer44",&m_yfit_layer44,"yfit_layer44/D");
m_tt->Branch("zfit_layer44",&m_zfit_layer44,"zfit_layer44/D");
m_tt->Branch("phi_fit_layer44",&m_phi_fit_layer44,"phi_fit_layer44/D");
m_tt->Branch("dphi_fit_layer44",&m_dphi_fit_layer44,"dphi_fit_layer44/D");
m_tt->Branch("dRphi_fit_layer44",&m_dRphi_fit_layer44,"dRphi_fit_layer44/D");
// additional charge bookkeeping
m_tt->Branch("hit_charge",&m_hit_charge);
m_tt->Branch("total_charge_layer",&m_total_charge_layer,"total_charge_layer/D");
m_tt->Branch("max_charge_layer",&m_max_charge_layer,"max_charge_layer/D");
// debug vectors (only filled with hits that pass selection)
m_tt->Branch("hitkey",&m_hitkeys);
m_tt->Branch("hitsetkey",&m_hitsetkeys);
m_tt->Branch("iphi",&m_iphi);
m_tt->Branch("tbin",&m_tbin);

  if(m_write_display_ntuple)
  {
    m_tt_display_intersections = new TTree("truth_intersections","truth track-cylinder intersections");
    m_tt_display_intersections->Branch("event",&m_display_intersection.event,"event/I");
    m_tt_display_intersections->Branch("trackid",&m_display_intersection.trackid,"trackid/I");
    m_tt_display_intersections->Branch("layer",&m_display_intersection.layer,"layer/I");
    m_tt_display_intersections->Branch("side",&m_display_intersection.side,"side/I");
    m_tt_display_intersections->Branch("gx",&m_display_intersection.gx,"gx/D");
    m_tt_display_intersections->Branch("gy",&m_display_intersection.gy,"gy/D");
    m_tt_display_intersections->Branch("gz",&m_display_intersection.gz,"gz/D");
    m_tt_display_intersections->Branch("r",&m_display_intersection.r,"r/D");
    m_tt_display_intersections->Branch("phi",&m_display_intersection.phi,"phi/D");
    m_tt_display_intersections->Branch("path",&m_display_intersection.path,"path/D");
    m_tt_display_intersections->Branch("used_in_seed",&m_display_intersection.used_in_seed,"used_in_seed/I");

    m_tt_display_g4hits = new TTree("truth_g4hits","sampled G4 hits for truth tracks");
    m_tt_display_g4hits->Branch("event",&m_display_g4hit.event,"event/I");
    m_tt_display_g4hits->Branch("trackid",&m_display_g4hit.trackid,"trackid/I");
    m_tt_display_g4hits->Branch("layer",&m_display_g4hit.layer,"layer/I");
    m_tt_display_g4hits->Branch("side",&m_display_g4hit.side,"side/I");
    m_tt_display_g4hits->Branch("gx",&m_display_g4hit.gx,"gx/D");
    m_tt_display_g4hits->Branch("gy",&m_display_g4hit.gy,"gy/D");
    m_tt_display_g4hits->Branch("gz",&m_display_g4hit.gz,"gz/D");
    m_tt_display_g4hits->Branch("r",&m_display_g4hit.r,"r/D");
    m_tt_display_g4hits->Branch("phi",&m_display_g4hit.phi,"phi/D");
    m_tt_display_g4hits->Branch("sample",&m_display_g4hit.sample,"sample/I");
    m_tt_display_g4hits->Branch("sample_frac",&m_display_g4hit.sample_frac,"sample_frac/D");
    m_tt_display_g4hits->Branch("step_path",&m_display_g4hit.step_path,"step_path/D");
    m_tt_display_g4hits->Branch("edep",&m_display_g4hit.edep,"edep/D");
    m_tt_display_g4hits->Branch("eion",&m_display_g4hit.eion,"eion/D");
    m_tt_display_g4hits->Branch("t0",&m_display_g4hit.t0,"t0/D");
    m_tt_display_g4hits->Branch("t1",&m_display_g4hit.t1,"t1/D");
    m_tt_display_g4hits->Branch("used_in_track",&m_display_g4hit.used_in_track,"used_in_track/I");
    m_tt_display_g4hits->Branch("hitid",&m_display_g4hit.hitid,"hitid/l");
  }
return 0;
}

int TpcLaserDNL::InitRun(PHCompositeNode* topNode)
{
  if(m_use_reco_seeds)
  {
    m_track_map = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  }
  m_hitsets = findNode::getClass<TrkrHitSetContainer>(topNode,"TRKR_HITSET");
  m_geom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  m_acts = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  m_truth_tracks = findNode::getClass<TrkrTruthTrackContainer>(topNode, "TRKR_TRUTHTRACKCONTAINER");
  m_g4hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  if(m_use_clusters)
  {
    m_clusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  const bool have_reco_tracks = (m_use_reco_seeds && m_track_map && !m_track_map->empty());
  const bool have_truth_tracks = (m_truth_tracks != nullptr && m_g4hits != nullptr);
  if(!have_reco_tracks && !have_truth_tracks)
  {
    std::cout << Name() << ": missing track sources. SvtxTrackMap="
              << (m_track_map!=nullptr)
              << " TRKR_TRUTHTRACKCONTAINER=" << (m_truth_tracks!=nullptr)
              << " G4HIT_TPC=" << (m_g4hits!=nullptr)
              << std::endl;
    return -1;
  }
  if(!m_geom || !m_acts)
  {
    std::cout << Name() << ": missing geometry/Acts nodes. Geom="
              << (m_geom!=nullptr)
              << " Acts=" << (m_acts!=nullptr) << std::endl;
    return -1;
  }
  if(!m_use_clusters && !m_hitsets)
  {
    std::cout << Name() << ": missing TRKR_HITSET node." << std::endl;
    return -1;
  }
  if(m_use_clusters && !m_clusters)
  {
    std::cout << Name() << ": missing TRKR_CLUSTER node." << std::endl;
    return -1;
  }
  std::cout << Name() << ": InitRun OK. mode="
            << (m_use_clusters? "clusters" : "hits")
            << ", has_reco=" << have_reco_tracks
            << ", has_truth=" << have_truth_tracks
            << ", vdrift=" << m_acts->get_drift_velocity() << " cm/ns" << std::endl;
  return 0;
}

static inline double sqr(double x){ return x*x; }

bool TpcLaserDNL::cylinder_intersection(double x0,double y0,double z0,
                                        double vx,double vy,double vz,
                                        double R, double& t_out,
                                        double& xi,double& yi,double& zi)
{
  const double a = sqr(vx)+sqr(vy);
  const double b = 2.0*(vx*x0 + vy*y0);
  const double c = sqr(x0)+sqr(y0) - sqr(R);
  const double disc = b*b - 4.0*a*c;
  if(a==0 || disc<0) return false;
  const double s = std::sqrt(disc);
  const double t1 = (-b + s)/(2.0*a);
  const double t2 = (-b - s)/(2.0*a);
  // choose smallest positive
  double t = 1e300;
  if(t1>0) t = std::min(t,t1);
  if(t2>0) t = std::min(t,t2);
  if(!std::isfinite(t) || t<=0) return false;
  t_out = t;
  xi = x0 + t*vx;
  yi = y0 + t*vy;
  zi = z0 + t*vz;
  return true;
}

double TpcLaserDNL::wrap_dphi(double d)
{
  while(d >  M_PI) d -= 2*M_PI;
  while(d < -M_PI) d += 2*M_PI;
  return d;
}

namespace
{
  struct FitPoint
  {
    Eigen::Vector3d pos{Eigen::Vector3d::Zero()};
    double weight{0.};
    unsigned int layer{0};
    int side{0};
  };

  bool weighted_line_fit(const std::vector<FitPoint>& points,
                         Eigen::Vector3d& origin,
                         Eigen::Vector3d& direction)
  {
    origin.setZero();
    direction.setZero();
    if(points.size() < 2) return false;

    double wsum = 0.;
    for(const auto& p : points)
    {
      if(p.weight <= 0) continue;
      origin += p.weight * p.pos;
      wsum += p.weight;
    }
    if(wsum <= 0) return false;
    origin /= wsum;

    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for(const auto& p : points)
    {
      if(p.weight <= 0) continue;
      const Eigen::Vector3d diff = p.pos - origin;
      cov += p.weight * (diff * diff.transpose());
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
    if(solver.info() != Eigen::Success) return false;
    direction = solver.eigenvectors().col(2);
    const double norm = direction.norm();
    if(norm < 1e-9) return false;
    direction /= norm;
    return true;
  }

  bool select_cylinder_intersection(const Eigen::Vector3d& origin,
                                    const Eigen::Vector3d& direction,
                                    double radius,
                                    int side_filter,
                                    const Eigen::Vector3d* target_point,
                                    Eigen::Vector3d& intersection)
  {
    const double a = direction.x()*direction.x() + direction.y()*direction.y();
    const double b = 2.0*(direction.x()*origin.x() + direction.y()*origin.y());
    const double c = origin.x()*origin.x() + origin.y()*origin.y() - radius*radius;

    if(std::abs(a) < 1e-12) return false;
    double disc = b*b - 4.0*a*c;
    if(disc < 0) return false;
    disc = std::max(0.0, disc);
    const double sdisc = std::sqrt(disc);

    const double t_candidates[2] = {
      (-b - sdisc)/(2.0*a),
      (-b + sdisc)/(2.0*a)
    };

    bool found = false;
    double best_metric = std::numeric_limits<double>::infinity();
    for(double t : t_candidates)
    {
      if(!std::isfinite(t)) continue;
      const Eigen::Vector3d candidate = origin + t*direction;
      if(side_filter == 0 && candidate.z() > 0) continue;
      if(side_filter == 1 && candidate.z() < 0) continue;
      double metric = std::abs(t);
      if(target_point)
      {
        metric = (candidate - *target_point).squaredNorm();
      }
      if(metric < best_metric)
      {
        best_metric = metric;
        intersection = candidate;
        found = true;
      }
    }
    return found;
  }
}

int TpcLaserDNL::process_event(PHCompositeNode* topNode)
{
  static int ievt = 0;

  if (auto* eh = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    m_event = eh->get_EvtSequence();
  }
  else
  {
    m_event = ievt++;
  }

  // quick visibility of container content each event (only meaningful when using hits)
  std::size_t tpc_hitset_count = 0;
  if(!m_use_clusters && m_hitsets)
  {
    auto hr = m_hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
    for(auto it = hr.first; it != hr.second; ++it) ++tpc_hitset_count;
  }
  std::cout << Name() << ": process_event evt=" << m_event
            << " tpc_hitsets=" << tpc_hitset_count
            << " use_clusters=" << (m_use_clusters?1:0)
            << std::endl;

  std::vector<TrackSeed> seeds;
  if(m_use_reco_seeds && m_track_map && !m_track_map->empty())
  {
    build_reco_seeds(seeds);
  }
  if(seeds.empty())
  {
    build_truth_seeds(seeds);
  }
  // fallback: if no seeds from truth and a track map exists, try reco seeds even if disabled
  if(seeds.empty() && m_track_map && !m_track_map->empty())
  {
    build_reco_seeds(seeds);
  }

  if(seeds.empty())
  {
    std::cout << Name() << ": evt " << m_event << " has no truth track seeds (truth container present="
              << (m_truth_tracks!=nullptr)
              << ", g4hits present=" << (m_g4hits!=nullptr)
              << ")" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  for(const auto& seed : seeds)
  {
    m_trkid = seed.id;

    struct LayerRecoResult
    {
      unsigned int layer{0};
      int side{0};
      double r{0.};
      double phi_true{0.};
      double phi_reco{std::numeric_limits<double>::quiet_NaN()};
      double dphi{std::numeric_limits<double>::quiet_NaN()};
      double dRphi{std::numeric_limits<double>::quiet_NaN()};
      int nused{0};
      int nhit_scanned{0};
      double adcsum{0.};
      double xtrue{0.};
      double ytrue{0.};
      double ztrue{0.};
      double xreco{std::numeric_limits<double>::quiet_NaN()};
      double yreco{std::numeric_limits<double>::quiet_NaN()};
      double zreco{std::numeric_limits<double>::quiet_NaN()};
      int npad_used{0};
      int ntbin_used{0};
      int nbins_used{0};
      double phi_pad_max{std::numeric_limits<double>::quiet_NaN()};
      double phase{std::numeric_limits<double>::quiet_NaN()};
      double phase_reco{std::numeric_limits<double>::quiet_NaN()};
      std::vector<double> pad_phi_centers;
      std::vector<ULong64_t> hitkeys;
      std::vector<ULong64_t> hitsetkeys;
      std::vector<unsigned int> iphi;
      std::vector<unsigned int> tbin;
      std::vector<double> hit_charge;
      double total_charge_layer{0.};
      double max_charge_layer{0.};
      double weight_sum{0.};
      double xfit{std::numeric_limits<double>::quiet_NaN()};
      double yfit{std::numeric_limits<double>::quiet_NaN()};
      double zfit{std::numeric_limits<double>::quiet_NaN()};
      double phi_fit{std::numeric_limits<double>::quiet_NaN()};
      double dphi_fit{std::numeric_limits<double>::quiet_NaN()};
      double dRphi_fit{std::numeric_limits<double>::quiet_NaN()};
    };

  std::vector<LayerRecoResult> layer_results;
  layer_results.reserve(seed.layers.size());
    const bool do_fit_layer44 = m_enable_fit_layer44_residuals;
    constexpr unsigned int kTargetLayer = 44;
    std::vector<FitPoint> fit_points;
    if(do_fit_layer44) fit_points.reserve(seed.layers.size());

    for(const auto& layerPoint : seed.layers)
    {
      auto* layergeom = m_geom->GetLayerCellGeom(static_cast<int>(layerPoint.layer));
      if(!layergeom) continue;

      const double x0 = layerPoint.x;
      const double y0 = layerPoint.y;
      const double z0 = layerPoint.z;
      const double vx = layerPoint.dirx;
      const double vy = layerPoint.diry;
      const double vz = layerPoint.dirz;
      const double v2 = vx*vx + vy*vy + vz*vz;
      if(v2 == 0) continue;

      const unsigned int layer = layerPoint.layer;
      const int side = layerPoint.side;
      m_layer = layer;
      if (layerPoint.radius <= 0)
      {
        std::cout << Name() << ": layer " << layer
                  << " track " << seed.id
                  << " has non-positive radius " << layerPoint.radius
                  << ", falling back to geometry radius." << std::endl;
      }
      const double radius = layerPoint.radius > 0 ? layerPoint.radius : layergeom->get_radius();
      const double phi_true = std::atan2(y0, x0);

      if(!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0))
      {
        continue;
      }

      if(
        m_write_display_ntuple &&
        m_tt_display_intersections &&
        (layerPoint.from_g4hit || m_include_fallback_intersections)
      )
      {
        m_display_intersection.event = m_event;
        m_display_intersection.trackid = seed.id;
        m_display_intersection.layer = static_cast<int>(layer);
        m_display_intersection.side = side;
        m_display_intersection.gx = x0;
        m_display_intersection.gy = y0;
        m_display_intersection.gz = z0;
        m_display_intersection.r = radius;
        m_display_intersection.phi = phi_true;
        m_display_intersection.path = layerPoint.path;
        m_display_intersection.used_in_seed = 1;
        m_tt_display_intersections->Fill();
      }

      const double AdcClockPeriod = layergeom->get_zstep();
      const unsigned short NTBins = static_cast<unsigned short>(layergeom->get_zbins());
      const double tdriftmax = AdcClockPeriod * NTBins / 2.0;
      const double vdrift = m_acts->get_drift_velocity();

      LayerRecoResult result;
      result.layer = layer;
      result.side = side;
      result.r = radius;
      result.phi_true = phi_true;
      result.xtrue = x0;
      result.ytrue = y0;
      result.ztrue = z0;

      double wx = 0;
      double wy = 0;
      double wz = 0;
      double wsum = 0;
      std::map<unsigned short, double> padWeights;
      std::vector<ULong64_t> hitkeys;
      std::vector<ULong64_t> hitsetkeys_vec;
      std::vector<unsigned int> iphi_vec;
      std::vector<unsigned int> tbin_vec;
      std::vector<double> hit_charge;
      std::set<unsigned int> unique_tbins;

      if(!m_use_clusters)
      {
        if(!m_hitsets) continue;
        TrkrHitSetContainer::ConstRange hr = m_hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
        for(auto hsit = hr.first; hsit != hr.second; ++hsit)
        {
          const TrkrDefs::hitsetkey& hsk = hsit->first;
          if(TrkrDefs::getLayer(hsk) != layer) continue;
          if(TpcDefs::getSide(hsk) != static_cast<unsigned int>(side)) continue;

          TrkrHitSet* hitset = hsit->second;
          if(!hitset) continue;

          TrkrHitSet::ConstRange hits = hitset->getHits();
          for(auto hitit = hits.first; hitit != hits.second; ++hitit)
          {
            const auto hitkey = hitit->first;
            TrkrHit* hit = hitit->second;
            if(!hit) continue;
            ++result.nhit_scanned;

            double weight = m_weight_by_adc ? static_cast<double>(hit->getAdc())
                                            : static_cast<double>(hit->getEnergy());
            if(m_weight_by_adc && m_use_pedestal) weight -= m_pedestal;
            if(weight < m_min_adc) continue;

            const unsigned short iphi = TpcDefs::getPad(hitkey);
            const unsigned short tbin = TpcDefs::getTBin(hitkey);

            const double phi_c = layergeom->get_phicenter(static_cast<int>(iphi), side);
            const double xh = radius*std::cos(phi_c);
            const double yh = radius*std::sin(phi_c);

            const double zcenter = layergeom->get_zcenter(tbin);
            double zdriftlen = zcenter * vdrift;
            double zh = tdriftmax * vdrift - zdriftlen;
            if(side == 0) zh = -zh;

            const double ox = xh - x0;
            const double oy = yh - y0;
            const double oz = zh - z0;
            const double tproj = (vx*ox + vy*oy + vz*oz) / v2;
            const double px = x0 + tproj*vx;
            const double py = y0 + tproj*vy;
            const double pz = z0 + tproj*vz;
            const double dca = std::sqrt(sqr(xh-px)+sqr(yh-py)+sqr(zh-pz));

            const double dzline = zh - z0;

            if(dca > m_max_dca) continue;
            if(std::abs(dzline) > m_max_dz) continue;

            wx += weight * xh;
            wy += weight * yh;
            wz += weight * zh;
            wsum += weight;
            result.adcsum += weight;
            result.nused++;

            hitkeys.push_back(static_cast<ULong64_t>(hitkey));
            hitsetkeys_vec.push_back(static_cast<ULong64_t>(hsk));
            iphi_vec.push_back(static_cast<unsigned int>(iphi));
            tbin_vec.push_back(static_cast<unsigned int>(tbin));
            hit_charge.push_back(weight);
            unique_tbins.insert(static_cast<unsigned int>(tbin));
            result.total_charge_layer += weight;
            if (weight > result.max_charge_layer) result.max_charge_layer = weight;
            padWeights[iphi] += weight;
          }
        }
      }
      else
      {
        if(!m_clusters) continue;
        auto hitsetkeys = m_clusters->getHitSetKeys(TrkrDefs::TrkrId::tpcId);
        for(const auto& hsk : hitsetkeys)
        {
          if(TrkrDefs::getLayer(hsk) != layer) continue;
          if(TpcDefs::getSide(hsk) != static_cast<unsigned int>(side)) continue;

          auto crange = m_clusters->getClusters(hsk);
          for(auto cit = crange.first; cit != crange.second; ++cit)
          {
            const auto ckey = cit->first;
            TrkrCluster* clus = cit->second;
            if(!clus) continue;

            double weight = static_cast<double>(clus->getAdc());
            if(weight < m_min_adc) continue;

            Acts::Vector3 g = m_acts->getGlobalPosition(ckey, clus);

            const double ox = g.x() - x0;
            const double oy = g.y() - y0;
            const double oz = g.z() - z0;
            const double tproj = (vx*ox + vy*oy + vz*oz) / v2;
            const double px = x0 + tproj*vx;
            const double py = y0 + tproj*vy;
            const double pz = z0 + tproj*vz;
            const double dca = std::sqrt(sqr(g.x()-px)+sqr(g.y()-py)+sqr(g.z()-pz));

            const double dzline = g.z() - z0;

            if(dca > m_max_dca) continue;
            if(std::abs(dzline) > m_max_dz) continue;

            wx += weight * g.x();
            wy += weight * g.y();
            wz += weight * g.z();
            wsum += weight;
            result.adcsum += weight;
            result.nused++;
            hit_charge.push_back(weight);
            result.total_charge_layer += weight;
            if (weight > result.max_charge_layer) result.max_charge_layer = weight;
          }
        }
      }

      if(wsum > 0)
      {
        result.xreco = wx/wsum;
        result.yreco = wy/wsum;
        result.zreco = wz/wsum;
        result.phi_reco = std::atan2(result.yreco, result.xreco);
        result.dphi = wrap_dphi(result.phi_reco - phi_true);
        result.dRphi = radius * result.dphi;
        result.weight_sum = wsum;

        if(!m_use_clusters && !padWeights.empty())
        {
          double maxWeight = -std::numeric_limits<double>::infinity();
          for(const auto& entry : padWeights)
          {
            const auto pad = entry.first;
            const auto weight = entry.second;
            const double phi_c = layergeom->get_phicenter(static_cast<int>(pad), side);
            result.pad_phi_centers.push_back(phi_c);
            if(weight > maxWeight)
            {
              maxWeight = weight;
              result.phi_pad_max = phi_c;
            }
          }
          result.npad_used = static_cast<int>(padWeights.size());

          const double phi_width = std::abs(layergeom->get_phistep());
          if(phi_width > 1e-12 && std::isfinite(result.phi_pad_max))
          {
            const double dphi_phase_true = wrap_dphi(phi_true - result.phi_pad_max);
            const double dphi_phase_reco = wrap_dphi(result.phi_reco - result.phi_pad_max);
            result.phase = dphi_phase_true / phi_width;
            result.phase_reco = dphi_phase_reco / phi_width;
          }
          else
          {
            result.phase = std::numeric_limits<double>::quiet_NaN();
            result.phase_reco = std::numeric_limits<double>::quiet_NaN();
          }
        }

        result.hitkeys = std::move(hitkeys);
        result.hitsetkeys = std::move(hitsetkeys_vec);
        result.iphi = std::move(iphi_vec);
        result.tbin = std::move(tbin_vec);
        result.hit_charge = std::move(hit_charge);
        result.ntbin_used = static_cast<int>(unique_tbins.size());
        result.nbins_used = static_cast<int>(result.hitkeys.size());
        layer_results.push_back(std::move(result));
        if(do_fit_layer44 && layer != kTargetLayer)
        {
          // use equal weighting across clusters for the fit
          fit_points.push_back({Eigen::Vector3d(layer_results.back().xreco, layer_results.back().yreco, layer_results.back().zreco), 1.0, layer, side});
        }
      }
    } // layer loop

    // perform straight-line fit to the reconstructed clusters (optional)
    Eigen::Vector3d intersection{Eigen::Vector3d::Zero()};
    bool have_intersection = false;
    if(do_fit_layer44)
    {
      Eigen::Vector3d fit_origin{Eigen::Vector3d::Zero()};
      Eigen::Vector3d fit_dir{Eigen::Vector3d::Zero()};
      const bool have_fit = weighted_line_fit(fit_points, fit_origin, fit_dir);
      auto* target_geom = m_geom->GetLayerCellGeom(static_cast<int>(kTargetLayer));
      Eigen::Vector3d target_cluster = Eigen::Vector3d::Zero();
      bool have_target_cluster = false;
      int target_side = -1;
      for(const auto& res : layer_results)
      {
        if(res.layer == kTargetLayer)
        {
          target_cluster = Eigen::Vector3d(res.xreco, res.yreco, res.zreco);
          have_target_cluster = true;
          target_side = res.side;
          break;
        }
      }

      have_intersection = have_fit && target_geom &&
        select_cylinder_intersection(
          fit_origin,
          fit_dir,
          target_geom->get_radius(),
          target_side,
          have_target_cluster ? &target_cluster : nullptr,
          intersection);
    }

    for(auto& res : layer_results)
    {
      if(res.layer == kTargetLayer && have_intersection)
      {
        res.xfit = intersection.x();
        res.yfit = intersection.y();
        res.zfit = intersection.z();
        res.phi_fit = std::atan2(intersection.y(), intersection.x());
        res.dphi_fit = wrap_dphi(res.phi_reco - res.phi_fit);
        res.dRphi_fit = res.r * res.dphi_fit;
      }
      else
      {
        res.xfit = std::numeric_limits<double>::quiet_NaN();
        res.yfit = std::numeric_limits<double>::quiet_NaN();
        res.zfit = std::numeric_limits<double>::quiet_NaN();
        res.phi_fit = std::numeric_limits<double>::quiet_NaN();
        res.dphi_fit = std::numeric_limits<double>::quiet_NaN();
        res.dRphi_fit = std::numeric_limits<double>::quiet_NaN();
      }

      m_layer = res.layer;
      m_side = res.side;
      m_r = res.r;
      m_phi_true = res.phi_true;
      m_phi_reco = res.phi_reco;
      m_dphi = res.dphi;
      m_dRphi = res.dRphi;
      m_nused = res.nused;
      m_nhit_scanned = res.nhit_scanned;
      m_adcsum = res.adcsum;
      m_xtrue = res.xtrue;
      m_ytrue = res.ytrue;
      m_ztrue = res.ztrue;
      m_xreco = res.xreco;
      m_yreco = res.yreco;
      m_zreco = res.zreco;
      m_npad_used = res.npad_used;
      m_ntbin_used = res.ntbin_used;
      m_nbins_used = res.nbins_used;
      m_phi_pad_max = res.phi_pad_max;
      m_phase = res.phase;
      m_phase_reco = res.phase_reco;
      m_pad_phi_centers = res.pad_phi_centers;
      m_hitkeys = res.hitkeys;
      m_hitsetkeys = res.hitsetkeys;
      m_iphi = res.iphi;
      m_tbin = res.tbin;
      m_hit_charge = res.hit_charge;
      m_total_charge_layer = res.total_charge_layer;
      m_max_charge_layer = res.max_charge_layer;
      m_xfit_layer44 = res.xfit;
      m_yfit_layer44 = res.yfit;
      m_zfit_layer44 = res.zfit;
      m_phi_fit_layer44 = res.phi_fit;
      m_dphi_fit_layer44 = res.dphi_fit;
      m_dRphi_fit_layer44 = res.dRphi_fit;
      m_tt->Fill();
    }
  } // track loop

  if(m_write_display_ntuple && m_tt_display_g4hits && m_g4hits)
  {
    std::unordered_set<int> truth_ids;
    truth_ids.reserve(seeds.size());
    for(const auto& seed : seeds)
    {
      truth_ids.insert(seed.id);
    }

    if(!truth_ids.empty())
    {
      const unsigned int subsamples = std::max(1u, m_display_hit_subsamples);
      const double denom = static_cast<double>(subsamples);
      PHG4HitContainer::ConstRange hitrange = m_g4hits->getHits();
      for(auto hitit = hitrange.first; hitit != hitrange.second; ++hitit)
      {
        PHG4Hit* hit = hitit->second;
        if(!hit) continue;

        const int trkid = hit->get_trkid();
        if(trkid < 0) continue;
        if(truth_ids.find(trkid) == truth_ids.end()) continue;

        const unsigned int layer = hit->get_layer();
        if(layer == std::numeric_limits<unsigned int>::max()) continue;

        const double x0 = hit->get_x(0);
        const double y0 = hit->get_y(0);
        const double z0 = hit->get_z(0);
        const double x1 = hit->get_x(1);
        const double y1 = hit->get_y(1);
        const double z1 = hit->get_z(1);

        const double dx = x1 - x0;
        const double dy = y1 - y0;
        const double dz = z1 - z0;
        const double path = std::sqrt(dx*dx + dy*dy + dz*dz);

        const double edep = hit->get_edep();
        const double eion = hit->get_eion();
        const double t0 = hit->get_t(0);
        const double t1 = hit->get_t(1);
        const ULong64_t hitid = static_cast<ULong64_t>(hit->get_hit_id());

        for(unsigned int isample = 0; isample <= subsamples; ++isample)
        {
          const double frac = (subsamples > 0) ? static_cast<double>(isample) / denom : 0.0;
          const double gx = x0 + frac * dx;
          const double gy = y0 + frac * dy;
          const double gz = z0 + frac * dz;
          if(!std::isfinite(gx) || !std::isfinite(gy) || !std::isfinite(gz)) continue;
          const double r = std::sqrt(gx*gx + gy*gy);

          m_display_g4hit.event = m_event;
          m_display_g4hit.trackid = trkid;
          m_display_g4hit.layer = static_cast<int>(layer);
          m_display_g4hit.side = (gz >= 0) ? 1 : 0;
          m_display_g4hit.gx = gx;
          m_display_g4hit.gy = gy;
          m_display_g4hit.gz = gz;
          m_display_g4hit.r = r;
          m_display_g4hit.phi = std::atan2(gy, gx);
          m_display_g4hit.sample = static_cast<int>(isample);
          m_display_g4hit.sample_frac = frac;
          m_display_g4hit.step_path = path * frac;
          m_display_g4hit.edep = edep;
          m_display_g4hit.eion = eion;
          m_display_g4hit.t0 = t0;
          m_display_g4hit.t1 = t1;
          m_display_g4hit.used_in_track = 1;
          m_display_g4hit.hitid = hitid;

          m_tt_display_g4hits->Fill();
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


void TpcLaserDNL::build_reco_seeds(std::vector<TrackSeed>& seeds) const
{
  if(!m_use_reco_seeds || !m_track_map || !m_geom) return;

  for(const auto& it : *m_track_map)
  {
    const auto* trk = it.second;
    if(!trk) continue;

    TrackSeed seed;
    seed.id = trk->get_id();
    seed.origin[0] = trk->get_x();
    seed.origin[1] = trk->get_y();
    seed.origin[2] = trk->get_z();

    const double base_x = trk->get_x();
    const double base_y = trk->get_y();
    const double base_z = trk->get_z();
    const double vx = trk->get_px();
    const double vy = trk->get_py();
    const double vz = trk->get_pz();
    const double v2 = vx*vx + vy*vy + vz*vz;
    if(v2 == 0) continue;
    const double vmag = std::sqrt(v2);
    if(vmag > 0)
    {
      seed.dir[0] = vx / vmag;
      seed.dir[1] = vy / vmag;
      seed.dir[2] = vz / vmag;
      seed.dir_valid = true;
    }

    auto lr = m_geom->get_begin_end();
    for(auto lit = lr.first; lit != lr.second; ++lit)
    {
      auto* layergeom = lit->second;
      if(!layergeom) continue;

      const double R = layergeom->get_radius();
      double tR{}, xi{}, yi{}, zi{};
      if(!cylinder_intersection(base_x, base_y, base_z, vx, vy, vz, R, tR, xi, yi, zi)) continue;

      LayerPoint lp;
      lp.layer = layergeom->get_layer();
      lp.radius = R;
      lp.x = xi;
      lp.y = yi;
      lp.z = zi;
      lp.dirx = vx;
      lp.diry = vy;
      lp.dirz = vz;
      lp.side = (zi > 0) ? 1 : 0;
      lp.path = std::sqrt(v2);
      lp.from_g4hit = true;
      seed.layers.push_back(lp);
    }

    if(!seed.layers.empty())
    {
      std::sort(seed.layers.begin(), seed.layers.end(),
                [](const LayerPoint& a, const LayerPoint& b){ return a.layer < b.layer; });
      seeds.push_back(std::move(seed));
    }
  }
}

void TpcLaserDNL::build_truth_seeds(std::vector<TrackSeed>& seeds) const
{
  if(!m_truth_tracks || !m_geom || !m_g4hits) return;

  TrkrTruthTrackContainer::ConstRange range = m_truth_tracks->getTruthTrackRange();
  std::unordered_map<unsigned int, std::size_t> index_by_id;
  seeds.reserve(std::distance(range.first, range.second));
  for(auto it = range.first; it != range.second; ++it)
  {
    TrkrTruthTrack* truth = it->second;
    if(!truth) continue;
    TrackSeed seed;
    seed.id = static_cast<int>(truth->getTrackid());
    seed.origin[0] = truth->getX0();
    seed.origin[1] = truth->getY0();
    seed.origin[2] = truth->getZ0();
    const double pt = truth->getPt();
    const double phi = truth->getPhi();
    const double eta = truth->getPseudoRapidity();
    const double px = pt * std::cos(phi);
    const double py = pt * std::sin(phi);
    const double pz = pt * std::sinh(eta);
    const double pmag = std::sqrt(px*px + py*py + pz*pz);
    if(pmag > 0)
    {
      seed.dir[0] = px / pmag;
      seed.dir[1] = py / pmag;
      seed.dir[2] = pz / pmag;
      seed.dir_valid = true;
    }
    index_by_id[truth->getTrackid()] = seeds.size();
    seeds.push_back(std::move(seed));
  }

  if(seeds.empty()) return;

  std::vector<unsigned int> geometry_layers;
  {
    auto lr = m_geom->get_begin_end();
    for(auto it = lr.first; it != lr.second; ++it)
    {
      if(it->second)
      {
        geometry_layers.push_back(static_cast<unsigned int>(it->second->get_layer()));
      }
    }
  }

  std::unordered_map<unsigned int, std::map<unsigned int, LayerPoint>> layer_cache;
  const double axial_threshold = 1e-9;
  const double param_tolerance = 1e-6;
  const double radial_tolerance = 5e-3;
  PHG4HitContainer::ConstRange hitrange = m_g4hits->getHits();
  for(auto hitit = hitrange.first; hitit != hitrange.second; ++hitit)
  {
    PHG4Hit* hit = hitit->second;
    if(!hit) continue;

    const int trkid = hit->get_trkid();
    if(trkid < 0) continue;

    auto idxIt = index_by_id.find(static_cast<unsigned int>(trkid));
    if(idxIt == index_by_id.end())
    {
      // fallback: create a seed directly from the first G4 hit we see for this trkid
      TrackSeed seed;
      seed.id = trkid;
      seed.origin[0] = hit->get_x(0);
      seed.origin[1] = hit->get_y(0);
      seed.origin[2] = hit->get_z(0);
      const double vx = hit->get_x(1) - hit->get_x(0);
      const double vy = hit->get_y(1) - hit->get_y(0);
      const double vz = hit->get_z(1) - hit->get_z(0);
      const double v2 = vx*vx + vy*vy + vz*vz;
      if(v2 > 0)
      {
        const double vmag = std::sqrt(v2);
        seed.dir[0] = vx / vmag;
        seed.dir[1] = vy / vmag;
        seed.dir[2] = vz / vmag;
        seed.dir_valid = true;
      }
      index_by_id[static_cast<unsigned int>(trkid)] = seeds.size();
      seeds.push_back(std::move(seed));
      idxIt = index_by_id.find(static_cast<unsigned int>(trkid));
    }

    unsigned int layer = hit->get_layer();
    if(layer == std::numeric_limits<unsigned int>::max()) continue;

    // Optimization: Check layers in the vicinity [layer-2, layer+2]
    // TPC layers are typically 7 to 54 (or similar range depending on geometry)
    // We clamp the range to avoid invalid lookups
    const int layer_center = static_cast<int>(layer);
    const int layer_min = std::max(0, layer_center - 2);
    const int layer_max = layer_center + 2;

    for(int ilayer = layer_min; ilayer <= layer_max; ++ilayer)
    {
      auto* layergeom = m_geom->GetLayerCellGeom(ilayer);
      if(!layergeom) continue;

      const double x0 = hit->get_x(0);
      const double y0 = hit->get_y(0);
      const double z0 = hit->get_z(0);
      const double x1 = hit->get_x(1);
      const double y1 = hit->get_y(1);
      const double z1 = hit->get_z(1);

      const double dx = x1 - x0;
      const double dy = y1 - y0;
      const double dz = z1 - z0;
      const double path = std::sqrt(dx*dx + dy*dy + dz*dz);
      if(path < 1e-6) continue;

      const double R = layergeom->get_radius();
      const double a = dx*dx + dy*dy;
      const double b = 2.0*(dx*x0 + dy*y0);
      const double c = x0*x0 + y0*y0 - R*R;

      LayerPoint candidate;
      bool have_candidate = false;
      double best_residual = std::numeric_limits<double>::max();

      auto try_record = [&](double t)
      {
        const double xi = x0 + t*dx;
        const double yi = y0 + t*dy;
        const double zi = z0 + t*dz;
        const double radial_residual = std::abs(std::sqrt(xi*xi + yi*yi) - R);
        if(radial_residual > radial_tolerance) return;
        if(!have_candidate || radial_residual < best_residual)
        {
          candidate.layer = static_cast<unsigned int>(ilayer);
          candidate.radius = R;
          candidate.x = xi;
          candidate.y = yi;
          candidate.z = zi;
          candidate.dirx = dx;
          candidate.diry = dy;
          candidate.dirz = dz;
          candidate.side = (zi > 0) ? 1 : 0;
          candidate.path = path;
          candidate.from_g4hit = true;
          best_residual = radial_residual;
          have_candidate = true;
        }
      };

      if(a < axial_threshold)
      {
        const double r0 = std::sqrt(x0*x0 + y0*y0);
        const double r1 = std::sqrt(x1*x1 + y1*y1);
        const double res0 = std::abs(r0 - R);
        const double res1 = std::abs(r1 - R);
        if(res0 <= radial_tolerance || res1 <= radial_tolerance)
        {
          const double t = (res0 <= res1) ? 0.0 : 1.0;
          try_record(t);
        }
      }
      else
      {
        double disc = b*b - 4.0*a*c;
        if(disc >= -1e-12)
        {
          disc = std::max(0.0, disc);
          const double sqrt_disc = std::sqrt(disc);
          const double t_candidates[2] = {
              (-b - sqrt_disc) / (2.0*a),
              (-b + sqrt_disc) / (2.0*a)};

          for(double t_candidate : t_candidates)
          {
            if(t_candidate < -param_tolerance || t_candidate > 1.0 + param_tolerance) continue;
            double t_clamped = t_candidate;
            if(t_clamped < 0.0) t_clamped = 0.0;
            if(t_clamped > 1.0) t_clamped = 1.0;
            try_record(t_clamped);
          }
        }
      }

      if(!have_candidate) continue;

      auto& layer_map = layer_cache[static_cast<unsigned int>(trkid)];
      auto layer_it = layer_map.find(static_cast<unsigned int>(ilayer));
      if(layer_it == layer_map.end() || candidate.path > layer_it->second.path)
      {
        layer_map[static_cast<unsigned int>(ilayer)] = candidate;
      }
    }
  }

  std::vector<TrackSeed> filtered;
  filtered.reserve(seeds.size());
  for(auto& seed : seeds)
  {
    auto cache_it = layer_cache.find(static_cast<unsigned int>(seed.id));
    if(cache_it == layer_cache.end()) continue;
    std::map<unsigned int, LayerPoint> completed = cache_it->second;

    if(seed.dir_valid)
    {
      const double ox = seed.origin[0];
      const double oy = seed.origin[1];
      const double oz = seed.origin[2];
      const double dx = seed.dir[0];
      const double dy = seed.dir[1];
      const double dz = seed.dir[2];
      const double dir_norm = std::sqrt(dx*dx + dy*dy + dz*dz);

      if(dir_norm > 1e-9)
      {
        const double a = dx*dx + dy*dy;

        for(const auto layer_id : geometry_layers)
        {
          if(completed.find(layer_id) != completed.end()) continue;
          auto* layergeom = m_geom->GetLayerCellGeom(static_cast<int>(layer_id));
          if(!layergeom) continue;
          const double R = layergeom->get_radius();

          const double b = 2.0 * (ox*dx + oy*dy);
          const double c = ox*ox + oy*oy - R*R;

          if(std::abs(a) < 1e-12)
          {
            continue;
          }
          const double disc = b*b - 4.0*a*c;
          if(disc < 0) continue;
          const double sqrt_disc = std::sqrt(disc);
          const double t_candidates[2] = {
              (-b - sqrt_disc) / (2.0*a),
              (-b + sqrt_disc) / (2.0*a)};

          double t_selected = std::numeric_limits<double>::infinity();
          for(double tval : t_candidates)
          {
            if(tval > 1e-6 && tval < t_selected)
            {
              t_selected = tval;
            }
          }
          if(!std::isfinite(t_selected)) continue;

          if(m_include_fallback_intersections)
          {
            LayerPoint fallback;
            fallback.layer = layer_id;
            fallback.radius = R;
            fallback.x = ox + dx * t_selected;
            fallback.y = oy + dy * t_selected;
            fallback.z = oz + dz * t_selected;
            fallback.dirx = dx;
            fallback.diry = dy;
            fallback.dirz = dz;
            fallback.side = (fallback.z > 0) ? 1 : 0;
            fallback.path = t_selected * dir_norm;
            fallback.from_g4hit = false;
            completed[layer_id] = fallback;
          }
        }
      }
    }

    if(completed.empty())
    {
      std::cout << Name() << ": track " << seed.id
                << " has no intersections recorded in layer cache" << std::endl;
      continue;
    }
    if(completed.size() < geometry_layers.size())
    {
      std::cout << Name() << ": track " << seed.id
                << " intersections " << completed.size()
                << " < expected " << geometry_layers.size()
                << ". Missing layers:";
      for(const auto layer_id : geometry_layers)
      {
        if(completed.find(layer_id) == completed.end())
        {
          std::cout << " " << layer_id;
        }
      }
      std::cout << std::endl;
    }

    seed.layers.clear();
    seed.layers.reserve(completed.size());
    for(const auto& kv : completed)
    {
      seed.layers.push_back(kv.second);
    }

    std::sort(seed.layers.begin(), seed.layers.end(),
              [](const LayerPoint& a, const LayerPoint& b){ return a.layer < b.layer; });
    filtered.push_back(std::move(seed));
  }
  seeds.swap(filtered);
}

int TpcLaserDNL::End(PHCompositeNode*)
{
if(m_tf)
{
m_tf->cd();
if(m_tt) m_tt->Write();
if(m_tt_display_intersections) m_tt_display_intersections->Write();
if(m_tt_display_g4hits) m_tt_display_g4hits->Write();
m_tf->Close();
}
return 0;
}
