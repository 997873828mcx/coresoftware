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

#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <cmath>
#include <iostream>

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
// debug vectors (only filled with hits that pass selection)
m_tt->Branch("hitkey",&m_hitkeys);
m_tt->Branch("hitsetkey",&m_hitsetkeys);
m_tt->Branch("iphi",&m_iphi);
m_tt->Branch("tbin",&m_tbin);
return 0;
}

int TpcLaserDNL::InitRun(PHCompositeNode* topNode)
{
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  m_hitsets = findNode::getClass<TrkrHitSetContainer>(topNode,"TRKR_HITSET");
  m_geom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  m_acts = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(m_use_clusters)
  {
    m_clusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  if(!(m_track_map && m_geom && m_acts && (m_use_clusters ? (m_clusters!=nullptr) : (m_hitsets!=nullptr))))
  {
    std::cout << Name() << ": missing nodes. SvtxTrackMap="
              << (m_track_map!=nullptr)
              << (m_use_clusters ? " TRKR_CLUSTER=" : " TRKR_HITSET=")
              << (m_use_clusters ? (m_clusters!=nullptr) : (m_hitsets!=nullptr))
              << " Geom=" << (m_geom!=nullptr) << " Acts=" << (m_acts!=nullptr) << std::endl;
    return -1;
  }
  std::cout << Name() << ": InitRun OK. mode="
            << (m_use_clusters? "clusters" : "hits")
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

int TpcLaserDNL::process_event(PHCompositeNode* topNode)
{
// optional EventHeader if desired; not required for tree
// m_event can be incremented locally if no EventHeader present
static int ievt=0;

// Prefer the EventHeader sequence if available; else fall back to local counter
if (auto* eh = findNode::getClass<EventHeader>(topNode, "EventHeader"))
{
  m_event = eh->get_EvtSequence();
}
else
{
  m_event = ievt++;
}

if(!m_track_map || !m_geom || !m_acts) 
{
  std::cout << "wrong" << std::endl;
  return 0;
}
if(!m_use_clusters && !m_hitsets) 
{
  std::cout << "wrong" << std::endl;
  return 0;
}
if(m_use_clusters && !m_clusters) 
{
  std::cout << "wrong" << std::endl;
  return 0;
}

// quick visibility of container content each event
std::size_t tpc_hitset_count = 0;
{
  auto hr = m_hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
  for(auto it = hr.first; it != hr.second; ++it) ++tpc_hitset_count;
}
std::cout << Name() << ": process_event evt=" << m_event
          << " tpc_hitsets=" << tpc_hitset_count
          << " use_clusters=" << (m_use_clusters?1:0)
          << std::endl;

int ntracks = 0;

for(const auto& it : *m_track_map)
{
  const auto* trk = it.second;
  if(!trk) continue;
  m_trkid = trk->get_id();
  ++ntracks;
m_trkid = trk->get_id();
const double x0 = trk->get_x();
const double y0 = trk->get_y();
const double z0 = trk->get_z();
const double vx = trk->get_px();
const double vy = trk->get_py();
const double vz = trk->get_pz();
const double v2 = vx*vx+vy*vy+vz*vz;
if(v2==0) continue;

// loop TPC layers
auto lr = m_geom->get_begin_end();
  for(auto lit = lr.first; lit != lr.second; ++lit)
  {
    PHG4TpcCylinderGeom* layergeom = lit->second;
    if(!layergeom) continue;

    const double R = layergeom->get_radius();   // layer center radius
    double tR{}, xi{}, yi{}, zi{};
    if(!cylinder_intersection(x0,y0,z0,vx,vy,vz,R,tR,xi,yi,zi)) continue;

    // truth ref at intersection
    m_layer = layergeom->get_layer();
    m_r = R;
    m_xtrue = xi; m_ytrue = yi; m_ztrue = zi;
    m_phi_true = std::atan2(yi, xi);
    m_side = (zi > 0) ? 1 : 0;

    // convert hits to global and select those near the line
    const double AdcClockPeriod = layergeom->get_zstep();
    const unsigned short NTBins = (unsigned short)layergeom->get_zbins();
    const double tdriftmax = AdcClockPeriod * NTBins / 2.0;
    const double vdrift = m_acts->get_drift_velocity();

    double wx=0, wy=0, wz=0, wsum=0;
    m_nused = 0;
    m_nhit_scanned = 0;
    m_adcsum = 0;
    m_hitkeys.clear();
    m_hitsetkeys.clear();
    m_iphi.clear();
    m_tbin.clear();

  if(!m_use_clusters)
  {
    // iterate TPC hitsets
    TrkrHitSetContainer::ConstRange hr = m_hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
    for(auto hsit = hr.first; hsit != hr.second; ++hsit)
    {
      const TrkrDefs::hitsetkey& hsk = hsit->first;
      if(TrkrDefs::getLayer(hsk) != m_layer) continue;
      if(TpcDefs::getSide(hsk) != (unsigned)m_side) continue;

      TrkrHitSet* hitset = hsit->second;
      if(!hitset) continue;

      TrkrHitSet::ConstRange hits = hitset->getHits();
      for(auto hitit = hits.first; hitit != hits.second; ++hitit)
      {
        const auto hitkey = hitit->first;
        TrkrHit* hit = hitit->second;
        if(!hit) continue;
        ++m_nhit_scanned;

        // choose weight: ADC (default) or charge (hit energy)
        double weight = m_weight_by_adc ? static_cast<double>(hit->getAdc())
                                        : static_cast<double>(hit->getEnergy());
        if(m_weight_by_adc && m_use_pedestal) weight -= m_pedestal;
        if(weight < m_min_adc) continue;

        const unsigned short iphi = TpcDefs::getPad(hitkey);
        const unsigned short tbin = TpcDefs::getTBin(hitkey);

        // Use geometry-provided pad-center calculation to avoid sector/orientation mismatches
        const double phi_c = layergeom->get_phicenter(static_cast<int>(iphi), m_side);
        const double xh = R*std::cos(phi_c);
        const double yh = R*std::sin(phi_c);

        const double zcenter = layergeom->get_zcenter(tbin);
        double zdriftlen = zcenter * vdrift;
        double zh = tdriftmax * vdrift - zdriftlen;
        if(m_side == 0) zh = -zh;

        // distance to laser line
        const double ox = xh - x0;
        const double oy = yh - y0;
        const double oz = zh - z0;
        const double tproj = (vx*ox + vy*oy + vz*oz) / v2;
        const double px = x0 + tproj*vx;
        const double py = y0 + tproj*vy;
        const double pz = z0 + tproj*vz;
        const double dca = std::sqrt(sqr(xh-px)+sqr(yh-py)+sqr(zh-pz));

        const double dzline = zh - zi; // compare to intersection z at this layer

        if(dca > m_max_dca) continue;
        if(std::abs(dzline) > m_max_dz) continue;

        wx += weight * xh;
        wy += weight * yh;
        wz += weight * zh;
        wsum += weight;
        m_adcsum += weight;
        m_nused++;
        // record used-hit debug info
        m_hitkeys.push_back(static_cast<ULong64_t>(hitkey));
        m_hitsetkeys.push_back(static_cast<ULong64_t>(hsk));
        m_iphi.push_back(static_cast<unsigned int>(iphi));
        m_tbin.push_back(static_cast<unsigned int>(tbin));
      }
    }
  }
  else
  {
    // iterate TPC clusters grouped by hitset
    auto hitsetkeys = m_clusters->getHitSetKeys(TrkrDefs::TrkrId::tpcId);
    for(const auto& hsk : hitsetkeys)
    {
      if(TrkrDefs::getLayer(hsk) != m_layer) continue;
      if(TpcDefs::getSide(hsk) != (unsigned)m_side) continue;

      auto crange = m_clusters->getClusters(hsk);
      for(auto cit = crange.first; cit != crange.second; ++cit)
      {
        const auto ckey = cit->first;
        TrkrCluster* clus = cit->second;
        if(!clus) continue;

        // Use cluster ADC as weight
        double weight = static_cast<double>(clus->getAdc());
        if(weight < m_min_adc) continue;

        // global cluster position
        Acts::Vector3 g = m_acts->getGlobalPosition(ckey, clus);

        // distance to laser line (3D)
        const double ox = g.x() - x0;
        const double oy = g.y() - y0;
        const double oz = g.z() - z0;
        const double tproj = (vx*ox + vy*oy + vz*oz) / v2;
        const double px = x0 + tproj*vx;
        const double py = y0 + tproj*vy;
        const double pz = z0 + tproj*vz;
        const double dca = std::sqrt(sqr(g.x()-px)+sqr(g.y()-py)+sqr(g.z()-pz));

        const double dzline = g.z() - zi; // compare to intersection z at this layer

        if(dca > m_max_dca) continue;
        if(std::abs(dzline) > m_max_dz) continue;

        wx += weight * g.x();
        wy += weight * g.y();
        wz += weight * g.z();
        wsum += weight;
        m_adcsum += weight;
        m_nused++;
      }
    }
  }

  // Always write one row per layer if any hits were scanned.
  // If no hits passed selection (wsum==0), keep reco equal to truth and adcsum=0 for diagnostics.
  //if (m_nhit_scanned > 0)
  //{
    if(wsum > 0)
    {
      m_xreco = wx/wsum;
      m_yreco = wy/wsum;
      m_zreco = wz/wsum;
      m_phi_reco = std::atan2(m_yreco, m_xreco);
      m_dphi = wrap_dphi(m_phi_reco - m_phi_true);
      m_dRphi = m_r * m_dphi;
      m_tt->Fill();
    }
    
   
    /* std::cout << Name() << ": evt " << m_event
              << " trk " << m_trkid
              << " layer " << m_layer
              << " side " << m_side
              << " scanned " << m_nhit_scanned
              << " used " << m_nused
              << " adcsum " << m_adcsum
              << " wsum " << wsum
              << std::endl; */
  //}
  //else
  //{
    // no hits seen for this layer for this track; print a minimal hint when containers had content
  //  if (tpc_hitset_count > 0)
   // {
   //   std::cout << Name() << ": evt " << m_event
    //            << " trk " << m_trkid
    //            << " layer " << m_layer
    //            << " side " << m_side
    //            << " scanned=0 (no matching hitsets/hits)"
    //            << std::endl;
   // }
  //}
} // layers
} // tracks

if (ntracks==0)
{
  std::cout << Name() << ": evt " << m_event << " has no tracks in SvtxTrackMap" << std::endl;
}

return 0;
}

int TpcLaserDNL::End(PHCompositeNode*)
{
if(m_tf)
{
m_tf->cd();
if(m_tt) m_tt->Write();
m_tf->Close();
}
return 0;
}
