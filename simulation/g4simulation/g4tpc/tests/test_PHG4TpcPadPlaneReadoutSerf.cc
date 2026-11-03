#include "PHG4TpcPadPlaneReadout.h"

#include <g4detectors/PHG4CylinderCellDefs.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <g4detectors/PHG4TpcGeomv2.h>
#include <g4main/PHG4Hitv1.h>
#include <phool/recoConsts.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
  using Key = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;
  using HitMap = std::map<Key, unsigned int>;
  constexpr unsigned long seed = 932771;
  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }
  HitMap hits(const TrkrHitSetContainer& container)
  {
    HitMap result;
    const auto sets = container.getHitSets();
    for (auto s = sets.first; s != sets.second; ++s)
    {
      const auto range = s->second->getHits();
      for (auto h = range.first; h != range.second; ++h)
        result[{s->first, h->first}] = h->second->getAdc();
    }
    return result;
  }
}

// Synthetic zigzag polygons make this regression independent of CDB and
// personal BRD paths. The real polygon integrator and MapToPadPlane are used.
struct PHG4TpcPadPlaneReadoutTestAccess
{
  static void geometry(PHG4TpcPadPlaneReadout& p, PHG4TpcGeomContainer& geometry)
  {
    p.GeomContainer = &geometry;
    p.sigmaT = .08;
    p.NTBins = 512;
    p.averageGEMGain = 2000.;
    p.MinRadius = {34.5, 50., 70.};
    p.MaxRadius = {37.5, 60., 80.};
    std::array<std::vector<double>, 2> lo, hi;
    const double step = (M_PI / 6.) / 94.;
    for (int side = 0; side < 2; ++side)
      for (int sector = 0; sector < 12; ++sector)
      {
        hi[side].push_back(11. * M_PI / 12. - sector * M_PI / 6.);
        lo[side].push_back(hi[side].back() - M_PI / 6.);
      }
    for (unsigned layer = 7; layer <= 9; ++layer)
    {
      auto* g = new PHG4TpcGeomv2;
      const double radius = 35. + layer - 7;
      g->set_binning(PHG4CylinderCellDefs::sizebinning);
      g->set_radius(radius); g->set_thickness(1.);
      g->set_phibins(1128); g->set_phistep(step); g->set_phimin(-M_PI);
      g->set_sector_min_phi(lo); g->set_sector_max_phi(hi);
      g->set_zbins(512); g->set_zmin(0.); g->set_zstep(53.3262);
      require(geometry.AddLayerCellGeom(layer, g) == 0, "geometry insertion failed");
      for (int pad = 0; pad < 94; ++pad)
      {
        const double phi = hi[0][2] - (pad + .5) * step;
        PHG4TpcPadPlaneReadout::PadInfo polygon{};
        polygon.pad_bin = pad;
        // Radial edges have a zigzag in the local tangential coordinate;
        // 20% of the tangential pitch remains a physical pad gap.
        for (const auto& xy : std::vector<std::pair<double,double>>{
               {-.48,-.4},{-.24,-.28},{0.,-.4},{.24,-.28},{.48,-.4},
               {.48,.4},{.24,.28},{0.,.4},{-.24,.28},{-.48,.4}})
        {
          const double r = radius + xy.first;
          const double tangential = xy.second * radius * step;
          polygon.vertices.push_back({r*std::cos(phi)-tangential*std::sin(phi),
                                      r*std::sin(phi)+tangential*std::cos(phi)});
        }
        p.Pads[layer].push_back(polygon);
      }
    }
  }

  static void seed_rng(PHG4TpcPadPlaneReadout& p) { gsl_rng_set(p.RandomGenerator, seed); }
  static void mask(PHG4TpcPadPlaneReadout& p, const Key& key)
  {
    p.m_maskDeadChannels = true;
    p.m_deadChannelMap.clear();
    require(p.m_deadChannelMap.add(key.first, TpcDefs::getPad(key.second)), "mask insertion failed");
  }

  // Physics reference: integrate each candidate pad directly, then apply the
  // documented gap normalization, SAMPA fractions and stored-hit quantization.
  // It intentionally does not use the new cache or a copied MapToPadPlane.
  static HitMap expected(PHG4TpcPadPlaneReadout& p, double radius, double phi,
                         unsigned side, bool preserve, double time,
                         std::ostream* dump, const std::string& label)
  {
    struct Contribution { unsigned layer; int pad; double mass; };
    std::vector<Contribution> contributions;
    double total = 0.;
    for (auto layer : p.layersInRadialWindow(radius, p.sigmaT, p._nsigmas))
    {
      p.LayerGeom = p.getGeomForLayer(layer);
      p.sector_min_Phi = p.LayerGeom->get_sector_min_phi();
      p.sector_max_Phi = p.LayerGeom->get_sector_max_phi();
      p.phi_bin_width = p.LayerGeom->get_phistep();
      std::vector<int> pads;
      std::vector<double> masses;
      p.SERF_zigzag_phibins(side, layer, p.check_phi(side, phi, radius), radius, p.sigmaT, pads, masses);
      require(pads.size() == masses.size(), "pad/mass vector mismatch");
      for (std::size_t i = 0; i < pads.size(); ++i)
      {
        total += masses[i]; contributions.push_back({layer,pads[i],masses[i]});
        if (dump) *dump << label << " overlap " << layer << ' ' << pads[i] << ' ' << std::hexfloat << masses[i] << '\n';
      }
    }
    const double normalization = preserve || total <= 1e-16 ? 1. : total;
    seed_rng(p);
    const double gain = p.getSingleEGEMAmplification();
    std::vector<int> bins;
    std::vector<double> shares;
    p.sampaTimeDistribution(time, bins, shares);
    // Peak-normalized voltage averages retain the physical analog pulse area;
    // neither a finite time window nor the eight-clock tail is renormalized.
    std::map<Key,TrkrHitv2> expected_hits;
    for (const auto& c : contributions)
    {
      if (c.pad < 0 || c.pad >= 1128 || c.mass/normalization <= 0.) continue;
      const auto set = TpcDefs::genHitSetKey(c.layer, c.pad/94, side);
      if (p.isChannelMasked(set,c.pad)) continue;
      for (std::size_t i = 0; i < bins.size(); ++i)
      {
        if (bins[i] < 0 || bins[i] >= 512) continue;
        const float energy = gain * (c.mass/normalization) * shares[i];
        if (energy < p.neffelectrons_threshold) continue;
        expected_hits[{set,TpcDefs::genHitKey(c.pad,bins[i])}].addEnergy(energy);
      }
    }
    HitMap result;
    for (const auto& h : expected_hits) result[h.first] = h.second.getAdc();
    return result;
  }

  static HitMap map(PHG4TpcPadPlaneReadout& p, double radius, double phi,
                    unsigned side, bool preserve, double time)
  {
    seed_rng(p);
    p.PreservePadGapLosses(preserve);
    PHG4HitContainer g4hits("G4HIT_TPC");
    g4hits.AddHit(7,new PHG4Hitv1);
    TpcClusterBuilder truth;
    TrkrHitSetContainerv1 single, accumulated;
    p.MapToPadPlane(truth,&single,&accumulated,nullptr,radius*std::cos(phi),radius*std::sin(phi),time,side,g4hits.getHits().first,nullptr,nullptr);
    require(hits(single)==hits(accumulated), "single/global hit maps disagree");
    return hits(accumulated);
  }

  static void empty_polygons(PHG4TpcPadPlaneReadout& p)
  {
    for (auto& layer : p.Pads)
      for (auto& pad : layer)
        for (auto& point : pad.vertices) { point.x += 1000.; point.y += 1000.; }
  }
};

int main(int argc, char** argv)
{
  try
  {
    recoConsts::instance()->set_IntFlag("RANDOMSEED",seed);
    std::ofstream dump;
    if (argc > 1) { dump.open(argv[1]); require(bool(dump),"cannot open dump"); }
    const int repeats = argc > 2 ? std::stoi(argv[2]) : 1;
    require(repeats > 0,"invalid repeat count");
    PHG4TpcGeomContainer geometry;
    PHG4TpcPadPlaneReadout p;
    PHG4TpcPadPlaneReadoutTestAccess::geometry(p,geometry);
    const double step = (M_PI/6.)/94.;
    const double center = M_PI/2.-.5*step;
    const std::vector<std::pair<std::string,std::pair<double,double>>> cases{
      {"pad_center",{36.,center}}, {"pad_edge",{36.,center+.4*step}},
      {"pad_gap",{36.,M_PI/2.}}, {"sector_edge",{36.,7.*M_PI/12.-1e-6}},
      {"radial_boundary",{36.5+1e-6,center}}, {"phi_wrap",{36.,-M_PI+1e-6}}};
    double runtime = 0.; unsigned count = 0;
    for (unsigned side : {0U,1U})
      for (bool preserve : {false,true})
        for (const auto& c : cases)
        {
          const std::string label = c.first+"_side"+std::to_string(side)+"_gap"+std::to_string(preserve);
          const auto [radius,phi] = c.second;
          const auto expected = PHG4TpcPadPlaneReadoutTestAccess::expected(p,radius,phi,side,preserve,500.,dump.is_open()?&dump:nullptr,label);
          require(!expected.empty(),label+": empty fixture");
          const auto actual = PHG4TpcPadPlaneReadoutTestAccess::map(p,radius,phi,side,preserve,500.);
          require(actual==expected,label+": stored pad/time/charge mismatch");
          if (dump.is_open())
            for (const auto& h : actual) dump << label << " hit " << h.first.first << ' ' << h.first.second << ' ' << h.second << '\n';
          const auto begin = std::chrono::steady_clock::now();
          for (int repeat=0; repeat<repeats; ++repeat)
            require(PHG4TpcPadPlaneReadoutTestAccess::map(p,radius,phi,side,preserve,500.)==expected,label+": nondeterministic repeated cloud");
          runtime += std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
          ++count;
        }
    const auto live = PHG4TpcPadPlaneReadoutTestAccess::map(p,36.,center,0,true,500.);
    auto target = live.begin(); while(target!=live.end() && target->second==0) ++target;
    require(target!=live.end(),"fixture has no positive ADC");
    PHG4TpcPadPlaneReadoutTestAccess::mask(p,target->first);
    for (bool preserve : {false,true})
    {
      const auto expected = PHG4TpcPadPlaneReadoutTestAccess::expected(p,36.,center,0,preserve,500.,nullptr,"");
      const auto actual = PHG4TpcPadPlaneReadoutTestAccess::map(p,36.,center,0,preserve,500.);
      require(actual==expected,"masked map differs from physical reference");
      for(const auto& h:actual)
        require(h.first.first!=target->first.first || TpcDefs::getPad(h.first.second)!=TpcDefs::getPad(target->first.second),"masked channel created a hit");
    }
    PHG4TpcPadPlaneReadoutTestAccess::empty_polygons(p);
    for(bool preserve : {false,true})
      require(PHG4TpcPadPlaneReadoutTestAccess::map(p,36.,center,0,preserve,500.).empty(),"zero recovered mass produced hits");
    std::cout << "PASS SERF: " << count << " boundary/side/normalization fixtures; mask and zero-mass checks; "
              << count*repeats << " timed clouds in " << runtime << " seconds\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "test_PHG4TpcPadPlaneReadoutSerf: " << error.what() << '\n';
    return 1;
  }
}
