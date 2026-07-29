#include "PHG4TpcPadPlaneReadout.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CellDefs.h>  // for genkey, keytype
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4HitContainer.h>

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
// Move to new storage containers
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, hitse...
#include <trackbase/TrkrHit.h>   // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <g4tracking/TrkrTruthTrackContainerv1.h>
#include <g4tracking/TrkrTruthTrackv1.h>

#include <phool/phool.h>  // for PHWHERE

#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMarker.h>
#include <TPad.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_alloc

#include <boost/format.hpp>

#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>  // for getenv
#include <format>
#include <iostream>
#include <map>      // for _Rb_tree_cons...
#include <utility>  // for pair
#include <fstream>  // for std::ifstream
#include <chrono>
#include <iomanip>

class PHCompositeNode;
class TrkrHitTruthAssoc;

namespace
{
  //! convenient square function
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }

  template <class T>
  inline T get_r(const T &x, const T &y)
  {
    return std::sqrt(square(x) + square(y));
  }

  //! return normalized gaussian centered on zero and of width sigma
  template <class T>
  inline T gaus(const T &x, const T &sigma)
  {
    return std::exp(-square(x / sigma) / 2) / (sigma * std::sqrt(2 * M_PI));
  }

  static constexpr int print_layer = 18;

}  // namespace

PHG4TpcPadPlaneReadout::PHG4TpcPadPlaneReadout(const std::string &name)
  : PHG4TpcPadPlane(name)
{
  InitializeParameters();
  // if(m_flagToUseGain==1)
  ReadGain();
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, PHRandomSeed());  // fixed seed is handled in this funtcion


  return;
}

PHG4TpcPadPlaneReadout::~PHG4TpcPadPlaneReadout()
{
  gsl_rng_free(RandomGenerator);
  for (auto his : h_gain)
  {
    delete his;
  }

  for(int i=0; i<2; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<12; ++k)
        delete flangau[i][j][k];
}

// pick all layers whose radial annulus intersects [rad - nsig*sigma, rad + nsig*sigma]
std::vector<unsigned int>
PHG4TpcPadPlaneReadout::layersInRadialWindow(double rad, double sigma, double nsig) const
{
  std::vector<unsigned int> out;
  const double rmin = rad - nsig*sigma;
  const double rmax = rad + nsig*sigma;

  PHG4TpcGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
  for (auto it = layerrange.first; it != layerrange.second; ++it)
  {
    const auto* g = it->second;
    const double rLow  = g->get_radius() - 0.5*g->get_thickness();
    const double rHigh = g->get_radius() + 0.5*g->get_thickness();

    const double lo = std::max(rLow,  rmin);
    const double hi = std::min(rHigh, rmax);
    if (hi > lo) out.push_back(g->get_layer());
  }
  return out;
}

// get geometry for a given layer (utility; linear scan is fine here)
PHG4TpcGeom*
PHG4TpcPadPlaneReadout::getGeomForLayer(unsigned int layer) const
{
  PHG4TpcGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
  for (auto it = layerrange.first; it != layerrange.second; ++it)
    if (static_cast<unsigned int>(it->second->get_layer()) == layer) return it->second;
  return nullptr;
}

double PHG4TpcPadPlaneReadout::GetAveragePadAreaForLayer(unsigned int layer) const
{
  if (layer >= m_average_pad_area_cm2.size() || m_pad_polygon_count_by_layer[layer] == 0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return m_average_pad_area_cm2[layer];
}

double PHG4TpcPadPlaneReadout::polygonArea(const std::vector<Point>& vertices) const
{
  if (vertices.size() < 3)
  {
    return 0.0;
  }

  double twice_area = 0.0;
  for (std::size_t i = 0; i < vertices.size(); ++i)
  {
    const auto& current = vertices[i];
    const auto& next = vertices[(i + 1) % vertices.size()];
    twice_area += current.x * next.y - next.x * current.y;
  }

  return 0.5 * std::abs(twice_area);
}

void PHG4TpcPadPlaneReadout::updatePadPolygonAreaSummary()
{
  m_average_pad_area_cm2.fill(0.0);
  m_pad_polygon_count_by_layer.fill(0);

  for (std::size_t layer = 0; layer < Pads.size(); ++layer)
  {
    double total_area_cm2 = 0.0;
    unsigned int polygon_count = 0;

    for (const auto& pad : Pads[layer])
    {
      if (pad.vertices.size() < 3)
      {
        continue;
      }

      total_area_cm2 += polygonArea(pad.vertices);
      ++polygon_count;
    }

    if (polygon_count == 0)
    {
      continue;
    }

    m_pad_polygon_count_by_layer[layer] = polygon_count;
    m_average_pad_area_cm2[layer] = total_area_cm2 / static_cast<double>(polygon_count);
  }
}

void PHG4TpcPadPlaneReadout::printPadPolygonAreaSummary() const
{
  std::size_t layers_with_polygons = 0;

  for (std::size_t layer = 0; layer < m_pad_polygon_count_by_layer.size(); ++layer)
  {
    const unsigned int polygon_count = m_pad_polygon_count_by_layer[layer];
    if (polygon_count == 0)
    {
      continue;
    }

    ++layers_with_polygons;

    const int module = (layer >= 7) ? static_cast<int>((layer - 7) / 16) : -1;
    const int local_layer = (layer >= 7) ? static_cast<int>((layer - 7) % 16) : -1;
    std::cout << Name() << ": pad polygon area layer=" << layer
              << " module=" << module
              << " local_layer=" << local_layer
              << " pad_count=" << polygon_count
              << " average_area_cm2=" << m_average_pad_area_cm2[layer]
              << std::endl;
  }

  std::cout << Name() << ": pad polygon area summary layers_with_polygons="
            << layers_with_polygons << std::endl;
}

void PHG4TpcPadPlaneReadout::plotPadPolygonAreaSummary() const
{
  if (m_pad_area_plot_output.empty())
  {
    std::cout << Name() << ": pad polygon area plot skipped (no output file configured)." << std::endl;
    return;
  }

  std::vector<double> layer_numbers;
  std::vector<double> average_areas;
  layer_numbers.reserve(m_pad_polygon_count_by_layer.size());
  average_areas.reserve(m_pad_polygon_count_by_layer.size());

  for (std::size_t layer = 0; layer < m_pad_polygon_count_by_layer.size(); ++layer)
  {
    if (m_pad_polygon_count_by_layer[layer] == 0)
    {
      continue;
    }

    layer_numbers.push_back(static_cast<double>(layer));
    average_areas.push_back(m_average_pad_area_cm2[layer]);
  }

  if (layer_numbers.empty())
  {
    std::cout << Name() << ": pad polygon area plot skipped (no polygon data)." << std::endl;
    return;
  }

  const bool was_batch = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  const std::string canvas_name = std::string("c_pad_area_by_layer_") + Name();
  TCanvas canvas(canvas_name.c_str(), "Average pad area by layer", 1000, 600);
  canvas.cd();
  canvas.SetGridx();
  canvas.SetGridy();

  TGraph graph(static_cast<int>(layer_numbers.size()), layer_numbers.data(), average_areas.data());
  graph.SetTitle("Average TPC Pad Polygon Area by Layer;Layer Number;Average Area [cm^{2}]");
  graph.SetMarkerStyle(20);
  graph.SetMarkerSize(1.0);
  graph.SetMarkerColor(kBlue + 2);
  graph.SetLineColor(kBlue + 2);
  graph.SetLineWidth(2);
  graph.SetMinimum(0.0);
  graph.Draw("ALP");

  canvas.Update();
  canvas.SaveAs(m_pad_area_plot_output.c_str());
  gROOT->SetBatch(was_batch);

  std::cout << Name() << ": saved pad polygon area plot to "
            << m_pad_area_plot_output << std::endl;
}

//_________________________________________________________
int PHG4TpcPadPlaneReadout::InitRun(PHCompositeNode *topNode)
{
  // base class
  const auto reply = PHG4TpcPadPlane::InitRun(topNode);
  if (reply != Fun4AllReturnCodes::EVENT_OK)
  {
    return reply;
  }
  const std::string seggeonodename = "TPCGEOMCONTAINER";
  GeomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, seggeonodename);
  assert(GeomContainer);

  // The geometry node is the authoritative source for TPC timing. This keeps
  // electron drift, pad time binning, and reconstruction on one configuration.
  auto *z_geometry = GeomContainer->GetLayerCellGeom(20);
  if (!z_geometry)
  {
    std::cerr << "PHG4TpcPadPlaneReadout: missing TPC layer geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const double tpc_adc_clock = z_geometry->get_adc_clock();
  const double max_drift_length = z_geometry->get_max_driftlength();
  const double drift_velocity_sim = z_geometry->get_drift_velocity_sim();
  const double extended_readout_time = z_geometry->get_extended_readout_time();

  if (tpc_adc_clock <= 0.0 || max_drift_length <= 0.0 || drift_velocity_sim <= 0.0)
  {
    std::cerr << "PHG4TpcPadPlaneReadout: invalid TPC timing geometry"
              << " (adc clock " << tpc_adc_clock
              << ", max drift length " << max_drift_length
              << ", drift velocity " << drift_velocity_sim << ")" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const double max_time =
      extended_readout_time + 2.0 * max_drift_length / drift_velocity_sim;
  NTBins = static_cast<int>(max_time / tpc_adc_clock) + 1;

  if(m_use_module_gain_weights)
    {
      int side, region, sector;
      double weight;
      std::ifstream weights_file(m_tpc_module_gain_weights_file);
      if(!weights_file.is_open())
	{
	  std::cout << ".In PHG4TpcPadPlaneReadout: Option to use module gain weights enabled, but weights file not found. Aborting." << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}

      for(int iside =0; iside < 2; ++iside)
	{
	  for(int isec = 0; isec < 12; ++isec)
	    {
	      for(int ir = 0; ir < 3; ++ir)
		{
		  weights_file >> side >> region >> sector >> weight;
		  m_module_gain_weight[side][region][sector] = weight;
		  std::cout << " iside " << iside << " side " << side << " ir " << ir
			    << " region " << region << " isec " << isec
			    << " sector " << sector << " weight " << weight << std::endl;
		}
	    }
	}
    }

  if(m_useLangau)
    {
      int side, region, sector;
      double par0; double par1; double par2; double par3;
      std::ifstream pars_file(m_tpc_langau_pars_file);
      if(!pars_file.is_open())
	{
	  std::cout << ".In PHG4TpcPadPlaneReadout: Option to use Langau parameters enabled, but parameter file not found. Aborting." << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}

      for(int iside =0; iside < 2; ++iside)
	{
	  for(int isec = 0; isec < 12; ++isec)
	    {
	      for(int ir = 0; ir < 3; ++ir)
		{
		  pars_file >> side >> region >> sector >> par0 >> par1 >> par2 >> par3;
		  flangau[side][region][sector] = new TF1((boost::format("flangau_%d_%d_%d") % side % region % sector).str().c_str(), [](double *x, double *par)
		  {
		    Double_t invsq2pi = 0.3989422804014;
		    Double_t mpshift  = -0.22278298;
		    Double_t np = 100.0;
		    Double_t sc =   5.0;
		    Double_t xx;
		    Double_t mpc;
		    Double_t fland;
		    Double_t sum = 0.0;
		    Double_t xlow,xupp;
		    Double_t step;
		    Double_t i;
		    mpc = par[1] - mpshift * par[0];
		    xlow = x[0] - sc * par[3];
		    xupp = x[0] + sc * par[3];
		    step = (xupp-xlow) / np;
		    for(i=1.0; i<=np/2; i++)
		    {
		      xx = xlow + (i-.5) * step;
		      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		      sum += fland * TMath::Gaus(x[0],xx,par[3]);

		      xx = xupp - (i-.5) * step;
		      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		      sum += fland * TMath::Gaus(x[0],xx,par[3]);
		    }

		    return (par[2] * step * sum * invsq2pi / par[3]);
		  }, 0, 5000, 4);


		  flangau[side][region][sector]->SetParameters(par0,par1,par2,par3);
		  //std::cout << " iside " << iside << " side " << side << " ir " << ir
		  //	    << " region " << region << " isec " << isec
		  //	    << " sector " << sector << " weight " << weight << std::endl;
		}
	    }
	}
    }
    if (m_maskDeadChannels)
  {
    makeChannelMask(m_deadChannelMap, m_deadChannelMapName, "TotalDeadChannels");
  }
  if (m_maskHotChannels)
  {
    makeChannelMask(m_hotChannelMap, m_hotChannelMapName, "TotalHotChannels");
  }

	  loadPadPlanes();
	  if (m_enable_pad_polygon_area_calculation)
	  {
	    updatePadPolygonAreaSummary();
	  }
	  else
	  {
	    m_average_pad_area_cm2.fill(0.0);
	    m_pad_polygon_count_by_layer.fill(0);
	  }
	  if (m_check_pad_geom_consistency)
	  {
	    runPadGeomConsistencyCheck();
	  }

	  // Summarize SERF polygon availability and optionally enforce requirement
  {
    std::size_t poly_layers = 0;
    std::vector<unsigned int> missing_layers;

    // We consider readout layers to start at 7 (Pads is sized as 7 + 16*3)
    constexpr unsigned int kFirstPolygonLayer = 7;

    // Check if any layer has polygons at all
    for (const auto &layerPads : Pads)
    {
      for (const auto &pad : layerPads)
      {
        if (pad.vertices.size() >= 3) { poly_layers++; break; }
      }
    }
    m_serf_polygons_present = (poly_layers > 0);

    // If SERF is desired, verify full coverage of all readout layers present in geometry
    if (m_use_serf_padsharing)
    {
      if (!GeomContainer)
      {
        std::cerr << "PHG4TpcPadPlaneReadout: Geometry container missing in InitRun." << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      PHG4TpcGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
      for (auto it = layerrange.first; it != layerrange.second; ++it)
      {
        unsigned int layer = static_cast<unsigned int>(it->second->get_layer());
        if (layer < kFirstPolygonLayer) continue; // allow early layers to miss polygons
        bool has_poly = false;
        if (layer < Pads.size())
        {
          const auto &layerPads = Pads[layer];
          for (const auto &pad : layerPads)
          {
            if (pad.vertices.size() >= 3) { has_poly = true; break; }
          }
        }
        if (!has_poly) missing_layers.push_back(layer);
      }

      if (!m_serf_polygons_present)
      {
        std::cout << "PHG4TpcPadPlaneReadout: SERF requested but pad polygons not available; using analytic sharing."
                  << " Call RequireSerfPadSharing(true) to abort instead." << std::endl;
        if (m_require_serf)
        {
          std::cerr << "PHG4TpcPadPlaneReadout: SERF pad polygons required but missing. Aborting run." << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
      }

      if (!missing_layers.empty())
      {
        std::cout << "PHG4TpcPadPlaneReadout: SERF polygons missing for "
                  << missing_layers.size() << " readout layer(s) starting at layer 7." << std::endl;
        if (m_require_serf_full)
        {
          std::cerr << "PHG4TpcPadPlaneReadout: RequireSerfFullCoverage enabled; aborting due to missing polygons in readout layers." << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
      }
	  }
  }

	  if (m_enable_pad_polygon_area_calculation)
	  {
	    printPadPolygonAreaSummary();
	    plotPadPolygonAreaSummary();
	  }

	  /* for (int j=0; j<7+16*3;j++){
  for(size_t i = 0; i < Pads[j].size(); i++)
  {
    std::cout<<"Module "<<(j-7)/16<<" layer = "<<j<<" pad_number "<<i<<" pad name "<<Pads[j][i].name<<" pad_bin "<<Pads[j][i].pad_bin<<" ( "<<ntpc_phibins_sector[(j-7)/16] - Pads[j][i].pad_bin -1 <<" ) "<<" cx "<<Pads[j][i].cx<<" cy "<<Pads[j][i].cy  <<" rad "<<Pads[j][i].rad<<" phi "<<Pads[j][i].phi<<" Number of verticies "<<Pads[j][i].vertices.size()<<std::endl;

  }
   }

*/
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcPadPlaneReadout::BeginEvent(unsigned int event)
{
  if (!m_enable_side_layer_debug)
  {
    return;
  }

  m_side_layer_debug_event = static_cast<int>(event);
  m_side_layer_no_layer = {{0, 0}};
  m_side_layer_debug_counters.clear();
}

void PHG4TpcPadPlaneReadout::EndEvent(unsigned int event)
{
  if (!m_enable_side_layer_debug)
  {
    return;
  }

  if (m_side_layer_debug_event < 0)
  {
    m_side_layer_debug_event = static_cast<int>(event);
  }

  std::uint64_t total_clouds = 0;
  std::uint64_t total_bins = 0;
  for (const auto &[key, cnt] : m_side_layer_debug_counters)
  {
    (void) key;
    total_clouds += cnt.cloud_calls;
    total_bins += cnt.bins_written;
  }
  if (total_clouds == 0 && m_side_layer_no_layer[0] == 0 && m_side_layer_no_layer[1] == 0)
  {
    return;
  }

  std::cout << Name() << ": PadPlane side/layer summary evt=" << m_side_layer_debug_event
            << " clouds=" << total_clouds
            << " bins_written=" << total_bins
            << " no_layer(side0/side1)=" << m_side_layer_no_layer[0] << "/" << m_side_layer_no_layer[1]
            << std::endl;

  for (const auto &[key, cnt] : m_side_layer_debug_counters)
  {
    const unsigned int side = key.first;
    const unsigned int layer = key.second;
    if (cnt.cloud_calls == 0)
    {
      continue;
    }
    std::cout << "  side=" << side
              << " layer=" << layer
              << " clouds=" << cnt.cloud_calls
              << " cloud_no_hits=" << cnt.cloud_no_hits
              << " pad_contributors=" << cnt.pad_contributors
              << " bins_written=" << cnt.bins_written
              << " neff_sum=" << cnt.neff_sum
              << std::endl;
  }
}


const std::vector<std::string>
  PHG4TpcPadPlaneReadout::brdMaps_ = {
    "/sphenix/user/mitrankova/Simulation/PadPlane/AutoPad-R1-RevA.brd",
    "/sphenix/user/mitrankova/Simulation/PadPlane/AutoPad-R2-RevA-Pads.brd",
    "/sphenix/user/mitrankova/Simulation/PadPlane/AutoPad-R3-RevA.brd"
};

//_________________________________________________________

void PHG4TpcPadPlaneReadout::loadPadPlanes() {
 for (size_t i = 0; i < brdMaps_.size(); ++i) {

    std::ifstream in(brdMaps_[i].c_str());
    if (!in) {
        std::cerr << "Cannot open " << brdMaps_[i] << "\n";
        return;
    }

    //std::cout<<"!!!!!getPadCoordinates filename "<<brdMaps_[i]<<std::endl;
    std::string line;
    bool        inSignal  = false;
    bool        inPolygon = false;
    bool        keepSignal = false;

    std::string pname_tmp;
    int iter = 0;
    double sumX = 0, sumY = 0;
    PadInfo p;
    int layer = -1, pad = -1;
    p.clear();
    p.isedge = false;
    for(int l=0;l<7;l++)
        Pads[l].push_back(p);

    while (std::getline(in, line)) {

        size_t pos = line.find_first_not_of(" \t");
        if (pos != std::string::npos) line = line.substr(pos);


        if (!inSignal && line.find("<signal ") == 0) {
            inSignal = true;

            size_t n1 = line.find("name=\"");
            if (n1 != std::string::npos) {
                n1 += 6;
                size_t n2 = line.find('"', n1);
                pname_tmp = line.substr(n1, n2 - n1);
                //std::cout<<"!!!!!pname_tmp "<<pname_tmp<<std::endl;

            }
            keepSignal = (pname_tmp.rfind("ZZ", 0) == 0);
            continue;
        }

        if ( inSignal && keepSignal && !inPolygon && line.find("<polygon") == 0) {
            p.clear();
            p.name = pname_tmp;
            layer = -1, pad = -1;
            if (std::sscanf(pname_tmp.c_str(), "ZZ.%2d.%3d", &layer, &pad) == 2) {
                p.pad_bin = pad;
                if (p.pad_bin==0 || p.pad_bin==ntpc_phibins_sector[i]-1)
                {
                    p.isedge = true; // edge pads
                }

            }
            inPolygon = true;
            sumX = 0, sumY = 0;
            continue;
        }

        if (inPolygon && line.find("<vertex") == 0) {
            size_t x1 = line.find("x=\"");
            size_t y1 = line.find("y=\"");
            if (x1!=std::string::npos && y1!=std::string::npos) {
                x1 += 3; size_t x2 = line.find('"', x1);
                y1 += 3; size_t y2 = line.find('"', y1);
                double x = std::atof(line.substr(x1, x2-x1).c_str())/10.0;
                double y = std::atof(line.substr(y1, y2-y1).c_str())/10.0;
                p.vertices.push_back(Point{x,y});
                sumX += x;
                sumY += y;
            }
            continue;
        }

        if (inPolygon && line.find("</polygon>") == 0) {
            inPolygon = false;
            p.cx = sumX / p.vertices.size();
            p.cy = sumY / p.vertices.size();
            p.rad = get_r(p.cx, p.cy);

            p.phi = std::atan2(p.cy, p.cx);

            p.pad_number=iter;

            iter++;
            continue;
        }


        if (inSignal && line.find("</signal>") == 0) {
            inSignal = false;
            keepSignal = false;
          //  std::cout<<"Module "<<i<<" pad_number "<<p.pad_number <<" pad name "<<p.name<<" pad layer "<<7 + i * 16 + layer<<" pad_bin "<<p.pad_bin<<" ( "<<ntpc_phibins_sector[i] - p.pad_bin -1 <<" ) "<<" cx "<<p.cx<<" cy "<<p.cy  <<" rad "<<p.rad<<" phi "<<p.phi<<" Number of verticies "<<p.vertices.size()<<std::endl;

            Pads[7 + i * 16 + layer].push_back(p);

            continue;
        }

    }
}
}

void PHG4TpcPadPlaneReadout::runPadGeomConsistencyCheck()
{
  if (!GeomContainer)
  {
    std::cout << Name() << ": pad-geom consistency check skipped (no GeomContainer)." << std::endl;
    return;
  }

  constexpr unsigned int kFirstReadoutLayer = 7;
  const bool detailed = (Verbosity() > 0);
  const bool dump_requested = (m_consistency_dump_layer >= 0 && m_consistency_dump_sector >= 0);
  const bool dump_sector_valid = (m_consistency_dump_sector >= 0 &&
                                  m_consistency_dump_sector < NSectors);
  const bool dump_enabled = dump_requested && dump_sector_valid;

  auto wrap = [](double dphi) {
    while (dphi > M_PI) dphi -= 2.0 * M_PI;
    while (dphi < -M_PI) dphi += 2.0 * M_PI;
    return dphi;
  };

  std::size_t layers_checked = 0;
  std::size_t layers_missing_polygons = 0;
  std::size_t layers_with_boundary_data = 0;
  std::size_t side_reports = 0;
  double global_max_abs_dphi = 0.0;
  double global_max_abs_dr = 0.0;
  double global_max_abs_dphi_low = 0.0;
  double global_max_abs_dphi_high = 0.0;
  double global_max_abs_dr_low = 0.0;
  double global_max_abs_dr_high = 0.0;

  std::cout << Name() << ": running LayerGeom/BRD pad consistency check..." << std::endl;
  if (dump_requested && !dump_sector_valid)
  {
    std::cout << Name()
              << ": per-pad consistency dump disabled due to invalid sector="
              << m_consistency_dump_sector
              << " (expected 0-" << (NSectors - 1) << ")." << std::endl;
  }
  else if (dump_enabled)
  {
    std::cout << Name() << ": per-pad consistency dump target layer="
              << m_consistency_dump_layer
              << " sector=" << m_consistency_dump_sector
              << " side=" << m_consistency_dump_side
              << " (side<0 means both sides)." << std::endl;
  }

  PHG4TpcGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
  for (auto layeriter = layerrange.first; layeriter != layerrange.second; ++layeriter)
  {
    auto* thisGeom = layeriter->second;
    if (!thisGeom) continue;

    const unsigned int layer = static_cast<unsigned int>(thisGeom->get_layer());
    if (layer < kFirstReadoutLayer) continue;
    ++layers_checked;

    const int tpc_module = static_cast<int>((layer - kFirstReadoutLayer) / 16);
    if (tpc_module < 0 || tpc_module > 2) continue;

    if (layer >= Pads.size() || Pads[layer].empty())
    {
      ++layers_missing_polygons;
      if (detailed)
      {
        std::cout << Name() << ": consistency layer=" << layer
                  << " skipped (no BRD polygons loaded)." << std::endl;
      }
      continue;
    }

    const int phibins = thisGeom->get_phibins();
    const int nlocal = ntpc_phibins_sector[tpc_module];
    if (phibins <= 0 || nlocal <= 0) continue;

    // Layer radial boundary check:
    // Compare GeomContainer annulus [r-thickness/2, r+thickness/2]
    // with BRD polygon radial envelope [min(vertex_r), max(vertex_r)].
    const double geom_r_low = thisGeom->get_radius() - 0.5 * thisGeom->get_thickness();
    const double geom_r_high = thisGeom->get_radius() + 0.5 * thisGeom->get_thickness();
    double brd_r_low = std::numeric_limits<double>::infinity();
    double brd_r_high = -std::numeric_limits<double>::infinity();
    int brd_vertex_count = 0;
    for (const auto& pad : Pads[layer])
    {
      for (const auto& v : pad.vertices)
      {
        const double rv = std::hypot(v.x, v.y);
        if (rv < brd_r_low) brd_r_low = rv;
        if (rv > brd_r_high) brd_r_high = rv;
        ++brd_vertex_count;
      }
    }

    if (brd_vertex_count > 0)
    {
      ++layers_with_boundary_data;
      const double dr_low = brd_r_low - geom_r_low;
      const double dr_high = brd_r_high - geom_r_high;
      const double abs_dr_low = std::abs(dr_low);
      const double abs_dr_high = std::abs(dr_high);
      if (abs_dr_low > global_max_abs_dr_low) global_max_abs_dr_low = abs_dr_low;
      if (abs_dr_high > global_max_abs_dr_high) global_max_abs_dr_high = abs_dr_high;

      if (detailed)
      {
        std::cout << Name() << ": boundary layer=" << layer
                  << " geom[r_low,r_high]=[" << geom_r_low << "," << geom_r_high << "]"
                  << " brd[r_low,r_high]=[" << brd_r_low << "," << brd_r_high << "]"
                  << " dR_low(brd-geom)=" << dr_low
                  << " dR_high(brd-geom)=" << dr_high
                  << std::endl;
      }
    }
    else if (detailed)
    {
      std::cout << Name() << ": boundary layer=" << layer
                << " skipped (no BRD vertices)." << std::endl;
    }

    for (unsigned int side = 0; side < static_cast<unsigned int>(NSides); ++side)
    {
      LayerGeom = thisGeom;
      sector_min_Phi = LayerGeom->get_sector_min_phi();
      sector_max_Phi = LayerGeom->get_sector_max_phi();

      if (detailed)
      {
        const std::streamsize old_precision = std::cout.precision();
        std::cout << std::setprecision(17);
        std::cout << Name() << ": sector bounds layer=" << layer
                  << " side=" << side << std::endl;
        for (int s = 0; s < NSectors; ++s)
        {
          if (s >= static_cast<int>(sector_min_Phi[side].size()) ||
              s >= static_cast<int>(sector_max_Phi[side].size()))
          {
            continue;
          }

          const double phi_min = sector_min_Phi[side][s];
          const double phi_max = sector_max_Phi[side][s];
          std::cout << Name()
                    << ": sector bounds layer=" << layer
                    << " side=" << side
                    << " sector=" << s
                    << " phi_min=" << phi_min
                    << " phi_max=" << phi_max
                    << " phi_span=" << (phi_max - phi_min)
                    << std::endl;
        }
        std::cout << std::setprecision(old_precision);
      }

      int matched = 0;
      int missing = 0;
      double sum_abs_dphi = 0.0;
      double max_abs_dphi = 0.0;
      double sum_abs_dphi_low = 0.0;
      double max_abs_dphi_low = 0.0;
      double sum_abs_dphi_high = 0.0;
      double max_abs_dphi_high = 0.0;
      double sum_abs_dr = 0.0;
      double max_abs_dr = 0.0;
      int worst_global_pad = -1;
      int worst_local_pad = -1;
      int worst_global_pad_low = -1;
      int worst_local_pad_low = -1;
      int worst_global_pad_high = -1;
      int worst_local_pad_high = -1;
      const double half_phi = 0.5 * std::abs(thisGeom->get_phistep());
      const double phistep_abs = std::abs(thisGeom->get_phistep());
      const int pads_per_sector = (NSectors > 0) ? (phibins / NSectors) : 0;
      std::array<double, NSectors> sector_max_abs_dphi_low{};
      std::array<double, NSectors> sector_max_abs_dphi_high{};
      std::array<int, NSectors> sector_matched{};
      const bool dump_this_side = dump_enabled &&
                                  static_cast<int>(layer) == m_consistency_dump_layer &&
                                  static_cast<int>(side) == m_consistency_dump_side;
      const bool dump_all_sides = dump_enabled &&
                                  static_cast<int>(layer) == m_consistency_dump_layer &&
                                  m_consistency_dump_side < 0;

      for (int pad_now = 0; pad_now < phibins; ++pad_now)
      {
        int look_pad = ((pad_now % nlocal) + nlocal) % nlocal;
        if (side == 1)
        {
          look_pad = nlocal - look_pad - 1;
        }

        if (look_pad < 0 || static_cast<std::size_t>(look_pad) >= Pads[layer].size())
        {
          ++missing;
          continue;
        }

        const auto& pad = Pads[layer][look_pad];
        if (pad.vertices.size() < 3)
        {
          ++missing;
          continue;
        }

        const double r_geom = thisGeom->get_radius();
        const double phi_geom = thisGeom->get_phicenter(pad_now, side);
        const double x = r_geom * std::cos(phi_geom);
        const double y = r_geom * std::sin(phi_geom);

        int sector_found = -1;
        double x_new = 0.0;
        double y_new = 0.0;
        rotatePointToSector(x, y, side, sector_found, x_new, y_new);
        (void) sector_found;

        const double phi_local = std::atan2(y_new, x_new);
        const double r_local = std::hypot(x_new, y_new);
        const double dphi = wrap(phi_local - pad.phi);
        const double dr = r_local - pad.rad;

        // Compare phi boundaries:
        // BRD boundaries are min/max vertex phi around BRD pad center.
        double brd_min_delta = std::numeric_limits<double>::infinity();
        double brd_max_delta = -std::numeric_limits<double>::infinity();
        for (const auto& v : pad.vertices)
        {
          const double vphi = std::atan2(v.y, v.x);
          const double d = wrap(vphi - pad.phi);
          if (d < brd_min_delta) brd_min_delta = d;
          if (d > brd_max_delta) brd_max_delta = d;
        }
        const double brd_phi_low = pad.phi + brd_min_delta;
        const double brd_phi_high = pad.phi + brd_max_delta;
        const double cdb_phi_low = phi_local - half_phi;
        const double cdb_phi_high = phi_local + half_phi;
        const double dphi_low = wrap(cdb_phi_low - brd_phi_low);
        const double dphi_high = wrap(cdb_phi_high - brd_phi_high);
        const double dphi_low_over_step = (phistep_abs > 0.0) ? (dphi_low / phistep_abs) : -1.0;
        const double dphi_high_over_step = (phistep_abs > 0.0) ? (dphi_high / phistep_abs) : -1.0;

        if (pads_per_sector > 0 && (dump_this_side || dump_all_sides))
        {
          const int sector_idx = pad_now / pads_per_sector;
          if (sector_idx == m_consistency_dump_sector)
          {
            const int local_pad = pad_now - sector_idx * pads_per_sector;
            std::cout << Name() << ": per-pad consistency layer=" << layer
                      << " side=" << side
                      << " sector=" << sector_idx
                      << " pad_global=" << pad_now
                      << " pad_local=" << local_pad
                      << " dphi_low/|phistep|=" << dphi_low_over_step
                      << " dphi_high/|phistep|=" << dphi_high_over_step
                      << std::endl;
          }
        }

        const double abs_dphi = std::abs(dphi);
        const double abs_dphi_low = std::abs(dphi_low);
        const double abs_dphi_high = std::abs(dphi_high);
        const double abs_dr = std::abs(dr);
        sum_abs_dphi += abs_dphi;
        sum_abs_dphi_low += abs_dphi_low;
        sum_abs_dphi_high += abs_dphi_high;
        sum_abs_dr += abs_dr;

        if (abs_dphi > max_abs_dphi)
        {
          max_abs_dphi = abs_dphi;
          worst_global_pad = pad_now;
          worst_local_pad = look_pad;
        }
        if (abs_dphi_low > max_abs_dphi_low)
        {
          max_abs_dphi_low = abs_dphi_low;
          worst_global_pad_low = pad_now;
          worst_local_pad_low = look_pad;
        }
        if (abs_dphi_high > max_abs_dphi_high)
        {
          max_abs_dphi_high = abs_dphi_high;
          worst_global_pad_high = pad_now;
          worst_local_pad_high = look_pad;
        }
        if (abs_dr > max_abs_dr)
        {
          max_abs_dr = abs_dr;
        }

        if (pads_per_sector > 0)
        {
          int sector_idx = pad_now / pads_per_sector;
          if (sector_idx < 0) sector_idx = 0;
          if (sector_idx >= NSectors) sector_idx = NSectors - 1;
          if (abs_dphi_low > sector_max_abs_dphi_low[sector_idx])
          {
            sector_max_abs_dphi_low[sector_idx] = abs_dphi_low;
          }
          if (abs_dphi_high > sector_max_abs_dphi_high[sector_idx])
          {
            sector_max_abs_dphi_high[sector_idx] = abs_dphi_high;
          }
          sector_matched[sector_idx]++;
        }

        ++matched;
      }

      if (matched <= 0)
      {
        if (detailed)
        {
          std::cout << Name() << ": consistency layer=" << layer
                    << " side=" << side
                    << " matched=0 missing=" << missing << std::endl;
        }
        continue;
      }

      const double mean_abs_dphi = sum_abs_dphi / static_cast<double>(matched);
      const double mean_abs_dphi_low = sum_abs_dphi_low / static_cast<double>(matched);
      const double mean_abs_dphi_high = sum_abs_dphi_high / static_cast<double>(matched);
      const double mean_abs_dr = sum_abs_dr / static_cast<double>(matched);
      const double phistep = phistep_abs;
      const double mean_over_phistep_low = (phistep > 0.0) ? (mean_abs_dphi_low / phistep) : -1.0;
      const double mean_over_phistep_high = (phistep > 0.0) ? (mean_abs_dphi_high / phistep) : -1.0;
      const double max_over_phistep = (phistep > 0.0) ? (max_abs_dphi / phistep) : -1.0;
      const double max_over_phistep_low = (phistep > 0.0) ? (max_abs_dphi_low / phistep) : -1.0;
      const double max_over_phistep_high = (phistep > 0.0) ? (max_abs_dphi_high / phistep) : -1.0;

      if (max_abs_dphi > global_max_abs_dphi) global_max_abs_dphi = max_abs_dphi;
      if (max_abs_dphi_low > global_max_abs_dphi_low) global_max_abs_dphi_low = max_abs_dphi_low;
      if (max_abs_dphi_high > global_max_abs_dphi_high) global_max_abs_dphi_high = max_abs_dphi_high;
      if (max_abs_dr > global_max_abs_dr) global_max_abs_dr = max_abs_dr;

      if (detailed)
      {
        std::cout << Name() << ": consistency layer=" << layer
                  << " side=" << side
                  << " matched=" << matched
                  << " missing=" << missing
                  << " mean|dphi|=" << mean_abs_dphi
                  << " max|dphi|=" << max_abs_dphi
                  << " max|dphi|/|phistep|=" << max_over_phistep
                  << " mean|dphi_low(cdb-brd)|=" << mean_abs_dphi_low
                  << " mean|dphi_low|/|phistep|=" << mean_over_phistep_low
                  << " max|dphi_low(cdb-brd)|=" << max_abs_dphi_low
                  << " max|dphi_low|/|phistep|=" << max_over_phistep_low
                  << " mean|dphi_high(cdb-brd)|=" << mean_abs_dphi_high
                  << " mean|dphi_high|/|phistep|=" << mean_over_phistep_high
                  << " max|dphi_high(cdb-brd)|=" << max_abs_dphi_high
                  << " max|dphi_high|/|phistep|=" << max_over_phistep_high
                  << " mean|dr|=" << mean_abs_dr
                  << " max|dr|=" << max_abs_dr
                  << " worst(global/local)=" << worst_global_pad << "/" << worst_local_pad
                  << " worst_low(global/local)=" << worst_global_pad_low << "/" << worst_local_pad_low
                  << " worst_high(global/local)=" << worst_global_pad_high << "/" << worst_local_pad_high
                  << std::endl;

        if (pads_per_sector > 0)
        {
          for (int s = 0; s < NSectors; ++s)
          {
            const double sec_low_over_step = (phistep > 0.0) ? (sector_max_abs_dphi_low[s] / phistep) : -1.0;
            const double sec_high_over_step = (phistep > 0.0) ? (sector_max_abs_dphi_high[s] / phistep) : -1.0;
            std::cout << Name() << ": consistency-by-sector layer=" << layer
                      << " side=" << side
                      << " sector=" << s
                      << " matched=" << sector_matched[s]
                      << " max|dphi_low(cdb-brd)|=" << sector_max_abs_dphi_low[s]
                      << " max|dphi_low|/|phistep|=" << sec_low_over_step
                      << " max|dphi_high(cdb-brd)|=" << sector_max_abs_dphi_high[s]
                      << " max|dphi_high|/|phistep|=" << sec_high_over_step
                      << std::endl;
          }
        }
      }
      ++side_reports;
    }
  }

  std::cout << Name() << ": consistency summary layers_checked=" << layers_checked
            << " layers_missing_polygons=" << layers_missing_polygons
            << " layers_with_boundary_data=" << layers_with_boundary_data
            << " side_reports=" << side_reports
            << " global_max|dphi|=" << global_max_abs_dphi
            << " global_max|dphi_low(cdb-brd)|=" << global_max_abs_dphi_low
            << " global_max|dphi_high(cdb-brd)|=" << global_max_abs_dphi_high
            << " global_max|dr|=" << global_max_abs_dr
            << " global_max|dR_low(brd-geom)|=" << global_max_abs_dr_low
            << " global_max|dR_high(brd-geom)|=" << global_max_abs_dr_high
            << std::endl;
}

//_________________________________________________________
double PHG4TpcPadPlaneReadout::getSingleEGEMAmplification()
{
  // Jin H.: For the GEM gain in sPHENIX TPC,
  //         Bob pointed out the PHENIX HBD measured it as the Polya function with theta parameter = 0.8.
  //         Just talked with Tom too, he suggest us to start the TPC modeling with simpler exponential function
  //         with lambda parameter of 1/2000, (i.e. Polya function with theta parameter = 0, q_bar = 2000). Please note, this gain variation need to be sampled for each initial electron individually.
  //         Summing over ~30 initial electrons, the distribution is pushed towards more Gauss like.
  // Bob A.: I like Tom's suggestion to use the exponential distribution as a first approximation
  //         for the single electron gain distribution -
  //         and yes, the parameter you're looking for is of course the slope, which is the inverse gain.
  if (averageGEMGain <= 0)
  {
    return 0.0;
  }

  double nelec = gsl_ran_exponential(RandomGenerator, averageGEMGain);
  if (m_usePolya)
  {
    double y;
    double xmax = 5000;
    double ymax = 0.376;
    while (true)
    {
      nelec = gsl_ran_flat(RandomGenerator, 0, xmax);
      y = gsl_rng_uniform(RandomGenerator) * ymax;
      if (y <= pow((1 + polyaTheta) * (nelec / averageGEMGain), polyaTheta) * exp(-(1 + polyaTheta) * (nelec / averageGEMGain)))
      {
        break;
      }
    }
  }
  // Put gain reading here

  return nelec;
}

//_________________________________________________________
double PHG4TpcPadPlaneReadout::getSingleEGEMAmplification(double weight)
{
  // Jin H.: For the GEM gain in sPHENIX TPC,
  //         Bob pointed out the PHENIX HBD measured it as the Polya function with theta parameter = 0.8.
  //         Just talked with Tom too, he suggest us to start the TPC modeling with simpler exponential function
  //         with lambda parameter of 1/2000, (i.e. Polya function with theta parameter = 0, q_bar = 2000). Please note, this gain variation need to be sampled for each initial electron individually.
  //         Summing over ~30 initial electrons, the distribution is pushed towards more Gauss like.
  // Bob A.: I like Tom's suggestion to use the exponential distribution as a first approximation
  //         for the single electron gain distribution -
  //         and yes, the parameter you're looking for is of course the slope, which is the inverse gain.
  const double q_bar = averageGEMGain * weight;


  double nelec = gsl_ran_exponential(RandomGenerator, q_bar);
  if (m_usePolya)
  {
    double y;
    double xmax = 5000;
    double ymax = 0.376;
    while (true)
    {
      nelec = gsl_ran_flat(RandomGenerator, 0, xmax);
      y = gsl_rng_uniform(RandomGenerator) * ymax;
      if (y <= pow((1 + polyaTheta) * (nelec / q_bar), polyaTheta) * exp(-(1 + polyaTheta) * (nelec / q_bar)))
      {
        break;
      }
    }
  }
  // Put gain reading here

  return nelec;
}

//_________________________________________________________
double PHG4TpcPadPlaneReadout::getSingleEGEMAmplification(TF1 *f)
{
  double nelec = f->GetRandom(0,5000);
  // Put gain reading here

  return nelec;
}



void PHG4TpcPadPlaneReadout::rotatePointToSector(
    double x, double y,
    unsigned int side,
    int&    sectorFound,
    double& xNew, double& yNew
)
{
  // ---- helpers -------------------------------------------------------------
  const double PI = std::acos(-1.0);
  const double TWOPI = 2.0*PI;

  auto wrap = [&](double a) {
    while (a <= -PI) a += TWOPI;
    while (a >   PI) a -= TWOPI;
    return a;
  };

  // check a ∈ [lo,hi) with wrap at ±π
  auto inInterval = [&](double a, double lo, double hi) {
    a  = wrap(a); lo = wrap(lo); hi = wrap(hi);
    if (lo <= hi) return (a >= lo && a < hi);
    // interval crosses the branch cut
    return (a >= lo || a < hi);
  };

  // midpoint on the circle between lo..hi (shorter arc)
  auto mid = [&](double lo, double hi) {
    double d = wrap(hi - lo);
    return wrap(lo + 0.5*d);
  };

  // ---- 1) angle & sector ---------------------------------------------------
  const double R   = std::hypot(x, y);
  const double phi = wrap(std::atan2(y, x));

  // sector_min_Phi / sector_max_Phi are class members filled from LayerGeom
  sectorFound = -1;
  for (int s = 0; s < 12; ++s)
  {
    if (inInterval(phi, sector_min_Phi[side][s], sector_max_Phi[side][s]))
    { sectorFound = s; break; }
  }

  // If exactly on a boundary, pick nearest sector center
  if (sectorFound < 0)
  {
    double best = 1e9; int bestS = 0;
    for (int s = 0; s < 12; ++s)
    {
      double c = mid(sector_min_Phi[side][s], sector_max_Phi[side][s]);
      double d = std::fabs(wrap(phi - c));
      if (d < best) { best = d; bestS = s; }
    }
    sectorFound = bestS;
  }

  // ---- 2) rotate to the reference sector (sector 2) ------------------------
  // The SERF reference polygons are defined in the frame of sector 2.
  // Rotate points from the found sector into the sector-2 frame by the
  // boundary difference:
  //   side 0 (South): dphi = sec_min[found] - sec_min[2]
  //   side 1 (North): dphi = sec_max[found] - sec_max[2]
  // We then rotate the point by -dphi so that boundaries align.
  const int refSector = 2;
  const double currentBoundary = (side == 0)
      ? sector_min_Phi[side][sectorFound]
      : sector_max_Phi[side][sectorFound];
  const double refBoundary = (side == 0)
      ? sector_min_Phi[side][refSector]
      : sector_max_Phi[side][refSector];
  const double dphi = wrap(currentBoundary - refBoundary);
  const double phiRot = wrap(phi - dphi);

  xNew = R * std::cos(phiRot);
  yNew = R * std::sin(phiRot);

  // ---- 3) mirror for South side to match reference polygon handedness ------
  // Looking from the South side outward (toward +z), the x-axis is to the left
  // while the reference sector has x to the right, so flip x for side 0.
  constexpr unsigned SOUTH_SIDE = 0;
  if (side == SOUTH_SIDE)
  {
    xNew = -xNew;
  }
}




/*
//_________________________________________________________
std::vector<double> computeShaperKernel() {
    int NT = static_cast<int>(DetectorParams::window_ns / DetectorParams::adc_dt);
    std::vector<double> h(NT);
    for(int i = 0; i < NT; ++i) {
        double t = (i + 0.5) * DetectorParams::adc_dt;
        h[i] = (t / std::pow(DetectorParams::tau_shaper,2)) * std::exp(-t / DetectorParams::tau_shaper);
    }
    // normalize
    double sum = 0;
    for(double v : h) sum += v * DetectorParams::adc_dt;
    for(double &v : h) v /= sum;
    return h;
}

//_________________________________________________________

double gaussianIntegral1D(double a, double b, double mu, double sigma) {
    if(sigma <= 0) return 0.0;
    return 0.5 * (TMath::Erf((b-mu)/(std::sqrt(2)*sigma)) - TMath::Erf((a-mu)/(std::sqrt(2)*sigma)));
}

//_________________________________________________________

*/

bool PHG4TpcPadPlaneReadout::pointInPolygon(
    double x, double y,
    const std::vector<Point>& poly
) {
    size_t n = poly.size();
    if (n < 3) return false;  // no area → never inside

    bool inside = false;
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        double xi = poly[i].x, yi = poly[i].y;
        double xj = poly[j].x, yj = poly[j].y;

        // does edge (j→i) straddle the horizontal ray at y?
        bool cond = ((yi > y) != (yj > y));
        if (cond) {
            // compute intersection's X coordinate on the edge at height y
            double x_intersect = xj + (xi - xj) * (y - yj) / (yi - yj);
            if (x < x_intersect) {
                inside = !inside;
            }
        }
    }
    return inside;
}

/*
//------------------------------------------------------------------------------
// 2) findPadForPoint: returns padNumber or −1 if none
int PHG4TpcPadPlaneReadout::findPadForPoint( double x, double y, int tpc_module) {

    for (size_t i = 0; i < Pads[tpc_module].size(); ++i) {
        if (pointInPolygon(x, y, Pads[tpc_module][i].vertices)) {
            return  Pads[tpc_module][i].pad_bin;
            std::cout<<"Center of the pad is ("<<Pads[tpc_module][i].cx<<", "
                     <<Pads[tpc_module][i].cy<<")"<<std::endl;
        }
    }
    return -1;
}
*/
//------------------------------------------------------------------------------
//, TH2* h_adc_ref , TH2* h_adc_serf
void PHG4TpcPadPlaneReadout::MapToPadPlane(
    TpcClusterBuilder &tpc_truth_clusterer,
    TrkrHitSetContainer *single_hitsetcontainer,
    TrkrHitSetContainer *hitsetcontainer,
    TrkrHitTruthAssoc * /*hittruthassoc*/,
    const double x_gem, const double y_gem, const double t_gem, const unsigned int side,
    PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/ )
{
  // One electron per call of this method
  // The x_gem and y_gem values have already been randomized within the transverse drift diffusion width
  // The t_gem value already reflects the drift time of the primary electron from the production point, and is randomized within the longitudinal diffusion witdth

  double phi = atan2(y_gem, x_gem);
  if (phi > +M_PI)
  {
    phi -= 2 * M_PI;
  }
  if (phi < -M_PI)
  {
    phi += 2 * M_PI;
  }

  double rad_gem = get_r(x_gem, y_gem);

  // Moving electrons from dead area to a closest pad
  for (int iregion = 0; iregion < 3; ++iregion)
  {
    double daR = 0;
    if (iregion == 0 || iregion == 2)
    {
      daR = 1.0;  // 1.0cm edge to collect electrons from
    }
    else
    {
      daR = MinRadius[iregion] - MaxRadius[iregion - 1];
    }
    if (rad_gem <= MinRadius[iregion] && rad_gem >= MinRadius[iregion] - daR)
    {
      if (rad_gem <= MinRadius[iregion] - daR / 2)
      {
        rad_gem = MinRadius[iregion] - (1.1 * daR);
      }
      else
      {
        rad_gem = MinRadius[iregion] + 0.1 * daR;
      }
    }
  }
//std::cout<<" phi = "<<phi<<" rad_gem = "<<rad_gem<<std::endl;

  unsigned int layernum = 0;
  /* TpcClusterBuilder pass_data {}; */

  // Find which readout layer this electron ends up in

  PHG4TpcGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
  for (PHG4TpcGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    double rad_low = layeriter->second->get_radius() - layeriter->second->get_thickness() / 2.0;
    double rad_high = layeriter->second->get_radius() + layeriter->second->get_thickness() / 2.0;
//std::cout<<" rad_low "<<rad_low<<" rad_high = "<<rad_high<<std::endl;
    if (rad_gem > rad_low && rad_gem < rad_high)
    {
      // capture the layer where this electron hits the gem stack
      LayerGeom = layeriter->second;

      layernum = LayerGeom->get_layer();
      /* pass_data.layerGeom = LayerGeom; */
      /* pass_data.layer = layernum; */
      if (Verbosity() > 1000)
      {
        std::cout << " g4hit id " << hiter->first << " rad_gem " << rad_gem << " rad_low " << rad_low << " rad_high " << rad_high
                  << " layer  " << hiter->second->get_layer() << " want to change to " << layernum << std::endl;
      }
      hiter->second->set_layer(layernum);  // have to set here, since the stepping action knows nothing about layers
    }
  }

  if (layernum == 0)
  {
    if (m_enable_side_layer_debug && side < static_cast<unsigned int>(NSides))
    {
      ++m_side_layer_no_layer[side];
    }
    return;
  }

  SideLayerDebugCounters *debug_counter = nullptr;
  std::uint64_t bins_written_this_cloud = 0;
  if (m_enable_side_layer_debug && side < static_cast<unsigned int>(NSides))
  {
    auto &counter = m_side_layer_debug_counters[{side, layernum}];
    counter.cloud_calls++;
    debug_counter = &counter;
  }

  // store phi bins and tbins upfront to avoid repetitive checks on the phi methods
  /* const auto phibins = LayerGeom->get_phibins(); */
  /* const auto tbins = LayerGeom->get_zbins(); */

  sector_min_Phi = LayerGeom->get_sector_min_phi();
  sector_max_Phi = LayerGeom->get_sector_max_phi();
  phi_bin_width = LayerGeom->get_phistep();

  phi = check_phi(side, phi, rad_gem);

  // Create the distribution function of charge on the pad plane around the electron position

  // The resolution due to pad readout includes the charge spread during GEM multiplication.
  // this now defaults to 400 microns during construction from Tom (see 8/11 email).
  // Use the setSigmaT(const double) method to update...
  // We use a double gaussian to represent the smearing due to the SAMPA chip shaping time - default values of fShapingLead and fShapingTail are for 80 ns SAMPA

  // amplify the single electron in the gem stack
  //===============================

  double nelec = getSingleEGEMAmplification();

  // consider all layers whose radial annulus intersects a 5σ window
  const double nsig_window = _nsigmas;

  std::vector<unsigned int> cand_layers = layersInRadialWindow(rad_gem, sigmaT, nsig_window);
  if (cand_layers.empty()) {
    cand_layers.push_back(layernum);
  }


  // Applying weight with respect to the rad_gem and phi after electrons are redistributed
  double phi_gain = phi;
  if (phi < 0)
  {
    phi_gain += 2 * M_PI;
  }
  double gain_weight = 1.0;
  if (m_flagToUseGain == 1)
  {
    gain_weight = h_gain[side]->GetBinContent(h_gain[side]->FindBin(rad_gem * 10, phi_gain));  // rad_gem in cm -> *10 to get mm
    nelec = nelec * gain_weight;
  }

  if(m_use_module_gain_weights)
    {
      double phistep = 30.0;
      int sector = 0;

      if( (phi_gain*180.0/M_PI) >=15 && (phi_gain*180.0 / M_PI) < 345)
	{
	  sector = 1 + (int) ( (phi_gain*180.0/M_PI - 15) / phistep);
	}
      else
	{
	  sector = 0;
	}

      int this_region = -1;
      for (int iregion = 0; iregion < 3; ++iregion)
	{
	  if (rad_gem < MaxRadius[iregion] && rad_gem > MinRadius[iregion])
	    {
	      this_region = iregion;
	    }
	}
      if(this_region > -1)
	{
	  gain_weight = m_module_gain_weight[side][this_region][sector];
	}
      // regenerate nelec with the new distribution
      //    double original_nelec = nelec;
      nelec = getSingleEGEMAmplification(gain_weight);
      //  std::cout << " side " << side << " this_region " << this_region
      //	<<  " sector " << sector << " original nelec "
      //	<< original_nelec << " new nelec " << nelec << std::endl;
    }

  if(m_useLangau)
  {
    double phistep = 30.0;
    int sector = 0;

    if( (phi_gain*180.0/M_PI) >=15 && (phi_gain*180.0 / M_PI) < 345)
    {
      sector = 1 + (int) ( (phi_gain*180.0/M_PI - 15) / phistep);
    }
    else
    {
      sector = 0;
    }

    int this_region = -1;
    for (int iregion = 0; iregion < 3; ++iregion)
    {
      if (rad_gem < MaxRadius[iregion] && rad_gem > MinRadius[iregion])
      {
	this_region = iregion;
      }
    }
    if(this_region > -1)
    {
      nelec = getSingleEGEMAmplification(flangau[side][this_region][sector]);
    }
    else
    {
      nelec = getSingleEGEMAmplification();
    }
  }

  // std::cout<<"PHG4TpcPadPlaneReadout::MapToPadPlane gain_weight = "<<gain_weight<<std::endl;
  /* pass_data.neff_electrons = nelec; */

  // Distribute the charge between the pads in phi
  //====================================

  if (Verbosity() > 200)
  {
    std::cout << "  populate phi bins for "
              << " layernum " << layernum
              << " phi " << phi
              << " rad_gem " << rad_gem
              << " sigmaT " << sigmaT
              //<< " zigzag_pads " << zigzag_pads
              << std::endl;
  }



  // NOTE: we no longer compute single-layer pad sharing here. Instead, we
  // compute per-layer pad shares inside the candidate-layers loop below
  // (using SERF_zigzag_phibins) to avoid double counting and to include
  // all pads within 5σ radially across layers.
/*
norm1 = 0.0;
  for (unsigned int ipad = 0; ipad < pad_phibin_ref.size(); ++ipad)
  {
    double pad_share = pad_phibin_share_ref[ipad];
    norm1 += pad_share;
  }
  for (unsigned int iphi = 0; iphi < pad_phibin_ref.size(); ++iphi)
  {
    pad_phibin_share_ref[iphi] /= norm1;
  }*/

  /*std::cout<<"--------------------------"<<std::endl;
    for (unsigned int iphi = 0; iphi < pad_phibin.size(); ++iphi)
  {
   //  std::cout<<" phibin ref "<<pad_phibin[iphi]<<" phibin serf "<<pad_phibin_serf[iphi]<<" charge ref "<< pad_phibin_share[iphi]<<" charge serf "<< pad_phibin_share_serf[iphi]<<std::endl;
    std::cout<<" phibin ref "<<pad_phibin_ref[iphi]<<" phibin serf "<<pad_phibin[iphi]<<" charge ref "<< pad_phibin_share_ref[iphi]<<" charge serf "<< pad_phibin_share[iphi]<<std::endl;

  }
  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;*/
  // Distribute the charge between the pads in t
  //====================================
  if (Verbosity() > 100 && layernum == print_layer)
  {
    std::cout << "  populate t bins for layernum " << layernum
              << " with t_gem " << t_gem << " sigmaL[0] " << sigmaL[0] << " sigmaL[1] " << sigmaL[1] << std::endl;
  }

  std::vector<int> adc_tbin;
  std::vector<double> adc_tbin_share;
  sampaTimeDistribution(t_gem, adc_tbin, adc_tbin_share);
  /* if (adc_tbin.size() == 0)  { */
  /* pass_data.neff_electrons = 0; */
  /* } else { */
  /* pass_data.fillTimeBins(adc_tbin); */
  /* } */

  // Normalize the shares so that they add up to 1
  double tnorm = 0.0;
  for (unsigned int it = 0; it < adc_tbin.size(); ++it)
  {
    double bin_share = adc_tbin_share[it];
    tnorm += bin_share;
  }
  for (unsigned int it = 0; it < adc_tbin.size(); ++it)
  {
    adc_tbin_share[it] /= tnorm;
  }
  // accumulate quick centroid info across all candidate layers
  double phi_integral = 0.0;
  double t_integral = 0.0;
  double weight = 0.0;

  // Determine if SERF pad sharing is possible (pad polygons loaded)
  bool have_polygons = true;
  for (unsigned int layer_cand : cand_layers)
  {
    if (layer_cand >= Pads.size() || Pads[layer_cand].empty())
    {
      have_polygons = false;
      break;
    }
  }

  if (m_use_serf_padsharing && have_polygons)
  {
    // Pass 1: compute total Gaussian mass across all candidate layers
    double total_mass = 0.0;
    for (unsigned int layer_cand : cand_layers)
    {
      PHG4TpcGeom* thisGeom = getGeomForLayer(layer_cand);
      if (!thisGeom) continue;
      LayerGeom = thisGeom;
      sector_min_Phi = LayerGeom->get_sector_min_phi();
      sector_max_Phi = LayerGeom->get_sector_max_phi();
      phi_bin_width  = LayerGeom->get_phistep();
      const double phi_for_layer = check_phi(side, phi, rad_gem);

      std::vector<int>    pad_phibin_tmp;
      std::vector<double> pad_mass_tmp; // unnormalized integrated mass per pad
      SERF_zigzag_phibins(side, layer_cand, phi_for_layer, rad_gem, sigmaT, pad_phibin_tmp, pad_mass_tmp);
      for (double v : pad_mass_tmp) total_mass += v;
    }

    // Normalization strategy:
    //  - legacy: renormalize to recovered pad mass (fills polygon gaps)
    //  - preserve gap losses: use raw SERF overlap integral directly
    //    (no renormalization), i.e. pad_fraction = pad_mass[i].
    double normalization_mass = total_mass;
    if (m_preserve_pad_gap_losses)
    {
      normalization_mass = 1.0;
    }
    else if (normalization_mass <= 1e-16)
    {
      normalization_mass = 1.0;  // avoid division by zero
    }

    if (Verbosity() > 1)
    {
      std::cout << Name() << ": SERF normalization side=" << side
                << " layer_center=" << layernum
                << " recovered_mass=" << total_mass
                << " norm_mass=" << normalization_mass
                << " preserve_gap_losses=" << static_cast<int>(m_preserve_pad_gap_losses)
                << std::endl;
    }

    // Pass 2: fill hits with correct per-layer mass fraction
    for (unsigned int layer_cand : cand_layers)
    {
      PHG4TpcGeom* thisGeom = getGeomForLayer(layer_cand);
      if (!thisGeom) continue;
      LayerGeom = thisGeom;
      sector_min_Phi = LayerGeom->get_sector_min_phi();
      sector_max_Phi = LayerGeom->get_sector_max_phi();
      phi_bin_width  = LayerGeom->get_phistep();
      const double phi_for_layer = check_phi(side, phi, rad_gem);

      std::vector<int>    pad_phibin;
      std::vector<double> pad_mass; // unnormalized integrated mass per pad
      SERF_zigzag_phibins(side, layer_cand, phi_for_layer, rad_gem, sigmaT, pad_phibin, pad_mass);

      const auto phibins = LayerGeom->get_phibins();
      const auto tbins   = LayerGeom->get_zbins();
      const unsigned int pads_per_sector = phibins / 12;

      for (size_t i = 0; i < pad_phibin.size(); ++i)
      {
        const int pad_num = pad_phibin[i];
        const double pad_fraction = pad_mass[i] / normalization_mass;
        if (pad_fraction <= 0) continue;
        if (debug_counter)
        {
          ++debug_counter->pad_contributors;
        }

        const unsigned int sector = (pad_num >= 0) ? (static_cast<unsigned>(pad_num) / pads_per_sector) : 0;
        TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layer_cand, sector, side);
        auto hitsetit        = hitsetcontainer->findOrAddHitSet(hitsetkey);
        auto single_hitsetit = single_hitsetcontainer->findOrAddHitSet(hitsetkey);

        for (unsigned int itb = 0; itb < adc_tbin.size(); ++itb)
        {
          const int tbin_num = adc_tbin[itb];
          if (tbin_num < 0 || tbin_num >= static_cast<int>(tbins)) continue;
          if (pad_num  < 0 || pad_num  >= static_cast<int>(phibins)) continue;

          const double tshare = adc_tbin_share[itb];
          const float neffelectrons_bin = nelec * pad_fraction * tshare;
          if (neffelectrons_bin < neffelectrons_threshold) continue;

          TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(static_cast<unsigned>(pad_num), static_cast<unsigned>(tbin_num));
          TrkrHit* hit = hitsetit->second->getHit(hitkey);
          if (!hit) { hit = new TrkrHitv2(); hitsetit->second->addHitSpecificKey(hitkey, hit); }
          hit->addEnergy(neffelectrons_bin);

          TrkrHit* single_hit = single_hitsetit->second->getHit(hitkey);
          if (!single_hit) { single_hit = new TrkrHitv2(); single_hitsetit->second->addHitSpecificKey(hitkey, single_hit); }
          single_hit->addEnergy(neffelectrons_bin);

          if (Verbosity() > 50)
          {
            const int layer_print = static_cast<int>(layer_cand);
            if (m_dbg_pad_hit_layer < 0 || layer_print == m_dbg_pad_hit_layer)
            {
              std::cout << "PadHit: layer=" << layer_cand
                        << " side=" << side
                        << " pad=" << pad_num
                        << " tbin=" << tbin_num
                        << " energy=" << neffelectrons_bin
                        << std::endl;
            }
          }

          tpc_truth_clusterer.addhitset(hitsetkey, hitkey, neffelectrons_bin);

          const double tcenter  = LayerGeom->get_zcenter(tbin_num);
          const double phicenter = LayerGeom->get_phicenter(pad_num, side);
          phi_integral += phicenter * neffelectrons_bin;
          t_integral   += tcenter   * neffelectrons_bin;
          weight       += neffelectrons_bin;
          if (debug_counter)
          {
            ++debug_counter->bins_written;
            debug_counter->neff_sum += neffelectrons_bin;
            ++bins_written_this_cloud;
          }
        }
      }
    }
  }
  else
  {
    if (m_use_serf_padsharing && !have_polygons && !m_warned_serf_fallback)
    {
      std::cout << "PHG4TpcPadPlaneReadout: SERF polygons incomplete for some layers; per-hit fallback to analytic sharing is active." << std::endl;
      m_warned_serf_fallback = true;
    }
    // Analytic triangular response across all layers within 5σ
    // Pass 1: compute total radial weight across candidate layers
    const double inv_sqrt2_sigma = 1.0 / (M_SQRT2 * sigmaT);
    double total_radw = 0.0;
    for (unsigned int layer_cand : cand_layers)
    {
      PHG4TpcGeom* thisGeom = getGeomForLayer(layer_cand);
      if (!thisGeom) continue;
      const double rcen = thisGeom->get_radius();
      const double thk  = thisGeom->get_thickness();
      const double rlow = rcen - 0.5*thk;
      const double rhigh= rcen + 0.5*thk;
      const double dx1  = (rhigh - rad_gem) * inv_sqrt2_sigma;
      const double dx0  = (rlow  - rad_gem) * inv_sqrt2_sigma;
      const double radw = 0.5 * (std::erf(dx1) - std::erf(dx0));
      if (radw > 0) total_radw += radw;
    }
    if (total_radw <= 1e-16) total_radw = 1.0;

    // Pass 2: per-layer phi sharing and fill
    for (unsigned int layer_cand : cand_layers)
    {
      PHG4TpcGeom* thisGeom = getGeomForLayer(layer_cand);
      if (!thisGeom) continue;
      LayerGeom = thisGeom;
      sector_min_Phi = LayerGeom->get_sector_min_phi();
      sector_max_Phi = LayerGeom->get_sector_max_phi();
      phi_bin_width  = LayerGeom->get_phistep();

      const double rcen = LayerGeom->get_radius();
      const double thk  = LayerGeom->get_thickness();
      const double rlow = rcen - 0.5*thk;
      const double rhigh= rcen + 0.5*thk;
      const double dx1  = (rhigh - rad_gem) * inv_sqrt2_sigma;
      const double dx0  = (rlow  - rad_gem) * inv_sqrt2_sigma;
      const double radw = std::max(0.0, 0.5 * (std::erf(dx1) - std::erf(dx0)));
      if (radw <= 0) continue;

      const double phi_for_layer = check_phi(side, phi, rad_gem);
      std::vector<int>    pad_phibin;
      std::vector<double> pad_share_phi;
      populate_zigzag_phibins(side, layer_cand, phi_for_layer, sigmaT, pad_phibin, pad_share_phi);

      double norm_phi = 0.0; for (double v : pad_share_phi) norm_phi += v;
      if (norm_phi <= 1e-16) continue;

      const auto phibins = LayerGeom->get_phibins();
      const auto tbins   = LayerGeom->get_zbins();
      const unsigned int pads_per_sector = phibins / 12;

      for (size_t i = 0; i < pad_phibin.size(); ++i)
      {
        const int pad_num = pad_phibin[i];
        const double pshare_phi = pad_share_phi[i] / norm_phi;
        const double pad_fraction = (radw / total_radw) * pshare_phi;
        if (pad_fraction <= 0) continue;
        if (debug_counter)
        {
          ++debug_counter->pad_contributors;
        }

        const unsigned int sector = (pad_num >= 0) ? (static_cast<unsigned>(pad_num) / pads_per_sector) : 0;
        TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layer_cand, sector, side);
        if (m_maskDeadChannels)
{
  TrkrDefs::hitkey checkkey = TpcDefs::genHitKey((unsigned int) pad_num, 0);
  if (m_deadChannelMap.contains(hitsetkey) &&
      std::find(m_deadChannelMap[hitsetkey].begin(), m_deadChannelMap[hitsetkey].end(), checkkey) != m_deadChannelMap[hitsetkey].end())
  {
    continue; // Skip this hit, the channel is dead
  }
}
if (m_maskHotChannels)
{
  TrkrDefs::hitkey checkkey = TpcDefs::genHitKey((unsigned int) pad_num, 0);
  if (m_hotChannelMap.contains(hitsetkey) &&
      std::find(m_hotChannelMap[hitsetkey].begin(), m_hotChannelMap[hitsetkey].end(), checkkey) != m_hotChannelMap[hitsetkey].end())
  {
    continue; // Skip this hit, the channel is hot
  }
}
        auto hitsetit        = hitsetcontainer->findOrAddHitSet(hitsetkey);
        auto single_hitsetit = single_hitsetcontainer->findOrAddHitSet(hitsetkey);

        for (unsigned int itb = 0; itb < adc_tbin.size(); ++itb)
        {
          const int tbin_num = adc_tbin[itb];
          if (tbin_num < 0 || tbin_num >= static_cast<int>(tbins)) continue;
          if (pad_num  < 0 || pad_num  >= static_cast<int>(phibins)) continue;

          const double tshare = adc_tbin_share[itb];
          const float neffelectrons_bin = nelec * pad_fraction * tshare;
          if (neffelectrons_bin < neffelectrons_threshold) continue;

          TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(static_cast<unsigned>(pad_num), static_cast<unsigned>(tbin_num));
          TrkrHit* hit = hitsetit->second->getHit(hitkey);
          if (!hit) { hit = new TrkrHitv2(); hitsetit->second->addHitSpecificKey(hitkey, hit); }
          hit->addEnergy(neffelectrons_bin);

          TrkrHit* single_hit = single_hitsetit->second->getHit(hitkey);
          if (!single_hit) { single_hit = new TrkrHitv2(); single_hitsetit->second->addHitSpecificKey(hitkey, single_hit); }
          single_hit->addEnergy(neffelectrons_bin);

          if (Verbosity() > 50)
          {
            const int layer_print = static_cast<int>(layer_cand);
            if (m_dbg_pad_hit_layer < 0 || layer_print == m_dbg_pad_hit_layer)
            {
              std::cout << "PadHit: layer=" << layer_cand
                        << " side=" << side
                        << " pad=" << pad_num
                        << " tbin=" << tbin_num
                        << " energy=" << neffelectrons_bin
                        << std::endl;
            }
          }

          tpc_truth_clusterer.addhitset(hitsetkey, hitkey, neffelectrons_bin);

          const double tcenter  = LayerGeom->get_zcenter(tbin_num);
          const double phicenter = LayerGeom->get_phicenter(pad_num, side);
          phi_integral += phicenter * neffelectrons_bin;
          t_integral   += tcenter   * neffelectrons_bin;
          weight       += neffelectrons_bin;
          if (debug_counter)
          {
            ++debug_counter->bins_written;
            debug_counter->neff_sum += neffelectrons_bin;
            ++bins_written_this_cloud;
          }
        }
      }
    }
  }
  // end of pad sharing and hit filling across candidate layers
  /* pass_data.phi_integral = phi_integral; */
  /* pass_data.time_integral = t_integral; */

  /*
  // Capture the input values at the gem stack and the quick clustering results, elecron-by-electron
  if (Verbosity() > 0)
  {
    assert(ntpad);
    ntpad->Fill(layernum, phi, phi_integral / weight, t_gem, t_integral / weight);
  }
  */

  if (Verbosity() > 100 && weight > 0)
  {
    std::cout << " hit " << m_NHits << " quick centroid for this electron " << std::endl;
    std::cout << "      phi centroid = " << phi_integral / weight << " phi in " << phi << " phi diff " << phi_integral / weight - phi << std::endl;
    std::cout << "      t centroid = " << t_integral / weight << " t in " << t_gem << " t diff " << t_integral / weight - t_gem << std::endl;
  }

  if (debug_counter && bins_written_this_cloud == 0)
  {
    ++debug_counter->cloud_no_hits;
  }

  m_NHits++;
  /* return pass_data; */
}
double PHG4TpcPadPlaneReadout::check_phi(const unsigned int side, const double phi, const double radius)
{
  double new_phi = phi;
  int p_region = -1;
  for (int iregion = 0; iregion < 3; ++iregion)
  {
    if (radius < MaxRadius[iregion] && radius > MinRadius[iregion])
    {
      p_region = iregion;
    }
  }

  if (p_region >= 0)
  {
    for (int s = 0; s < 12; s++)
    {
      double daPhi = 0;
      if (s == 0)
      {
        daPhi = std::fabs(sector_min_Phi[side][11] + 2 * M_PI - sector_max_Phi[side][s]);
      }
      else
      {
        daPhi = std::fabs(sector_min_Phi[side][s - 1] - sector_max_Phi[side][s]);
      }
      double min_phi = sector_max_Phi[side][s];
      double max_phi = sector_max_Phi[side][s] + daPhi;
      if (new_phi <= max_phi && new_phi >= min_phi)
      {
        if (std::fabs(max_phi - new_phi) > std::fabs(new_phi - min_phi))
        {
          new_phi = min_phi - phi_bin_width / 5;
        }
        else
        {
          new_phi = max_phi + phi_bin_width / 5;
        }
      }
    }
    if (new_phi < sector_min_Phi[side][11] && new_phi >= -M_PI)
    {
      new_phi += 2 * M_PI;
    }
  }

  return new_phi;
}

void PHG4TpcPadPlaneReadout::EnableSingleCloudVisualization(bool enable,
                                                            const std::string &output_file,
                                                            int target_side,
                                                            int target_layer,
                                                            double grid_step)
{
  m_visualize_single_cloud = enable;
  m_visualization_output = output_file;
  m_visualization_target_side = target_side;
  m_visualization_target_layer = target_layer;
  m_visualization_grid_step = grid_step;
  m_visualization_done = false;
  if (enable)
  {
    m_visualization_dump_index = 0;
    m_visualization_cloud_counter = 0;
    m_visualization_aggregate_samples.clear();
    m_visualization_circles.clear();
    m_visualization_pad_union.clear();
  }
}


void PHG4TpcPadPlaneReadout::SetVisualizationDumpFile(const std::string &file)
{
  m_visualization_dump_file = file;
  m_visualization_dump_index = 0;
}

void PHG4TpcPadPlaneReadout::SetVisualizeAllClouds(bool enable)
{
  m_visualize_all_matches = enable;
  m_visualization_cloud_counter = 0;
  m_visualization_aggregate_samples.clear();
  m_visualization_circles.clear();
  m_visualization_pad_union.clear();
}
/*
void PHG4TpcPadPlaneReadout::build_serf_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &phibin_pad, std::vector<double> &phibin_pad_share)
{
    const double radius = LayerGeom->get_radius();
    const double phistepsize = LayerGeom->get_phistep();
    const auto phibins = LayerGeom->get_phibins();

    double rphi = phi * radius;
    if (Verbosity() > 100)
    {
      if (LayerGeom->get_layer() == print_layer)
      {
        std::cout << " populate_zigzag_phibins for layer " << layernum << " with radius " << radius << " phi " << phi
                  << " rphi " << rphi << " phistepsize " << phistepsize << std::endl;
        std::cout << " fcharge created: radius " << radius << " rphi " << rphi << " cloud_sig_rp " << cloud_sig_rp << std::endl;
      }
    }


}
*/
double PHG4TpcPadPlaneReadout::integratedDensityOfCircleAndPad(
    double hitX,
    double hitY,
    double sigma,
    const std::vector<Point>& pad,
    double gridStep,
    std::vector<DebugSample>* debug_samples
) {
    // sanity checks
    if (sigma <= 0.0 || pad.empty()) return 0.0;

    // constants
    const double    nSigma       = _nsigmas;
    double R                  = nSigma * sigma;
    double gaussConst         = 1.0 / (2.0 * M_PI * sigma * sigma);
    double expDenominator     = 2.0 * sigma * sigma;

    // pick a default grid step if none provided
    if (gridStep <= 0.0) {
        gridStep = sigma / 50.0;
    }

    // 1) Compute pad bounding box
    // Initialize min/max from the same first vertex to ensure both extremas consider index 0
    double minx = pad[0].x, maxx = pad[0].x;
    double miny = pad[0].y, maxy = pad[0].y;
    for (size_t i = 1; i < pad.size(); ++i) {
        minx = std::min(minx, pad[i].x);
        maxx = std::max(maxx, pad[i].x);
        miny = std::min(miny, pad[i].y);
        maxy = std::max(maxy, pad[i].y);
     //  std::cout<<"Pad point "<<i<<" x = "<<pad[i].x<<", y = "<<pad[i].y<<std::endl;
    }
//std::cout<<"Pad 0: ( "<<pad[0].x<<" ; "<<pad[0].y<<" ) - ( "<<pad.back().x<<" ; "<<pad.back().y<<" );  Pad box "<<" minx = "<<minx<<", maxx = "<<maxx<<", miny = "<<miny<<", maxy = "<<maxy<<std::endl;
    // 2) Intersect that box with the nσ circle’s box (n = _nsigmas)
    double x0 = std::max(minx, hitX - R);
    double x1 = std::min(maxx, hitX + R);
    double y0 = std::max(miny, hitY - R);
    double y1 = std::min(maxy, hitY + R);
//std::cout<<"Circle box "<<" x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl;

    if (x1 <= x0 || y1 <= y0) {
        // no overlap
        return 0.0;
    }

    // 3) Determine grid dimensions
    int nx = static_cast<int>(std::ceil((x1 - x0) / gridStep));
    int ny = static_cast<int>(std::ceil((y1 - y0) / gridStep));
    //std::cout<<"Grid dimensions: nx = "<<nx<<", ny = "<<ny<<" gridStep = "<<gridStep<<" gridStep^2 = "<<gridStep*gridStep<<" gaussConst = "<<gaussConst<<" expDenominator = "<<expDenominator<<std::endl;
    double total = 0.0;

    // 4) Loop over cell centers
    for (int ix = 0; ix < nx; ++ix) {
        double x = x0 + (ix + 0.5) * gridStep;
        for (int iy = 0; iy < ny; ++iy) {
            double y = y0 + (iy + 0.5) * gridStep;
            // skip outside the circle
            double dx = x - hitX, dy = y - hitY;

            if (dx*dx + dy*dy > R*R) continue;
            // skip outside the pad polygon
            if (!pointInPolygon(x, y, pad)) continue;
            // accumulate Gaussian density
            const double density = gaussConst * std::exp(-(dx*dx + dy*dy) / expDenominator);
            total += density;
            if (debug_samples)
            {
              debug_samples->push_back({x, y, density});
            }
           // std::cout<<" dx = "<<dx<<" dy = "<<dy<<" exp = "<<std::exp(-(dx*dx + dy*dy))<<std::endl;
            // std::cout<<"Cell center ("<<x<<", "<<y<<") contributes to integral with density = "
            //    <<gaussConst * std::exp(-(dx*dx + dy*dy) / expDenominator)<<std::endl;
        }
    }
   //std::cout<<"  total = "<<total<<std::endl;
    // 5) multiply by cell area for the approximate integral
    return total * (gridStep * gridStep);
}

void PHG4TpcPadPlaneReadout::maybeVisualizeAvalanche(
    unsigned int side,
    unsigned int layernum,
    double phi,
    double rad_gem,
    double cloud_sig_rp,
    double x_center,
    double y_center,
    const std::vector<DebugPadContribution> &contribs,
    const std::vector<DebugSample> &samples,
    const std::vector<VisualizationCircle> &circles)
{
  if (!m_visualize_single_cloud) return;
  if (m_visualization_done && !m_visualize_all_matches) return;
  if (contribs.empty()) return;
  if (cloud_sig_rp <= 0.0) return;
  if (m_visualization_target_layer >= 0 &&
      static_cast<int>(layernum) != m_visualization_target_layer) return;
  if (m_visualization_target_side >= 0 &&
      static_cast<int>(side) != m_visualization_target_side) return;

  double xmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::lowest();
  double ymin = std::numeric_limits<double>::max();
  double ymax = std::numeric_limits<double>::lowest();
  bool have_bounds = false;
  for (const auto &sample : samples)
  {
    xmin = std::min(xmin, sample.x);
    xmax = std::max(xmax, sample.x);
    ymin = std::min(ymin, sample.y);
    ymax = std::max(ymax, sample.y);
    have_bounds = true;
  }

  double max_circle_radius = 0.0;
  for (const auto &circle : circles) max_circle_radius = std::max(max_circle_radius, circle.radius);
  for (const auto &pad : contribs)
  {
    for (const auto &v : pad.polygon)
    {
      xmin = std::min(xmin, v.x);
      xmax = std::max(xmax, v.x);
      ymin = std::min(ymin, v.y);
      ymax = std::max(ymax, v.y);
      have_bounds = true;
    }
  }
  for (const auto &circle : circles)
  {
    xmin = std::min(xmin, circle.x - circle.radius);
    xmax = std::max(xmax, circle.x + circle.radius);
    ymin = std::min(ymin, circle.y - circle.radius);
    ymax = std::max(ymax, circle.y + circle.radius);
    have_bounds = true;
  }

  if (!have_bounds)
  {
    const double fallback = std::max(3.0 * cloud_sig_rp, 0.2);
    xmin = x_center - fallback;
    xmax = x_center + fallback;
    ymin = y_center - fallback;
    ymax = y_center + fallback;
  }

  const double span = std::max(xmax - xmin, ymax - ymin);
  const double margin = std::max({0.12 * span, 0.6 * max_circle_radius, 0.15});
  xmin -= margin;
  xmax += margin;
  ymin -= margin;
  ymax += margin;

  double grid_step = m_visualization_grid_step;
  if (grid_step <= 0.0)
  {
    grid_step = cloud_sig_rp / 50.0;
  }
  grid_step = std::max(grid_step, cloud_sig_rp / 200.0);

  const int raw_nx = std::max(20, static_cast<int>(std::ceil((xmax - xmin) / grid_step)));
  const int raw_ny = std::max(20, static_cast<int>(std::ceil((ymax - ymin) / grid_step)));
  const int max_bins = 600;
  const int nx = std::min(max_bins, raw_nx);
  const int ny = std::min(max_bins, raw_ny);

  const std::string hname = "h_cloud_side" + std::to_string(side) + "_layer" + std::to_string(layernum);
  TH2F *hist = new TH2F(hname.c_str(), "", nx, xmin, xmax, ny, ymin, ymax);
  hist->SetTitle(";sector frame x [cm];sector frame y [cm];overlap density");
  hist->SetStats(false);

  for (const auto &sample : samples)
  {
    hist->Fill(sample.x, sample.y, sample.density);
  }

  double total_mass = 0.0;
  double max_bin = 0.0;
  double min_positive = std::numeric_limits<double>::max();
  int nonzero_bins = 0;
  for (int ix = 1; ix <= nx; ++ix)
  {
    for (int iy = 1; iy <= ny; ++iy)
    {
      const double val = hist->GetBinContent(ix, iy);
      total_mass += val;
      if (val > 0.0)
      {
        ++nonzero_bins;
        if (val > max_bin) max_bin = val;
        if (val < min_positive) min_positive = val;
      }
    }
  }
  double floor_fill = 0.0;
  if (Verbosity() > 0)
  {
    std::cout << "test_1_PHG4TpcPadPlaneReadout: histogram summary for visualization "
              << "(side " << side << ", layer " << layernum << "): total mass = "
              << total_mass << ", non-zero bins = " << nonzero_bins
              << ", sampled_points = " << samples.size();
    if (nonzero_bins > 0)
    {
      std::cout << ", min(bin) = " << min_positive << ", max(bin) = " << max_bin;
    }
    std::cout << std::endl;
  }

  const std::string cname = "c_cloud_side" + std::to_string(side) + "_layer" + std::to_string(layernum);
  TCanvas *canvas = new TCanvas(cname.c_str(), "Avalanche cloud vs zigzag pads", 900, 800);
  canvas->SetFillColor(kWhite);
  canvas->cd();
  if (gStyle)
  {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
  }
  gPad->SetFillColor(kWhite);
  gPad->SetTicks(1, 1);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.13);
  gPad->SetRightMargin(0.18);
  hist->SetDirectory(nullptr);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetTitleSize(0.050);
  hist->GetYaxis()->SetTitleSize(0.050);
  hist->GetZaxis()->SetTitleSize(0.032);
  hist->GetXaxis()->SetLabelSize(0.040);
  hist->GetYaxis()->SetLabelSize(0.040);
  hist->GetZaxis()->SetLabelSize(0.0);
  hist->GetXaxis()->SetTitleOffset(1.05);
  hist->GetYaxis()->SetTitleOffset(1.25);
  hist->GetZaxis()->SetTitleOffset(1.45);
  hist->GetZaxis()->SetTitle("Charge-cloud density on readout pads [arb. units]");
  hist->GetZaxis()->SetTickLength(0.0);
  if (nonzero_bins > 0)
  {
    const double min_disp = std::max(0.02 * max_bin, 1e-9);
    hist->SetMinimum(min_disp);
    hist->SetMaximum(1.1 * max_bin);
  }
  else
  {
    hist->SetMinimum(1e-9);
    hist->SetMaximum(1.0);
  }
  hist->SetContour(255);
  hist->Draw("COLZ");

  const std::vector<DebugPadContribution>* pads_to_draw = &contribs;
  std::vector<DebugPadContribution> pad_aggregated;
  if (m_visualize_all_matches)
  {
    pad_aggregated.reserve(m_visualization_pad_union.size());
    for (const auto &kv : m_visualization_pad_union) pad_aggregated.push_back(kv.second);
    pads_to_draw = &pad_aggregated;
  }

  std::vector<TGraph *> pad_graphs;
  pad_graphs.reserve(pads_to_draw->size());
  const int zigzag_colors[3] = {kRed + 1, kGray + 2, kMagenta + 2};
  std::size_t zigzag_color_index = 0;

  for (const auto &pad : *pads_to_draw)
  {
    if (pad.polygon.empty())
    {
      continue;
    }
    const size_t npts = pad.polygon.size();
    std::vector<double> xs(npts + 1);
    std::vector<double> ys(npts + 1);
    for (size_t i = 0; i < npts; ++i)
    {
      xs[i] = pad.polygon[i].x;
      ys[i] = pad.polygon[i].y;
    }
    xs[npts] = pad.polygon.front().x;
    ys[npts] = pad.polygon.front().y;

    TGraph *outline = new TGraph(static_cast<int>(xs.size()), xs.data(), ys.data());
    const int outline_color = zigzag_colors[zigzag_color_index % 3];
    ++zigzag_color_index;
    outline->SetLineColor(outline_color);
    outline->SetLineWidth(1);
    outline->SetFillStyle(0);
    outline->Draw("L SAME");
    pad_graphs.push_back(outline);
  }

  std::vector<TGraph *> center_markers;
  center_markers.reserve(circles.size());
  for (size_t ic = 0; ic < circles.size(); ++ic)
  {
    const auto &circle = circles[ic];
    double center_x[1] = {circle.x};
    double center_y[1] = {circle.y};
    TGraph *center = new TGraph(1, center_x, center_y);
    center->SetMarkerStyle(29);
    //center->SetMarkerSize(ic + 1 == circles.size() ? 1.8 : 1.2);
    center->SetMarkerSize(0.55);
    center->SetMarkerColor(kBlack);
    center->Draw("P SAME");
    center_markers.push_back(center);
  }

  TMarker *electron_label_marker = new TMarker(0.62, 0.88, 29);
  electron_label_marker->SetNDC(kTRUE);
  electron_label_marker->SetMarkerColor(kBlack);
  electron_label_marker->SetMarkerSize(0.65);
  electron_label_marker->Draw("SAME");

  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextColor(kBlack);
  label.SetTextAlign(13);
  label.SetTextSize(0.038);
  label.DrawLatex(0.18, 0.88, "#bf{#it{sPHENIX}} TPC Simulation");
  label.SetTextSize(0.032);
  label.DrawLatex(0.18, 0.825, "single #pi^{+}, p_{T}=10 GeV");
  label.SetTextAlign(13);
  label.SetTextSize(0.030);
  label.DrawLatex(0.64, 0.88, "ionizing electrons");

  if (!m_visualization_dump_file.empty())
  {
    TFile dumpFile(m_visualization_dump_file.c_str(), "UPDATE");
    if (!dumpFile.IsOpen())
    {
      std::cout << "PHG4TpcPadPlaneReadout: unable to open "
                << m_visualization_dump_file
                << " for visualization dump." << std::endl;
    }
    else
    {
      const std::string hname_save =
          "cloudHist_" + std::to_string(m_visualization_dump_index) +
          "_side" + std::to_string(side) +
          "_layer" + std::to_string(layernum);
      TH2F histOut(*hist);
      histOut.SetName(hname_save.c_str());
      histOut.SetTitle(hist->GetTitle());
      histOut.SetDirectory(&dumpFile);
      histOut.Write(hname_save.c_str(), TObject::kOverwrite);

      TNtuple *meta = dynamic_cast<TNtuple *>(dumpFile.Get("cloud_meta"));
      if (!meta)
      {
        meta = new TNtuple("cloud_meta",
                           "SERF cloud metadata",
                           "index:side:layer:centerX:centerY:phi:radius:sigma:floor:minBin:maxBin:totalMass:nonZero");
      }
      meta->SetDirectory(&dumpFile);
      const float minBinOut = (nonzero_bins > 0) ? static_cast<float>(min_positive) : 0.F;
      const float floorOut = static_cast<float>(floor_fill);
      meta->Fill(static_cast<float>(m_visualization_dump_index),
                 static_cast<float>(side),
                 static_cast<float>(layernum),
                 static_cast<float>(x_center),
                 static_cast<float>(y_center),
                 static_cast<float>(phi),
                 static_cast<float>(rad_gem),
                 static_cast<float>(cloud_sig_rp),
                 floorOut,
                 minBinOut,
                 static_cast<float>(max_bin),
                 static_cast<float>(total_mass),
                 static_cast<float>(nonzero_bins));
      meta->Write("", TObject::kOverwrite);
      dumpFile.Close();
      ++m_visualization_dump_index;
    }
  }

  canvas->Update();
  canvas->SaveAs(m_visualization_output.c_str());

  std::string canvas_output = m_visualization_output;
  const std::string::size_type dot_pos = canvas_output.find_last_of('.');
  if (dot_pos != std::string::npos)
  {
    canvas_output.replace(dot_pos, std::string::npos, ".root");
  }
  else
  {
    canvas_output += ".root";
  }
  TFile canvasFile(canvas_output.c_str(), "RECREATE");
  if (canvasFile.IsOpen())
  {
    canvasFile.cd();
    canvas->Write();
    canvasFile.Close();
  }
  else
  {
    std::cout << "PHG4TpcPadPlaneReadout: unable to write canvas ROOT file "
              << canvas_output << std::endl;
  }

  std::cout << "PHG4TpcPadPlaneReadout: saved single cloud visualization to "
            << m_visualization_output
            << " (side " << side << ", layer " << layernum
            << ", phi " << phi << ", radius " << rad_gem << ")" << std::endl;

  m_visualization_done = !m_visualize_all_matches;

  for (TGraph *c : center_markers) delete c;
  for (TGraph *g : pad_graphs) delete g;
  delete electron_label_marker;
  delete canvas;
  delete hist;
}

void PHG4TpcPadPlaneReadout::SERF_zigzag_phibins(const unsigned int side, const unsigned int layernum,  const double phi, const double rad_gem,const double cloud_sig_rp, std::vector<int> &phibin_pad, std::vector<double> &phibin_pad_share)
{
  const double radius = LayerGeom->get_radius();
  const double phistepsize = LayerGeom->get_phistep();
  const auto phibins = LayerGeom->get_phibins();

  int sector1 = 0;
  double x = rad_gem * cos(phi);
  double y = rad_gem * sin(phi);
  double xNew = 0.0;
  double yNew = 0.0;

  rotatePointToSector( x,  y,  side, sector1, xNew, yNew);
  //  double phiNew = std::atan2(yNew,xNew);
  //double rad_gem_new = std::sqrt(xNew*xNew+ yNew*yNew);
 // std::cout<<"PHG4TpcPadPlaneReadout::SERF_zigzag_phibins: x = "<<x<<", y = "<<y<<", xNew = "<<xNew<<", yNew = "<<yNew<<" phiNew = "<<phiNew<<" rad_gem_new = "<<rad_gem_new<<std::endl;
  int tpc_module = (int)(layernum - 7)/16;

  const bool matches_target =
      (m_visualization_target_layer < 0 || static_cast<int>(layernum) == m_visualization_target_layer) &&
      (m_visualization_target_side < 0 || static_cast<int>(side) == m_visualization_target_side);
  const bool capture_debug = m_visualize_single_cloud &&
                             matches_target &&
                             (!m_visualization_done || m_visualize_all_matches);
  std::vector<DebugPadContribution> debug_contribs;
  std::vector<DebugSample> debug_samples_storage;
  std::vector<DebugSample>* debug_samples = nullptr;
  const bool need_samples = capture_debug;
  if (need_samples)
  {
    debug_samples_storage.reserve(7000);
    debug_samples = &debug_samples_storage;
  }
  if (capture_debug && Verbosity() > 0)
  {
    std::cout << "PHG4TpcPadPlaneReadout::SERF debug cloud: side=" << side
              << " layer=" << layernum
              << " sector=" << sector1
              << " phi=" << phi
              << " r=" << rad_gem
              << " xy=(" << x << "," << y << ")"
              << " rotated=(" << xNew << "," << yNew << ")"
              << std::endl;
  }
/* int phi_bin = LayerGeom->get_phibin(phi, side);
 int sector = 0;
 for (int i=0;i<12;i++)
  {
    if (phi >= sector_min_Phi[side][i] && phi < sector_max_Phi[side][i])
      {
        sector = i;
        break;
      }
  }
*/
//std::cout<<"!!!!!! phi = "<<phi<<" radius = "<<rad_gem<<" layer radius = "<<radius<<" sector = "<<sector<<" layer = "<<layernum<<" tpc_module = "<< tpc_module <<" phibins = "<<phibins<<" layer in module "<<((int)(layernum - 7)%16)<<" phi_bin = "<<phi_bin<<" cx = "<<Pads[layernum][ntpc_phibins_sector[tpc_module] -1 - phi_bin].cx<<" cy = "<<Pads[layernum][ntpc_phibins_sector[tpc_module] -1 - phi_bin].cy<<" Sector phi min "<<sector_min_Phi[side][sector]<<" sector phi max = "<<sector_max_Phi[side][sector]<<std::endl;

  // make the charge distribution gaussian
  double rphi = phi * radius;
  if (Verbosity() > 100)
  {
    if (LayerGeom->get_layer() == print_layer)
    {
      std::cout << " !!!!!!! SERF for layer " << layernum << " with radius " << radius << " phi " << phi
                << " rphi " << rphi << " phistepsize " << phistepsize << std::endl;
      std::cout << " fcharge created: radius " << radius << " rphi " << rphi << " cloud_sig_rp " << cloud_sig_rp << std::endl;
    }
  }

  const double philim_low_calc = phi - (_nsigmas * cloud_sig_rp / rad_gem) - phistepsize;
  const double philim_high_calc = phi + (_nsigmas * cloud_sig_rp / rad_gem) + phistepsize;

  // Find the pad range that covers this phi range
  const double philim_low = check_phi(side, philim_low_calc, rad_gem);
  const double philim_high = check_phi(side, philim_high_calc, rad_gem);

//std::cout<<"   SERF    zigzags: phi " << phi << " philim_low_calc " << philim_low_calc << " philim_low " << philim_low
              //  << " philim_high_calc " << philim_high_calc << " philim_high " << philim_high << " phiNew = "<<phiNew<<" rad_gem_new = "<<rad_gem_new<< std::endl;
  int phibin_low = LayerGeom->get_phibin(philim_high, side);
  int phibin_high = LayerGeom->get_phibin(philim_low, side);
  int npads = phibin_high - phibin_low;
/*
  const double radlim_low_calc = rad_gem - (_nsigmas * cloud_sig_rp / radius) - phistepsize;
  const double radlim_high_calc = rad_gem + (_nsigmas * cloud_sig_rp / radius) + phistepsize;

  // Find the pad range that covers this phi range
  const double radlim_low = check_phi(side, phi, radlim_low_calc);
  const double radlim_high = check_phi(side, phi, radlim_high_calc);
  if (radlim_low < (LayerGeom->get_radius() - LayerGeom->get_thickness() / 2.0))
  {
    std::cout << " radlim_low " << radlim_low << " is in the previous layer " <<layernum<<" with upper radius "<< LayerGeom->get_radius() - LayerGeom->get_thickness() / 2.0 << std::endl;
  }
  if (radlim_high > (LayerGeom->get_radius() + LayerGeom->get_thickness() / 2.0))
  {
    std::cout << " radlim_high " << radlim_high << " is in the next layer " <<layernum<<" with lower radius "<< LayerGeom->get_radius() + LayerGeom->get_thickness() / 2.0 << std::endl;
}
*/
    //  std::cout << "   SERF    zigzags: phi " << phi << " philim_low " << philim_low << " phibin_low " << phibin_low << " philim_high " << philim_high << " phibin_high " << phibin_high << " npads " << npads << std::endl;
  if (capture_debug)
  {
    const int reserve_count = (npads >= 0) ? (npads + 2) : 2;
    debug_contribs.reserve(static_cast<size_t>(reserve_count));
  }

  for (int ipad = 0; ipad <= npads; ipad++)
  {
    int pad_now = phibin_low + ipad;
    // if(phibin_low<0 && phibin_high<0) pad_now = phibin_high + ipad;
    //  check that we do not exceed the maximum number of pads, wrap if necessary
    if (pad_now >= phibins)
    {
      pad_now -= phibins;
    }

    int look_pad = pad_now;
    const int nlocal = ntpc_phibins_sector[tpc_module];
    // Convert global pad index to sector-local index.
    look_pad = ((look_pad % nlocal) + nlocal) % nlocal;
    // Side 0 uses the same local ordering as polygon pad indices.
    // Side 1 keeps the historical reversed lookup.
    if (side == 1)
    {
      look_pad = nlocal - look_pad - 1;
    }

   // std::cout<<"pad now = "<<pad_now<<" look_pad = "<<look_pad<<" ntpc_phibins_sector[tpc_module] = "<<ntpc_phibins_sector[tpc_module];

   //std::cout<<"layernum = "<<layernum<<" n = ntpc_phibins_sector[tpc_module] = "<<ntpc_phibins_sector[tpc_module]<<" pad now "<<pad_now<<"  look up pad = "<<look_pad<<" Pads[layernum] "<<Pads[layernum].size()<<std::endl;


    //std::cout<<"   SERF    zigzags: ipad " << ipad << " pad_now " << pad_now << " phibin_low " << phibin_low
         //     << " phibin_high " << phibin_high << " npads " << npads << "ntpc_phibins_sector[tpc_module] = "<<ntpc_phibins_sector[tpc_module]<<" pad look "<<ntpc_phibins_sector[tpc_module] - (pad_now - ntpc_phibins_sector[tpc_module]*sector) << std::endl;
    // Build polygon: either use loaded zigzag vertices (default) or synthetic rectangular pad
    std::vector<Point> poly;
    if (m_use_rectangular_pad_response)
    {
      // rectangular pad approximated by radial and phi bounds
      const double pad_phi_center = LayerGeom->get_phicenter(pad_now, side);
      const double half_phi = 0.5 * LayerGeom->get_phistep();
      const double r_center = LayerGeom->get_radius();
      const double half_thickness = 0.5 * LayerGeom->get_thickness();
      const double r_low = r_center - half_thickness;
      const double r_high = r_center + half_thickness;

      const double phi_edges[2] = {pad_phi_center - half_phi, pad_phi_center + half_phi};
      const double r_edges[2] = {r_low, r_high};

      poly.reserve(4);
      // ordered rectangle: (r_low,phi_low) -> (r_low,phi_high) -> (r_high,phi_high) -> (r_high,phi_low)
      poly.push_back({r_edges[0] * std::cos(phi_edges[0]), r_edges[0] * std::sin(phi_edges[0])});
      poly.push_back({r_edges[0] * std::cos(phi_edges[1]), r_edges[0] * std::sin(phi_edges[1])});
      poly.push_back({r_edges[1] * std::cos(phi_edges[1]), r_edges[1] * std::sin(phi_edges[1])});
      poly.push_back({r_edges[1] * std::cos(phi_edges[0]), r_edges[1] * std::sin(phi_edges[0])});
    }
    else
    {
      // Guard against missing pad polygons
      if (layernum >= Pads.size() || look_pad < 0 || static_cast<size_t>(look_pad) >= Pads[layernum].size())
      {
        // No polygon data for this pad/layer; skip contribution
        continue;
      }
      auto  padinfo = Pads[layernum][look_pad];
      poly = padinfo.vertices;
    }
  // std::cout<<"Calculate charge for pad with cx = "<<padinfo.cx<<", cy = "<<padinfo.cy<<" phi from Pads = "<<padinfo.phi<<", phi center(get pad now) = "
 //  <<LayerGeom->get_phicenter(pad_now, side)<<", phi center(get look_pad) = "<<LayerGeom->get_phicenter(look_pad, side)<< "sigma = "<<cloud_sig_rp<< "sigma/r = "<<cloud_sig_rp/rad_gem <<" pad look "<<look_pad<<" number of vert = "<<padinfo.vertices.size()<<std::endl;
  //std::cout<<"n = ntpc_phibins_sector[tpc_module]-1 = "<<n<<"  look up pad = "<<look_pad<<" Pads[layernum] "<<Pads[layernum].size()<<" xNew = "<<xNew<<" yNew = "<<yNew<<std::endl;

    // pick coordinate frame consistent with polygon choice
    const double hit_x = m_use_rectangular_pad_response ? x : xNew;
    const double hit_y = m_use_rectangular_pad_response ? y : yNew;
    double charge = integratedDensityOfCircleAndPad(hit_x, hit_y, cloud_sig_rp , poly, 0.0, debug_samples);
    if (capture_debug)
    {
      DebugPadContribution dbg;
      dbg.pad_bin = pad_now;
      dbg.charge = charge;
      dbg.pad_phi = LayerGeom->get_phicenter(pad_now, side);
      dbg.pad_r = LayerGeom->get_radius();
      dbg.polygon = poly;
      debug_contribs.push_back(std::move(dbg));
    }

    /* std::cout<<"   SERF    zigzags: ipad " << ipad << " pad_now " << pad_now << " charge " << charge<<" look_pad "<<look_pad<<" ntpc_phibins_sector[tpc_module] "<<ntpc_phibins_sector[tpc_module]
              << " pad cx " << padinfo.cx << " cy " << padinfo.cy
              << " phi center " << LayerGeom->get_phicenter(pad_now, side) <<" pad phi center  "<<padinfo.phi<< std::endl; */
    phibin_pad.push_back(pad_now);
    phibin_pad_share.push_back(charge);
 /*   std::cout<<"   SERF    zigzags: ipad " << ipad << " pad_now " << pad_now << " charge " << charge<<" look_pad "<<look_pad<<" ntpc_phibins_sector[tpc_module] "<<ntpc_phibins_sector[tpc_module]
              << " pad cx " << padinfo.cx << " cy " << padinfo.cy
              << " phi center " << LayerGeom->get_phicenter(pad_now, side) <<" pad phi center  "<<padinfo.phi<< std::endl;
*/
  }

  if (capture_debug && debug_samples && !debug_contribs.empty())
  {
    VisualizationCircle this_circle{xNew, yNew, _nsigmas * cloud_sig_rp};
    if (m_visualize_all_matches)
    {
      m_visualization_aggregate_samples.insert(m_visualization_aggregate_samples.end(),
                                               debug_samples->begin(), debug_samples->end());
      m_visualization_circles.push_back(this_circle);
      for (const auto &pad : debug_contribs) m_visualization_pad_union[pad.pad_bin] = pad;

      const double vis_x = m_use_rectangular_pad_response ? x : xNew;
      const double vis_y = m_use_rectangular_pad_response ? y : yNew;
      maybeVisualizeAvalanche(side, layernum, phi, rad_gem, cloud_sig_rp, vis_x, vis_y, debug_contribs,
                              m_visualization_aggregate_samples, m_visualization_circles);
    }
    else
    {
      std::vector<VisualizationCircle> single_circle{this_circle};
      const double vis_x = m_use_rectangular_pad_response ? x : xNew;
      const double vis_y = m_use_rectangular_pad_response ? y : yNew;
      maybeVisualizeAvalanche(side, layernum, phi, rad_gem, cloud_sig_rp, vis_x, vis_y, debug_contribs,
                              *debug_samples, single_circle);
      m_visualization_done = true;
    }
    ++m_visualization_cloud_counter;
  }

  return;
}


void PHG4TpcPadPlaneReadout::populate_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &phibin_pad, std::vector<double> &phibin_pad_share)
{
  const double radius = LayerGeom->get_radius();
  const double phistepsize = LayerGeom->get_phistep();
  const auto phibins = LayerGeom->get_phibins();
  // make the charge distribution gaussian
  double rphi = phi * radius;
  if (Verbosity() > 100)
  {
    if (LayerGeom->get_layer() == print_layer)
    {
      std::cout << " populate_zigzag_phibins for layer " << layernum << " with radius " << radius << " phi " << phi
                << " rphi " << rphi << " phistepsize " << phistepsize << std::endl;
      std::cout << " fcharge created: radius " << radius << " rphi " << rphi << " cloud_sig_rp " << cloud_sig_rp << std::endl;
    }
  }

  // Get the range of phi values that completely contains all pads  that touch the charge distribution - (nsigmas + 1/2 pad width) in each direction
  const double philim_low_calc = phi - (_nsigmas * cloud_sig_rp / radius) - phistepsize;
  const double philim_high_calc = phi + (_nsigmas * cloud_sig_rp / radius) + phistepsize;

  // Find the pad range that covers this phi range
  const double philim_low = check_phi(side, philim_low_calc, radius);
  const double philim_high = check_phi(side, philim_high_calc, radius);

  int phibin_low = LayerGeom->get_phibin(philim_high, side);
  int phibin_high = LayerGeom->get_phibin(philim_low, side);
  int npads = phibin_high - phibin_low;
 int sector = 0;
 for (int i=0;i<12;i++)
  {
    if (phi >= sector_min_Phi[side][i] && phi < sector_max_Phi[side][i])
      {
        sector = i;
        break;
      }
  }

  if (Verbosity() > 1000)
  {
    if (layernum == print_layer)
    {
      std::cout << "    REF       zigzags: phi " << phi << " philim_low " << philim_low << " phibin_low " << phibin_low
                << " philim_high " << philim_high << " phibin_high " << phibin_high << " npads " << npads <<" sector "<<sector<< std::endl;
    }
  }

  if (npads < 0 || npads > 9)
  {
    npads = 9;  // can happen if phibin_high wraps around. If so, limit to 10 pads and fix below
  }

  // Calculate the maximum extent in r-phi of pads in this layer. Pads are assumed to touch the center of the next phi bin on both sides.
  const double pad_rphi = 2.0 * LayerGeom->get_phistep() * radius;

  // Make a TF1 for each pad in the phi range
  using PadParameterSet = std::array<double, 2>;
  std::array<PadParameterSet, 10> pad_parameters{};
  std::array<int, 10> pad_keep{};

  // Now make a loop that steps through the charge distribution and evaluates the response at that point on each pad
  std::array<double, 10> overlap = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  double pads_phi[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double sum_of_pads_phi = 0.;
  double sum_of_pads_absphi = 0.;
  for (int ipad = 0; ipad <= npads; ipad++)
  {
    int pad_now = phibin_low + ipad;
    if (pad_now >= phibins)
    {
      pad_now -= phibins;
    }
    pads_phi[ipad] = LayerGeom->get_phicenter(pad_now, side);
    sum_of_pads_phi += pads_phi[ipad];
    sum_of_pads_absphi += std::fabs(pads_phi[ipad]);
  }

  for (int ipad = 0; ipad <= npads; ipad++)
  {
    int pad_now = phibin_low + ipad;
    // if(phibin_low<0 && phibin_high<0) pad_now = phibin_high + ipad;
    //  check that we do not exceed the maximum number of pads, wrap if necessary
    if (pad_now >= phibins)
    {
      pad_now -= phibins;
    }

    pad_keep[ipad] = pad_now;
    double phi_now = pads_phi[ipad];
    const double rphi_pad_now = phi_now * radius;
    pad_parameters[ipad] = {{pad_rphi / 2.0, rphi_pad_now}};

    if (Verbosity() > 1000)
    {
      if (layernum == print_layer)
      {
        std::cout << " zigzags: make fpad for ipad " << ipad << " pad_now " << pad_now << " pad_rphi/2 " << pad_rphi / 2.0
                  << " rphi_pad_now " << rphi_pad_now << std::endl;
      }
    }
    //}

    // use analytic integral
    // for (int ipad = 0; ipad <= npads; ipad++)
    //{
    // const double pitch = pad_parameters[ipad][0];
    // const double x_loc = pad_parameters[ipad][1] - rphi;
    // const double sigma = cloud_sig_rp;

    const double pitch = pad_rphi / 2.0;     // eshulga
    double x_loc_tmp = rphi_pad_now - rphi;  // eshulga
    const double sigma = cloud_sig_rp;       // eshulga

    // Checking if the pads are on the same side of the TPC in phi
    if (std::fabs(sum_of_pads_phi) != sum_of_pads_absphi)
    {
      if (phi < -M_PI / 2 && phi_now > 0)
      {
        x_loc_tmp = (phi_now - 2 * M_PI) * radius - rphi;
      }
      if (phi > M_PI / 2 && phi_now < 0)
      {
        x_loc_tmp = (phi_now + 2 * M_PI) * radius - rphi;
      }
      if (phi < 0 && phi_now > 0)
      {
        x_loc_tmp = (phi_now + std::fabs(phi)) * radius;
      }
      if (phi > 0 && phi_now < 0)
      {
        x_loc_tmp = (2 * M_PI - phi_now + phi) * radius;
      }
    }

    const double x_loc = x_loc_tmp;
    // calculate fraction of the total charge on this pad
    if (m_use_rectangular_pad_response)
    {
      // flat response across the pad width (rectangular pad)
      const double half_width = pad_rphi / 2.0;
      const double upper = (x_loc + half_width) / (M_SQRT2 * sigma);
      const double lower = (x_loc - half_width) / (M_SQRT2 * sigma);
      overlap[ipad] = 0.5 * (std::erf(upper) - std::erf(lower));
    }
    else
    {
      /*
      this corresponds to integrating the charge distribution Gaussian function (centered on rphi and of width cloud_sig_rp),
      convoluted with a strip response function, which is triangular from -pitch to +pitch, with a maximum of 1. at strip center
      */
      overlap[ipad] =
          (pitch - x_loc) * (std::erf(x_loc / (M_SQRT2 * sigma)) - std::erf((x_loc - pitch) / (M_SQRT2 * sigma))) / (pitch * 2) + (pitch + x_loc) * (std::erf((x_loc + pitch) / (M_SQRT2 * sigma)) - std::erf(x_loc / (M_SQRT2 * sigma))) / (pitch * 2) + (gaus(x_loc - pitch, sigma) - gaus(x_loc, sigma)) * square(sigma) / pitch + (gaus(x_loc + pitch, sigma) - gaus(x_loc, sigma)) * square(sigma) / pitch;
    }
  }

  // now we have the overlap for each pad
  for (int ipad = 0; ipad <= npads; ipad++)
  {
    phibin_pad.push_back(pad_keep[ipad]);
    phibin_pad_share.push_back(overlap[ipad]);
        std::cout<<"   REF    zigzags: ipad " << ipad << " pad_now " <<  pad_keep[ipad] << " charge " << overlap[ipad]
              << " phi center " << LayerGeom->get_phicenter(pad_keep[ipad], side) << std::endl;
  }

  return;
}



void PHG4TpcPadPlaneReadout::UseGain(const int flagToUseGain)
{
  m_flagToUseGain = flagToUseGain;
  if (m_flagToUseGain == 1 && Verbosity() > 0)
  {
    std::cout << "PHG4TpcPadPlaneReadout: UseGain: TRUE " << std::endl;
  }
}

void PHG4TpcPadPlaneReadout::ReadGain()
{
  // Reading TPC Gain Maps from the file
  if (m_flagToUseGain == 1)
  {
    char *calibrationsroot = getenv("CALIBRATIONROOT");
    if (calibrationsroot != nullptr)
    {
      std::string gain_maps_filename = std::string(calibrationsroot) + std::string("/TPC/GainMaps/TPCGainMaps.root");
      TFile *fileGain = TFile::Open(gain_maps_filename.c_str(), "READ");
      h_gain[0] = (TH2 *) fileGain->Get("RadPhiPlot0")->Clone();
      h_gain[1] = (TH2 *) fileGain->Get("RadPhiPlot1")->Clone();
      h_gain[0]->SetDirectory(nullptr);
      h_gain[1]->SetDirectory(nullptr);
      fileGain->Close();
    }
  }
}
void PHG4TpcPadPlaneReadout::SetDefaultParameters()
{
  set_default_double_param("tpc_minradius_inner", 31.105);  // 30.0);  // cm
  set_default_double_param("tpc_minradius_mid", 41.153);    // 40.0);
  set_default_double_param("tpc_minradius_outer", 58.367);  // 60.0);

  set_default_double_param("tpc_maxradius_inner", 40.249);  // 40.0);  // cm
  set_default_double_param("tpc_maxradius_mid", 57.475);    // 60.0);
  set_default_double_param("tpc_maxradius_outer", 75.911);  // 77.0);  // from Tom
  set_default_double_param("tpc_sampa_peaking_time", 80.0); // ns

  // Minimum effective electrons per (pad,tbin) to create a hit
  // Set to 0.0 by default to include all contributions
  set_default_double_param("neffelectrons_threshold", 0.0);
  set_default_double_param("gem_cloud_sigma", 0.04);     // cm = 400 microns
  set_default_double_param("sampa_shaping_lead", 32.0);  // ns, for 80 ns SAMPA
  set_default_double_param("sampa_shaping_tail", 48.0);  // ns, for 80 ns SAMPA

  set_default_double_param("tpc_sector_phi_inner", 0.5024);  // 2 * M_PI / 12 );//sector size in phi for R1 sector
  set_default_double_param("tpc_sector_phi_mid", 0.5087);    // 2 * M_PI / 12 );//sector size in phi for R2 sector
  set_default_double_param("tpc_sector_phi_outer", 0.5097);  // 2 * M_PI / 12 );//sector size in phi for R3 sector

  set_default_int_param("ntpc_phibins_inner", 1128);  // 94 * 12
  set_default_int_param("ntpc_phibins_mid", 1536);    // 128 * 12
  set_default_int_param("ntpc_phibins_outer", 2304);  // 192 * 12

  // GEM Gain
  /*
  hp (2020/09/04): gain changed from 2000 to 1400, to accomodate gas mixture change
  from Ne/CF4 90/10 to Ne/CF4 50/50, and keep the average charge per particle per pad constant
  */
  set_default_double_param("gem_amplification", 1400);
  set_default_double_param("polya_theta", 0.8);
  return;
}

void PHG4TpcPadPlaneReadout::UpdateInternalParameters()
{
  neffelectrons_threshold = get_double_param("neffelectrons_threshold");

  MinRadius =
      {{get_double_param("tpc_minradius_inner"),
        get_double_param("tpc_minradius_mid"),
        get_double_param("tpc_minradius_outer")}};

  MaxRadius =
      {{get_double_param("tpc_maxradius_inner"),
        get_double_param("tpc_maxradius_mid"),
        get_double_param("tpc_maxradius_outer")}};

  sigmaT = get_double_param("gem_cloud_sigma");
  sigmaL = {{get_double_param("sampa_shaping_lead"),
             get_double_param("sampa_shaping_tail")}};

  averageGEMGain = get_double_param("gem_amplification");
  polyaTheta = get_double_param("polya_theta");
  Ts = get_double_param("tpc_sampa_peaking_time");

}
void PHG4TpcPadPlaneReadout::makeChannelMask(hitMaskTpc &aMask, const std::string &dbName, const std::string &totalChannelsToMask)
{
  CDBTTree *cdbttree;
  if (m_maskFromFile)
  {
    cdbttree = new CDBTTree(dbName);
  }
  else // mask using CDB TTree, default
  {
    std::string database = CDBInterface::instance()->getUrl(dbName);
    cdbttree = new CDBTTree(database);
  }

  std::cout << "Masking TPC Channel Map: " << dbName << std::endl;

  int NChan = -1;
  NChan = cdbttree->GetSingleIntValue(totalChannelsToMask);

  for (int i = 0; i < NChan; i++)
  {
    int Layer = cdbttree->GetIntValue(i, "layer");
    int Sector = cdbttree->GetIntValue(i, "sector");
    int Side = cdbttree->GetIntValue(i, "side");
    int Pad = cdbttree->GetIntValue(i, "pad");
    if (Verbosity() > VERBOSITY_A_LOT)
    {
      std::cout << dbName << ": Will mask layer: " << Layer << ", sector: " << Sector << ", side: " << Side << ", Pad: " << Pad << std::endl;
    }

    TrkrDefs::hitsetkey DeadChannelHitKey = TpcDefs::genHitSetKey(Layer, Sector, Side);
    TrkrDefs::hitkey DeadHitKey = TpcDefs::genHitKey((unsigned int) Pad, 0);
    aMask[DeadChannelHitKey].push_back(DeadHitKey);
  }

  delete cdbttree;

}
// -------------------------------------------------------------------------
// REPLACEMENT FUNCTIONS (From Code A)
// -------------------------------------------------------------------------

void PHG4TpcPadPlaneReadout::sampaTimeDistribution(double tzero, std::vector<int> &adc_tbin, std::vector<double> &adc_tbin_share)
{
  // tzero is the arrival time of the electron at the GEM
  // Ts is the sampa peaking time
  // Assume the response is over after 8 clock cycles (400 ns)
  int nclocks = 8;

  double tstepsize = LayerGeom->get_zstep();
  int tbinzero = LayerGeom->get_zbin(tzero);

  // the first clock bin is a special case
  double tfirst_end = LayerGeom->get_zcenter(tbinzero) + tstepsize/2.0;
  double vfirst_end = sampaShapingResponseFunction(tzero, tfirst_end);
  double first_integral = (vfirst_end / 2.0) * (tfirst_end - tzero);

  adc_tbin.push_back(tbinzero);
  adc_tbin_share.push_back(first_integral);

  for(int iclock = 1; iclock < nclocks; ++iclock)
  {
    int tbin = tbinzero + iclock;
    if (tbin < 0 || tbin > LayerGeom->get_zbins())
    {
      if (Verbosity() > 0)
	    {
	      std::cout << " t bin " << tbin << " is outside range of " << LayerGeom->get_zbins() << " so skip it" << std::endl;
	    }
      continue;
    }

    // get the beginning and end of this clock bin
    double tcenter = LayerGeom->get_zcenter(tbin);
    double tlow = tcenter - tstepsize/2.0;

    // sample the voltage in this bin at nsamples locations
    int nsamples = 6;
    double sample_step = tstepsize / (double) nsamples;
    double sintegral = 0;
    for(int isample = 0; isample < nsamples; ++isample)
    {
      double tnow = tlow + (double) isample * sample_step + sample_step / 2.0;
      double vnow = sampaShapingResponseFunction(tzero, tnow);
      sintegral += vnow * sample_step;
    }

    adc_tbin.push_back(tbin);
    adc_tbin_share.push_back(sintegral);
  }
}

double PHG4TpcPadPlaneReadout::sampaShapingResponseFunction(double tzero, double t) const
{
  // The specific SAMPA response function: V ~ (t/tau)^4 * exp(-4t/tau)
  //if (t < tzero) return 0.0;
  double v = exp(-4.0 * (t - tzero) / Ts) * pow((t - tzero) / Ts, 4.0);
  return v;
}
