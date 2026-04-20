/*!
 * \file TpcDnlCorrection.cc
 * \brief applies per-cluster TPC DNL phi corrections from a calibration ROOT file
 */

#include "TpcDnlCorrection.h"

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <Acts/Definitions/Algebra.hpp>

#include <TFile.h>
#include <TF1.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <memory>

namespace
{
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kPhaseTolerance = 1e-9;
}  // namespace

//____________________________________________________________________________..
TpcDnlCorrection::TpcDnlCorrection(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
double TpcDnlCorrection::wrap_phi(double phi)
{
  while (phi <= -kPi) phi += 2.0 * kPi;
  while (phi > kPi) phi -= 2.0 * kPi;
  return phi;
}

//____________________________________________________________________________..
bool TpcDnlCorrection::load_corrections()
{
  m_corrections.clear();

  if (m_correction_filename.empty())
  {
    std::cout << Name() << ": no DNL correction filename configured." << std::endl;
    return false;
  }

  std::unique_ptr<TFile> inputfile(TFile::Open(m_correction_filename.c_str(), "READ"));
  if (!inputfile || inputfile->IsZombie())
  {
    std::cout << Name() << ": cannot open DNL correction file " << m_correction_filename << std::endl;
    return false;
  }

  auto* correction_tree = dynamic_cast<TTree*>(inputfile->Get("dnl_corrections"));
  if (!correction_tree)
  {
    std::cout << Name() << ": missing dnl_corrections tree in " << m_correction_filename << std::endl;
    return false;
  }

  int cal_layer = -1;
  int cal_side = -1;
  int cal_npads = -1;
  int cal_fit_status = -999;
  double cal_fit_phase_min = std::numeric_limits<double>::quiet_NaN();
  double cal_fit_phase_max = std::numeric_limits<double>::quiet_NaN();
  char cal_func_name[256] = "";
  char cal_group_name[64] = "";
  char cal_fit_region[32] = "";

  correction_tree->SetBranchAddress("layer", &cal_layer);
  correction_tree->SetBranchAddress("side", &cal_side);
  correction_tree->SetBranchAddress("npads", &cal_npads);
  correction_tree->SetBranchAddress("fit_status", &cal_fit_status);
  correction_tree->SetBranchAddress("fit_phase_min", &cal_fit_phase_min);
  correction_tree->SetBranchAddress("fit_phase_max", &cal_fit_phase_max);
  correction_tree->SetBranchAddress("func_name", cal_func_name);
  correction_tree->SetBranchAddress("group_name", cal_group_name);
  correction_tree->SetBranchAddress("fit_region", cal_fit_region);

  int loaded_functions = 0;
  int missing_functions = 0;

  const auto nentries = correction_tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    correction_tree->GetEntry(i);
    if (cal_fit_status != 0) continue;
    if (std::strlen(cal_func_name) == 0) continue;

    const TString func_path = TString::Format("correction_functions/%s/%s", cal_group_name, cal_func_name);
    auto* source_function = dynamic_cast<TF1*>(inputfile->Get(func_path));
    if (!source_function)
    {
      ++missing_functions;
      continue;
    }

    CorrectionEntry entry;
    entry.source_side = cal_side;
    entry.phase_min = cal_fit_phase_min;
    entry.phase_max = cal_fit_phase_max;
    entry.fit_region = cal_fit_region;
    entry.func_name = cal_func_name;
    entry.function.reset(static_cast<TF1*>(source_function->Clone()));

    m_corrections[CorrectionKey(static_cast<unsigned int>(cal_layer), cal_side, cal_npads)].push_back(std::move(entry));
    ++loaded_functions;
  }

  for (auto& [key, entries] : m_corrections)
  {
    std::sort(entries.begin(), entries.end(),
        [](const CorrectionEntry& lhs, const CorrectionEntry& rhs)
        { return (lhs.phase_max - lhs.phase_min) < (rhs.phase_max - rhs.phase_min); });
  }

  std::cout << Name() << ": loaded " << loaded_functions << " DNL correction functions";
  if (missing_functions > 0)
  {
    std::cout << " (" << missing_functions << " missing from calibration file)";
  }
  std::cout << std::endl;

  if (loaded_functions == 0)
  {
    std::cout << Name() << ": no usable DNL corrections were loaded from " << m_correction_filename << std::endl;
    return false;
  }

  return true;
}

//____________________________________________________________________________..
int TpcDnlCorrection::InitRun(PHCompositeNode* topNode)
{
  m_clusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  m_cluster_hit_assoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  m_hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  m_geom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");

  if (!m_clusters || !m_cluster_hit_assoc || !m_hitsets || !m_tGeometry || !m_geom)
  {
    std::cout << Name() << ": missing required nodes."
              << " TRKR_CLUSTER=" << static_cast<bool>(m_clusters)
              << " TRKR_CLUSTERHITASSOC=" << static_cast<bool>(m_cluster_hit_assoc)
              << " TRKR_HITSET=" << static_cast<bool>(m_hitsets)
              << " ActsGeometry=" << static_cast<bool>(m_tGeometry)
              << " TPCGEOMCONTAINER=" << static_cast<bool>(m_geom)
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!load_corrections())
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
TpcDnlCorrection::MatchResult TpcDnlCorrection::find_correction(
    unsigned int layer,
    int side,
    int npads,
    double phase_reco) const
{
  MatchResult result;

  auto it = m_corrections.find(CorrectionKey(layer, side, npads));
  if (it != m_corrections.end())
  {
    const auto* best = [&]()
    {
      const CorrectionEntry* match = nullptr;
      double best_width = std::numeric_limits<double>::infinity();
      for (const auto& entry : it->second)
      {
        if (!(phase_reco >= entry.phase_min - kPhaseTolerance &&
              phase_reco <= entry.phase_max + kPhaseTolerance))
        {
          continue;
        }

        const double width = entry.phase_max - entry.phase_min;
        if (!match || width < best_width)
        {
          match = &entry;
          best_width = width;
        }
      }
      return match;
    }();

    result.have_any_for_key = true;
    result.entry = best;
    if (result.entry) return result;
  }

  if (!m_allow_combined_side_fallback) return result;

  it = m_corrections.find(CorrectionKey(layer, -1, npads));
  if (it != m_corrections.end())
  {
    const auto* best = [&]()
    {
      const CorrectionEntry* match = nullptr;
      double best_width = std::numeric_limits<double>::infinity();
      for (const auto& entry : it->second)
      {
        if (!(phase_reco >= entry.phase_min - kPhaseTolerance &&
              phase_reco <= entry.phase_max + kPhaseTolerance))
        {
          continue;
        }

        const double width = entry.phase_max - entry.phase_min;
        if (!match || width < best_width)
        {
          match = &entry;
          best_width = width;
        }
      }
      return match;
    }();

    result.have_any_for_key = true;
    result.entry = best;
    result.used_combined_side = (result.entry != nullptr);
  }

  return result;
}

//____________________________________________________________________________..
TpcDnlCorrection::ClusterPhaseInputs TpcDnlCorrection::build_cluster_phase_inputs(
    TrkrDefs::cluskey ckey,
    TrkrDefs::hitsetkey hsk,
    TrkrCluster* cluster) const
{
  ClusterPhaseInputs result;
  if (!cluster || !m_cluster_hit_assoc || !m_hitsets || !m_tGeometry || !m_geom)
  {
    return result;
  }

  auto* hitset = m_hitsets->findHitSet(hsk);
  if (!hitset) return result;

  const auto layer = TrkrDefs::getLayer(ckey);
  const int side = static_cast<int>(TpcDefs::getSide(hsk));
  auto* layergeom = m_geom->GetLayerCellGeom(layer);
  if (!layergeom) return result;

  std::map<unsigned short, double> pad_adc_sums;
  const auto hit_range = m_cluster_hit_assoc->getHits(ckey);
  for (auto hitit = hit_range.first; hitit != hit_range.second; ++hitit)
  {
    const auto hitkey = static_cast<TrkrDefs::hitkey>(hitit->second);
    auto* hit = hitset->getHit(hitkey);
    if (!hit) continue;

    const double adc = static_cast<double>(hit->getAdc());
    if (!(adc > 0.0)) continue;

    const unsigned short iphi = TpcDefs::getPad(hitkey);
    pad_adc_sums[iphi] += adc;
  }

  if (pad_adc_sums.empty()) return result;

  unsigned short best_iphi = 0;
  double best_adc_sum = -std::numeric_limits<double>::infinity();
  for (const auto& [iphi, adc_sum] : pad_adc_sums)
  {
    if (adc_sum > best_adc_sum)
    {
      best_iphi = iphi;
      best_adc_sum = adc_sum;
    }
  }

  const Acts::Vector3 global = m_tGeometry->getGlobalPosition(ckey, cluster);
  if (!std::isfinite(global.x()) || !std::isfinite(global.y()) || !std::isfinite(global.z()))
  {
    return result;
  }

  const double r = std::hypot(global.x(), global.y());
  const double phi_width = std::abs(layergeom->get_phistep());
  if (!(r > 0.0) || !(phi_width > 0.0))
  {
    return result;
  }

  const double phi_reco = std::atan2(global.y(), global.x());
  const double phi_pad_max = layergeom->get_phicenter(static_cast<int>(best_iphi), side);
  const double phase_reco = wrap_phi(phi_reco - phi_pad_max) / phi_width;
  if (!std::isfinite(phase_reco))
  {
    return result;
  }

  result.valid = true;
  result.npads = static_cast<int>(pad_adc_sums.size());
  result.r = r;
  result.z = global.z();
  result.phi_reco = phi_reco;
  result.phi_pad_max = phi_pad_max;
  result.phi_width = phi_width;
  result.phase_reco = phase_reco;
  return result;
}

//____________________________________________________________________________..
int TpcDnlCorrection::process_event(PHCompositeNode* /*topNode*/)
{
  if (!m_clusters || !m_cluster_hit_assoc || !m_hitsets || !m_tGeometry || !m_geom)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (const auto& hsk : m_clusters->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto cluster_range = m_clusters->getClusters(hsk);
    for (auto clusit = cluster_range.first; clusit != cluster_range.second; ++clusit)
    {
      const auto ckey = clusit->first;
      auto* cluster = clusit->second;
      if (!cluster) continue;

      ++m_nclusters_total;

      const auto inputs = build_cluster_phase_inputs(ckey, hsk, cluster);
      if (!inputs.valid)
      {
        ++m_ninvalid_inputs;
        continue;
      }

      const auto layer = TrkrDefs::getLayer(ckey);
      const int side = static_cast<int>(TpcDefs::getSide(hsk));
      const auto match = find_correction(layer, side, inputs.npads, inputs.phase_reco);
      if (!match.entry)
      {
        if (match.have_any_for_key)
        {
          ++m_nphase_gap;
        }
        else
        {
          ++m_nmissing_key;
        }
        continue;
      }

      const double dRphi_correction = match.entry->function->Eval(inputs.phase_reco);
      if (!std::isfinite(dRphi_correction))
      {
        ++m_ninvalid_inputs;
        continue;
      }

      const double phi_corr = wrap_phi(inputs.phi_reco - dRphi_correction / inputs.r);
      const Acts::Vector3 corrected_global(
          inputs.r * std::cos(phi_corr),
          inputs.r * std::sin(phi_corr),
          inputs.z);

      auto surface = m_tGeometry->maps().getSurface(ckey, cluster);
      if (!surface)
      {
        ++m_nsurface_failures;
        continue;
      }

      Acts::Vector3 corrected_local =
          surface->transform(m_tGeometry->geometry().getGeoContext()).inverse() *
          (corrected_global * Acts::UnitConstants::cm);
      corrected_local /= Acts::UnitConstants::cm;

      if (!std::isfinite(corrected_local(0)))
      {
        ++m_nsurface_failures;
        continue;
      }

      cluster->setLocalX(corrected_local(0));
      ++m_ncorrected;
      if (match.used_combined_side) ++m_ncorrected_combined;

      if (Verbosity() > 3)
      {
        std::cout << Name()
                  << ": layer=" << layer
                  << " side=" << side
                  << " npads=" << inputs.npads
                  << " phase_reco=" << inputs.phase_reco
                  << " dRphi=" << dRphi_correction
                  << (match.used_combined_side ? " combined-side" : " exact-side")
                  << std::endl;
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcDnlCorrection::End(PHCompositeNode* /*topNode*/)
{
  std::cout << Name()
            << ": processed " << m_nclusters_total
            << " TPC clusters, corrected " << m_ncorrected
            << " (combined-side fallback " << m_ncorrected_combined
            << "), invalid inputs " << m_ninvalid_inputs
            << ", missing key " << m_nmissing_key
            << ", phase outside fit " << m_nphase_gap
            << ", surface failures " << m_nsurface_failures
            << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
