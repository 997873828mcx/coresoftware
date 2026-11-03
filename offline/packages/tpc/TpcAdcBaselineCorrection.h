#ifndef TPC_TPCADCBASELINECORRECTION_H
#define TPC_TPCADCBASELINECORRECTION_H

#include <pdbcalbase/PdbParameterMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

/** Shared event-level ADC normalization for reconstruction consumers.
 * Header-only to avoid a dependency from the standard TPC clusterizer on the
 * polytracking package. No truth information is consulted or changed.
 */
namespace TpcAdcBaseline
{
  struct Options
  {
    double declared_pedestal_adc {std::numeric_limits<double>::quiet_NaN()};
    bool input_already_subtracted {false};
    bool allow_missing_metadata {false};
    std::string correction_name {"TPC_ADC_BASELINE"};
  };

  struct Result
  {
    bool ok {false};
    bool has_metadata {false};
    bool corrected {false};
    bool already_subtracted {false};
    double pedestal_adc {std::numeric_limits<double>::quiet_NaN()};
    double effective_integer_shift_adc {std::numeric_limits<double>::quiet_NaN()};
    std::uint64_t hits_corrected {0};
    std::uint64_t hits_clamped {0};
    std::uint64_t hits_zero {0};
    std::string error;
  };

  namespace detail
  {
    inline bool getInt(const PdbParameterMap& map, const std::string& name, int& value)
    {
      const auto range = map.get_iparam_iters();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        if (iter->first == name)
        {
          value = iter->second;
          return true;
        }
      }
      return false;
    }

    inline bool getDouble(const PdbParameterMap& map, const std::string& name, double& value)
    {
      const auto range = map.get_dparam_iters();
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        if (iter->first == name)
        {
          value = iter->second;
          return true;
        }
      }
      return false;
    }

    inline bool samePedestal(const double a, const double b)
    {
      return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= 1.e-5;
    }
  }

  inline Result normalize(PHCompositeNode* topNode, const Options& options)
  {
    Result result;
    auto fail = [&result](const std::string& message)
    {
      result.error = message;
      return result;
    };
    if (!std::isnan(options.declared_pedestal_adc) &&
        (!std::isfinite(options.declared_pedestal_adc) ||
         options.declared_pedestal_adc < 0.0 || options.declared_pedestal_adc > 1023.0))
    {
      return fail("declared pedestal is outside [0,1023] ADC");
    }
    auto* baseline = findNode::getClass<PdbParameterMap>(topNode, "TPC_ADC_BASELINE");
    int schema = 0;
    int state = -1;
    double pedestal = std::numeric_limits<double>::quiet_NaN();
    double removed = std::numeric_limits<double>::quiet_NaN();
    const bool hasMetadata = baseline && detail::getInt(*baseline, "schema_version", schema);
    if (baseline && !hasMetadata &&
        (baseline->get_iparam_iters().first != baseline->get_iparam_iters().second ||
         baseline->get_dparam_iters().first != baseline->get_dparam_iters().second ||
         baseline->get_cparam_iters().first != baseline->get_cparam_iters().second))
    {
      return fail("nonempty TPC_ADC_BASELINE has no schema version");
    }
    if (hasMetadata)
    {
      if (schema != 1 || !detail::getInt(*baseline, "state", state) ||
          !detail::getDouble(*baseline, "pedestal_adc", pedestal) ||
          !detail::getDouble(*baseline, "removed_pedestal_adc", removed) ||
          !std::isfinite(pedestal) || pedestal < 0.0 || pedestal > 1023.0 ||
          (state != 0 && state != 1) ||
          !detail::samePedestal(removed, state == 0 ? 0.0 : pedestal))
      {
        return fail("invalid TPC_ADC_BASELINE event contract");
      }
      if ((std::isfinite(options.declared_pedestal_adc) &&
           !detail::samePedestal(pedestal, options.declared_pedestal_adc)) ||
          (options.input_already_subtracted && state != 1))
      {
        return fail("declared input baseline disagrees with event metadata");
      }
    }
    else
    {
      if (!options.input_already_subtracted && !std::isfinite(options.declared_pedestal_adc))
      {
        if (options.allow_missing_metadata)
        {
          // Real data and unmarked historical inputs retain their explicitly
          // configured consumer convention; never infer a pedestal from ADC.
          result.ok = true;
          return result;
        }
        return fail("missing TPC_ADC_BASELINE; explicitly declare legacy pedestal or already-subtracted input");
      }
      pedestal = std::isfinite(options.declared_pedestal_adc) ? options.declared_pedestal_adc : 0.0;
      state = options.input_already_subtracted ? 1 : 0;
      if (!baseline)
      {
        PHNodeIterator iter(topNode);
        auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
        if (!dstNode)
        {
          return fail("DST node is missing");
        }
        baseline = new PdbParameterMap;
        dstNode->addNode(new PHIODataNode<PHObject>(baseline, "TPC_ADC_BASELINE", "PHObject"));
      }
      baseline->Reset();
      baseline->set_int_param("schema_version", 1);
      baseline->set_double_param("pedestal_adc", pedestal);
      baseline->set_string_param("source", options.input_already_subtracted ? "declared_already_subtracted" : "declared_legacy_simulation");
    }
    result.has_metadata = true;
    result.pedestal_adc = pedestal;
    result.effective_integer_shift_adc = std::ceil(pedestal - 0.5);
    auto* hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    if (!hits)
    {
      return fail("TRKR_HITSET is missing");
    }
    if (state == 1)
    {
      result.already_subtracted = true;
      result.ok = true;
      baseline->set_int_param("state", 1);
      baseline->set_double_param("removed_pedestal_adc", pedestal);
      return result;
    }
    const auto hitsets = hits->getHitSets(TrkrDefs::TrkrId::tpcId);
    for (auto hitset = hitsets.first; hitset != hitsets.second; ++hitset)
    {
      const auto hitRange = hitset->second->getHits();
      for (auto hit = hitRange.first; hit != hitRange.second; ++hit)
      {
        const double corrected = static_cast<double>(hit->second->getAdc()) - pedestal;
        result.hits_clamped += corrected < 0.0;
        hit->second->setAdc(corrected > 0.0 ? static_cast<unsigned int>(std::floor(corrected + 0.5)) : 0U);
        result.hits_zero += hit->second->getAdc() == 0;
        ++result.hits_corrected;
      }
    }
    baseline->set_int_param("state", 1);
    baseline->set_double_param("removed_pedestal_adc", pedestal);
    baseline->set_double_param("effective_integer_shift_adc", result.effective_integer_shift_adc);
    baseline->set_double_param("hits_corrected", static_cast<double>(result.hits_corrected));
    baseline->set_double_param("hits_clamped", static_cast<double>(result.hits_clamped));
    baseline->set_double_param("hits_zero", static_cast<double>(result.hits_zero));
    baseline->set_string_param("quantization", "max(0,floor(raw_adc-pedestal_adc+0.5))");
    baseline->set_string_param("correction", options.correction_name);
    result.corrected = true;
    result.ok = true;
    return result;
  }

  inline Result normalizeIfMetadata(PHCompositeNode* topNode,
                                    const std::string& name = "automatic_reconstruction_consumer")
  {
    Options options;
    options.allow_missing_metadata = true;
    options.correction_name = name;
    return normalize(topNode, options);
  }
}

#endif
