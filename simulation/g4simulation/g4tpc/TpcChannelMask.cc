#include "TpcChannelMask.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <trackbase/TpcDefs.h>

#include <iostream>
#include <memory>

namespace
{
  int padsPerSector(const int layer)
  {
    if (layer >= 7 && layer <= 22)
    {
      return 94;
    }
    if (layer >= 23 && layer <= 38)
    {
      return 128;
    }
    if (layer >= 39 && layer <= 54)
    {
      return 192;
    }
    return 0;
  }
}  // namespace

void TpcChannelMask::clear()
{
  m_mask.clear();
  m_size = 0;
}

bool TpcChannelMask::add(const TrkrDefs::hitsetkey hitsetkey,
                         const unsigned int pad)
{
  const auto channel_key = TpcDefs::genHitKey(pad, 0);
  const bool inserted = m_mask[hitsetkey].insert(channel_key).second;
  if (inserted)
  {
    ++m_size;
  }
  return inserted;
}

bool TpcChannelMask::contains(const TrkrDefs::hitsetkey hitsetkey,
                              const unsigned int pad) const
{
  const auto hitset_iter = m_mask.find(hitsetkey);
  if (hitset_iter == m_mask.end())
  {
    return false;
  }

  const auto channel_key = TpcDefs::genHitKey(pad, 0);
  return hitset_iter->second.contains(channel_key);
}

bool TpcChannelMask::load(const std::string& map_name,
                          const std::string& total_channels_key,
                          const bool from_file,
                          const int verbosity)
{
  clear();

  if (map_name.empty())
  {
    std::cerr << "TpcChannelMask::load - empty map name" << std::endl;
    return false;
  }

  const std::string source =
      from_file ? map_name : CDBInterface::instance()->getUrl(map_name);
  if (source.empty())
  {
    std::cerr << "TpcChannelMask::load - could not resolve " << map_name
              << std::endl;
    return false;
  }

  auto cdbttree = std::make_unique<CDBTTree>(source);
  const int channel_count =
      cdbttree->GetSingleIntValue(total_channels_key);
  if (channel_count <= 0)
  {
    std::cerr << "TpcChannelMask::load - invalid or empty "
              << total_channels_key << " in " << source << std::endl;
    return false;
  }

  for (int index = 0; index < channel_count; ++index)
  {
    const int layer = cdbttree->GetIntValue(index, "layer");
    const int sector = cdbttree->GetIntValue(index, "sector");
    const int side = cdbttree->GetIntValue(index, "side");
    const int pad = cdbttree->GetIntValue(index, "pad");

    const int pads_per_sector = padsPerSector(layer);
    if (pads_per_sector == 0 ||
        sector < 0 || sector >= 12 ||
        side < 0 || side >= 2 ||
        pad < 0 || pad >= 12 * pads_per_sector ||
        pad / pads_per_sector != sector)
    {
      std::cerr << "TpcChannelMask::load - invalid channel at row "
                << index << " in " << source
                << ": layer=" << layer
                << " sector=" << sector
                << " side=" << side
                << " pad=" << pad << std::endl;
      clear();
      return false;
    }

    const auto hitsetkey = TpcDefs::genHitSetKey(layer, sector, side);
    if (!add(hitsetkey, static_cast<unsigned int>(pad)))
    {
      std::cerr << "TpcChannelMask::load - duplicate channel at row "
                << index << " in " << source << std::endl;
      clear();
      return false;
    }

    if (verbosity > 5)
    {
      std::cout << map_name << ": masking layer=" << layer
                << " sector=" << sector
                << " side=" << side
                << " pad=" << pad << std::endl;
    }
  }

  if (m_size != static_cast<std::size_t>(channel_count))
  {
    std::cerr << "TpcChannelMask::load - expected " << channel_count
              << " unique channels but loaded " << m_size
              << " from " << source << std::endl;
    clear();
    return false;
  }

  std::cout << "Loaded TPC channel mask " << source
            << " with " << m_size << " channels" << std::endl;
  return true;
}
