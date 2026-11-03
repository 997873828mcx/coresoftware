#ifndef G4TPC_TPCCHANNELMASK_H
#define G4TPC_TPCCHANNELMASK_H

#include <trackbase/TrkrDefs.h>

#include <cstddef>
#include <map>
#include <string>
#include <unordered_set>

/**
 * Compact, time-bin-independent TPC channel mask.
 *
 * A TPC electronics channel is identified by (layer, sector, side, pad).
 * The time bin is deliberately not part of the mask, so one entry disables
 * every time bin on the corresponding pad.
 */
class TpcChannelMask
{
 public:
  using ChannelSet = std::unordered_set<TrkrDefs::hitkey>;
  using Mask = std::map<TrkrDefs::hitsetkey, ChannelSet>;

  bool load(const std::string& map_name,
            const std::string& total_channels_key,
            bool from_file,
            int verbosity = 0);

  void clear();
  bool add(TrkrDefs::hitsetkey hitsetkey, unsigned int pad);
  bool contains(TrkrDefs::hitsetkey hitsetkey, unsigned int pad) const;

  std::size_t size() const { return m_size; }
  bool empty() const { return m_size == 0; }

 private:
  Mask m_mask;
  std::size_t m_size = 0;
};

#endif
