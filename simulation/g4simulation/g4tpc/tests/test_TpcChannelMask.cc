#include "TpcChannelMask.h"

#include <trackbase/TpcDefs.h>

#include <cstdlib>
#include <iostream>

namespace
{
  void require(const bool condition, const char* message)
  {
    if (!condition)
    {
      std::cerr << "test_TpcChannelMask: " << message << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}  // namespace

int main(const int argc, char* argv[])
{
  TpcChannelMask mask;

  const auto target = TpcDefs::genHitSetKey(23, 7, 1);
  const auto wrong_layer = TpcDefs::genHitSetKey(24, 7, 1);
  const auto wrong_sector = TpcDefs::genHitSetKey(23, 8, 1);
  const auto wrong_side = TpcDefs::genHitSetKey(23, 7, 0);

  require(mask.empty(), "new mask is not empty");
  require(mask.add(target, 731), "could not add target channel");
  require(!mask.add(target, 731), "duplicate channel was accepted");
  require(mask.size() == 1, "duplicate changed the mask size");

  require(mask.contains(target, 731), "target channel is not masked");
  require(!mask.contains(target, 730), "lower neighboring pad is masked");
  require(!mask.contains(target, 732), "upper neighboring pad is masked");
  require(!mask.contains(wrong_layer, 731), "wrong layer is masked");
  require(!mask.contains(wrong_sector, 731), "wrong sector is masked");
  require(!mask.contains(wrong_side, 731), "wrong side is masked");

  // A channel mask is deliberately independent of the time bin. The same
  // pad number is extracted from every time-bin hit key on this channel.
  for (unsigned int tbin : {0U, 1U, 137U, 511U})
  {
    const auto hitkey = TpcDefs::genHitKey(731, tbin);
    require(
        mask.contains(target, TpcDefs::getPad(hitkey)),
        "mask unexpectedly depends on the time bin");
  }

  mask.clear();
  require(mask.empty(), "clear did not empty the mask");
  require(!mask.contains(target, 731), "cleared channel remains masked");

  // An explicitly requested but unresolved map must fail closed.
  require(
      !mask.load("", "TotalDeadChannels", true),
      "empty map name did not fail closed");
  require(mask.empty(), "failed load left a partially populated mask");

  // An optional real CDB payload can be supplied when running this test
  // manually; the normal test suite remains self-contained.
  if (argc == 2)
  {
    require(
        mask.load(argv[1], "TotalDeadChannels", true),
        "could not load supplied real map");
    require(!mask.empty(), "supplied real map loaded as empty");
  }

  std::cout << "TpcChannelMask tests passed" << std::endl;
  return 0;
}
