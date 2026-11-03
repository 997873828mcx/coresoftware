#include "../CDBInterface.h"

#include <ffaobjects/CdbUrlSave.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition) { throw std::runtime_error(message); }
  }
}

int main()
{
  char filename[] = "/tmp/CDBInterfaceCache.XXXXXX";
  const int descriptor = mkstemp(filename);
  require(descriptor >= 0, "cannot create cache fixture");
  close(descriptor);
  try
  {
    {
      std::ofstream fixture(filename);
      fixture << "used cached-payload.root\nunused unused-payload.root\n";
      require(bool(fixture), "cannot write cache fixture");
    }
    auto* server = Fun4AllServer::instance();
    auto* cache = CDBInterface::instance();
    cache->ReadCalibrationsFromFile(filename);
    auto* constants = recoConsts::instance();
    require(!constants->FlagExist("TIMESTAMP"), "fixture must begin without timestamp");
    require(cache->getUrl("used") == "cached-payload.root", "cached URL changed");
    require(cache->getUrl("missing", "must-not-fall-back.root").empty(), "cache miss used a live/default URL");
    require(cache->UpdateRunNode(server->topNode()) == 0, "cannot save unknown-time request");
    constants->set_uint64Flag("TIMESTAMP", 76905);
    require(cache->getUrl("used") == "cached-payload.root", "timestamped cached URL changed");
    require(cache->getUrl("used") == "cached-payload.root", "repeated cached URL changed");
    cache->Disable();
    require(cache->getUrl("unused").empty(), "disabled interface returned a URL");
    cache->Enable();
    require(cache->UpdateRunNode(server->topNode()) == 0, "cannot save timestamped request");
    require(cache->UpdateRunNode(server->topNode()) == 0, "repeat RUN update failed");
    const auto* saved = findNode::getClass<CdbUrlSave>(server->topNode(), "CdbUrl");
    require(saved != nullptr, "RUN payload provenance missing");
    std::set<std::tuple<std::string, std::string, std::uint64_t>> observed;
    unsigned int count = 0;
    for (const auto& row : *saved) { observed.insert(row); ++count; }
    const std::set<std::tuple<std::string, std::string, std::uint64_t>> expected{
        {"used", "cached-payload.root", 0}, {"used", "cached-payload.root", 76905}};
    require(observed == expected && count == expected.size(), "cached RUN provenance changed, duplicated, or includes unused URLs");
    delete server;
    std::remove(filename);
    std::cout << "PASS: cached payloads retain URL semantics and exact RUN access provenance\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    std::remove(filename);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
