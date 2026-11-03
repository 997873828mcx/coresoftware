#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

namespace
{
#include "production_helpers.inc"

  void require(bool condition,const char *message)
  { if(!condition) throw std::runtime_error(message); }
  template<typename Function> void rejected(Function operation)
  {
    bool caught=false;
    try { operation(); } catch(const std::runtime_error &) { caught=true; }
    require(caught,"invalid geometry/ID was not rejected");
  }
}

int main()
{
  try
  {
    const MvtxSensorLayout physical{{0,12,2.6},{1,16,3.4},{2,20,4.2}};
    std::vector<std::array<unsigned int,3>> expected;
    for(const auto &layer:physical)
      for(unsigned int stave=0;stave<layer.staves;++stave)
        for(unsigned int chip=0;chip<9;++chip) expected.push_back({layer.layer,stave,chip});
    for(const bool misalignment:{false,true})
    {
      // Identity decoding deliberately has no misalignment parameter. The
      // caller may shift centers, while the same physical address must result.
      const double displacement_corrected_radius=misalignment?3.42:3.4;
      const auto merged=mvtx_group_layouts(physical,{432},{displacement_corrected_radius},9);
      require(merged.size()==1,"merged group was split");
      std::set<std::array<unsigned int,3>> addresses;
      for(std::size_t index=0;index<expected.size();++index)
      {
        const auto result=mvtx_sensor_address(merged.front(),index+1,9);
        require(result==expected[index],"merged ID has wrong physical address");
        require(addresses.insert(result).second,"merged duplicate sensor address");
      }
      require(addresses.size()==432,"merged coverage incomplete");
      require(mvtx_sensor_address(merged.front(),109,9)==std::array<unsigned int,3>{1,0,0},"first layer1 sensor regressed");
      require(mvtx_sensor_address(merged.front(),253,9)==std::array<unsigned int,3>{2,0,0},"first layer2 sensor regressed");
      const double shift=misalignment?.02:0.;
      const auto legacy=mvtx_group_layouts(physical,{180,108,144},{4.2+shift,2.6-shift,3.4+shift},9);
      addresses.clear();
      for(const auto &group:legacy)
      {
        const auto &layer=group.front();
        unsigned int sensitive=1;
        for(unsigned int stave=0;stave<layer.staves;++stave)
          for(unsigned int chip=0;chip<9;++chip,++sensitive)
          {
            const auto result=mvtx_sensor_address(group,sensitive,9);
            require(result==std::array<unsigned int,3>{layer.layer,stave,chip},"legacy reset ID has wrong address");
            require(addresses.insert(result).second,"legacy duplicate physical sensor");
          }
      }
      require(addresses.size()==432,"legacy reset coverage incomplete");
    }
    const MvtxSensorLayout equal_counts{{0,12,2.5},{1,12,3.5},{2,12,4.5}};
    const auto resolved=mvtx_group_layouts(equal_counts,{108,108,108},{4.48,2.52,3.51},9);
    require(resolved[0][0].layer==2&&resolved[1][0].layer==0&&resolved[2][0].layer==1,"radius did not resolve equal counts");
    rejected([&]{mvtx_sensor_address(physical,0,9);});
    rejected([&]{mvtx_sensor_address(physical,433,9);});
    rejected([&]{mvtx_sensor_address({physical[1]},145,9);});
    rejected([&]{mvtx_sensor_address(physical,1,0);});
    rejected([&]{mvtx_group_layouts(physical,{431},{3.4},9);});
    rejected([&]{mvtx_group_layouts(physical,{108,108,180},{2.6,2.6,4.2},9);});
    rejected([&]{mvtx_group_layouts(equal_counts,{108,108,108},{3.0,3.5,4.5},9);});
    rejected([&]{mvtx_group_layouts(physical,{432},{std::numeric_limits<double>::quiet_NaN()},9);});
    rejected([&]{mvtx_group_layouts({{0,0,2.6}},{1},{2.6},9);});
    rejected([&]{mvtx_group_layouts({{0,12,2.6},{0,16,3.4}},{252},{3.0},9);});
    std::cout<<"PASS: all432 merged and legacy-reset MVTX addresses, nominal/misaligned identity independence, layer boundaries, equal-count radius resolution, invalid/duplicate/ambiguous/missing coverage rejection\n";
    return 0;
  }
  catch(const std::exception &error) { std::cerr<<"FAIL: "<<error.what()<<'\n';return 1; }
}
