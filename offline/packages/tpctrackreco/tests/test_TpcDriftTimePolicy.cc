#include "TpcDriftTimePolicy.h"

#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
  void require(bool condition,const char *message)
  { if (!condition) { throw std::runtime_error(message); } }
}

int main()
{
  try
  {
    // Independently reproduce the previous two callers' accepted-time
    // expressions across exact endpoints and a dense range of normal/late
    // samples. CM classification must have no effect when the option is zero.
    for (const double edge:{0.0,56.881262})
    {
      for (const bool cm:{false,true})
      {
        for (int index=0;index<20000;++index)
        {
          const double end=317*56.881262;
          const double time=end-100+index*.01;
          const bool old_accepted=time<=end+edge;
          const auto current=resolve_tpc_drift_time(time,end,edge,0,cm);
          require(current.accepted==old_accepted,"zero option changed upstream acceptance");
          if (old_accepted)
          { require(current.time_ns==std::min(time,end),"zero option changed upstream sampled time"); }
          require(!current.central_membrane_clamped,"zero option reported a CM extension");
        }
      }
    }
    // Three points separated by 50 ns end at 100 ns. ReverseDrift omits its
    // terminal partial step; a CM status does not imply a stored point at CM.
    // The configured tolerance is measured from this LAST STORED point.
    auto decision=resolve_tpc_drift_time(115,100,0,15,true);
    require(decision.accepted&&decision.central_membrane_clamped&&decision.time_ns==100,
            "bounded CM endpoint clamp failed");
    require(!resolve_tpc_drift_time(115.00001,100,0,15,true).accepted,"CM extension exceeded its bound");
    for (const double time:{99.0,100.0})
    {
      decision=resolve_tpc_drift_time(time,100,0,500,true);
      require(decision.accepted&&!decision.central_membrane_clamped&&decision.time_ns==time,
              "CM option shifted normal drift times");
    }
    // Every non-CM terminal status (radial, z, stuck, invalid velocity,
    // maxsteps or unknown) supplies false and cannot gain extra acceptance.
    require(!resolve_tpc_drift_time(101,100,0,10000,false).accepted,"non-CM endpoint was extended");
    decision=resolve_tpc_drift_time(105,100,10,50,false);
    require(decision.accepted&&!decision.central_membrane_clamped&&decision.time_ns==100,
            "crossing's existing edge policy changed");
    require(!resolve_tpc_drift_time(111,100,10,50,false).accepted,"CM option enlarged radial crossing edge");
    decision=resolve_tpc_drift_time(150,100,10,50,true);
    require(decision.accepted&&decision.central_membrane_clamped&&decision.time_ns==100,
            "independent CM crossing extension failed");
    require(!resolve_tpc_drift_time(150.001,100,10,50,true).accepted,"CM and legacy tolerances were incorrectly added");
    require(!resolve_tpc_drift_time(std::numeric_limits<double>::quiet_NaN(),100,0,50,true).accepted,
            "invalid time was accepted");
    std::cout<<"PASS: 80,000 default-policy comparisons; bounded CM-only clamping; normal times unchanged; "
                "radial/invalid endpoints excluded; independent legacy crossing edge; truncated-step semantics\n";
    return 0;
  }
  catch(const std::exception &error) { std::cerr<<"FAIL: "<<error.what()<<'\n';return 1; }
}
