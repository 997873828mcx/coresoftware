#include "TpcSingleTrackDiagnostics.h"
#include "TpcPolyTrackTruthEvaluator.h"

#include <tpctrackreco/Tpc_PolyClusterContainerv1.h>
#include <tpctrackreco/Tpc_PolyClusterv1.h>
#include <tpctrackreco/Tpc_PolyTrackContainerv1.h>
#include <tpctrackreco/Tpc_PolyTrackv1.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <g4detectors/PHG4TpcGeomv2.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPointv1.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssocv1.h>

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unistd.h>

namespace
{
  void require(const bool condition,const char *message)
  { if (!condition) { throw std::runtime_error(message); } }
  void node(PHCompositeNode *parent,PHObject *object,const std::string &name)
  { parent->addNode(new PHIODataNode<PHObject>(object,name,"PHObject")); }
  TTree *tree(TFile &file,const char *name)
  { auto *value=dynamic_cast<TTree *>(file.Get(name)); require(value,"missing output tree"); return value; }
}

int main()
{
  try
  {
    const auto scratch=std::filesystem::temp_directory_path()/ ("tpc_single_track_diagnostics_"+std::to_string(getpid()));
    require(std::filesystem::create_directory(scratch),"test scratch directory already exists");
    auto top=std::make_unique<PHCompositeNode>("TOP");
    auto *clusters=new Tpc_PolyClusterContainerv1();
    auto *tracks=new Tpc_PolyTrackContainerv1();
    auto *geometry=new PHG4TpcGeomContainer();
    node(top.get(),clusters,"TPC_POLYCLUSTERS"); node(top.get(),tracks,"TPC_POLYTRACKS");
    node(top.get(),geometry,"TPCGEOMCONTAINER");
    auto *good=new Tpc_PolyTrackv1(); good->set_track_id(0); good->set_fit_status(1);
    good->set_source_assembled_track_id(17); good->set_charge(1);
    good->set_x(0);good->set_y(0);good->set_z(2);good->set_px(0.8*std::cos(0.7));good->set_py(0.8*std::sin(0.7));good->set_pz(.2);
    auto config=TpcTrackDiagnosticFitter::yuri_config(nullptr);config.analytic_uniform_propagation=true;
    const std::array<double,6> truth{{0,0,2,.7,-1/.8,.25}};
    for (unsigned int index=0;index<31;++index)
    {
      const auto state=TpcTrackKalmanFitter::propagate_state(truth,32+1.5*index,config);
      auto *cluster=new Tpc_PolyClusterv1();
      const unsigned int layer=7+index;
      const auto clusterkey=TpcDefs::genClusKey(layer,2,1,100+index);
      cluster->set_cluster_id(100+index*7);cluster->set_trkr_cluster_key(clusterkey);
      cluster->set_source_assembled_track_id(index<30?17:999);cluster->set_side(1);
      cluster->set_centroid_x(state[0]);cluster->set_centroid_y(state[1]);cluster->set_centroid_z(state[2]);
      cluster->set_adc(100);cluster->set_phi_width(2);cluster->set_time_width(3);
      cluster->add_hit(TpcDefs::genHitSetKey(layer,2,1),index,state[0],state[1],state[2]);
      clusters->add_cluster(cluster);
      auto *geom=new PHG4TpcGeomv2();geom->set_radius(std::hypot(state[0],state[1]));geom->set_thickness(1);
      geometry->AddLayerCellGeom(layer,geom);
      if (index<30) { good->add_cluster_key(clusterkey); }
    }
    tracks->add_track(good);
    auto *bad=new Tpc_PolyTrackv1();bad->set_track_id(1);bad->set_source_assembled_track_id(19);
    bad->add_cluster_key(good->get_cluster_key(0));bad->add_cluster_key(good->get_cluster_key(1));tracks->add_track(bad);
    const auto enabled=scratch/"enabled.root",disabled=scratch/"disabled.root",evaluated=scratch/"truth.root";
    for (const bool loo:{false,true})
    {
      TpcSingleTrackDiagnostics diagnostics("synthetic",(loo?enabled:disabled).string());
      diagnostics.set_sample_id(123456);diagnostics.set_require_field_map(false);
      diagnostics.set_write_leave_one_out(loo);diagnostics.set_provenance("synthetic_fixture");
      require(diagnostics.Init(top.get())==0,"diagnostic Init failed");
      require(diagnostics.process_event(top.get())==0,"diagnostic event failed");
      require(diagnostics.End(top.get())==0,"diagnostic End failed");
    }
    TFile on(enabled.c_str(),"READ"),off(disabled.c_str(),"READ");
    require(!on.IsZombie()&&!off.IsZombie()&&!on.TestBit(TFile::kRecovered)&&!off.TestBit(TFile::kRecovered),"output ROOT recovered/unreadable");
    require(tree(on,"eventTree")->GetEntries()==1,"event identity missing");
    require(tree(on,"trackTree")->GetEntries()==2,"failed input track was omitted");
    require(tree(on,"clusterTree")->GetEntries()==31,"unassociated cluster was omitted");
    require(tree(on,"residualTree")->GetEntries()==64,"failed/rejected cluster residual row was omitted");
    unsigned int failed_excluded=0;std::vector<std::uint64_t> *failed_retained=nullptr;
    tree(on,"residualTree")->SetBranchAddress("excluded_measurements",&failed_excluded);
    tree(on,"residualTree")->SetBranchAddress("retained_cluster_ids",&failed_retained);
    tree(on,"residualTree")->GetEntry(63);
    require(failed_excluded==1&&failed_retained&&*failed_retained==std::vector<std::uint64_t>{100},
            "unavailable-input LOO omitted intended exclusion/retained identities");
    int success=0;auto *track_tree=tree(on,"trackTree");track_tree->SetBranchAddress("kalman_success",&success);
    track_tree->GetEntry(0);require(success==1,"nominal synthetic track failed");
    track_tree->GetEntry(1);require(success==0,"invalid polytrack unexpectedly fitted");
    std::array<double,6> on_state{},off_state{};
    track_tree->SetBranchAddress("kalman_pca_state",on_state.data());
    tree(off,"trackTree")->SetBranchAddress("kalman_pca_state",off_state.data());
    track_tree->GetEntry(0);tree(off,"trackTree")->GetEntry(0);
    require(on_state==off_state,"enabling LOO changed nominal momentum/state");
    require(good->get_px()==.8*std::cos(.7)&&clusters->size()==31,"diagnostics mutated pattern output");
    // No-hit generated pion must survive the evaluator's denominator catalog.
    auto *truth_info=new PHG4TruthInfoContainer();
    auto *particle=new PHG4Particlev2("pi+",211,.8,0,.2);
    particle->set_track_id(1);particle->set_primary_id(1);particle->set_parent_id(0);particle->set_vtx_id(1);particle->set_e(1);
    truth_info->AddParticle(1,particle);truth_info->AddVertex(1,new PHG4VtxPointv1(0,0,0,0,1));
    auto *truth_hits=new PHG4HitContainer();
    node(top.get(),truth_info,"G4TruthInfo");node(top.get(),truth_hits,"G4HIT_TPC");
    node(top.get(),new TrkrHitSetContainerv1(),"TRKR_HITSET");node(top.get(),new TrkrHitTruthAssocv1(),"TRKR_HITTRUTHASSOC");
    // Empty pattern containers avoid dangling raw-hit references in this
    // deliberately hit-free denominator fixture.
    tracks->Reset();clusters->Reset();
    const double reference_radius=geometry->GetFirstLayerCellGeom()->get_radius();
    for (int passage=0;passage<3;++passage)
    {
      auto *hit=new PHG4Hitv1();hit->set_trkid(1);hit->set_edep(.001);hit->set_eion(.001);
      const bool outward=passage!=1;
      for (int endpoint=0;endpoint<2;++endpoint)
      {
        const double offset=(endpoint==0?-1:1)*(outward?1:-1);
        hit->set_x(endpoint,reference_radius+offset);hit->set_y(endpoint,0);hit->set_z(endpoint,0);
        hit->set_px(endpoint,.8);hit->set_py(endpoint,0);hit->set_pz(endpoint,.2);
        hit->set_t(endpoint,passage*10+endpoint);
      }
      truth_hits->AddHit(7,hit);
    }
    const auto reference_output=scratch/"truth_reference.root";
    TpcSingleTrackDiagnostics references("truth_surface_fixture",reference_output.string());
    references.set_sample_id(123456);references.set_require_field_map(false);
    references.set_write_leave_one_out(false);references.set_write_truth_diagnostics(true);
    require(references.Init(top.get())==0&&references.process_event(top.get())==0&&references.End(top.get())==0,
            "truth surface diagnostics failed");
    TFile reference_file(reference_output.c_str(),"READ");
    auto *reference_tree=tree(reference_file,"truthReferenceTree");
    require(reference_tree->GetEntries()==1,"primary truth reference missing");
    double reference_time=0;unsigned int reference_intersections=0;int reference_success=0;
    reference_tree->SetBranchAddress("time_ns",&reference_time);
    reference_tree->SetBranchAddress("intersections",&reference_intersections);
    reference_tree->SetBranchAddress("success",&reference_success);
    reference_tree->GetEntry(0);
    require(reference_success==1&&reference_intersections==2&&std::abs(reference_time-.5)<1e-6,
            "truth reference did not select earliest outward crossing independently");
    require(!reference_file.IsZombie()&&!reference_file.TestBit(TFile::kRecovered),"truth-reference ROOT unreadable");
    reference_file.Close();truth_hits->Reset();
    TpcPolyTrackTruthEvaluator evaluator("primary_pion_fixture",evaluated.string());
    evaluator.set_primary_pion_mode(true);evaluator.set_sample_id(123456);
    require(evaluator.Init(top.get())==0&&evaluator.process_event(top.get())==0&&evaluator.End(top.get())==0,"zero-hit evaluator failed");
    TFile truth_file(evaluated.c_str(),"READ");
    require(!truth_file.IsZombie()&&!truth_file.TestBit(TFile::kRecovered),"truth ROOT output unreadable");
    require(tree(truth_file,"truthTree")->GetEntries()==1,"zero-hit primary pion was omitted");
    int selected=0;tree(truth_file,"truthTree")->SetBranchAddress("is_selected_primary_pion",&selected);
    tree(truth_file,"truthTree")->GetEntry(0);require(selected==1,"zero-hit primary selection missing");
    on.Close();off.Close();truth_file.Close();std::filesystem::remove_all(scratch);
    std::cout<<"PASS: data-only diagnostics, failed tracks/residuals, unassociated clusters, stable identities, LOO noninterference, earliest outward truth crossing, zero-hit primary truth, ROOT reopen\n";
    return 0;
  }
  catch (const std::exception &error) { std::cerr<<"FAIL: "<<error.what()<<'\n';return 1; }
}
