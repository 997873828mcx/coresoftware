#include "PHGarfieldCompactMap.h"
#include "PHGarfieldRootOutput.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TNamed.h>
#include <TParameter.h>
#include <algorithm>
#include <cerrno>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace
{
using Map=PHGarfieldCompactMap;
void need(bool value,const std::string& text){if(!value)throw std::runtime_error(text);}
template<class T>T parameter(TDirectory* directory,const char* name)
{auto* object=dynamic_cast<TParameter<T>*>(directory->Get(name));need(object,"missing/incorrectly typed source parameter "+std::string(name));return object->GetVal();}
std::string named(TDirectory* directory,const char* name)
{auto* object=dynamic_cast<TNamed*>(directory->Get(name));need(object,"missing source string "+std::string(name));return object->GetTitle();}
std::uint64_t unsigned_option(const std::string& value)
{
  need(!value.empty()&&std::all_of(value.begin(),value.end(),[](char c){return c>='0'&&c<='9';}),"resource caps require unsigned decimal integers");
  std::size_t used=0;const auto parsed=std::stoull(value,&used);need(used==value.size(),"invalid resource cap suffix");return parsed;
}
std::string quote(const std::string& value)
{
  std::ostringstream out;out<<'"';
  for(const unsigned char c:value)
  {
    if(c=='"'||c=='\\')out<<'\\'<<static_cast<char>(c);
    else if(c<32)out<<"\\u"<<std::hex<<std::setfill('0')<<std::setw(4)<<static_cast<unsigned int>(c)<<std::dec;
    else out<<static_cast<char>(c);
  }
  out<<'"';return out.str();
}
void copy_metadata(TDirectory* input,TDirectory* output)
{
  TDirectory::TContext context(output);TIter next(input->GetListOfKeys());
  while(auto* key=dynamic_cast<TKey*>(next()))
  {
    std::unique_ptr<TObject> object(key->ReadObj());need(object!=nullptr,"unreadable source metadata");
    const std::string type=object->ClassName();need(type=="TNamed"||type.rfind("TParameter<",0)==0,"unexpected source metadata object "+type);
    need(object->Write(key->GetName())>0,"cannot copy source metadata");
  }
}
void exclusive_text(const std::string& path,const std::string& text)
{
  const int fd=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);need(fd>=0,"cannot create exclusive proof temporary");
  std::size_t done=0;
  while(done<text.size())
  {
    const auto count=::write(fd,text.data()+done,text.size()-done);
    if(count<0&&errno==EINTR)continue;
    if(count<=0){::close(fd);throw std::runtime_error("proof write failed; partial evidence preserved");}
    done+=static_cast<std::size_t>(count);
  }
  const bool synced=::fsync(fd)==0;const bool closed=::close(fd)==0;need(synced&&closed,"proof flush/close failed");
}
}

int main(int argc,char**argv)
{
  try
  {
    std::map<std::string,std::string> options;
    for(int i=1;i<argc;i+=2)need(i+1<argc&&options.emplace(argv[i],argv[i+1]).second,"options require unique key/value pairs");
    const auto required=[&](const char* key)->const std::string&{need(options.count(key),std::string("missing ")+key);return options.at(key);};
    const std::vector<std::string> allowed{"--input","--output","--proof-output","--expect-source-sha256","--expect-source-configuration","--expect-source-content",
      "--max-source-points","--max-runtime-nodes","--max-runtime-values","--max-runtime-corners","--max-runtime-bytes"};
    for(const auto& [key,value]:options){(void)value;need(std::find(allowed.begin(),allowed.end(),key)!=allowed.end(),"unsupported option "+key);}
    const auto input_path=std::filesystem::canonical(required("--input"));
    const auto output_path=std::filesystem::absolute(required("--output")).lexically_normal();
    const auto proof_path=std::filesystem::absolute(required("--proof-output")).lexically_normal();
    need(!std::filesystem::exists(output_path)&&!std::filesystem::exists(proof_path),"preserve existing output/proof");
    need(output_path!=proof_path,"ROOT and proof outputs must differ");
    Map::Limits limits;
    for(const auto& [key,target]:std::vector<std::pair<std::string,std::uint64_t*>>{
          {"--max-source-points",&limits.max_source_points},{"--max-runtime-nodes",&limits.max_nodes},
          {"--max-runtime-values",&limits.max_values},{"--max-runtime-corners",&limits.max_corners},{"--max-runtime-bytes",&limits.max_bytes}})
      if(options.count(key))*target=unsigned_option(options.at(key));
    std::ostringstream errors;
    const bool valid_limits=Map::ValidateLimits(limits,errors);need(valid_limits,errors.str());
    Map::ExpectedSource expected{required("--expect-source-sha256"),required("--expect-source-configuration"),required("--expect-source-content")};
    need(Map::FileSha256(input_path)==expected.file_sha256,"input does not match externally pinned source SHA256");
    TFile input(input_path.c_str(),"READ");need(!input.IsZombie()&&!input.TestBit(TFile::kRecovered),"training ROOT file unreadable/recovered");
    auto* metadata=input.GetDirectory("GarfieldFullTransport/Metadata");need(metadata,"source metadata missing");
    need(parameter<int>(metadata,"schema_version")==4&&named(metadata,"schema_name")==PHGarfieldForwardMap::SchemaName,"converter requires complete schema4 training input");
    need(named(metadata,"configuration_fingerprint_fnv1a64")==expected.configuration_fingerprint_fnv1a64,"source configuration differs from external expectation");
    PHGarfieldAdaptiveMap::Geometry geometry;
    geometry.phi_nodes=parameter<int>(metadata,"phi_nodes");geometry.radial_nodes=parameter<int>(metadata,"radial_nodes");geometry.longitudinal_nodes=parameter<int>(metadata,"longitudinal_nodes");
    geometry.inner_radius_cm=parameter<double>(metadata,"gas_inner_radius_cm");geometry.outer_radius_cm=parameter<double>(metadata,"gas_outer_radius_cm");
    geometry.central_membrane_half_width_cm=parameter<double>(metadata,"central_membrane_half_width_cm");geometry.readout_abs_z_cm=parameter<double>(metadata,"readout_abs_z_cm");
    geometry.nominal_velocity_cm_ns=parameter<double>(metadata,"nominal_phg4_drift_velocity_cm_per_ns");
    PHGarfieldAdaptiveMap::Data source;
    need(PHGarfieldAdaptiveMap::Read(input.GetDirectory("GarfieldFullTransport/Adaptive"),geometry,source,errors,true),errors.str());
    need(PHGarfieldAdaptiveMap::ContentFingerprint(source)==expected.content_fingerprint_fnv1a64,"source content differs from external expectation");
    Map::SourceBinding binding{expected.file_sha256,expected.configuration_fingerprint_fnv1a64,expected.content_fingerprint_fnv1a64,
      PHGarfieldAdaptiveMap::BuildCanonicalPolicy(source),Map::FileSha256("/proc/self/exe")};
    Map::Data data;need(Map::Convert(source,binding,data,limits,errors,true),errors.str());
    const auto temp_path=output_path.string()+".temporary-"+std::to_string(getpid());
    {
      TFile output(temp_path.c_str(),"CREATE");need(!output.IsZombie(),"cannot create unique temporary ROOT output");
      auto* top=output.mkdir("GarfieldFullTransport");need(top,"cannot create transport directory");
      auto* identity=top->mkdir("Metadata");need(identity,"cannot create compact identity metadata");
      {
        TDirectory::TContext context(identity);
        need(TParameter<int>("schema_version",Map::SchemaVersion).Write()>0&&TNamed("schema_name",Map::SchemaName).Write()>0&&
          TParameter<int>("source_schema_version",4).Write()>0&&TNamed("compact_layout",Map::Layout).Write()>0&&
          TNamed("source_file_sha256",binding.file_sha256.c_str()).Write()>0&&
          TNamed("source_configuration_fingerprint_fnv1a64",binding.configuration_fingerprint_fnv1a64.c_str()).Write()>0&&
          TNamed("source_content_fingerprint_fnv1a64",binding.content_fingerprint_fnv1a64.c_str()).Write()>0,"cannot write compact identity metadata");
      }
      auto* copied=top->mkdir("SourceMetadata");need(copied,"cannot create source metadata directory");copy_metadata(metadata,copied);
      need(Map::Write(top->mkdir("Runtime"),data,limits,errors),errors.str());
      need(PHGarfieldFlushAndCloseRootOutput(output,errors),errors.str());
    }
    const auto content=Map::ContentFingerprint(data);const auto proof=data.proof;const auto nodes=data.nodes.size(),corners=data.corners.size(),values=data.values.size();
    data={}; // Release converted vectors before independent ROOT readback.
    {
      TFile check(temp_path.c_str(),"READ");need(!check.IsZombie()&&!check.TestBit(TFile::kRecovered),"compact ROOT output reopened recovered/unreadable");
      need(Map::Read(check.GetDirectory("GarfieldFullTransport/Runtime"),geometry,expected,data,limits,errors,true),errors.str());
      need(Map::ContentFingerprint(data)==content&&Map::VerifyConversion(source,data,errors),errors.str());
    }
    input.Close();need(Map::FileSha256(input_path)==expected.file_sha256,"source changed during conversion");
    const auto output_sha=Map::FileSha256(temp_path);
    need(::link(temp_path.c_str(),output_path.c_str())==0,"atomic ROOT publication failed; preserved temporary file");
    need(::unlink(temp_path.c_str())==0,"published ROOT but temporary unlink failed");
    std::ostringstream json;
    json<<"{\n  \"schema\":\"phgarfield_compact_conversion_manifest_v1\",\n  \"status\":\"PASS_CONVERSION_ONLY\",\n"
      <<"  \"source\":{\"path\":"<<quote(input_path)<<",\"sha256\":"<<quote(expected.file_sha256)<<"},\n"
      <<"  \"compact\":{\"path\":"<<quote(output_path)<<",\"sha256\":"<<quote(output_sha)<<"},\n"
      <<"  \"source_configuration_fingerprint_fnv1a64\":"<<quote(expected.configuration_fingerprint_fnv1a64)<<",\n"
      <<"  \"source_content_fingerprint_fnv1a64\":"<<quote(expected.content_fingerprint_fnv1a64)<<",\n"
      <<"  \"compact_content_fingerprint_fnv1a64\":"<<quote(content)<<",\n  \"conversion_proof_sha256\":"<<quote(proof.sha256)
      <<",\n  \"converter_sha256\":"<<quote(binding.converter_sha256)<<",\n  \"nodes\":"<<nodes<<",\n  \"certified_corner_records\":"<<corners<<",\n  \"unique_corner_values\":"<<values<<",\n"
      <<"  \"runtime_limits\":{\"max_source_points\":"<<limits.max_source_points<<",\"max_nodes\":"<<limits.max_nodes
      <<",\"max_corners\":"<<limits.max_corners<<",\"max_values\":"<<limits.max_values<<",\"max_roots\":"<<limits.max_roots<<",\"max_bytes\":"<<limits.max_bytes<<"},\n"
      <<"  \"checks\":{\"all_source_bounds_routes_corner_bits_exact\":true,\"independent_ROOT_readback_unrecovered\":true,\"source_immutable\":true},\n"
      <<"  \"scope\":\"Conversion and structure only; source accuracy/coverage and detector-production approval remain external prerequisites.\"\n}\n";
    const auto temporary_proof=proof_path.string()+".temporary-"+std::to_string(getpid());
    exclusive_text(temporary_proof,json.str());
    need(::link(temporary_proof.c_str(),proof_path.c_str())==0,"atomic proof publication failed; ROOT/proof evidence preserved");
    need(::unlink(temporary_proof.c_str())==0,"proof temporary unlink failed");
    std::cout<<json.str();return 0;
  }
  catch(const std::exception& error){std::cerr<<"Compact conversion failed: "<<error.what()<<'\n';return 1;}
}
