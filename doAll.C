#include "TString.h"
#include "ScanChain.C"
#include "DefineDatasets.C"
#include "ConfigParser.C"
#include "modes.C"

void runScanChain(ConfigParser* conf){
  cout<<"Using config:"<<endl;
  conf->print();

  bool hadoop, skimmed;
  hadoop = (conf->get("hadoop") == "true") ? true : false;
  skimmed = (conf->get("skimmed") == "true") ? true : false;

  ScanChain(getTChain(conf->get("data_set"), hadoop, skimmed), conf);
}

void doAll ( TString config_name, TString config_file="configs/run_modes.conf" ) {

  ConfigParser *conf = new ConfigParser(config_file.Data());

  if (config_name == "all") {
    while ( conf->loadNextConfig() )
    {
      runScanChain(conf);
    }
  }

  else if ( conf->loadConfig(config_name.Data()) ){
    runScanChain(conf);
  }

  else{
    cout<<"Configuration "<<config_name<<" could not be loaded, please check the name and config file, "<<config_file<<"and try again."<<endl;
    exit(1);
  }
}

int main(int argC, char** argV) {
  TString config_name="GammaZtoNuNu";
  TString config_file="configs/gamma_mu_studies/SingleMuon_ttg/run_modes.conf";
  doAll(config_name, config_file);

  return 0;
}