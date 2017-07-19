//Helper functions for Config parsing.
#ifndef INCLUDED_CONFIG_HELPER
#define INCLUDED_CONFIG_HELPER

#include "ConfigParser.C"

TString PLOT_OUTPUT_LOCATION="/home/users/bhashemi/public_html/WWWCrossSection/";
TString HIST_OUTPUT_LOCATION="/nfs-7/userdata/bobak/WWWCrossSection_Hists/";

TString parseConfDir(TString conf_path){
  /* Replace *configs/.../FNAME.conf with .../ */
  conf_path = TString(conf_path(conf_path.Index("configs/")+8, conf_path.Last('/')-conf_path.Index("configs/")-7)); 
  return conf_path;
}

TString getOutputDir(ConfigParser *conf, TString type){
  /*Determines the proper output locations of files by parsing the option conf_path to get the directory structure above the level 'configs/'*/
	//cout<<__LINE__<<endl;
  if (type == "hist"){
    //cout<<__LINE__<<endl;
		if (conf->get("histo_output_dir") != ""){
      //cout<<__LINE__<<endl;
			return TString(conf->get("histo_output_dir"));
		}
		else{
			TString output_dir = parseConfDir(conf->get("conf_path"));
      //cout<<__LINE__<<endl;
			return output_dir.Prepend(HIST_OUTPUT_LOCATION);
		}
    //cout<<__LINE__<<endl;
	}
	else if (type == "plot" )
	{
		if (conf->get("save_dir") != ""){
			return TString(conf->get("save_dir"));
		}
		else{
      TString conf_path = conf->get("conf_path");
      TString output_dir = parseConfDir(conf_path);
      output_dir.Prepend(PLOT_OUTPUT_LOCATION);
      conf_path = TString(conf_path(conf_path.Last('/')+1, conf_path.Last('.')-conf_path.Last('/')-1)); //get name of config file
      output_dir+=conf_path+"/";
      //also add filename for conf script
      return output_dir;
		}	
	}
  else{
    TString error = type.Prepend("In ConfigHelper::getOutputDir -- Did not recieve valid type, either hist or plot... got: ");
    throw std::invalid_argument(error.Data());
    return TString("");
  }
}

TString getDefaultHistDir(ConfigParser *conf){
	/*Returns the default hist output location + the conf_path defined by the location of the config file*/
	return HIST_OUTPUT_LOCATION+parseConfDir(conf->get("conf_path"));
}

TString parseLatex(TString opt){
  /* Replaces \ with # to be in line with ROOT's latex syntax */
  opt.ReplaceAll("\\","#");
  return opt;
}

vector<double> parseVector(TString opt){
	/* Parses options in the config files which are formatted to be vectors [double,double,double,...]*/
	vector<double> ret;
  TString token;
  Ssiz_t from=0;
  //cout<<"got vector in string form: "<<opt<<endl;
  while(opt.Tokenize(token, from, "[,]")){
    token.ReplaceAll("[", "");
    token.ReplaceAll("]", "");
    token.ReplaceAll(",", "");
    token.ReplaceAll(" ", "");
    //cout<<"token: "<<token<<endl;
    ret.push_back(stod(token.Data()));
  }
  return ret;
}

vector<int> iparseVector(TString opt){
  /* Parses options in the config files which are formatted to be vectors [int,int,int,...]*/
  vector<int> ret;
  TString token;
  Ssiz_t from=0;
  //cout<<"got vector in string form: "<<opt<<endl;
  while(opt.Tokenize(token, from, "[,]")){
    token.ReplaceAll("[", "");
    token.ReplaceAll("]", "");
    token.ReplaceAll(",", "");
    token.ReplaceAll(" ", "");
    //cout<<"token: "<<token<<endl;
    ret.push_back(stoi(token.Data()));
  }
  return ret;
}

vector<TString> sParseVector(TString opt){
  /* Parses options in the config files which are formatted to be vectors [string,string,string,...]*/
  vector<TString> ret;
  TString token;
  Ssiz_t from=0;
  //cout<<"got vector in string form: "<<opt<<endl;
  while(opt.Tokenize(token, from, "[,]")){
    token.ReplaceAll("[", "");
    token.ReplaceAll("]", "");
    token.ReplaceAll(",", "");
    token.ReplaceAll(" ", "");
    //cout<<"token: "<<token<<endl;
    ret.push_back(TString(token.Data()));
  }
  return ret;
}

/*TString combineTH1D(TString output_dir, TString sample_names, hist_name){
  // This method will eventually be used to combine my files....

  vector<TH1D*> hists;
  vector<TString> fnames = sParseVector(sample_names);
  for (int i = 0; i<(int)samples.size(); i++){
    TFile *f = new TFile(output_dir+fnames[i]+".root", "r");
    hists.push_back((TH1D*) fget(hist_name)->Clone("h_"+to_string(i)));
  }
  output_name = sample_names;
  output_name.ReplaceAll("[", "");
  output_name.ReplaceAll("]", "");
  output_name.ReplaceAll(",", " ");
  output_name.ReplaceAll(" ", "_");

  TString new_file_name = TString(output_dir+'/'+output_name+".root");

  TFile *f_new = new TFile(new_file_name, "w");
  f_new->cd();
  
  TH1D* h_new = hists[0]->Clone(hist_name);
  h_new->SetDirectory(rootdir);

  for (int i = 1; i <(int) samples.size(); i++) h_new->Add(hists[i]);

  h_new->Write();
  f_new.close();

  return new_file_name;
}*/

#endif