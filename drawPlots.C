#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>

#include "TChain.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "THStack.h"
#include "TGaxis.h"
#include "TCut.h"
#include "TH1D.h"
#include "TColor.h"
#include "TString.h"

//#include "External/CTable.cpp"
#include "External/tdrstyle_SUSY.C"

#include "ConfigParser.C"
#include "ConfigHelper.C"
#include "HistTools.C"

using namespace std;

std::array<int, 14> ROOT_COLOR_PALATE = {46,8,9,38,40,2,30,6,28,42,3,5,7,41};

TColor rc_light_blue   = TColor(2001,91/255.,187/255.,241/255.);
TColor rc_blue         = TColor(2002,60/255.,144/255.,196/255.);
TColor rc_orange       = TColor(2003,230/255.,159/255.,0/255.);
TColor rc_brown        = TColor(2004,180/255.,117/255.,0/255.);
TColor rc_yellow       = TColor(2005,245/255.,236/255.,69/255.);
TColor rc_dark_yellow  = TColor(2006,215/255.,200/255.,0/255.);
TColor rc_blue_violet  = TColor(2007,70/255.,109/255.,171/255.);
TColor rc_violet       = TColor(2008,70/255.,90/255.,134/255.);
TColor rc_dark_violet  = TColor(2009,55/255.,65/255.,100/255.);
TColor rc_light_green  = TColor(2010,120/255.,160/255.,0/255.);
TColor rc_green        = TColor(2011,0/255.,158/255.,115/255.);
TColor rc_pink         = TColor(2012,204/255.,121/255.,167/255.);

std::map<TString, Int_t> s_colors;

/*std::pair<double,double> getLegendLocation(TH1D* bg_sum){
  // Returns the best (xmin, ymin) pair for the location of the legend 

}*/

double errMult(double A, double B, double errA, double errB, double C) {
  return sqrt(C*C*(pow(errA/A,2) + pow(errB/B,2)));
}

float err_binomial(float A, float B, float errA, float errB) {
  /* 
  returns the error on C = A/(A+B) 
  note that if A and B are integers, simplifies to sqrt((C * (1-C)) / (A+B))
  or thinking of an efficiency, sqrt((eff * (1-eff)) / N) 
  */
 return (1/pow(A+B,2)) * sqrt(pow(B*errA,2) + pow(A*errB,2));
}

void drawLatexFromTString(TString text, double x_low, double y_low){
  TLatex *lumitex = NULL;
  
  lumitex = new TLatex(x_low, y_low , text );      
  lumitex->SetNDC();    
  lumitex->SetTextSize(30);    
  lumitex->SetLineWidth(2);
  lumitex->SetTextFont(43);    
  lumitex->Draw();

  return;
}

bool TH1DIntegralSort(TH1D* hist_1, TH1D* hist_2){
  return (hist_1->Integral() < hist_2->Integral()) ;
}

bool LabeledHistSort(pair<TH1D*,TString> hist_1, pair<TH1D*,TString> hist_2){
  return (hist_1.first->Integral() < hist_2.first->Integral()) ;
}

void drawCMSLatex(double luminosity){
  TLatex *lumitex = NULL;
  double height=1.01-gPad->GetTopMargin();
  float left_margin = gPad->GetLeftMargin();

  // lumitex = new TLatex(0.66,0.955, Form("%.1f fb^{-1} (13 TeV)", luminosity) );
  if (luminosity < 1){  
    lumitex = new TLatex(.9-left_margin, height , Form("%.1f pb^{-1} (13 TeV)", luminosity*1000) );    
  }
  else{  
    lumitex = new TLatex(.9-left_margin, height , Form("%.1f fb^{-1} (13 TeV)", luminosity) );    
  }    
  // lumitex = new TLatex(0.66,0.955, Form("few pb^{-1} (13 TeV)") );    
  lumitex->SetNDC();    
  lumitex->SetTextSize(0.03);    
  lumitex->SetLineWidth(2);
  lumitex->SetTextFont(42);    
  lumitex->Draw();

  TLatex *cmstex = NULL;
  cmstex = new TLatex(left_margin, height, "#it{CMS #bf{Preliminary}}" );    
  cmstex->SetNDC();    
  cmstex->SetTextSize(0.03);    
  cmstex->SetLineWidth(2);
  cmstex->SetTextFont(62);    
  cmstex->Draw();

  return;
}

void drawSRText(TString SR, double high_y, double low_x){
  cout<<"Drawing SR Text"<<endl;
  TString text;
  float left_margin = gPad->GetLeftMargin();

  if(SR == "Strong_Bveto_2j"){
    text="#splitline{2-3 jets; No b-tags}{H_{T} > 500 GeV; M_{T2} > 80 GeV}";
  }
  else if(SR == "Strong_Bveto_4j"){
    text="#splitline{4-5 jets; No b-tags}{H_{T} > 500 GeV; M_{T2} > 80 GeV}";
  }
  else if(SR == "Strong_Bveto_6j"){
    text="#splitline{#geq 6 jets; No b-tags}{M_{T2} > 80 GeV}";
  }
  else if(SR == "Strong_Btag_2j"){
    text="#splitline{2-3 jets; #geq 1 b-tags}{H_{T} > 200 GeV; M_{T2} > 100 GeV}";
  }
  else if(SR == "Strong_Btag_4j"){
    text="#splitline{4-5 jets; #geq 1 b-tags}{H_{T} > 200 GeV; M_{T2} > 100 GeV}";
  }
  else if(SR == "Strong_Btag_6j"){
    text="#splitline{#geq 6 jets; #geq 1 b-tags}{M_{T2} > 100 GeV}";
  }
  else if(SR == "TChiWZ"){
    text="EWK WZ/ZZ Region";
  }
  else if(SR == "TChiHZ"){
    text="EWK HZ Region";
  }
  else{
    cout<<"Coult not match SR: "<<SR<<endl;
  }

  cout<<"y_max: "<<high_y<<endl;

  //TLatex *SRText = new TLatex(low_x, high_y, text.Data());   //Doesn't work because splitline was coded by an asshole
  TLatex *SRText = new TLatex(left_margin+.05, .85, text.Data());   //just awful style, thanks ROOT.
  SRText->SetNDC();
  SRText->SetTextSize(0.03);    
  SRText->SetLineWidth(2);
  SRText->SetTextFont(62);    
  SRText->Draw();

  return;
}

TString drawArbitraryNumberWithResidual(ConfigParser *conf){
  // This method expects conf to have a plot config loaded in already. 
  //In the conf, we expect there to be hist names of the form file_N_path,
  //hist_n_name, starting with 0 for the primary histogram, which is normally 
  //going to be the data events in our signal region. The rest of the hists, starting
  //from 1, are added to a THStack which is normalized to hist_0 in the bin 0-50. 
  //num_hists should be the number of the number of histograms in the plot.

  TString errors="";
  TGraphAsymmErrors *prediction_errors;

  int num_hists=0;
  int BG_sum_from=(conf->get("BG_sum_from") != "") ? stoi(conf->get("BG_sum_from")) : 1;

  //Add files from which to obtain histos
  TString default_hist_dir = getDefaultHistDir(conf);
  vector<TFile*> hist_files;
  //cout<<"num_hists: "<<num_hists<<endl;
  //cout<<"file_"<<num_hists<<"_path: "<<conf->get("file_"+to_string(num_hists)+"_path")<<endl;
  //cout<<"sample_"<<num_hists<<": "<<conf->get("sample_"+to_string(num_hists))<<endl;
  while ((conf->get("file_"+to_string(num_hists)+"_path") != "") || (conf->get("sample_"+to_string(num_hists)) != "")){
    //cout<<"Getting hist "<<num_hists<<endl;
    TString sample_loc = "";
    if (conf->get("file_"+to_string(num_hists)+"_path") != ""){
      sample_loc = TString(conf->get("file_"+to_string(num_hists)+"_path"));
    }
    else{
      sample_loc = default_hist_dir+conf->get("sample_"+to_string(num_hists))+".root";
    }
    hist_files.push_back(new TFile(sample_loc));
    num_hists++;
  }

  num_hists = (conf->get("num_hists") != "") ? stoi(conf->get("num_hists")) : num_hists;

  if (num_hists < 2){
    return TString("Less than Two hists can not be turned into a residual plot, please call drawSingleTH1");
  } 
  cout << "Found files for "<<num_hists<<" hists"<<endl;

  TString plot_name = conf->get("plot_name");
  double bin_size;
  if (conf->get("bin_size") != ""){
    bin_size = stod(conf->get("bin_size"));
  }
  else{
    bin_size=1;
  }

  cout<<"/////////////////////////////////////////////////"<<endl;
  cout << "Making Plots for: "<<plot_name<<endl;
  cout<<"/////////////////////////////////////////////////"<<endl;

  //Get name of hist to read from file
  vector<TString> hist_names (num_hists);
  for (int i = 0; i<num_hists; i++){
    if (conf->get("hist_"+to_string(i)+"_name") != ""){
      hist_names[i]=conf->get("hist_"+to_string(i)+"_name");    
    }
    else{
      hist_names[i]= conf->get("hist_0_name");
    }
  }

  //Get labels for TLegend
  vector<TString> hist_labels (num_hists);
  for (int i = 0; i<num_hists; i++){
    hist_labels[i]=parseLatex(conf->get("hist_"+to_string(i)+"_label"));    
  } 

  //Get sample names for Colors
  vector<TString> sample_names (num_hists);
  for (int i = 0; i<num_hists; i++){
    sample_names[i]=conf->get("sample_"+to_string(i));    
  }  


  cout<<"Hist names set"<<endl;
  TString xlabel=parseLatex(conf->get("xlabel"));
  TString ylabel=parseLatex(conf->get("ylabel"));
  TString save_dir=(conf->get("save_dir") != "") ? conf->get("save_dir") : getOutputDir(conf, "plot");
  TString plot_title=parseLatex(conf->get("title"));

  vector<TH1D*> hists (num_hists);
  for (int i = 0; i<num_hists; i++){
    if (sParseVector(conf->get("sample_"+to_string(i))).size() > 1){
      /*//parse names;
      vector<TString> snames = sParseVector(conf->get("sample_"+to_string(i)));
      //get vector of files for all the names
      vector<TFile> sfiles;
      for (int i = 0; i<(int) snames.size(); i++) sfiles.push_back(TFile(default_hist_dir+snames[i]+".root","r"));
      hists[i] = (TH1D*) ((TH1D*) sfiles[i].Get(hist_names[i]))->Clone("hist_"+to_string(i)+"_"+plot_name);
      for (int i = 0; i<(int) snames.size(); i++) hists[i]->Add((TH1D*) ((TH1D*) sfiles[i].Get(hist_names[i]))->Clone("hist_"+to_string(i)+"_"+plot_name));
      */
      cout<<"can't do multi-hists yet"<<endl;
      exit(1);
    }
    else{
      cout<<__LINE__<<endl;
      cout<<"Hist name: "<<hist_names[i]<<endl;
      hists[i] = (TH1D*) ((TH1D*) hist_files[i]->Get(hist_names[i]))->Clone("hist_"+to_string(i)+"_"+plot_name);
      cout<<__LINE__<<endl;
      cout<<hist_names[i]<<" found in "<<hist_files[i]->GetName()<<endl;
    }
    cout<<hist_labels[i]<<" Integral bin 0 to bin 100 Content: "<<hists[i]->Integral(0,100)<<endl;
  }  
  cout << "Histograms pulled from files, adding draw options"<<endl;
  
  //cout<<__LINE__<<endl;
  if (conf->get("relabel_x_axis") == "true"){
    hists[0]->LabelsDeflate();
  }

  double xmax = (conf->get("xmax") != "" ) ? stod(conf->get("xmax")) : hists[0]->GetNbinsX();
  double xmin = (conf->get("xmin") != "" ) ? stod(conf->get("xmin")) : 0;


  //============================================
  // Draw Data-MC Plots
  //============================================
  
  //cout<<__LINE__<<endl;
  TCanvas * c = new TCanvas("c","",2000,2000);
  c->cd();
  gPad->SetRightMargin(0.05);
  gPad->Modified();
  gStyle->SetOptStat(kFALSE);
  TPad *fullpad = new TPad("fullpad", "fullpad", 0,0,1,1);
  //cout<<__LINE__<<endl;
  fullpad->Draw();
  fullpad->cd();
  //cout<<__LINE__<<endl;
  TPad *plotpad = new TPad("plotpad", "plotpad",0,0.2,1.0,0.99);
  
  plotpad->SetRightMargin(0.05);
  if (conf->get("ExtraRightMargin") == "true")
  {
    plotpad->SetRightMargin(0.08);
  }
  plotpad->SetBottomMargin(0.12);
  //cout<<__LINE__<<endl;
  plotpad->Draw();
  plotpad->cd();
  //cout<<__LINE__<<endl;
  if (conf->get("logy") == "true")
  {
    cout<<"Plot tagged for log y-axis"<<endl;
    plotpad->SetLogy();
  }
  //cout<<__LINE__<<endl;
  if (conf->get("bin_size") != ""){
    for (int i = 0; i<num_hists; i++){
      hists[i]->Rebin(bin_size);
    }
  }
  else if (conf->get("binning") != ""){
    vector<double> binning = parseVector(conf->get("binning"));
    for (int i = 0; i<num_hists; i++){
      hists[i] = (TH1D*) hists[i]->Rebin(binning.size()-1, TString(hist_names[i]+"_rebin"), &binning[0]);
    }
  }

  cout<<"===========================================\nRebin\n===========================================\n";
  for (int i = 0; i<num_hists; i++){
    cout<<hist_labels[i]<<" Integral bin 0 to bin 100 Content: "<<hists[i]->Integral(hists[i]->FindBin(0),hists[i]->FindBin(100))<<endl;
  } 

  //===========================
  // Normalize
  //===========================
  //cout<<__LINE__<<endl;
  TH1D* clonedBG_norm = NULL;
  TH1D* clonedPrimary_norm = (TH1D*) hists[0]->Clone("clonedPrimary_forNorm_"+plot_name);
  double numEventsData, numEventsMC, errEventsMC;
  int norm_bin;
  //Add scale factors like RSFOF
  cout<<"Scaling Hists"<<endl;
  for (int i=0; i < num_hists; i++){
    if (conf->get("hist_"+to_string(i)+"_scale") != ""){
      hists[i]->Scale(stod(conf->get("hist_"+to_string(i)+"_scale")));
    }
  }

  cout<<"Zeroing Negative Numbers"<<endl;
  for (int i=0; i < num_hists; i++){
    zeroNegatives(hists[i]);
  }

  if (conf->get("normalize") == "true"){
    TString hist_nums_for_norm = conf->get("normalize_hist_nums");
    cout<<"Normalizing hists: "<<hist_nums_for_norm<<endl;
    //cout<<__LINE__<<endl;
    
    //This if statement is used when there is 
    //special normalization being done. This part of
    //the code makes a proxy for the BG sum and for the
    //primary hist (hists[0]). The primary proxy is just a clone,
    //the BG proxy is the sum of all the hists that are named in 
    //hist_nums_for_norm. If subrtact_non_normed is also true then
    //the hists not marked to be normalized will be subtracted from the 
    //primary before the scale factors are calculated.
    
    if (hist_nums_for_norm != ""){ //if not all hists to be normed.
      //cout<<__LINE__<<endl;
      for (int i=1; i<num_hists; i++){ //for each hist
        if (hist_nums_for_norm.Contains(to_string(i))){ //if hist is in list to be normed
          cout<<"Adding hist "<<i<<" to normalization BG."<<endl;
          //cout<<__LINE__<<endl;
          if (clonedBG_norm == NULL){ //Make new clonedBG
            //cout<<__LINE__<<endl;
            clonedBG_norm = (TH1D*) hists[i]->Clone("clonedBG_forNorm_"+plot_name);
            cout<<hist_labels[i]<<": original-num "<<hists[i]->GetBinContent(1)<<endl;
          }
          else{
            //cout<<__LINE__<<endl;
            clonedBG_norm->Add(hists[i]); //Add to clonedBG
            cout<<hist_labels[i]<<": original-num "<<hists[i]->GetBinContent(1)<<endl;
          }
          //cout<<__LINE__<<endl;
        }
      }
    }
    else{ //otherwise put all hists in clonedBG
      clonedBG_norm = (TH1D*) hists[1]->Clone("clonedBG_forNorm_"+plot_name);
      for (int i=2; i<num_hists; i++){
        clonedBG_norm->Add(hists[i]);
      }
    }
    //cout<<__LINE__<<endl;

    //check to make sure bg hist is not empty
    if (clonedBG_norm == NULL){
      return TString("Check the normalize_hist_nums opt, no hists in range labeled for normalization");
    }
    //cout<<__LINE__<<endl;

    if (conf->get("subtract_non_normed")=="true"){
      //cout<<__LINE__<<endl;
      for (int i=1; i<num_hists; i++){
        //cout<<__LINE__<<endl;
        if( ! hist_nums_for_norm.Contains(to_string(i))){
          //cout<<__LINE__<<endl;
          cout<<hist_labels[i]<<": subtracting "<<hists[i]->GetBinContent(1)<<endl;
          clonedPrimary_norm->Add(hists[i], -1); //subtract
        }
      }
      zeroNegatives(clonedPrimary_norm);
    }

    //cout<<__LINE__<<endl;
    //double numEventsData; -- Now made more global, defined above
    double scaleFactor;
    if (conf->get("norm_0_50") == "true")
    {
      numEventsData = clonedPrimary_norm->Integral(clonedPrimary_norm->FindBin(0),clonedPrimary_norm->FindBin(49.9));
      numEventsMC = clonedBG_norm->IntegralAndError(clonedBG_norm->FindBin(0),clonedBG_norm->FindBin(49.9), errEventsMC);
      norm_bin = -1;
    }
    else if (conf->get("norm_50_100") == "true"){
      numEventsData = clonedPrimary_norm->Integral(clonedPrimary_norm->FindBin(50),clonedPrimary_norm->FindBin(99.9));
      numEventsMC = clonedBG_norm->IntegralAndError(clonedBG_norm->FindBin(50),clonedBG_norm->FindBin(99.9), errEventsMC);
      norm_bin = 0;
    }
    else{
      numEventsData = clonedPrimary_norm->Integral(0,-1);
      numEventsMC = clonedBG_norm->IntegralAndError(0,-1, errEventsMC);
      norm_bin = -1;
    }
    //cout<<__LINE__<<endl;
    cout<<"Num Events Primary: "<<numEventsData<<endl;
    cout<<"Num Events BG: "<<numEventsMC<<endl;
    //rescale everything to scale factor
    scaleFactor = ((double) numEventsData/numEventsMC);
    cout<<"Applying scale factor "<<scaleFactor<<" to MET templates hist."<<endl;
    for (int i = 1; i<num_hists; i++){
      if (hist_nums_for_norm.Contains(to_string(i)) || hist_nums_for_norm == "" ){
        cout<<i<<" "<<hist_labels[i]<<" SF: "<<scaleFactor<<endl;
        hists[i]->Scale(scaleFactor);  //if hist is marked for norm or no hists marked for norms.
        cout<<hist_labels[i]<<" count: "<<hists[i]->GetBinContent(1)<<endl;
      }
    }
  }

  cout<<"===========================================\nNorm\n===========================================\n";
  for (int i = 0; i<num_hists; i++){
    cout<<hist_labels[i]<<" Integral bin 0 to bin 100 Content: "<<hists[i]->Integral(hists[i]->FindBin(0),hists[i]->FindBin(100))<<endl;
  } 

  //Create sum of background samples
  TH1D *bg_sum = (TH1D*) hists[BG_sum_from]->Clone("bg_sum_"+plot_name);
  bg_sum->SetTitle("Sum of background samples");
  
  //cout<<__LINE__<<endl;
  for (int i=BG_sum_from+1; i<num_hists; i++){
    bg_sum->Add(hists[i]);
  }
  //cout<<__LINE__<<endl;

  for (int i = 0; i<num_hists; i++){
    cout<<hist_labels[i]<<": after-reweight "<<hists[i]->GetBinContent(1)<<endl;
  }

  delete clonedPrimary_norm;
  delete clonedBG_norm;
  //cout<<__LINE__<<endl;

  //===========================
  // SET MC COLORS
  //===========================
  //cout<<__LINE__<<endl;
  if(conf->get("WWW_colors") == "true"){
    for (int i = 1; i<num_hists; i++){
      hists[i]->SetFillColor(s_colors[sample_names[i]]);
      //cout<<__LINE__<<endl;
      hists[i]->SetFillStyle(1001);
    }
  }
  else{
    for (int i = 1; i<num_hists; i++){
      //cout<<__LINE__<<endl;
      hists[i]->SetFillColor(ROOT_COLOR_PALATE[(i-1) % ROOT_COLOR_PALATE.size()]);
      //cout<<__LINE__<<endl;
      hists[i]->SetFillStyle(1001);
    }
  }
  //cout<<__LINE__<<endl;
  hists[0]->SetMarkerStyle(20);
  //cout<<__LINE__<<endl;

  //===========================
  // Find Plot Maxima
  //===========================
  //cout<<__LINE__<<endl;
  double ymax = 0;
  double ymin = 0;
  TH1D* clonedBG = (TH1D*) bg_sum->Clone("clonedBG_forReweight_"+plot_name);
  TH1D* clonedPrimary = (TH1D*) hists[0]->Clone("clonedPrimary_forReweight_"+plot_name);
  //cout<<__LINE__<<endl;
  clonedPrimary->GetXaxis()->SetRangeUser(xmin, xmax);
  clonedBG->GetXaxis()->SetRangeUser(xmin,xmax);
  //cout<<__LINE__<<endl;
  if (conf->get("ymax") != ""){
    ymax = stod(conf->get("ymax"));
  }
  else{
    if (clonedBG->GetMaximum() < clonedPrimary->GetMaximum()){
        ymax = 1.2*clonedPrimary->GetMaximum();
    }
    else {
        ymax = 1.2*clonedBG->GetMaximum();   
    }
  }

  if (conf->get("logy") == "true"){
    ymax *= 10;
    ymin = 0.1;
  }

  if (conf->get("ymin") != ""){
    ymin = stod(conf->get("ymin"));
  }
  //cout<<__LINE__<<endl;
  cout<<"Primary Max: "<< clonedPrimary->GetMaximum() << " Secondary Max: "<< clonedBG->GetMaximum() <<endl;
  cout<<"Proper plot maximum set to "<<ymax<<endl;
  
  delete clonedBG;
  delete clonedPrimary;
  //cout<<__LINE__<<endl;
  
  TH2F* h_axes = new TH2F(Form("%s_axes",plot_name.Data()),plot_title,hists[0]->GetNbinsX(),xmin,xmax,1000,ymin,ymax);
  
  if (conf->get("relabel_x_axis") == "true"){
    TString bin_label;
    for (int i = (int) xmin; i <= (int) xmax; i++)
    {
      bin_label=hists[0]->GetXaxis()->GetBinLabel(hists[0]->FindBin(i));
      h_axes->GetXaxis()->SetBinLabel(h_axes->FindBin(i), bin_label);
    }  
    //h_axes->GetXaxis()->LabelsOption("v");
    h_axes->GetXaxis()->SetLabelSize(.015);
  }
  //-----------------------
  // AXES FIX
  //-----------------------
  
  cout<<"Setting axis names"<<endl;
  h_axes->GetXaxis()->SetTitle(xlabel);
  h_axes->GetYaxis()->SetTitle(ylabel);
  //cout<<__LINE__<<endl;  

  gStyle->SetTitleW(0.6);
  gStyle->SetTitleH(0.06);
  gStyle->SetTitleFont(12);

  //===========================
  // Print Closure Stats
  //===========================
  cout<<"===========================================\nPrint Stats\n===========================================\n";
  for (int i = 0; i<num_hists; i++){
    cout<<hist_labels[i]<<" Integral bin 0 to bin 100 Content: "<<hists[i]->Integral(hists[i]->FindBin(0),hists[i]->FindBin(100))<<endl;
  } 

  /*if (conf->get("print_stats") == "true"){
    vector<pair<double,double>> stats_bins;
    int j = 0;
    
    //cout<<__LINE__<<endl;

    while (conf->get("stats_"+to_string(j)+"_low_val") != "" ){
      stats_bins.push_back(make_pair(stod(conf->get("stats_"+to_string(j)+"_low_val")),stod(conf->get("stats_"+to_string(j)+"_high_val"))));
      j++;
    }

    //cout<<__LINE__<<endl;

    if(conf->get("simple_errors") == "true"){
      vector<vector<pair<double, double>>> stats; //holds a pair of count error for each sample, and the bg sum
      double count, error;
      vector<pair<double,double>> stat_row;
      
      //cout<<__LINE__<<endl;
      //Loop over the stats bins
      // Build Table ========================================================
      for(int i = 0 ; i < (int)hists.size(); i++){
        for(int st_bin=0; st_bin < (int) stats_bins.size(); st_bin++){
          //cout<<__LINE__<<endl;
          count = hists[i]->IntegralAndError(hists[i]->FindBin(stats_bins[st_bin].first), hists[i]->FindBin(stats_bins[st_bin].second - 0.000001), error);
          stat_row.push_back(make_pair(count,error));
        }
        stats.push_back(stat_row);
        stat_row.clear();
      }
      
      if (conf->get("templates_closure") == "true"){
        //Tack on Ratio row
        double zjets_count, zjets_err;
        for(int st_bin=0; st_bin < (int) stats_bins.size(); st_bin++){
          count = bg_sum->IntegralAndError(bg_sum->FindBin(stats_bins[st_bin].first), bg_sum->FindBin(stats_bins[st_bin].second - 0.000001), error);
          zjets_count = stats[0][st_bin].first;
          zjets_err = stats[0][st_bin].second;
          cout<<"zjets_count: "<<zjets_count<<" zjets_err: "<<zjets_err<<endl;
          cout<<"count: "<<count<<"err: "<<error<<endl;
          error = (1./count)*(sqrt( pow(zjets_err, 2) + pow(zjets_count * error / count, 2) ) );
          count =  zjets_count / count;
          stat_row.push_back(make_pair(count,error)); 
        } 
        stats.push_back(stat_row);
      }
      else {
        //Tack on BG sum row
        for(int st_bin=0; st_bin < (int) stats_bins.size(); st_bin++){
          count = bg_sum->IntegralAndError(bg_sum->FindBin(stats_bins[st_bin].first), bg_sum->FindBin(stats_bins[st_bin].second - 0.000001  ), error);
          stat_row.push_back(make_pair(count,error)); 
        } 
        stats.push_back(stat_row);
      }
      // End Table Building ==================================================

      // Print Table =========================================================
      CTable table;
      table.setPrecision(2);
      //Set Column Labels
      //cout<<__LINE__<<endl;
      table.setTitle(Form("Efficiencies for %s",plot_name.Data()));
      table.useTitle();
      for (int st_bin=0; st_bin < (int) stats_bins.size(); st_bin++){
        //cout<<__LINE__<<endl;
        table.setColLabel(Form("%.3f-%.3f",stats_bins[st_bin].first, stats_bins[st_bin].second), st_bin);
      }
      //cout<<__LINE__<<endl;

      
      if (conf->get("templates_closure") == "true"){
        //Output Rows for samples
        for(int row = 0; row <= (int) hists.size(); row++ ){
          if (row == hists.size()){
            table.setRowLabel("Ratio", hists.size());  
          }
          else{
            table.setRowLabel(hist_labels[row], row);
          }
          for(int col=0; col < (int) stats_bins.size(); col++){
            //cout<<__LINE__<<endl;
            table.setCell(Form("%.2f+/-%.2f", stats[row][col].first, stats[row][col].second), row, col);
          }
        }
      }
      else{
        //Output Rows for samples
        for(int row = 0; row <= (int) hists.size(); row++ ){
          if (row == hists.size()){
            table.setRowLabel("Sum of BG", hists.size());  
          }
          else{
            table.setRowLabel(hist_labels[row], row);
          }
          for(int col=0; col < (int) stats_bins.size(); col++){
            //cout<<__LINE__<<endl;
            table.setCell(Form("%.2f+/-%.2f; Eff: %.2f", stats[row][col].first, stats[row][col].second, stats[row][col].first/stats[row][0].first), row, col);
          }
        }
      }
      //cout<<__LINE__<<endl;

      table.print();
      TString conf_id = conf->get("conf_path");
      conf_id.ReplaceAll("//","/");
      conf_id.ReplaceAll("/","_");
      conf_id = conf_id(8, conf_id.Last('.')-8);

      table.saveTex(Form("outputs/efficiency_table_%s.tex", conf_id.Data()));

    }
  }*/
  

  cout<<"===========================================\nUpdate Overflow\n===========================================\n";
  //Done stats after so that we don't double count these since they are not deleted from final bin.
  for (int i = 0; i<num_hists; i++){
    cout<<hist_labels[i]<<" Integral bin 0 to bin 100 Content: "<<hists[i]->Integral(hists[i]->FindBin(0),hists[i]->FindBin(100))<<endl;
  } 

  //----------------------
  // ADD OVERFLOW BIN
  //----------------------
  //cout<<__LINE__<<endl;
  if (conf->get("overflow")=="true"){
    cout<<"Plot tagged for overflow bin, building..."<<endl;
    updateOverUnderflow(bg_sum, xmax);
    for (int i = 0; i<num_hists; i++){
      updateOverUnderflow(hists[i], xmax);
    }
    /*double n_bins = hists[0]->GetNbinsX();
    double overflow, max;
    //cout<<__LINE__<<endl;
    for (int i = 0; i<num_hists; i++){
      overflow = hists[i]->GetBinContent(n_bins + 1);
      max = hists[i]->Integral(hists[i]->FindBin(xmax-.001), n_bins);
      hists[i]->SetBinContent(hists[i]->FindBin(xmax-.001), max+overflow);
    }
    overflow = bg_sum->GetBinContent(n_bins + 1);
    max = bg_sum->Integral(bg_sum->FindBin(xmax-.001), n_bins);
    bg_sum->SetBinContent(bg_sum->FindBin(xmax-.001), max+overflow);*/
    //cout<<__LINE__<<endl;
  }
  
      
  
  plotpad->SetLeftMargin(0.15);
  h_axes->GetYaxis()->SetTitleOffset(1.3);
  h_axes->GetYaxis()->SetTitleSize(0.05);
  h_axes->GetYaxis()->SetLabelSize(0.04);
  h_axes->GetXaxis()->SetTitleSize(0.05);
  h_axes->GetXaxis()->SetLabelSize(0.04);
  //cout<<__LINE__<<endl;
  cout<<"Drawing histograms"<<endl;
  h_axes->Draw();

   //==========================
  // Hist Groups
  //==========================

  //cout<<"---------------"<<endl;
  //for (int i = 0; i<num_hists; i++){
  //  cout<<"hist "<<i<<" "<<hist_labels[i]<<endl;
  //}
  //cout<<"---------------"<<endl;

  //Hold list of hists to be removed from list since they will be added to the first hist
  //in the group.
  vector<int> remove_after_groups;

  //loop through the conf group_i's
  int num_group = 0;
  while (conf->get("group_"+to_string(num_group)) != ""){
    //get list of samples in the group
    vector<int> hnums = iparseVector(conf->get("group_"+to_string(num_group)));

    cout<<"Making Group "<<num_group<<" elements are:"<<hist_labels[hnums[0]];

    //change the label of the first hist to the name of the group
    hist_labels[hnums[0]] = TString(conf->get("group_"+to_string(num_group)+"_label"));
    for (vector<int>::iterator n = hnums.begin()+1; n != hnums.end(); n++){
      //Add other hists in the group to the first one
      hists[hnums[0]]->Add(hists[*n]);
      //Store a list of the hists added so they can be removed.
      remove_after_groups.push_back(*n);
      cout<<", "<<hist_labels[*n];
    }
    cout<<"."<<endl;
    num_group++;
  }

  //Sort the list of hists to be removed so that they can be traversed backwards and we don't need to keep track of 
  //hist numbers changing
  sort(remove_after_groups.begin(), remove_after_groups.end());

  //Remove hists from the list of hists that have been grouped into another spot. Going through backwards ensures the 
  //indexes don't change.
  for (int i = remove_after_groups.size()-1; i>=0; i--){
    cout<<"Removing "<<hist_labels[remove_after_groups[i]]<<" since it was part of a group"<<endl;
    hists.erase(hists.begin()+remove_after_groups[i]);
    hist_labels.erase(hist_labels.begin()+remove_after_groups[i]);
    num_hists--;
  }


  //cout<<"---------------"<<endl;
  //for (int i = 0; i<num_hists; i++){
  //  cout<<"hist "<<i<<" "<<hist_labels[i]<<endl;
  //}
  //cout<<"---------------"<<endl;

  //===========================
  // MAKE STACK
  //===========================
  //Add all the background hists to a stack.
  THStack * stack = new THStack(("stack_"+conf->get("Name")).c_str(), conf->get("title").c_str());
  //cout<<__LINE__<<endl;

  
  vector<pair<TH1D*, TString>> hists_labeled; 

  for (int i=1; i<num_hists; i++)
  {
    hists_labeled.push_back(make_pair(hists[i], hist_labels[i]));
  } 

  cout<<"Sorting by total count"<<endl;

  sort(hists_labeled.begin(), hists_labeled.end(), LabeledHistSort);

  cout<<"Sorted, Adding to Stack"<<endl;

  //sort(hists.begin()+1, hists.end(), TH1DIntegralSort); 
  for (int i=0; i< (int)hists_labeled.size(); i++)
  {
    stack->Add(hists_labeled[i].first);
  } 
  stack->Draw("HIST SAME");
  cout<<"Stack Drawn"<<endl;
  if (conf->get("blindAfter") != ""){
    blindAfter(hists[0], stod(conf->get("blindAfter")));
  }
  TGraphAsymmErrors *bg_err = new TGraphAsymmErrors(bg_sum);

  if ((conf->get("print_stats") == "true") && (conf->get("simple_errors") != "true")){
    prediction_errors->SetFillStyle(3244);
    prediction_errors->SetFillColor(kGray+3);
    prediction_errors->Draw("SAME 2");
  }
  else if(conf->get("draw_bg_errs") != "false"){
    cout<<"Drawing BG hatch bands"<<endl;
    bg_err->SetFillStyle(3244);
    bg_err->SetFillColor(kGray+3);
    bg_err->Draw("SAME 2");
  }
  if (conf->get("print_stats") == "true" && conf->get("simple_errors") != "true"){
    hists[0]->SetMarkerSize(2.5);
  }
  else{
    hists[0]->SetMarkerSize(1.5);
  }
  hists[0]->Draw("same e0 x0 e1 p0");

  plotpad->RedrawAxis();
  //cout<<__LINE__<<endl;
  //===========================
  // BUILD LEGEND
  //===========================

  cout<<"Building Legend"<<endl;

  TLegend *l1;
  if (conf->get("small_legend") == "true"){
    l1 = new TLegend(0.78, 0.78, 0.93, 0.93);
  }
  else if (conf->get("horizontal_legend") == "true"){
    short n_cols = short(ceil(num_hists/3.0));
    //cout<<"n_cols "<<n_cols<<endl;
    l1 = new TLegend(0.151, 0.8, 0.95, 0.949); 
    l1->SetNColumns(n_cols);
  }
  else{
    /*cout<<"UtoPixel(0.65): "<<gPad->UtoPixel(.65)<<endl;
    cout<<"PixelToX(..): "<<gPad->PixeltoX(gPad->UtoPixel(.65))<<endl;
    double x_under_legend = gPad->PixeltoX(gPad->UtoPixel(.65));
    double max_count_under_legend = bg_sum->GetBinContent(bg_sum->FindBin(x_under_legend));
    
    cout<<"x under legend: "<<x_under_legend<<endl;
    cout<<"max count under legend: "<<max_count_under_legend<<endl;*/
    l1 = new TLegend(0.65, 0.6, 0.93, 0.93);
  }
  
  l1->SetLineColor(kWhite);  
  l1->SetShadowColor(kWhite);
  l1->SetFillColor(kWhite);
  l1->SetTextSize(.03);
  cout<<__LINE__<<endl;
  l1->AddEntry(hists[0], hist_labels[0], "pe");
  /* //Put objects in legend with the same order as the they go into the stack
  for (int i = hists_labeled.size()-1; i>=0; i--){
    l1->AddEntry(hists_labeled[i].first, hists_labeled[i].second, "f");
  }*/
  //Put objects in legend in the order they are written in the config
  for (int i = hists.size()-1; i>=1; i--){
    l1->AddEntry(hists[i], hist_labels[i], "f");
  }

  l1->Draw("same");
  
  //--------------------------
  // Fill in Residual Plot
  //--------------------------
  
  cout<<"Getting ready for residual plots"<<endl;
  fullpad->cd();
  TPad *ratiopad = new TPad("ratiopad", "ratiopad" ,0.,0.,1,0.21);
  ratiopad->SetTopMargin(0.05);
  ratiopad->SetLeftMargin(0.15);
  ratiopad->SetBottomMargin(0.1);
  ratiopad->SetRightMargin(0.05);
  ratiopad->SetGridy();  // doesn't actually appear for some reason..
  ratiopad->Draw();
  ratiopad->cd();
  //cout<<__LINE__<<endl;
  
  TH1D* residual = (TH1D*) hists[0]->Clone("residual");
  residual->Divide(bg_sum);
  //cout<<__LINE__<<endl;
  //cout<<"Fixing error bars"<<endl;
  //for (int count=1; count<=mc_sum->GetNbinsX(); count++){ 
  //  double relative_error = (mc_sum->GetBinError(count))/ (mc_sum->GetBinContent(count));
  //  residual->SetBinError(count, residual->GetBinContent(count)*relative_error);
  //}
  
  cout<<"Building axes"<<endl;
  TH1D* h_axis_ratio = new TH1D(Form("%s_residual_axes",plot_name.Data()),"",residual->GetNbinsX(),xmin,xmax);
  //cout<<__LINE__<<endl;
  h_axis_ratio->GetYaxis()->SetTitleOffset(0.33);
  h_axis_ratio->GetYaxis()->SetTitleSize(0.16);
  //h_axis_ratio->GetYaxis()->SetTitleFont(12);
  h_axis_ratio->GetYaxis()->SetNdivisions(5);
  h_axis_ratio->GetYaxis()->SetLabelSize(0.15);
  //h_axis_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  
  if(residual->GetMaximum() > 3 && residual->GetMaximum() <= 4 ){
    h_axis_ratio->GetYaxis()->SetRangeUser(0.001,4.0);
  }
  else if(residual->GetMaximum() > 2 && residual->GetMaximum() <= 3 ){
    h_axis_ratio->GetYaxis()->SetRangeUser(0.001,3.0);
  }
  else{
    h_axis_ratio->GetYaxis()->SetRangeUser(0.001,2.0);
  }


  if(conf->get("ratio_yaxis") != ""){
    h_axis_ratio->GetYaxis()->SetTitle(parseLatex(conf->get("ratio_yaxis")));  
  }
  else{
    h_axis_ratio->GetYaxis()->SetTitle("#frac{Signal}{Prediction}");
  }

  h_axis_ratio->GetXaxis()->SetTickLength(0.07);
  h_axis_ratio->GetXaxis()->SetTitleSize(0.005);
  h_axis_ratio->GetXaxis()->SetLabelSize(0.);

  h_axis_ratio->GetYaxis()->CenterTitle();
  
  //cout<<__LINE__<<endl;
  TLine* line1 = new TLine(xmin,1,xmax,1);
  line1->SetLineStyle(2);

  if (conf->get("print_stats") == "true" && conf->get("simple_errors") != "true"){
    residual->SetMarkerSize(2.5);
  }
  else{
    residual->SetMarkerSize(1.5);
  }
  
  
  cout<<"Drawing ratio plot"<<endl;
  h_axis_ratio->Draw("axis");
  //cout<<__LINE__<<endl;
  line1->Draw("same");
  //cout<<__LINE__<<endl;
  residual->Draw("same e0 x0 e1 p0");
  //cout<<__LINE__<<endl;
  c->Update();
  //cout<<__LINE__<<endl;
  c->cd();
  //cout<<__LINE__<<endl;
  
  //Draw luminosity and CMS tag
  if (conf->get("luminosity_fb") != ""){
    plotpad->cd();
    drawCMSLatex(stod(conf->get("luminosity_fb")));
  }
  if (conf->get("SR") != ""){
    plotpad->cd();
    drawSRText(conf->get("SR"), bg_sum->GetMaximum(), xmin);
  }
  //cout<<__LINE__<<endl;

  cout<<"Saving..."<<endl;
  c->SaveAs(save_dir+plot_name+TString(".pdf"));
  c->SaveAs(save_dir+plot_name+TString(".png"));
  //c->SaveAs(save_dir+plot_name+TString(".root"));
  //c->SaveAs(save_dir+plot_name+TString(".C"));
  //cout<<__LINE__<<endl;
  cout<<"Cleaning up plot variables"<<endl;
  delete l1;
  hists_labeled.clear();
  hists.clear();
  hist_names.clear();
  hist_labels.clear();
  delete residual;
  delete ratiopad;
  delete plotpad;
  delete fullpad;
  delete bg_err;
  delete c;
  //cout<<__LINE__<<endl;
  for (int i = 0; i<num_hists; i++){
    hist_files[i]->Close();
  }
  hist_files.clear();
  //cout<<__LINE__<<endl;
  return errors;
}

TString drawArbitraryNumber(ConfigParser *conf){
  // This method expects conf to have a plot config loaded in already. 
  //In the conf, we expect there to be hist names of the form file_N_path,
  //hist_n_name, starting with 0 for the primary histogram, which is normally 
  //going to be the data events in our signal region. The rest of the hists, starting
  //from 1, are added to a THStack which is normalized to hist_0 in the bin 0-50. 
  //num_hists should be the number of the number of histograms in the plot.
  TString errors="";

  int num_hists=stoi(conf->get("num_hists"));
  int BG_sum_from=(conf->get("BG_sum_from") != "") ? stoi(conf->get("BG_sum_from")) : 0;

  if (num_hists < 2){
    return TString("Less than Two hists can not be turned into a stack plot, please call drawSingleTH1 (replace config PLOT_TYPE with Single)");
  } 

  //Add files from which to obtain histos
  TString default_hist_dir = getDefaultHistDir(conf);
  vector<TFile*> hist_files (num_hists);
  for (int i = 0; i<num_hists; i++){
    TString sample_loc = "";
    if (conf->get("file_"+to_string(i)+"_path") != ""){
      sample_loc = TString(conf->get("file_"+to_string(i)+"_path"));
    }
    else{
      sample_loc = default_hist_dir+conf->get("sample_"+to_string(i))+".root";
    }

    hist_files[i]=new TFile(sample_loc);
  }
  cout << "Found files "<<endl;

  TString plot_name = conf->get("plot_name");
  double bin_size;
  if (conf->get("bin_size") != ""){
    bin_size = stod(conf->get("bin_size"));
  }
  else{
    bin_size=1;
  }
  //Get name of hist to read from file
  vector<TString> hist_names (num_hists);
  for (int i = 0; i<num_hists; i++){
    if (conf->get("hist_"+to_string(i)+"_name") != ""){
      hist_names[i]=conf->get("hist_"+to_string(i)+"_name");    
    }
    else{
      hist_names[i]= conf->get("hist_0_name");
    }
  }

  //Get labels for TLegend
  vector<TString> hist_labels (num_hists);
  for (int i = 0; i<num_hists; i++){
    hist_labels[i]=parseLatex(conf->get("hist_"+to_string(i)+"_label"));    
  }  

  //Get sample names for Colors
  vector<TString> sample_names (num_hists);
  for (int i = 0; i<num_hists; i++){
    sample_names[i]=conf->get("sample_"+to_string(i));    
  }  


  cout<<"Hist names set"<<endl;
  TString xlabel=parseLatex(conf->get("xlabel"));
  TString ylabel=parseLatex(conf->get("ylabel"));
  TString save_dir=(conf->get("save_dir") != "") ? conf->get("save_dir") : getOutputDir(conf, "plot");
  TString plot_title=parseLatex(conf->get("title"));


  cout << "Making Plots for: "<<plot_name<<endl;

  vector<TH1D*> hists (num_hists);
  for (int i = 0; i<num_hists; i++){
    hists[i] = (TH1D*) ((TH1D*) hist_files[i]->Get(hist_names[i]))->Clone("hist_"+to_string(i)+"_"+plot_name);
    cout<<hist_names[i]<<" found in "<<hist_files[i]->GetName()<<endl;
    cout<<hist_labels[i]<<" Integral bin 0 to bin 100 Content: "<<hists[i]->Integral(0,100)<<endl;
  }  

  cout << "Histograms pulled from files, adding draw options"<<endl;
  
  //cout<<__LINE__<<endl;

  if (conf->get("relabel_x_axis") == "true"){
    hists[0]->LabelsDeflate();
  }

  double xmax = (conf->get("xmax") != "" ) ? stod(conf->get("xmax")) : hists[0]->GetNbinsX();
  double xmin = (conf->get("xmin") != "" ) ? stod(conf->get("xmin")) : 0;


  //============================================
  // Draw Data-MC Plots
  //============================================
  
  //cout<<__LINE__<<endl;
  TCanvas * c = new TCanvas("c","",2000,2000);
  c->cd();
  gPad->SetRightMargin(0.05);
  gPad->Modified();
  gStyle->SetOptStat(kFALSE);
  TPad *fullpad = new TPad("fullpad", "fullpad", 0,0,1,1);
  //cout<<__LINE__<<endl;
  fullpad->Draw();
  fullpad->cd();
  //cout<<__LINE__<<endl;

  fullpad->SetRightMargin(0.05);
  if (conf->get("ExtraRightMargin") == "true")
  {
    fullpad->SetRightMargin(0.08);
  }
  fullpad->SetBottomMargin(0.12);
  //cout<<__LINE__<<endl;
  fullpad->Draw();
  fullpad->cd();
  //cout<<__LINE__<<endl;
  if (conf->get("logy") == "true")
  {
    cout<<"Plot tagged for log y-axis"<<endl;
    fullpad->SetLogy();
  }
  //cout<<__LINE__<<endl;
  if (conf->get("bin_size") != ""){
    for (int i = 0; i<num_hists; i++){
      hists[i]->Rebin(bin_size);
    }
  }
  else if (conf->get("binning") != ""){
    cout<<"Rebinning plots with variable ranges"<<endl;
    vector<double> binning = parseVector(conf->get("binning"));
    for (int i = 0; i<num_hists; i++){
      hists[i] = (TH1D*) hists[i]->Rebin(binning.size()-1, TString(hist_names[i]+"_rebin"), &binning[0]);
    }
  }

  //===========================
  // Normalize
  //===========================
  //cout<<__LINE__<<endl;
  
  //Add scale factors like RSFOF
  for (int i=0; i < num_hists; i++){
    if (conf->get("hist_"+to_string(i)+"_scale") != ""){
      hists[i]->Scale(stod(conf->get("hist_"+to_string(i)+"_scale")));
    }
  }

  //Create sum of background samples
  TH1D *bg_sum = (TH1D*) hists[BG_sum_from]->Clone("bg_sum_"+plot_name);
  bg_sum->SetTitle("Sum of background samples");
  
  //cout<<__LINE__<<endl;
  for (int i=BG_sum_from+1; i<num_hists; i++){
    bg_sum->Add(hists[i]);
  }
  //cout<<__LINE__<<endl;

  for (int i = 0; i<num_hists; i++){
    cout<<hist_labels[i]<<": after-reweight "<<hists[i]->GetBinContent(1)<<endl;
  }

  cout<<"===========================================\n2\n===========================================\n";
  for (int i = 0; i<num_hists; i++){
    cout<<hist_labels[i]<<" Integral bin 0 to bin 100 Content: "<<hists[i]->Integral(0,100)<<endl;
  } 

  //===========================
  // SET MC COLORS
  //===========================
  //cout<<__LINE__<<endl;
  if(conf->get("WWW_colors") == "true"){
    for (int i = 0; i<num_hists; i++){
      hists[i]->SetFillColor(s_colors[sample_names[i]]);
      //cout<<__LINE__<<endl;
      hists[i]->SetFillStyle(1001);
    }
  }
  else{
    for (int i = 0; i<num_hists; i++){
      //cout<<__LINE__<<endl;
      hists[i]->SetFillColor(ROOT_COLOR_PALATE[(i-1) % ROOT_COLOR_PALATE.size()]);
      //cout<<__LINE__<<endl;
      hists[i]->SetFillStyle(1001);
    }
  }
  //cout<<__LINE__<<endl;

  //===========================
  // Find Plot Maxima
  //===========================
  //cout<<__LINE__<<endl;
  double ymax = 0;
  double ymin = 0;
  TH1D* clonedBG = (TH1D*) bg_sum->Clone("clonedBG_forReweight_"+plot_name);
  //cout<<__LINE__<<endl;
  clonedBG->GetXaxis()->SetRangeUser(xmin,xmax);
  //cout<<__LINE__<<endl;
  if (conf->get("ymax") != ""){
    ymax = stod(conf->get("ymax"));
  }
  else{
      ymax = 1.2*clonedBG->GetMaximum();   
  }
  if (conf->get("logy") == "true"){
      ymax *= 10;
      ymin = 0.1;
  }

  if (conf->get("ymin") != ""){
    ymin = stod(conf->get("ymin"));
  }

  cout<<"Proper plot maximum set to "<<ymax<<endl;
  
  delete clonedBG;
  //cout<<__LINE__<<endl;
  
  TH2F* h_axes = new TH2F(Form("%s_axes",plot_name.Data()),plot_title,hists[0]->GetNbinsX(),xmin,xmax,1000,ymin,ymax);
  
  if (conf->get("relabel_x_axis") == "true"){
    TString bin_label;
    for (int i = (int) xmin; i <= (int) xmax; i++)
    {
      bin_label=hists[0]->GetXaxis()->GetBinLabel(hists[0]->FindBin(i));
      h_axes->GetXaxis()->SetBinLabel(h_axes->FindBin(i), bin_label);
    }  
    //h_axes->GetXaxis()->LabelsOption("v");
    h_axes->GetXaxis()->SetLabelSize(.015);
  }
  //-----------------------
  // AXES FIX
  //-----------------------
  
  cout<<"Setting axis names"<<endl;
  h_axes->GetXaxis()->SetTitle(xlabel);
  h_axes->GetYaxis()->SetTitle(ylabel);
  //cout<<__LINE__<<endl;  

  gStyle->SetTitleW(0.6);
  gStyle->SetTitleH(0.06);
  gStyle->SetTitleFont(12);

  //===========================
  // Print Closure Stats
  //===========================


  /*if (conf->get("print_stats") == "true"){
    vector<pair<double,double>> stats_bins;
    int j = 0;
    
    while (conf->get("stats_"+to_string(j)+"_low_val") != "" ){
      stats_bins.push_back(make_pair(stod(conf->get("stats_"+to_string(j)+"_low_val")),stod(conf->get("stats_"+to_string(j)+"_high_val"))));
      j++;
    }

    if(conf->get("simple_errors") == "true"){
      cout<<"Computing simple errors"<<endl;
      vector<vector<pair<double, double>>> stats; //holds a pair of count error for each sample, and the bg sum
      double count, error;
      vector<pair<double,double>> stat_row;
      
      //cout<<__LINE__<<endl;
      //Loop over the stats bins
      // Build Table ========================================================
      for(int i = 0 ; i < (int)hists.size(); i++){
        for(int st_bin=0; st_bin < (int) stats_bins.size(); st_bin++){
          //cout<<__LINE__<<endl;
          count = hists[i]->IntegralAndError(hists[i]->FindBin(stats_bins[st_bin].first), hists[i]->FindBin(stats_bins[st_bin].second - 0.01), error);
          stat_row.push_back(make_pair(count,error));
        }
        stats.push_back(stat_row);
        stat_row.clear();
      }
      //Tack on BG sum row
      for(int st_bin=0; st_bin < (int) stats_bins.size(); st_bin++){
        count = bg_sum->IntegralAndError(bg_sum->FindBin(stats_bins[st_bin].first), bg_sum->FindBin(stats_bins[st_bin].second - 0.01), error);
        stat_row.push_back(make_pair(count,error)); 
      } 
      stats.push_back(stat_row);
      // End Table Building ==================================================

      // Print Table =========================================================
      CTable table;
      table.setPrecision(2);
      //Set Column Labels
      //cout<<__LINE__<<endl;
      table.setTitle(Form("Efficiencies for %s",plot_name.Data()));
      table.useTitle();
      for (int st_bin=0; st_bin < (int) stats_bins.size(); st_bin++){
        //cout<<__LINE__<<endl;
        table.setColLabel(Form("%.2f-%.2f",stats_bins[st_bin].first, stats_bins[st_bin].second), st_bin);
      }
      //cout<<__LINE__<<endl;

      //Output Rows for samples
      for(int row = 0; row <= (int) hists.size(); row++ ){
        if (row == hists.size()){
          table.setRowLabel("Sum of BG", hists.size());  
        }
        else{
          table.setRowLabel(hist_labels[row], row);
        }
        for(int col=0; col < (int) stats_bins.size(); col++){
          //cout<<__LINE__<<endl;
          table.setCell(Form("%.2f+/-%.2f; Eff: %.2f", stats[row][col].first, stats[row][col].second, stats[row][col].first/stats[row][0].first), row, col);
        }
      }
      //cout<<__LINE__<<endl;

      table.print();
      TString conf_id = conf->get("conf_path");
      conf_id.ReplaceAll("//","/");
      conf_id.ReplaceAll("/","_");
      conf_id = conf_id(8, conf_id.Last('.')-8);

      table.saveTex(Form("outputs/efficiency_table_%s.tex", conf_id.Data()));

    }
  }*/
  
  cout<<"===========================================\n3\n===========================================\n";
  for (int i = 0; i<num_hists; i++){
    cout<<hist_labels[i]<<" Integral bin 0 to bin 100 Content: "<<hists[i]->Integral(0,100)<<endl;
  } 

  //----------------------
  // ADD OVERFLOW BIN
  //----------------------
  //cout<<__LINE__<<endl;
  if (conf->get("overflow")=="true"){
    cout<<"Plot tagged for overflow bin, building..."<<endl;

    updateOverUnderflow(bg_sum, xmax);
    for (int i = 0; i<num_hists; i++){
      updateOverUnderflow(hists[i], xmax);
    }
    
    /*double n_bins = hists[0]->GetNbinsX();
    double overflow, max;
    //cout<<__LINE__<<endl;
    for (int i = 0; i<num_hists; i++){
      overflow = hists[i]->GetBinContent(n_bins + 1);
      max = hists[i]->Integral(hists[i]->FindBin(xmax-.001), n_bins);
      hists[i]->SetBinContent(hists[i]->FindBin(xmax-.001), max+overflow);
    }
    overflow = bg_sum->GetBinContent(n_bins + 1);
    max = bg_sum->Integral(bg_sum->FindBin(xmax-.001), n_bins);
    bg_sum->SetBinContent(bg_sum->FindBin(xmax-.001), max+overflow);*/
    //cout<<__LINE__<<endl;
  }
  
      
  
  fullpad->SetLeftMargin(0.15);
  h_axes->GetYaxis()->SetTitleOffset(1.3);
  h_axes->GetYaxis()->SetTitleSize(0.05);
  h_axes->GetYaxis()->SetLabelSize(0.04);
  //cout<<__LINE__<<endl;
  cout<<"Drawing histograms"<<endl;
  h_axes->Draw();

  //==========================
  // Hist Groups
  //==========================

  //cout<<"---------------"<<endl;
  //for (int i = 0; i<num_hists; i++){
  //  cout<<"hist "<<i<<" "<<hist_labels[i]<<endl;
  //}
  //cout<<"---------------"<<endl;

  //Hold list of hists to be removed from list since they will be added to the first hist
  //in the group.
  vector<int> remove_after_groups;

  //loop through the conf group_i's
  int num_group = 0;
  while (conf->get("group_"+to_string(num_group)) != ""){
    //get list of samples in the group
    vector<int> hnums = iparseVector(conf->get("group_"+to_string(num_group)));

    cout<<"Making Group "<<num_group<<" elements are:"<<hist_labels[hnums[0]];

    //change the label of the first hist to the name of the group
    hist_labels[hnums[0]] = TString(conf->get("group_"+to_string(num_group)+"_label"));
    for (vector<int>::iterator n = hnums.begin()+1; n != hnums.end(); n++){
      //Add other hists in the group to the first one
      hists[hnums[0]]->Add(hists[*n]);
      //Store a list of the hists added so they can be removed.
      remove_after_groups.push_back(*n);
      cout<<", "<<hist_labels[*n];
    }
    cout<<"."<<endl;
    num_group++;
  }

  //Sort the list of hists to be removed so that they can be traversed backwards and we don't need to keep track of 
  //hist numbers changing
  sort(remove_after_groups.begin(), remove_after_groups.end());

  //Remove hists from the list of hists that have been grouped into another spot. Going through backwards ensures the 
  //indexes don't change.
  for (int i = remove_after_groups.size()-1; i>=0; i--){
    cout<<"Removing "<<hist_labels[remove_after_groups[i]]<<" since it was part of a group"<<endl;
    hists.erase(hists.begin()+remove_after_groups[i]);
    hist_labels.erase(hist_labels.begin()+remove_after_groups[i]);
    num_hists--;
  }


  //cout<<"---------------"<<endl;
  //for (int i = 0; i<num_hists; i++){
  //  cout<<"hist "<<i<<" "<<hist_labels[i]<<endl;
  //}
  //cout<<"---------------"<<endl;

  //===========================
  // MAKE STACK
  //===========================
  //Add all the background hists to a stack.
  THStack * stack = new THStack(("stack_"+conf->get("Name")).c_str(), conf->get("title").c_str());
  //cout<<__LINE__<<endl;

  vector<pair<TH1D*, TString>> hists_labeled; 

  for (int i=0; i<num_hists; i++)
  {
    hists_labeled.push_back(make_pair(hists[i], hist_labels[i]));
  } 

  sort(hists_labeled.begin(), hists_labeled.end(), LabeledHistSort);
  
  for (int i=0; i<num_hists; i++)
  {
    stack->Add(hists_labeled[i].first);
  } 
  stack->Draw("HIST");
  
  TGraphAsymmErrors *bg_err = new TGraphAsymmErrors(bg_sum);
  
  if(conf->get("draw_bg_errs") != "false"){
    cout<<"Drawing BG hatch bands"<<endl;
    bg_err->SetFillStyle(3244);
    bg_err->SetFillColor(kGray+3);
    bg_err->Draw("SAME 2");
  }

  fullpad->RedrawAxis();
  //cout<<__LINE__<<endl;

  //===========================
  // BUILD LEGEND
  //===========================

  cout<<"Building Legend"<<endl;
  TLegend *l1;
  if (conf->get("small_legend") == "true"){
    l1 = new TLegend(0.78, 0.78, 0.93, 0.93);
  }
  else if (conf->get("horizontal_legend") == "true"){
    short n_cols = short(ceil(num_hists/3.0));
    //cout<<"n_cols "<<n_cols<<endl;
    l1 = new TLegend(0.151, 0.8, 0.95, 0.949); 
    l1->SetNColumns(n_cols);
  }
  else{
    /*cout<<"UtoPixel(0.65): "<<gPad->UtoPixel(.65)<<endl;
    cout<<"PixelToX(..): "<<gPad->PixeltoX(gPad->UtoPixel(.65))<<endl;
    double x_under_legend = gPad->PixeltoX(gPad->UtoPixel(.65));
    double max_count_under_legend = bg_sum->GetBinContent(bg_sum->FindBin(x_under_legend));
    
    cout<<"x under legend: "<<x_under_legend<<endl;
    cout<<"max count under legend: "<<max_count_under_legend<<endl;*/
    l1 = new TLegend(0.65, 0.6, 0.93, 0.93);
  }
  
  
  l1->SetLineColor(kWhite);  
  l1->SetShadowColor(kWhite);
  l1->SetFillColor(kWhite);
  l1->SetTextSize(.03);
  //cout<<__LINE__<<endl;
  //Put objects in legend in the order they are written in the config
  for (int i = hists.size(); i>=0; i--){
    l1->AddEntry(hists[i], hist_labels[i], "f");
  }

  l1->Draw("same");
  cout<<"Drawing Latex String"<<endl;
 
  fullpad->cd();
  //Draw luminosity and CMS tag
  if (conf->get("luminosity_fb") != ""){
    drawCMSLatex(stod(conf->get("luminosity_fb")));
  }
  //cout<<__LINE__<<endl;

  cout<<"Saving..."<<endl;
  c->SaveAs(save_dir+plot_name+TString(".pdf"));
  c->SaveAs(save_dir+plot_name+TString(".png"));
  //c->SaveAs(save_dir+plot_name+TString(".root"));
  //c->SaveAs(save_dir+plot_name+TString(".C"));
  //cout<<__LINE__<<endl;
  cout<<"Cleaning up plot variables"<<endl;
  delete l1;
  hists_labeled.clear();
  hists.clear();
  hist_names.clear();
  hist_labels.clear();
  delete bg_err;
  delete fullpad;
  delete c;
  //cout<<__LINE__<<endl;
  for (int i = 0; i<num_hists; i++){
    hist_files[i]->Close();
  }
  hist_files.clear();
  //cout<<__LINE__<<endl;
  return errors;
}

TString drawSingleTH1(ConfigParser *conf){
  /* This method expects conf to have a plot config loaded in already. */
  TString errors="";

  TString sample_loc;
  if (conf->get("file_path") != ""){
      sample_loc = TString(conf->get("file_path"));
  }
  else{
    sample_loc = TString(getDefaultHistDir(conf)+conf->get("sample")+".root");
  }

  TFile *f_primary =new TFile(sample_loc);

  cout << "Found files "<<endl;

  TString plot_name = conf->get("plot_name");
  TString plot_title = parseLatex(conf->get("title"));
  TString hist_name=conf->get("hist_name");
  TString xlabel=parseLatex(conf->get("xlabel"));
  TString ylabel=parseLatex(conf->get("ylabel"));
  TString save_dir=(conf->get("save_dir") != "") ? conf->get("save_dir") : getOutputDir(conf, "plot");

  double bin_size;
  if (conf->get("bin_size") != ""){
    bin_size = stod(conf->get("bin_size"));
  }
  else{
    bin_size=1;
  }

  cout << "Making Plots for: "<<plot_name<<endl;

  TH1D* p_hist = (TH1D*) ((TH1D*) f_primary->Get(hist_name))->Clone("phist_"+plot_name);
  cout<<hist_name<<" found in "<<f_primary->GetName()<<endl;

  cout << "Histograms pulled from files, adding draw options"<<endl;
  
  if (conf->get("relabel_x_axis") == "true"){
    p_hist->LabelsDeflate();
  }

  double xmax = (conf->get("xmax") != "" ) ? stod(conf->get("xmax")) : p_hist->GetNbinsX();
  double xmin = (conf->get("xmin") != "" ) ? stod(conf->get("xmin")) : 0;

  //============================================
  // Draw Data-MC Plots
  //============================================
  
  TCanvas * c = new TCanvas("c","",2000,2000);
  c->cd();
  gPad->SetRightMargin(0.05);
  gPad->Modified();
  gStyle->SetOptStat(kFALSE);
  TPad *fullpad = new TPad("fullpad", "fullpad", 0,0,1,1);
  
  fullpad->Draw();
  fullpad->cd();
    
  fullpad->SetRightMargin(0.05);
  if (conf->get("ExtraRightMargin") == "true")
  {
    fullpad->SetRightMargin(0.08);
  }
  fullpad->SetBottomMargin(0.12);
  
  fullpad->Draw();
  fullpad->cd();
  
  if (conf->get("logy") == "true")
  {
    cout<<"Plot tagged for log y-axis"<<endl;
    fullpad->SetLogy();
  }
  

  if (conf->get("bin_size") != ""){
    p_hist->Rebin(bin_size);
  }
  else if (conf->get("binning") != ""){
    cout<<"Rebinning plots with variable ranges"<<endl;
    vector<double> binning = parseVector(conf->get("binning"));
    p_hist = (TH1D*) p_hist->Rebin(binning.size()-1, TString(hist_name+"_rebin"), &binning[0]);
  }
  
  //===========================
  // SET MC COLORS
  //===========================
  
  p_hist->SetFillColor(kAzure+5);
  p_hist->SetFillStyle(1001);

  //===========================
  // Find Plot Maxima
  //===========================
  
  double ymax = 0;
  double ymin = 0;
  
  if (conf->get("ymax") != ""){
    ymax = stod(conf->get("ymax"));
  }
  else{
    ymax = 1.2*p_hist->GetMaximum();
  }

  if (conf->get("logy") == "true"){
    ymax *= 10;
    ymin = 0.1;
  }

  if (conf->get("ymin") != ""){
    ymin = stod(conf->get("ymin"));
  }

  
  cout<<"Proper plot maximum set to "<<ymax<<endl;
  
  TH2F* h_axes = new TH2F(Form("%s_axes",plot_name.Data()),plot_title,p_hist->GetNbinsX(),xmin,xmax,1000,ymin,ymax);
  
  if (conf->get("relabel_x_axis") == "true"){
    TString bin_label;
    for (int i = (int) xmin; i <= (int) xmax; i++)
    {
      bin_label=p_hist->GetXaxis()->GetBinLabel(p_hist->FindBin(i));
      h_axes->GetXaxis()->SetBinLabel(h_axes->FindBin(i), bin_label);
    }  
    //h_axes->GetXaxis()->LabelsOption("v");
    h_axes->GetXaxis()->SetLabelSize(.015);
  }
  //-----------------------
  // AXES FIX
  //-----------------------
  
  cout<<"Setting axis names"<<endl;
  h_axes->GetXaxis()->SetTitle(xlabel);
  h_axes->GetYaxis()->SetTitle(ylabel);
  
  TString stat_string_1;

  //===========================
  // Print Closure Stats
  //===========================
  if (conf->get("print_stats") == "true")
  {
    int low_val = stoi(conf->get("stats_low_val"));
    int high_val = stoi(conf->get("stats_high_val"));

    Double_t p_evts_gtr150_err;
    double p_evts_gtr150 = p_hist->IntegralAndError(p_hist->FindBin(low_val), p_hist->FindBin(high_val-.001), p_evts_gtr150_err);
    
    stat_string_1 = "Number of Events in "+hist_name+" from "+conf->get("stats_low_val")+" to "+conf->get("stats_high_val")+" : "+to_string(p_evts_gtr150)+" Error: "+to_string(p_evts_gtr150_err);


    cout<<f_primary->GetName()<<" STATS: "<<stat_string_1<<endl;
  }
  
  //----------------------
  // ADD OVERFLOW BIN
  //----------------------
  if (conf->get("overflow")=="true"){
    cout<<"Plot tagged for overflow bin, building..."<<endl;
    
    updateOverUnderflow(p_hist, xmax);
    /*double n_bins = p_hist->GetNbinsX();
    
    double overflow_primary = p_hist->GetBinContent(n_bins + 1);

    double max_primary = p_hist->Integral(p_hist->FindBin(xmax) - 1, n_bins);

    p_hist->SetBinContent(p_hist->FindBin(xmax) - 1, max_primary+overflow_primary);*/
  }      
  
  fullpad->SetLeftMargin(0.15);
  h_axes->GetYaxis()->SetTitleOffset(1.3);
  h_axes->GetYaxis()->SetTitleSize(0.05);
  h_axes->GetYaxis()->SetLabelSize(0.03);
  p_hist->SetLineWidth(3);
  
  cout<<"Drawing histograms"<<endl;
  TString plot_opts = (conf->get("plot_opts") != "") ? conf->get("plot_opts") : "HIST E";
  h_axes->Draw();
  p_hist->Draw(plot_opts+" SAME");
  
  fullpad->RedrawAxis();
  

  if (conf->get("luminosity_fb") != ""){
    fullpad->cd();
    drawCMSLatex(stod(conf->get("luminosity_fb")));
  }
  
  drawLatexFromTString(stat_string_1, .52,.5);

  cout<<"Saving..."<<endl;
  c->SaveAs(save_dir+plot_name+TString(".pdf"));
  c->SaveAs(save_dir+plot_name+TString(".png"));
  //c->SaveAs(save_dir+plot_name+TString(".root"));
  //c->SaveAs(save_dir+plot_name+TString(".C"));
  
  cout<<"Cleaning up plot variables"<<endl;
  delete p_hist;
  delete fullpad;
  delete c;

  f_primary->Close();
  delete f_primary;

  return errors;
}

TString drawCutDebug(ConfigParser *conf){
  /* This method expects conf to have a plot config loaded in already. */
  TString errors="";

  TString sample_loc;
  if (conf->get("file_path") != ""){
      sample_loc = TString(conf->get("file_path"));
  }
  else{
    sample_loc = TString(getDefaultHistDir(conf)+conf->get("sample")+".root");
  }

  TFile *f_primary =new TFile(sample_loc);

  cout << "Found files "<<endl;

  TString plot_name = conf->get("plot_name");
  TString plot_title = parseLatex(conf->get("title"));
  double xmax = stod(conf->get("xmax"));
  double xmin = stod(conf->get("xmin"));
  double bin_size = stod(conf->get("bin_size"));
  TString hist_name=conf->get("hist_name");
  TString xlabel=parseLatex(conf->get("xlabel"));
  TString ylabel=parseLatex(conf->get("ylabel"));
  TString save_dir=(conf->get("save_dir") != "") ? conf->get("save_dir") : getOutputDir(conf, "plot");


  cout << "Making Plots for: "<<plot_name<<endl;

  TH1D* p_hist = (TH1D*) ((TH1D*) f_primary->Get(hist_name))->Clone("phist_"+plot_name);
  cout<<hist_name<<" found in "<<f_primary->GetName()<<endl;


  cout << "Histograms pulled from files, adding draw options"<<endl;
  
  //============================================
  // Draw Data-MC Plots
  //============================================
  
  TCanvas * c = new TCanvas("c","",2000,2000);
  c->cd();
  gPad->SetRightMargin(0.05);
  gPad->Modified();
  gStyle->SetOptStat(kFALSE);
  TPad *fullpad = new TPad("fullpad", "fullpad", 0,0,1,1);
  
  fullpad->Draw();
  fullpad->cd();
    
  fullpad->SetRightMargin(0.05);
  if (conf->get("ExtraRightMargin") == "true")
  {
    fullpad->SetRightMargin(0.08);
  }
  fullpad->SetBottomMargin(0.3);
  
  fullpad->Draw();
  fullpad->cd();
  
  if (conf->get("logy") == "true")
  {
    cout<<"Plot tagged for log y-axis"<<endl;
    fullpad->SetLogy();
  }
  
  p_hist->Rebin(bin_size);
  
  //===========================
  // SET MC COLORS
  //===========================
  
  p_hist->SetFillColor(kAzure+5);
  p_hist->SetFillStyle(1001);

  //===========================
  // Find Plot Maxima
  //===========================
  
  double ymax = 0;

  ymax = 1.2*p_hist->GetMaximum();

  cout<<"Proper plot maximum set to "<<ymax<<endl;
  
  TH2F* h_axes = new TH2F(Form("%s_axes",plot_name.Data()),plot_title,p_hist->GetNbinsX(),xmin,xmax,1000,0.001,ymax);
  
  
  //-----------------------
  // AXES FIX
  //-----------------------
  
  cout<<"Setting axis names"<<endl;
  h_axes->GetXaxis()->SetTitle(xlabel);
  h_axes->GetYaxis()->SetTitle(ylabel);
  
  
  //----------------------
  // SET AXIS LABELS
  //----------------------
  ConfigParser label_conf(conf->get("labels_file"));
  label_conf.loadConfig("Error Labels");
  TString bin_label;
  for (int i = xmin; i<xmax; i++)
  {
    bin_label=parseLatex(label_conf[to_string(i)]);
    bin_label+=" ("+to_string((int) p_hist->GetBinContent(h_axes->FindBin(i)))+")";
    h_axes->GetXaxis()->SetBinLabel(h_axes->FindBin(i), bin_label);
  }  
  h_axes->GetXaxis()->LabelsOption("v");
  h_axes->GetXaxis()->SetLabelSize(.015);
  
  fullpad->SetLeftMargin(0.15);
  h_axes->GetYaxis()->SetTitleOffset(1.3);
  h_axes->GetYaxis()->SetTitleSize(0.05);
  h_axes->GetYaxis()->SetLabelSize(0.04);
  
  cout<<"Drawing histograms"<<endl;
  h_axes->Draw();
  p_hist->Draw("HIST SAME");
  
  fullpad->RedrawAxis();
  
  cout<<"Saving..."<<endl;
  c->SaveAs(save_dir+plot_name+TString(".pdf"));
  c->SaveAs(save_dir+plot_name+TString(".png"));
  //c->SaveAs(save_dir+plot_name+TString(".root"));
  //c->SaveAs(save_dir+plot_name+TString(".C"));
  
  cout<<"Cleaning up plot variables"<<endl;
  delete p_hist;
  delete fullpad;
  delete c;

  f_primary->Close();
  delete f_primary;

  return errors;
}

TString drawCutDebug(TString sample_name, TString sample_loc, TString save_dir){
  TString errors="";

  TFile *f_primary =new TFile(sample_loc);

  cout<<"Found files for Debug"<<endl;

  TString plot_name = TString("cuts_")+sample_name;
  TString plot_title = TString("Event Debug For ")+sample_name;
  TString hist_name="numEvents";
  
  cout << "Making Debug Plots for: "<<sample_name<<endl;

  TH1I* p_hist = (TH1I*) ((TH1I*) f_primary->Get(hist_name))->Clone("phist_"+plot_name);
  cout<<hist_name<<" found in "<<f_primary->GetName()<<endl;


  cout << "NumEvents pulled from files, adding draw options"<<endl;
  
  //============================================
  // Draw Data-MC Plots
  //============================================
  
  TCanvas * c = new TCanvas("c","",2000,2000);
  c->cd();
  gPad->SetRightMargin(0.05);
  gPad->Modified();
  gStyle->SetOptStat(kFALSE);
  TPad *fullpad = new TPad("fullpad", "fullpad", 0,0,1,1);
  
  fullpad->Draw();
  fullpad->cd();
    
  fullpad->SetRightMargin(0.05);
  fullpad->SetBottomMargin(0.3);
  
  fullpad->Draw();
  fullpad->cd();
  
  fullpad->SetLogy();
  //===========================
  // SET MC COLORS
  //===========================
  
  p_hist->SetFillColor(kAzure+5);
  p_hist->SetFillStyle(1001);

  //===========================
  // Find Plot Maxima
  //===========================
  
  double xmax = (double) p_hist->GetNbinsX();
  double xmin = 0;

  cout<<"Debug hist has "<<xmax<<" bins"<<endl;

  double ymax = 0;

  ymax = 1.2*p_hist->GetMaximum();

  cout<<"Proper plot maximum set to "<<ymax<<endl;
  
  TH2F* h_axes = new TH2F(Form("%s_axes",plot_name.Data()),plot_title,p_hist->GetNbinsX(),xmin,xmax,1000,0.001,ymax);
  
  
  //-----------------------
  // AXES FIX
  //-----------------------
  
  cout<<"Setting axis names"<<endl;
  h_axes->GetYaxis()->SetTitle("Count");
  
  
  //----------------------
  // SET AXIS LABELS
  //----------------------
  ConfigParser label_conf("configs/cutlabels.conf");
  label_conf.loadConfig("Error Labels");
  TString bin_label;
  for (int i = xmin; i<xmax; i++)
  {
    bin_label=parseLatex(label_conf[to_string(i)]);
    bin_label+=" ("+to_string((int) p_hist->GetBinContent(h_axes->FindBin(i)))+")";
    h_axes->GetXaxis()->SetBinLabel(h_axes->FindBin(i), bin_label);
  }  
  h_axes->GetXaxis()->LabelsOption("v");
  h_axes->GetXaxis()->SetLabelSize(.015);
  
  fullpad->SetLeftMargin(0.15);
  h_axes->GetYaxis()->SetTitleOffset(1.3);
  h_axes->GetYaxis()->SetTitleSize(0.05);
  h_axes->GetYaxis()->SetLabelSize(0.04);
  
  cout<<"Drawing histograms"<<endl;
  h_axes->Draw();
  p_hist->Draw("HIST SAME");
  
  fullpad->RedrawAxis();
  
  cout<<"Saving..."<<endl;
  c->SaveAs(save_dir+"Debug/"+plot_name+TString(".pdf"));
  c->SaveAs(save_dir+"Debug/"+plot_name+TString(".png"));
  //c->SaveAs(save_dir+plot_name+TString(".root"));
  //c->SaveAs(save_dir+plot_name+TString(".C"));
  
  cout<<"Cleaning up plot variables"<<endl;
  delete p_hist;
  delete fullpad;
  delete c;

  f_primary->Close();
  delete f_primary;

  return errors;
}

TString drawWeightDebug(TString sample_name, TString sample_loc, TString save_dir, TString hist_name){
  TString errors="";

  TFile *f_primary =new TFile(sample_loc);

  cout<<"Found files for Debug"<<endl;

  TString plot_name = TString(hist_name+"_"+sample_name);
  TString plot_title = TString("Event Debug For ")+sample_name;
  
  cout << "Making Debug Plots for: "<<sample_name<<endl;

  TH1D* p_hist = (TH1D*) ((TH1D*) f_primary->Get(hist_name))->Clone("phist_"+plot_name);
  cout<<hist_name<<" found in "<<f_primary->GetName()<<endl;


  cout << "Weight Log pulled from files, adding draw options"<<endl;
  
  //============================================
  // Draw Data-MC Plots
  //============================================
  
  TCanvas * c = new TCanvas("c","",2000,2000);
  c->cd();
  gPad->SetRightMargin(0.05);
  gPad->Modified();
  gStyle->SetOptStat(kFALSE);
  TPad *fullpad = new TPad("fullpad", "fullpad", 0,0,1,1);
  
  fullpad->Draw();
  fullpad->cd();
    
  fullpad->SetRightMargin(0.05);
  fullpad->SetBottomMargin(0.3);
  fullpad->SetLeftMargin(0.15);

  
  fullpad->Draw();
  fullpad->cd();
  
  fullpad->SetLogy();
  
  //===========================
  // Find Plot Maxima
  //===========================
  
  double xmax = (double) p_hist->GetNbinsX();
  double xmin = 0;

  updateOverUnderflow(p_hist, p_hist->GetBinLowEdge(xmax));

  cout<<"Debug hist has "<<xmax<<" bins"<<endl;

  double ymax = 0;

  ymax = 1.2*p_hist->GetMaximum();

  cout<<"Proper plot maximum set to "<<ymax<<endl;
  
  TH2F* h_axes = new TH2F(Form("%s_axes",plot_name.Data()),plot_title, p_hist->GetNbinsX(), p_hist->GetBinLowEdge(1), p_hist->GetBinLowEdge(xmax+1),1000,0.5,ymax);
  
  
  //-----------------------
  // AXES FIX
  //-----------------------
  
  cout<<"Setting axis names"<<endl;
  h_axes->GetXaxis()->SetTitle("Event weight");
  h_axes->GetYaxis()->SetTitle("Count");
  
  
  //----------------------
  // SET AXIS LABELS
  //----------------------
  cout<<"Setting axis options"<<endl;
  h_axes->GetYaxis()->SetTitleOffset(1.3);
  h_axes->GetYaxis()->SetTitleSize(0.05);

  //Ensure Flat Binning
  TH1D* flat_hist = new TH1D("flat_hist", p_hist->GetTitle(), p_hist->GetNbinsX(), p_hist->GetBinLowEdge(1), p_hist->GetBinLowEdge(xmax+1));
  for(int i = 0; i<=xmax+1; i++) //0 is underflow, NbinsX is last non-overflow bin.
  {
    flat_hist->SetBinContent(i, p_hist->GetBinContent(i));
  }  

  cout<<"Setting axis labels"<<endl;
  TString bin_label;
  for (int i = 1; i<=xmax; i++)
  {
    if (plot_name.Contains("flat")){
      bin_label=to_string((double) p_hist->GetBinLowEdge(i));
    }
    else{
      bin_label=to_string(pow(10,(double) p_hist->GetBinLowEdge(i))); 
    }
    bin_label+=" ("+to_string((int) p_hist->GetBinContent(i))+")";
    h_axes->GetXaxis()->SetBinLabel(i, bin_label);
  } 
  
  h_axes->GetXaxis()->LabelsOption("v");
  h_axes->GetXaxis()->SetNdivisions(xmax+2);
  h_axes->GetXaxis()->SetLabelSize(.015);
  h_axes->GetXaxis()->SetTitleOffset(2);
  

  //===========================
  // SET PLOT COLORS
  //===========================
  
  flat_hist->SetFillColor(kAzure+5);
  flat_hist->SetFillStyle(1001);

  fullpad->cd();
  cout<<"Drawing histogram"<<endl;
  h_axes->Draw();
  flat_hist->Draw("HIST SAME");
  
  fullpad->RedrawAxis();
  
  cout<<"Saving..."<<endl;
  c->SaveAs(save_dir+"Debug/"+plot_name+TString(".pdf"));
  c->SaveAs(save_dir+"Debug/"+plot_name+TString(".png"));
  //c->SaveAs(save_dir+"Debug/"+plot_name+TString(".root"));
  //c->SaveAs(save_dir+plot_name+TString(".C"));

  TFile *f = new TFile(save_dir+"Debug/"+plot_name+TString(".root"), "recreate");
  f->cd();
  flat_hist->Write();
  h_axes->Write();
  f->Close();
  
  cout<<"Cleaning up plot variables"<<endl;
  delete p_hist;
  delete flat_hist;
  delete fullpad;
  delete c;

  f_primary->Close();
  delete f_primary;

  return errors;
}

TString drawDebugPlots(ConfigParser *conf){
  
  TString save_dir=(conf->get("save_dir") != "") ? conf->get("save_dir") : getOutputDir(conf, "plot");
  TString sample_name, sample_loc;
  TString default_hist_dir = getDefaultHistDir(conf);

  if (conf->get("PLOT_TYPE") == "ratio" || conf->get("PLOT_TYPE") == "stack"){
    int i = 0;
    //Add files from which to obtain histos
    while (conf->get("file_"+to_string(i)+"_path") != "" || conf->get("sample_"+to_string(i)) != ""){
      if (conf->get("file_"+to_string(i)+"_path") != ""){
        sample_loc = TString(conf->get("file_"+to_string(i)+"_path"));
        sample_name=sample_loc(sample_loc.Last('/')+1, sample_loc.Index(".root")-sample_loc.Last('/')-1); //gets the name of the .root file
      }
      else{
        sample_loc = TString(default_hist_dir+conf->get("sample_"+to_string(i))+".root");
        sample_name = TString(conf->get("sample_"+to_string(i)));
      }

      drawCutDebug(sample_name, sample_loc, save_dir);
      drawWeightDebug(sample_name, sample_loc, save_dir, "weight_log");
      drawWeightDebug(sample_name, sample_loc, save_dir, "weight_log_flat");
      i++;
    }
    if (i<1){
      cout<<"NO DEBUG PLOTS MADE!!! Did you start counting at 1? Not how to do it."<<endl;
      return TString("ERROR: Could not build debug plots, no sample_0 or file_0_path specified.\n");
    }
  }
  
  else if (conf->get("PLOT_TYPE") == "single" || conf->get("PLOT_TYPE") == "single2D"){
    if (conf->get("file_path") != ""){
      sample_loc = TString(conf->get("file_path"));
      sample_name=sample_loc(sample_loc.Last('/')+1, sample_loc.Index(".root")-sample_loc.Last('/')-1); //gets the name of the .root file
    }
    else{
      sample_loc=TString(default_hist_dir+conf->get("sample")+".root");
      sample_name=TString(conf->get("sample"));
    }
    drawCutDebug(sample_name, sample_loc, save_dir);
    drawWeightDebug(sample_name, sample_loc, save_dir, "weight_log");
    drawWeightDebug(sample_name, sample_loc, save_dir, "weight_log_flat");
  }
  else{
    return TString("ERROR: Could not build debug plots, unknown Plot Type: "+conf->get("PLOT_TYPE")+"\n");
  }

  return TString("Debug Plots Made.\n");
}

TString drawSingleTH2(ConfigParser *conf){
  TString errors="";
  TString plot_name = conf->get("plot_name");
  TString plot_title = parseLatex(conf->get("title"));
  
  double xmax = (conf->get("xmax") != "") ? stod(conf->get("xmax")) : 500;
  double xmin = (conf->get("xmin") != "") ? stod(conf->get("xmin")) : 500;
  double ymax = (conf->get("ymax") != "") ? stod(conf->get("ymax")) : 500;
  double ymin = (conf->get("ymin") != "") ? stod(conf->get("ymin")) : 500;

  double bin_size_x = (conf->get("bin_size_x") != "") ? stod(conf->get("bin_size_x")) : 1;
  double bin_size_y = (conf->get("bin_size_y") != "") ? stod(conf->get("bin_size_y")) : 1;
  
  TString hist_name=conf->get("hist_name");
  TString xlabel=parseLatex(conf->get("xlabel"));
  TString ylabel=parseLatex(conf->get("ylabel"));
  TString save_dir=(conf->get("save_dir") != "") ? conf->get("save_dir") : getOutputDir(conf, "plot");

  cout<<"Options set"<<endl;

  TString sample_loc;
  if (conf->get("file_path") != ""){
      sample_loc = TString(conf->get("file_path"));
  }
  else{
    sample_loc = TString(getDefaultHistDir(conf)+conf->get("sample")+".root");
  }

  TFile *f_primary =new TFile(sample_loc);

  cout << "Found files "<<endl;

  TH2D *h = (TH2D*) ((TH2D*) f_primary->Get(hist_name))->Clone("hist_"+plot_name);

  cout<<"Found histogram "<<hist_name<<" with "<<h->GetEntries()<<" entries"<<endl;

  TCanvas * c = new TCanvas("c","",2000,2000);
  c->cd();
  gPad->Modified();
  gStyle->SetOptStat(kFALSE);
  //gStyle->SetTitleW(0.6);
  //gStyle->SetTitleH(0.06);
  //gStyle->SetTitleFont(12);
  TPad *fullpad = new TPad("fullpad", "fullpad", 0,0,1,1);
  
  fullpad->Draw();
  fullpad->cd();
  
  TH2F* h_axes = new TH2F(Form("%s_axes",plot_name.Data()),plot_title,1,xmin,xmax,1,ymin,ymax);
  h_axes->GetXaxis()->SetTitle(xlabel);
  h_axes->GetYaxis()->SetTitle(ylabel);
  h_axes->GetYaxis()->SetTitleOffset(1.3);
  h_axes->GetYaxis()->SetTitleSize(0.05);
  h_axes->GetYaxis()->SetLabelSize(0.04);
  h_axes->GetXaxis()->SetTitleSize(0.05);
  h_axes->GetXaxis()->SetLabelSize(0.04);
  //cout<<__LINE__<<endl;
  cout<<"Drawing axes"<<endl;
  h_axes->Draw();

  h->Rebin2D(bin_size_x, bin_size_y);

  h->GetXaxis()->SetRangeUser(xmin, xmax);
  h->GetYaxis()->SetRangeUser(ymin, ymax);

  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  if (ymax >= 1000){
    fullpad->SetLeftMargin(0.15);
    h->GetYaxis()->SetTitleOffset(1.8);
  }
  else{
    fullpad->SetLeftMargin(0.1);
    h->GetYaxis()->SetTitleOffset(1.3); 
  }
  if (h->GetMaximum() > 100){
    fullpad->SetRightMargin(0.15);
  }
  else{
    fullpad->SetRightMargin(0.1);
  }

  h->SetTitle(plot_title);
  gStyle->SetTitleW(.9); //title width 
  //gStyle->SetTitleSize(1.5); //title size 
  gStyle->SetTitleH(.1); //title height 
  //h->SetTitleSize(2);

  h->Draw("colz same");

  cout<<"Saving..."<<endl;
  c->SaveAs(save_dir+plot_name+TString(".pdf"));
  c->SaveAs(save_dir+plot_name+TString(".png"));
  //c->SaveAs(save_dir+plot_name+TString(".root"));
  //c->SaveAs(save_dir+plot_name+TString(".C"));
  
  cout<<"Cleaning up plot variables"<<endl;
  delete h;
  delete fullpad;
  delete c;

  f_primary->Close();
  delete f_primary;

  return errors;
}

void drawPlots(TString config_file, bool draw_debugs = true){
  TString errors="";

  ConfigParser *configs=new ConfigParser(config_file.Data());

  setTDRStyle();

  s_colors[TString("WZ")]         =  rc_orange.GetNumber();
  s_colors[TString("WW")]         =  rc_blue_violet.GetNumber();
  s_colors[TString("WJets")]      =  rc_dark_yellow.GetNumber();
  s_colors[TString("ZZ")]         =  rc_green.GetNumber();
  s_colors[TString("TTBar2l")]    =  rc_violet.GetNumber();
  s_colors[TString("TTBar1l")]    =  rc_yellow.GetNumber();
  s_colors[TString("DY")]         =  rc_light_green.GetNumber();
  s_colors[TString("TTV")]        =  rc_brown.GetNumber();
  s_colors[TString("SingleTop")]  =  rc_dark_violet.GetNumber();
  s_colors[TString("Other")]      =  rc_pink.GetNumber();
  s_colors[TString("ZH")]         =  rc_pink.GetNumber();
  s_colors[TString("Wh")]         =  46;
  s_colors[TString("VVV")]        =  kMagenta;
  s_colors[TString("WWW")]        =  kRed;
  s_colors[TString("QCD")]        =  11;
  
  TGaxis::SetExponentOffset(-0.07, 0, "y"); // X and Y offset for Y axis
  TGaxis::SetExponentOffset(-.8, -0.07, "x"); // X and Y offset for X axis

  while(configs->loadNextConfig()) {
    if (configs->get("PLOT_TYPE") == "ratio"){
      errors+=drawArbitraryNumberWithResidual(configs);
    }
    else if (configs->get("PLOT_TYPE") == "single"){
      errors+=drawSingleTH1(configs);
    }
    else if (configs->get("PLOT_TYPE") == "stack"){
      errors+=drawArbitraryNumber(configs);
    }
    else if (configs->get("PLOT_TYPE") == "debug"){
      errors+=drawCutDebug(configs);
    }
    else if (configs->get("PLOT_TYPE") == "single2D"){
      errors+=drawSingleTH2(configs);
    }
  }
  
  if (draw_debugs){
    errors+=drawDebugPlots(configs);
  }

  cout<<errors<<endl;
  return;
}

