//Helper functions for Config parsing.
#ifndef INCLUDED_HIST_TOOLS
#define INCLUDED_HIST_TOOLS


void updateOverUnderflow( TH1D * &hist, double xmax, double xmin = -100000 ){
  /* updates bins at the edges of xmax (xmin) with everything above (below) including over(under)flow */
  int overflowbin = hist->FindBin(xmax-0.0001);
  for( int bini = overflowbin; bini <= hist->GetNbinsX(); bini++ ){
    hist->SetBinContent( overflowbin, hist->GetBinContent( overflowbin ) + hist->GetBinContent( bini + 1 ) ); 
    hist->SetBinError  ( overflowbin, sqrt( pow(hist->GetBinError  ( overflowbin ), 2 ) + pow( hist->GetBinError( bini + 1 ), 2 ) ) );  
    hist->SetBinContent( bini + 1, 0 );
    hist->SetBinError  ( bini + 1, 0 );
  }

  if (xmin > -100000){
    int underflowbin = hist->FindBin(xmin+0.0001);
    for(int bini = overflowbin; bini > 0; bini-- ){
      hist->SetBinContent( underflowbin, hist->GetBinContent( underflowbin ) + hist->GetBinContent( bini - 1 ) ); 
      hist->SetBinError  ( underflowbin, sqrt( pow(hist->GetBinError  ( underflowbin ), 2 ) + pow( hist->GetBinError( bini - 1 ), 2 ) ) );  
      hist->SetBinContent( bini - 1, 0 );
      hist->SetBinError  ( bini - 1, 0 );
    }
  }
}

void blindAfter(TH1D * &hist, double xmax){
  /* Sets all the bins starting from xmax to 0, count and error */
  int max_bin = hist->FindBin(xmax);
  for( int bini = max_bin; bini <= hist->GetNbinsX()+1; bini++ ){
    hist->SetBinContent( bini, 0 ); 
    hist->SetBinError  ( bini, 0 );  
  }
}

void zeroNegatives(TH1D* &hist){
  /* Zeros all the bins with negative weights */
  for( int i = 0; i <= hist->GetNbinsX()+1; i++ ){
    if ( ((double) hist->GetBinContent(i)) < 0 ){
      hist->SetBinContent( i, 0 );
    }
  }
}

#endif