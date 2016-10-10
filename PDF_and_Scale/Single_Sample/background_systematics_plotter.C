#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
using namespace std;

void plot(string filename, Float_t nevents, Float_t xsec)
{
  // Change these values as needed
  Float_t int_lumi = 12900.0;
  TString xaxis_title = TString("Photon #it{p}_{T} [GeV]");
  TString yaxis_title = TString("Events / bin");
  const int nBins = 5;
  const int nPDFreplicas = 100;

  TString histname_facUp = TString("Photon_Et_range_102");
  TString histname_facDown = TString("Photon_Et_range_103");
  TString histname_renUp = TString("Photon_Et_range_104");
  TString histname_renDown = TString("Photon_Et_range_107");

  std::vector<TString> histnames_pdf;
  histnames_pdf.clear();
  for(int i = 0; i < nPDFreplicas+1; i++)
  {
    char histindex_char[100];
    sprintf(histindex_char, "_%d", i);
    std::string histindex(histindex_char);
    TString histname = TString("Photon_Et_range"+histindex);
    histnames_pdf.push_back(histname);
  }

  TFile* inputfile = new TFile(TString(filename+".root"));

  // Scale variations
  
  TH1F* histo_facUp = (TH1F*)((TH1F*)inputfile->Get(histname_facUp))->Clone(TString("h_"+filename+"_facUp"));
  histo_facUp->SetStats(0);
  histo_facUp->Scale(xsec*int_lumi/nevents);
  histo_facUp->SetTitle("");
  histo_facUp->GetXaxis()->SetTitle(xaxis_title);
  histo_facUp->GetYaxis()->SetTitle(yaxis_title);
  TH1F* histo_facDown = (TH1F*)((TH1F*)inputfile->Get(histname_facDown))->Clone(TString("h_"+filename+"_facDown"));
  histo_facDown->SetStats(0);
  histo_facDown->Scale(xsec*int_lumi/nevents);
  histo_facDown->SetTitle("");
  histo_facDown->GetXaxis()->SetTitle(xaxis_title);
  histo_facDown->GetYaxis()->SetTitle(yaxis_title);
  
  TH1F* histo_renUp = (TH1F*)((TH1F*)inputfile->Get(histname_renUp))->Clone(TString("h_"+filename+"_renUp"));
  histo_renUp->SetStats(0);
  histo_renUp->Scale(xsec*int_lumi/nevents);
  histo_renUp->SetTitle("");
  histo_renUp->GetXaxis()->SetTitle(xaxis_title);
  histo_renUp->GetYaxis()->SetTitle(yaxis_title);
  TH1F* histo_renDown = (TH1F*)((TH1F*)inputfile->Get(histname_renDown))->Clone(TString("h_"+filename+"_renDown"));
  histo_renDown->SetStats(0);
  histo_renDown->Scale(xsec*int_lumi/nevents);
  histo_renDown->SetTitle("");
  histo_renDown->GetXaxis()->SetTitle(xaxis_title);
  histo_renDown->GetYaxis()->SetTitle(yaxis_title);

  // PDF variations

  std::vector<Float_t> pdf_sum_bin1;
  pdf_sum_bin1.clear();
  std::vector<Float_t> pdf_sum_bin2;
  pdf_sum_bin2.clear();
  std::vector<Float_t> pdf_sum_bin3;
  pdf_sum_bin3.clear();
  std::vector<Float_t> pdf_sum_bin4;
  pdf_sum_bin4.clear();
  std::vector<Float_t> pdf_sum_bin5;
  pdf_sum_bin5.clear();
  
  for(int i = 0; i <= nPDFreplicas; i++)
  {
    TH1F* hist_pdf = (TH1F*)inputfile->Get(histnames_pdf[i]);
    pdf_sum_bin1.push_back(hist_pdf->GetBinContent(1));
    pdf_sum_bin2.push_back(hist_pdf->GetBinContent(2));
    pdf_sum_bin3.push_back(hist_pdf->GetBinContent(3));
    pdf_sum_bin4.push_back(hist_pdf->GetBinContent(4));
    pdf_sum_bin5.push_back(hist_pdf->GetBinContent(5));
  }
  
  double sum_of_sums_bin1 = 0.0;
  double sum_of_sums_bin2 = 0.0;
  double sum_of_sums_bin3 = 0.0;
  double sum_of_sums_bin4 = 0.0;
  double sum_of_sums_bin5 = 0.0;
  for(int i = 1; i <= nPDFreplicas; i++)
  {
    sum_of_sums_bin1 += pdf_sum_bin1[i];
    sum_of_sums_bin2 += pdf_sum_bin2[i];
    sum_of_sums_bin3 += pdf_sum_bin3[i];
    sum_of_sums_bin4 += pdf_sum_bin4[i];
    sum_of_sums_bin5 += pdf_sum_bin5[i];
  }
  double mean_Npassing_bin1 = sum_of_sums_bin1/nPDFreplicas;
  double mean_Npassing_bin2 = sum_of_sums_bin2/nPDFreplicas;
  double mean_Npassing_bin3 = sum_of_sums_bin3/nPDFreplicas;
  double mean_Npassing_bin4 = sum_of_sums_bin4/nPDFreplicas;
  double mean_Npassing_bin5 = sum_of_sums_bin5/nPDFreplicas;
  
  double sum_of_squared_residuals_bin1 = 0.0;
  double sum_of_squared_residuals_bin2 = 0.0;
  double sum_of_squared_residuals_bin3 = 0.0;
  double sum_of_squared_residuals_bin4 = 0.0;
  double sum_of_squared_residuals_bin5 = 0.0;
  for(int i = 1; i <= nPDFreplicas; i++)
  {
    sum_of_squared_residuals_bin1 += pow((pdf_sum_bin1[i] - mean_Npassing_bin1),2.0);
    sum_of_squared_residuals_bin2 += pow((pdf_sum_bin2[i] - mean_Npassing_bin2),2.0);
    sum_of_squared_residuals_bin3 += pow((pdf_sum_bin3[i] - mean_Npassing_bin3),2.0);
    sum_of_squared_residuals_bin4 += pow((pdf_sum_bin4[i] - mean_Npassing_bin4),2.0);
    sum_of_squared_residuals_bin5 += pow((pdf_sum_bin5[i] - mean_Npassing_bin5),2.0);
  }
  Float_t rms_error_bin1 = sqrt(sum_of_squared_residuals_bin1/(nPDFreplicas-1)); // RMS error in event weight sum for this bin
  Float_t rms_error_bin2 = sqrt(sum_of_squared_residuals_bin2/(nPDFreplicas-1));
  Float_t rms_error_bin3 = sqrt(sum_of_squared_residuals_bin3/(nPDFreplicas-1));
  Float_t rms_error_bin4 = sqrt(sum_of_squared_residuals_bin4/(nPDFreplicas-1));
  Float_t rms_error_bin5 = sqrt(sum_of_squared_residuals_bin5/(nPDFreplicas-1));
  
  TH1F* histo_pdfUp = (TH1F*)((TH1F*)inputfile->Get(histnames_pdf[0]))->Clone(TString("h_"+filename+"_pdfUp"));
  histo_pdfUp->SetBinContent(1,histo_pdfUp->GetBinContent(1) + rms_error_bin1);
  histo_pdfUp->SetBinContent(2,histo_pdfUp->GetBinContent(2) + rms_error_bin2);
  histo_pdfUp->SetBinContent(3,histo_pdfUp->GetBinContent(3) + rms_error_bin3);
  histo_pdfUp->SetBinContent(4,histo_pdfUp->GetBinContent(4) + rms_error_bin4);
  histo_pdfUp->SetBinContent(5,histo_pdfUp->GetBinContent(5) + rms_error_bin5);
  histo_pdfUp->SetStats(0);
  histo_pdfUp->Scale(xsec*int_lumi/nevents);
  histo_pdfUp->SetTitle("");
  histo_pdfUp->GetXaxis()->SetTitle(xaxis_title);
  histo_pdfUp->GetYaxis()->SetTitle(yaxis_title);
  TH1F* histo_pdfDown = (TH1F*)((TH1F*)inputfile->Get(histnames_pdf[0]))->Clone(TString("h_"+filename+"_pdfDown"));
  histo_pdfDown->SetBinContent(1,TMath::Max(histo_pdfDown->GetBinContent(1) - rms_error_bin1,0.0));
  histo_pdfDown->SetBinContent(2,TMath::Max(histo_pdfDown->GetBinContent(2) - rms_error_bin2,0.0));
  histo_pdfDown->SetBinContent(3,TMath::Max(histo_pdfDown->GetBinContent(3) - rms_error_bin3,0.0));
  histo_pdfDown->SetBinContent(4,TMath::Max(histo_pdfDown->GetBinContent(4) - rms_error_bin4,0.0));
  histo_pdfDown->SetBinContent(5,TMath::Max(histo_pdfDown->GetBinContent(5) - rms_error_bin5,0.0));
  histo_pdfDown->SetStats(0);
  histo_pdfDown->Scale(xsec*int_lumi/nevents);
  histo_pdfDown->SetTitle("");
  histo_pdfDown->GetXaxis()->SetTitle(xaxis_title);
  histo_pdfDown->GetYaxis()->SetTitle(yaxis_title);
  
  // Set zero bins to a small value to avoid fitting issues
  // Different values (0.97e-6, 1.03e-6) are meant to preserve the fact that pdfUp and pdfDown refer
  // to up and down shifts in the cross section itself, and that fac/ren up & down refer to shifts in the
  // QCD scales (resulting in oppositely-directed shifts in the cross section).
  for(int i = 1; i <= nBins; i++)
  {
    if(histo_facUp->GetBinContent(i) == 0.0)
      histo_facUp->SetBinContent(i,0.97e-6);
    if(histo_facDown->GetBinContent(i) == 0.0)
      histo_facDown->SetBinContent(i,1.03e-6);
    if(histo_renUp->GetBinContent(i) == 0.0)
      histo_renUp->SetBinContent(i,0.97e-6);
    if(histo_renDown->GetBinContent(i) == 0.0)
      histo_renDown->SetBinContent(i,1.03e-6);
    if(histo_pdfUp->GetBinContent(i) == 0.0)
      histo_pdfUp->SetBinContent(i,1.03e-6);
    if(histo_pdfDown->GetBinContent(i) == 0.0)
      histo_pdfDown->SetBinContent(i,0.97e-6);
  }
  
  TFile *outputFile = new TFile(TString("histos_"+filename+".root"),"RECREATE");
  outputFile->cd();
  histo_facUp->Write();
  histo_facDown->Write();
  histo_renUp->Write();
  histo_renDown->Write();
  histo_pdfUp->Write();
  histo_pdfDown->Write();

  inputfile->Close();
  outputFile->Close();
}

void background_systematics_plotter()
{
  std::vector<string> filenames;
  filenames.clear();
  std::vector<Float_t> nevents;
  nevents.clear();
  std::vector<Float_t> xsecs;
  xsecs.clear();

  // Make plots from as many files as desired
  filenames.push_back("ZnnG_pdfscale_ZLLGJets");
  nevents.push_back(440320.0);
  xsecs.push_back(0.143);

  for(int i = 0; i < filenames.size(); i++)
  {
    plot(filenames[i],nevents[i],xsecs[i]);
  }
}