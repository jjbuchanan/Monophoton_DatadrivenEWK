#include <fstream>
#include <vector>
#include <iomanip>
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

void plot(TString histname, Double_t yaxis_compress, Double_t leg_xoffset, Double_t leg_yoffset, TString xaxis_title, TString plotname)
{
  const int nBins = 5;
  
  TCanvas *c = new TCanvas("c", "canvas",700,640);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  c->SetLeftMargin(0.15);
  c->cd();

  TFile* f_ZnnG_ZNuNuGJets = new TFile("ZnnG_ZNuNuGJets.root");
  TH1F* histo_ZnnG_ZNuNuGJets = (TH1F*)f_ZnnG_ZNuNuGJets->Get(histname);
  histo_ZnnG_ZNuNuGJets->SetStats(0);
  histo_ZnnG_ZNuNuGJets->Scale(1.0/375920.0);
  histo_ZnnG_ZNuNuGJets->Scale(0.1903);
  histo_ZnnG_ZNuNuGJets->SetTitle("");
  histo_ZnnG_ZNuNuGJets->GetXaxis()->SetTitle("");
  histo_ZnnG_ZNuNuGJets->GetXaxis()->SetTickLength(0);
  histo_ZnnG_ZNuNuGJets->GetXaxis()->SetLabelOffset(999);
  histo_ZnnG_ZNuNuGJets->GetYaxis()->SetTitle("");
  histo_ZnnG_ZNuNuGJets->GetYaxis()->SetTickLength(0);
  histo_ZnnG_ZNuNuGJets->GetYaxis()->SetLabelOffset(999);
  
  TFile* f_ZnnG_WGJets = new TFile("ZnnG_WGJets.root");
  TH1F* histo_ZnnG_WGJets = (TH1F*)f_ZnnG_WGJets->Get(histname);
  histo_ZnnG_WGJets->SetStats(0);
  histo_ZnnG_WGJets->Scale(1.0/483496.0);
  histo_ZnnG_WGJets->Scale(0.6637);
  histo_ZnnG_WGJets->SetTitle("");
  histo_ZnnG_WGJets->GetXaxis()->SetTitle("");
  histo_ZnnG_WGJets->GetXaxis()->SetTickLength(0);
  histo_ZnnG_WGJets->GetXaxis()->SetLabelOffset(999);
  histo_ZnnG_WGJets->GetYaxis()->SetTitle("");
  histo_ZnnG_WGJets->GetYaxis()->SetTickLength(0);
  histo_ZnnG_WGJets->GetYaxis()->SetLabelOffset(999);
  
  TGraphAsymmErrors *transfer_factor = new TGraphAsymmErrors();
  for(int i = 1; i <= nBins; i++)
  {
    Double_t x_value = histo_ZnnG_ZNuNuGJets->GetBinCenter(i);
    Double_t x_halfwidth = (histo_ZnnG_ZNuNuGJets->GetBinWidth(i))/2.0;
    Double_t y_value = 0.0;
    Double_t y_error_up = 0.0;
    Double_t y_error_down = 0.0;
    if(histo_ZnnG_ZNuNuGJets->GetBinContent(i) > 0.0){
      Double_t N_bin_ZnnGinZnnG = histo_ZnnG_ZNuNuGJets->GetBinContent(i);
      Double_t N_bin_WGinZnnG = histo_ZnnG_WGJets->GetBinContent(i);
      y_value = N_bin_WGinZnnG / N_bin_ZnnGinZnnG;
      y_error_up = y_value*sqrt((N_bin_WGinZnnG*((0.6637/483496.0)-N_bin_WGinZnnG/483496.0))/(N_bin_WGinZnnG*N_bin_WGinZnnG) + (N_bin_ZnnGinZnnG*((0.1903/375920.0)-N_bin_ZnnGinZnnG/375920.0))/(N_bin_ZnnGinZnnG*N_bin_ZnnGinZnnG));
      y_error_down = y_error_up; //Symmetric, for now
    }
    transfer_factor->SetPoint(i-1,x_value,y_value);
    transfer_factor->SetPointError(i-1,x_halfwidth,x_halfwidth,y_error_down,y_error_up);
  }
  transfer_factor->SetTitle("");
  transfer_factor->GetXaxis()->SetTitle("");
  transfer_factor->GetXaxis()->SetTickLength(0);
  transfer_factor->GetXaxis()->SetLabelOffset(999);
  transfer_factor->GetYaxis()->SetTitle("");
  transfer_factor->GetYaxis()->SetTickLength(0);
  transfer_factor->GetYaxis()->SetLabelOffset(999);
  transfer_factor->SetMarkerStyle(21);
  transfer_factor->SetFillColor(kGreen);
  transfer_factor->GetYaxis()->SetRangeUser(0.0,0.4);
  transfer_factor->Draw("AP5");
  TGraphAsymmErrors* transfer_factor_clone = (TGraphAsymmErrors*)transfer_factor->Clone();
  for(int i = 1; i <= nBins; i++)
  {
    Double_t x_halfwidth = (histo_ZnnG_WGJets->GetBinWidth(i))/2.0;
    Double_t y_error_up = 0.0;
    Double_t y_error_down = 0.0;
    transfer_factor_clone->SetPointError(i-1,x_halfwidth,x_halfwidth,y_error_down,y_error_up);
  }
  transfer_factor_clone->Draw("P SAME");

//   //Central location of leg defined to be location of leg in phoPt plot
//   TLegend* leg = new TLegend(0.5+leg_xoffset,0.58075+leg_yoffset,0.885387+leg_xoffset,0.862969+leg_yoffset,"");
//   leg->AddEntry(histo_data,"Data");
//   leg->AddEntry(histo_WG,"W#gamma#rightarrowl#nu#gamma","F");
// //  leg->AddEntry(histo_diphoton,"#gamma#gamma","F");
//   leg->AddEntry(histo_TTG,"tt#gamma","F");
//   leg->AddEntry(histo_GJets_40to100,"#gamma+jet","F");
//   leg->AddEntry(histo_ZllG,"Z(ll)#gamma","F");
//   leg->AddEntry(histo_elefake,"e#rightarrow#gamma MisID","F");
//   leg->AddEntry(histo_WWG,"jet#rightarrow#gamma MisID,WW#gamma, t#gamma","F");
//   leg->SetFillColor(kWhite);
//   leg->SetFillStyle(0);
//   leg->SetTextSize(0.035);
//   leg->Draw();

  TLatex *texS = new TLatex(0.54023,0.907173,"#sqrt{s} = 13 TeV");
  texS->SetNDC();
  texS->SetTextFont(42);
  texS->SetTextSize(0.045);
  texS->Draw();
  TLatex *texS1 = new TLatex(0.16092,0.907173,"#bf{CMS} #it{Preliminary}");
  texS1->SetNDC();
  texS1->SetTextFont(42);
  texS1->SetTextSize(0.045);
  texS1->Draw();

  c->Update();

  double xmin = c->GetUxmin();
  double ymin = c->GetUymin();
  double xmax = c->GetUxmax();
  double ymax = c->GetUymax();

  TGaxis *xaxis = new TGaxis(xmin,ymin,xmax,ymin,xmin,xmax,510);
  xaxis->SetTitle(xaxis_title);
  xaxis->SetLabelFont(42);
  xaxis->SetLabelSize(0.040);
  xaxis->SetTitleFont(42);
  xaxis->SetTitleSize(0.045);
  xaxis->Draw("SAME");

  TGaxis *xaxis_top = new TGaxis(xmin,ymax,xmax,ymax,xmin,xmax,510,"-");
  xaxis_top->SetTitle("");
  xaxis_top->SetLabelOffset(999);
  xaxis_top->Draw("SAME");

  TGaxis *yaxis = new TGaxis(xmin,ymin,xmin,ymax,ymin,ymax,510);
  yaxis->SetTitle("f");
  yaxis->SetLabelFont(42);
  yaxis->SetLabelSize(0.040);
  yaxis->SetTitleFont(42);
  yaxis->SetTitleSize(0.045);
  yaxis->SetTitleOffset(1.3);
  yaxis->Draw("SAME");

  TGaxis *yaxis_right = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+");
  yaxis_right->SetTitle("");
  yaxis_right->SetLabelOffset(999);
  yaxis_right->Draw("SAME");

  c->SaveAs(TString("transfer_znngZNuNuGJets_to_znngWGJets_"+plotname+".png"));
  c->SaveAs(TString("transfer_znngZNuNuGJets_to_znngWGJets_"+plotname+".pdf"));
  delete(c);
}

void znng_transfer_factor_plotter()
{
  std::vector<TString> histnames;
  histnames.clear();
  std::vector<Double_t> yaxis_compressions;
  yaxis_compressions.clear();
  std::vector<Double_t> leg_xoffsets;
  leg_xoffsets.clear();
  std::vector<Double_t> leg_yoffsets;
  leg_yoffsets.clear();
  std::vector<TString> xaxis_titles;
  xaxis_titles.clear();
  std::vector<TString> plotnames;
  plotnames.clear();

//  histnames.push_back(TString("_4"));
//  yaxis_compressions.push_back(1.);
//  leg_xoffsets.push_back(0.);
//  leg_yoffsets.push_back(0.);
//  xaxis_titles.push_back(TString(""));
//  plotnames.push_back(TString(""));

  histnames.push_back(TString("Photon_Et_range_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("phoPt"));

  for(int i = 0; i < histnames.size(); i++)
  {
    plot(histnames[i],yaxis_compressions[i],leg_xoffsets[i],leg_yoffsets[i],xaxis_titles[i],plotnames[i]);
  }
}