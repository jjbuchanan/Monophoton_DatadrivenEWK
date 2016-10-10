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

void plot(TString histname, Double_t yaxis_compression, Double_t leg_xoffset, Double_t leg_yoffset, TString xaxis_title, TString plotname)
{
  const int nBins = 5;
  double maxLikelihoodFactors[nBins] = { 1.47902,1.08201,0.894284,0.503651,1.46734 };
  double maxLikelihoodPlus[nBins] = { 0.294795,0.194683,0.217482,0.337619,1.35248 };
  double maxLikelihoodMinus[nBins] = { -0.261142,-0.174514,-0.191378,-0.277054,-0.892669 };

  TCanvas *c = new TCanvas("c", "canvas",700,640);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  c->SetLeftMargin(0.15);
  c->cd();

  TFile *f_data_0000 = new TFile("ZllG_data_all.root");
  TH1F* histo_dilep_mass_data0000 = (TH1F*)f_data_0000->Get(histname);
  cout<<"Data: "<<(histo_dilep_mass_data0000->Integral()+histo_dilep_mass_data0000->GetBinContent(0)+histo_dilep_mass_data0000->GetBinContent(nBins+1))<<endl;
  histo_dilep_mass_data0000->SetStats(0);
  histo_dilep_mass_data0000->SetLineWidth(3);
  histo_dilep_mass_data0000->SetLineColor(kWhite);
  histo_dilep_mass_data0000->SetMarkerStyle(kFullSquare);
  histo_dilep_mass_data0000->SetMarkerColor(kWhite);
  histo_dilep_mass_data0000->SetTitle("");
  histo_dilep_mass_data0000->GetXaxis()->SetTitle("");
  histo_dilep_mass_data0000->GetXaxis()->SetTickLength(0);
  histo_dilep_mass_data0000->GetXaxis()->SetLabelOffset(999);
  histo_dilep_mass_data0000->GetYaxis()->SetTitle("");
  histo_dilep_mass_data0000->GetYaxis()->SetTickLength(0);
  histo_dilep_mass_data0000->GetYaxis()->SetLabelOffset(999);
  //Draw data first in white to establish axis boundaries
  histo_dilep_mass_data0000->Draw();

  c->Update();
  double ymax_data = c->GetUymax();

  TFile *f_wenu = new TFile("ZllG_wenu_all.root");
  TH1F* histo_elefake = (TH1F*)f_wenu->Get(histname);
  histo_elefake->Scale(0.0245); //R_e
  cout<<"Electron fakes: "<<(histo_elefake->Integral()+histo_elefake->GetBinContent(0)+histo_elefake->GetBinContent(nBins+1))<<endl;
  histo_elefake->SetStats(0);

  TFile* f_ZllG = new TFile("ZllG_ZLLGJets.root");
  TH1F* histo_dilep_mass_ZllG = (TH1F*)f_ZllG->Get(histname);
  histo_dilep_mass_ZllG->SetStats(0);
  histo_dilep_mass_ZllG->Scale(1.0/440320.0);
  histo_dilep_mass_ZllG->Scale(12900.0);
  histo_dilep_mass_ZllG->Scale(0.143);
  cout<<"Z(ll)+Gamma: "<<(histo_dilep_mass_ZllG->Integral()+histo_dilep_mass_ZllG->GetBinContent(0)+histo_dilep_mass_ZllG->GetBinContent(nBins+1))<<endl;
  histo_dilep_mass_ZllG->SetTitle("");
  histo_dilep_mass_ZllG->GetXaxis()->SetTitle("");
  histo_dilep_mass_ZllG->GetXaxis()->SetTickLength(0);
  histo_dilep_mass_ZllG->GetXaxis()->SetLabelOffset(999);
  histo_dilep_mass_ZllG->GetYaxis()->SetTitle("");
  histo_dilep_mass_ZllG->GetYaxis()->SetTickLength(0);
  histo_dilep_mass_ZllG->GetYaxis()->SetLabelOffset(999);
  histo_dilep_mass_ZllG->SetFillColor(kOrange);

//  for(int i = 1; i <= nBins; i++)
//  {
//    double bincontent = histo_dilep_mass_ZllG->GetBinContent(i);
//    double scaleFactor = maxLikelihoodFactors[i-1];
//    bincontent = bincontent*scaleFactor;
//    histo_dilep_mass_ZllG->SetBinContent(i,bincontent);
//  }

  TFile *f_TTG = new TFile("ZllG_TTGJets.root");
  TH1F* histo_dilep_mass_TTG = (TH1F*)f_TTG->Get(histname);
  histo_dilep_mass_TTG->SetStats(0);
  histo_dilep_mass_TTG->Scale(12900.0);
  histo_dilep_mass_TTG->Scale(1.0/4874091.0);
  histo_dilep_mass_TTG->Scale(3.697);
  cout<<"tt+Gamma: "<<(histo_dilep_mass_TTG->Integral()+histo_dilep_mass_TTG->GetBinContent(0)+histo_dilep_mass_TTG->GetBinContent(nBins+1))<<endl;
  histo_dilep_mass_TTG->SetTitle("");
  histo_dilep_mass_TTG->GetXaxis()->SetTitle("");
  histo_dilep_mass_TTG->GetXaxis()->SetTickLength(0);
  histo_dilep_mass_TTG->GetXaxis()->SetLabelOffset(999);
  histo_dilep_mass_TTG->GetYaxis()->SetTitle("");
  histo_dilep_mass_TTG->GetYaxis()->SetTickLength(0);
  histo_dilep_mass_TTG->GetYaxis()->SetLabelOffset(999);
  histo_dilep_mass_TTG->SetFillColor(kGreen);

  TFile *f_Zll_100to200 = new TFile("ZllG_DYJetsToLL_HT-100to200.root");
  TFile *f_Zll_200to400 = new TFile("ZllG_DYJetsToLL_HT-200to400.root");
  TFile *f_Zll_400to600 = new TFile("ZllG_DYJetsToLL_HT-400to600.root");
  TFile *f_Zll_600toInf = new TFile("ZllG_DYJetsToLL_HT-600toInf.root");
  TH1F* histo_dilep_mass_100to200 = (TH1F*)f_Zll_100to200->Get(histname);
  TH1F* histo_dilep_mass_200to400 = (TH1F*)f_Zll_200to400->Get(histname);
  TH1F* histo_dilep_mass_400to600 = (TH1F*)f_Zll_400to600->Get(histname);
  TH1F* histo_dilep_mass_600toInf = (TH1F*)f_Zll_600toInf->Get(histname);
  histo_dilep_mass_100to200->SetStats(0);
  histo_dilep_mass_200to400->SetStats(0);
  histo_dilep_mass_400to600->SetStats(0);
  histo_dilep_mass_600toInf->SetStats(0);
  histo_dilep_mass_100to200->Scale(12900.0);
  histo_dilep_mass_200to400->Scale(12900.0);
  histo_dilep_mass_400to600->Scale(12900.0);
  histo_dilep_mass_600toInf->Scale(12900.0);
  histo_dilep_mass_100to200->Scale(1.0/8149435.0);
  histo_dilep_mass_200to400->Scale(1.0/8637733.0);
  histo_dilep_mass_400to600->Scale(1.0/2490961.0);
  histo_dilep_mass_600toInf->Scale(1.0/3989990.0);
  histo_dilep_mass_100to200->Scale(147.4);
  histo_dilep_mass_200to400->Scale(40.99);
  histo_dilep_mass_400to600->Scale(5.678);
  histo_dilep_mass_600toInf->Scale(2.198);
  cout<<"Z(ll) HT-100to200: "<<(histo_dilep_mass_100to200->Integral()+histo_dilep_mass_100to200->GetBinContent(0)+histo_dilep_mass_100to200->GetBinContent(nBins+1))<<endl;
  cout<<"Z(ll) HT-200to400: "<<(histo_dilep_mass_200to400->Integral()+histo_dilep_mass_200to400->GetBinContent(0)+histo_dilep_mass_200to400->GetBinContent(nBins+1))<<endl;
  cout<<"Z(ll) HT-400to600: "<<(histo_dilep_mass_400to600->Integral()+histo_dilep_mass_400to600->GetBinContent(0)+histo_dilep_mass_400to600->GetBinContent(nBins+1))<<endl;
  cout<<"Z(ll) HT-600toInf: "<<(histo_dilep_mass_600toInf->Integral()+histo_dilep_mass_600toInf->GetBinContent(0)+histo_dilep_mass_600toInf->GetBinContent(nBins+1))<<endl;
  histo_dilep_mass_200to400->Add(histo_dilep_mass_100to200);
  histo_dilep_mass_200to400->Add(histo_dilep_mass_400to600);
  histo_dilep_mass_200to400->Add(histo_dilep_mass_600toInf);
  cout<<"Z(ll) combined: "<<(histo_dilep_mass_200to400->Integral()+histo_dilep_mass_200to400->GetBinContent(0)+histo_dilep_mass_200to400->GetBinContent(nBins+1))<<endl;
  histo_dilep_mass_200to400->SetTitle("");
  histo_dilep_mass_200to400->GetXaxis()->SetTitle("");
  histo_dilep_mass_200to400->GetXaxis()->SetTickLength(0);
  histo_dilep_mass_200to400->GetXaxis()->SetLabelOffset(999);
  histo_dilep_mass_200to400->GetYaxis()->SetTitle("");
  histo_dilep_mass_200to400->GetYaxis()->SetTickLength(0);
  histo_dilep_mass_200to400->GetYaxis()->SetLabelOffset(999);
  histo_dilep_mass_200to400->SetFillColor(kPink+1);

  TFile *f_WG = new TFile("ZllG_WGJets.root");
  TH1F* histo_dilep_mass_WG = (TH1F*)f_WG->Get(histname);
  histo_dilep_mass_WG->SetStats(0);
  histo_dilep_mass_WG->Scale(12900.0);
  histo_dilep_mass_WG->Scale(1.0/483496.0);
  histo_dilep_mass_WG->Scale(0.6637);
  cout<<"WG: "<<(histo_dilep_mass_WG->Integral()+histo_dilep_mass_WG->GetBinContent(0)+histo_dilep_mass_WG->GetBinContent(nBins+1))<<endl;

  TFile* f_WWG = new TFile("ZllG_WWG.root");
  TH1F* histo_dilep_mass_WWG = (TH1F*)f_WWG->Get(histname);
  histo_dilep_mass_WWG->SetStats(0);
  histo_dilep_mass_WWG->Scale(12900.0);
  histo_dilep_mass_WWG->Scale(1.0/999400.0);
  histo_dilep_mass_WWG->Scale(0.2147);
  cout<<"WWG: "<<(histo_dilep_mass_WWG->Integral()+histo_dilep_mass_WWG->GetBinContent(0)+histo_dilep_mass_WWG->GetBinContent(nBins+1))<<endl;

  TFile* f_GJets_40to100 = new TFile("ZllG_GJets_HT-40To100.root");
  TFile* f_GJets_100to200 = new TFile("ZllG_GJets_HT-100To200.root");
  TFile* f_GJets_200to400 = new TFile("ZllG_GJets_HT-200To400.root");
  TFile* f_GJets_400to600 = new TFile("ZllG_GJets_HT-400To600.root");  
  TFile* f_GJets = new TFile("ZllG_GJets_HT-600ToInf.root");
  TH1F* histo_GJets_40to100 = (TH1F*)f_GJets_40to100->Get(histname);
  TH1F* histo_GJets_100to200 = (TH1F*)f_GJets_100to200->Get(histname);
  TH1F* histo_GJets_200to400 = (TH1F*)f_GJets_200to400->Get(histname);
  TH1F* histo_GJets_400to600 = (TH1F*)f_GJets_400to600->Get(histname);
  TH1F* histo_dilep_mass_GJets = (TH1F*)f_GJets->Get(histname);
  histo_GJets_40to100->SetStats(0);
  histo_GJets_100to200->SetStats(0);
  histo_GJets_200to400->SetStats(0);
  histo_GJets_400to600->SetStats(0);
  histo_dilep_mass_GJets->SetStats(0);
  histo_GJets_40to100->Scale(12900.0);
  histo_GJets_100to200->Scale(12900.0);
  histo_GJets_200to400->Scale(12900.0);
  histo_GJets_400to600->Scale(12900.0);
  histo_dilep_mass_GJets->Scale(12900.0);
  histo_GJets_40to100->Scale(1.0/4468724.0);
  histo_GJets_100to200->Scale(1.0/5142782.0);
  histo_GJets_200to400->Scale(1.0/10281790.0);
  histo_GJets_400to600->Scale(1.0/2480186.0);
  histo_dilep_mass_GJets->Scale(1.0/2416755.0);
  histo_GJets_40to100->Scale(20790.0);
  histo_GJets_100to200->Scale(9238.0);
  histo_GJets_200to400->Scale(2305.0);
  histo_GJets_400to600->Scale(274.4);
  histo_dilep_mass_GJets->Scale(93.46);
  cout<<"Gamma+jets HT-40To100: "<<(histo_GJets_40to100->Integral()+histo_GJets_40to100->GetBinContent(0)+histo_GJets_40to100->GetBinContent(nBins+1))<<endl;
  cout<<"Gamma+jets HT-100To200: "<<(histo_GJets_100to200->Integral()+histo_GJets_100to200->GetBinContent(0)+histo_GJets_100to200->GetBinContent(nBins+1))<<endl;
  cout<<"Gamma+jets HT-200To400: "<<(histo_GJets_200to400->Integral()+histo_GJets_200to400->GetBinContent(0)+histo_GJets_200to400->GetBinContent(nBins+1))<<endl;
  cout<<"Gamma+jets HT-400To600: "<<(histo_GJets_400to600->Integral()+histo_GJets_400to600->GetBinContent(0)+histo_GJets_400to600->GetBinContent(nBins+1))<<endl;
  cout<<"Gamma+jets HT-600ToInf: "<<(histo_dilep_mass_GJets->Integral()+histo_dilep_mass_GJets->GetBinContent(0)+histo_dilep_mass_GJets->GetBinContent(nBins+1))<<endl;
  histo_dilep_mass_GJets->Add(histo_GJets_40to100);
  histo_dilep_mass_GJets->Add(histo_GJets_100to200);
  histo_dilep_mass_GJets->Add(histo_GJets_200to400);
  histo_dilep_mass_GJets->Add(histo_GJets_400to600);
  cout<<"Gamma+jets combined: "<<(histo_dilep_mass_GJets->Integral()+histo_dilep_mass_GJets->GetBinContent(0)+histo_dilep_mass_GJets->GetBinContent(nBins+1))<<endl;

  histo_dilep_mass_WG->Add(histo_elefake);
  histo_dilep_mass_WG->Add(histo_dilep_mass_WWG);
  // histo_dilep_mass_WG->Add(histo_dilep_mass_GJets);
  histo_dilep_mass_WG->SetTitle("");
  histo_dilep_mass_WG->GetXaxis()->SetTitle("");
  histo_dilep_mass_WG->GetXaxis()->SetTickLength(0);
  histo_dilep_mass_WG->GetXaxis()->SetLabelOffset(999);
  histo_dilep_mass_WG->GetYaxis()->SetTitle("");
  histo_dilep_mass_WG->GetYaxis()->SetTickLength(0);
  histo_dilep_mass_WG->GetYaxis()->SetLabelOffset(999);
  histo_dilep_mass_WG->SetFillColor(kBlue);

  THStack *stackHisto = new THStack("stackHisto","Title");
  stackHisto->Add(histo_dilep_mass_WG);
  stackHisto->Add(histo_dilep_mass_200to400);
  stackHisto->Add(histo_dilep_mass_TTG);
  stackHisto->Add(histo_dilep_mass_ZllG);
  stackHisto->SetTitle("");
  stackHisto->Draw("HIST");

  c->Update();
  double ymax_background = c->GetUymax();

  //Accommodate both the data and background plots
  if(ymax_background > ymax_data)
    ymax_data = ymax_background;
  ymax_data /= yaxis_compression;
  histo_dilep_mass_data0000->GetYaxis()->SetRangeUser(0.,ymax_data);
  histo_dilep_mass_data0000->Draw();
  stackHisto->Draw("HIST SAME");
  histo_dilep_mass_data0000->SetLineColor(kBlack);
  histo_dilep_mass_data0000->SetMarkerColor(kBlack);
  histo_dilep_mass_data0000->Draw("SAME");

  TLegend* leg = new TLegend(0.368195+leg_xoffset,0.690468+leg_yoffset,0.755014+leg_xoffset,0.872757+leg_yoffset,"");
  leg->AddEntry(histo_dilep_mass_data0000,"Data");
  leg->AddEntry(histo_dilep_mass_ZllG,"Z(ll)#gamma","F");
  leg->AddEntry(histo_dilep_mass_TTG,"tt#gamma","F");
  leg->AddEntry(histo_dilep_mass_200to400,"Z(ll)","F");
  leg->AddEntry(histo_dilep_mass_WG,"Electron#rightarrow#gamma MisID, WW#gamma","F");
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->Draw();

  TLatex *texS = new TLatex(0.54023,0.907173,"#sqrt{s} = 13 TeV, 12.9 fb^{-1}");
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
  yaxis->SetTitle("Events / bin");
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

  c->SaveAs(TString("zllg_"+plotname+".png"));
  c->SaveAs(TString("zllg_"+plotname+".pdf"));
  delete(c);
}

void zllg_plotter()
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

  histnames.push_back(TString("Photon_SCeta_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #eta"));
  plotnames.push_back(TString("phoEta"));

  histnames.push_back(TString("Photon_SCphi_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #phi"));
  plotnames.push_back(TString("phoPhi"));

  histnames.push_back(TString("pfMET_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("pfMET [GeV]"));
  plotnames.push_back(TString("pfMET"));

  histnames.push_back(TString("nJet_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Number of Jets"));
  plotnames.push_back(TString("nJet"));

  histnames.push_back(TString("h_photonic_recoil_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photonic Recoil [GeV]"));
  plotnames.push_back(TString("recoil"));

  histnames.push_back(TString("h_dPhi_phoRecoil_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("#Delta#phi(Photon,Recoil)"));
  plotnames.push_back(TString("dPhiPhoRecoil"));

  histnames.push_back(TString("h_leadingLeptonPt_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("elePt_leading"));

  histnames.push_back(TString("h_leadingLeptonEta_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #eta"));
  plotnames.push_back(TString("leadingEleEta"));

  histnames.push_back(TString("h_leadingLeptonPhi_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #phi"));
  plotnames.push_back(TString("leadingElePhi"));

  histnames.push_back(TString("h_subleadingLeptonPt_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("elePt_subleading"));

  histnames.push_back(TString("h_subleadingLeptonEta_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #eta"));
  plotnames.push_back(TString("subleadingEleEta"));

  histnames.push_back(TString("h_subleadingLeptonPhi_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #phi"));
  plotnames.push_back(TString("subleadingElePhi"));

  histnames.push_back(TString("h_dileptonPt_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("dielePt"));

  histnames.push_back(TString("h_dileptonM_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron Invariant Mass [GeV]"));
  plotnames.push_back(TString("dieleM"));

  histnames.push_back(TString("h_phoPT_over_dileptonPt_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} / Dielectron #it{p}_{T}"));
  plotnames.push_back(TString("phoPtOverDielectronPt"));

  histnames.push_back(TString("h_phoPT_over_photonicRecoil_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} / Recoil"));
  plotnames.push_back(TString("phoPtOverRecoil"));

  histnames.push_back(TString("h_dileptonPt_over_pfMET_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron #it{p}_{T} / pfMET"));
  plotnames.push_back(TString("dielePtOverMET"));

  histnames.push_back(TString("h_photonicRecoil_over_pfMET_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Recoil / pfMET"));
  plotnames.push_back(TString("recoilOverMET"));

  histnames.push_back(TString("h_photonicRecoil_over_dileptonPt_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Recoil / Dielectron #it{p}_{T}"));
  plotnames.push_back(TString("recoilOverDielePt"));

  for(int i = 0; i < histnames.size(); i++)
  {
    plot(histnames[i],yaxis_compressions[i],leg_xoffsets[i],leg_yoffsets[i],xaxis_titles[i],plotnames[i]);
  }
}