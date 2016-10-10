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

void plot(TString histname, Double_t yaxis_compress, Double_t leg_xoffset, Double_t leg_yoffset, TString xaxis_title, TString plotname)
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

  TFile *f_data = new TFile("ZnnG_data_all.root");
  TH1F* histo_data = (TH1F*)f_data->Get(histname);
  cout<<"Data: "<<(histo_data->Integral()+histo_data->GetBinContent(0)+histo_data->GetBinContent(nBins+1))<<endl;
  histo_data->SetStats(0);
  histo_data->SetLineWidth(3);
  histo_data->SetLineColor(kWhite);
  histo_data->SetMarkerStyle(kFullSquare);
  histo_data->SetMarkerColor(kWhite);
  histo_data->SetTitle("");
  histo_data->GetXaxis()->SetTitle("");
  histo_data->GetXaxis()->SetTickLength(0);
  histo_data->GetXaxis()->SetLabelOffset(999);
  histo_data->GetYaxis()->SetTitle("");
  histo_data->GetYaxis()->SetTickLength(0);
  histo_data->GetYaxis()->SetLabelOffset(999);
  //Draw data first in white to establish axis boundaries
  histo_data->Draw();

  c->Update();
  double ymax_data = c->GetUymax();

  TFile *f_jetfake = new TFile("ZnnG_qcd_all.root");
  TH1F* histo_jetfake = (TH1F*)f_jetfake->Get(histname);
  cout<<"QCD fakes: "<<(histo_jetfake->Integral()+histo_jetfake->GetBinContent(0)+histo_jetfake->GetBinContent(nBins+1))<<endl;
  histo_jetfake->SetStats(0);
  histo_jetfake->SetTitle("");
  histo_jetfake->GetXaxis()->SetTitle("");
  histo_jetfake->GetXaxis()->SetTickLength(0);
  histo_jetfake->GetXaxis()->SetLabelOffset(999);
  histo_jetfake->GetYaxis()->SetTitle("");
  histo_jetfake->GetYaxis()->SetTickLength(0);
  histo_jetfake->GetYaxis()->SetLabelOffset(999);
  histo_jetfake->SetFillColor(kBlue-4);

  TFile *f_elefake = new TFile("ZnnG_wenu_all.root");
  TH1F* histo_elefake = (TH1F*)f_elefake->Get(histname);
  histo_elefake->Scale(0.0245); //R_e
  cout<<"Electron fakes: "<<(histo_elefake->Integral()+histo_elefake->GetBinContent(0)+histo_elefake->GetBinContent(nBins+1))<<endl;
  histo_elefake->SetStats(0);
  histo_elefake->SetTitle("");
  histo_elefake->GetXaxis()->SetTitle("");
  histo_elefake->GetXaxis()->SetTickLength(0);
  histo_elefake->GetXaxis()->SetLabelOffset(999);
  histo_elefake->GetYaxis()->SetTitle("");
  histo_elefake->GetYaxis()->SetTickLength(0);
  histo_elefake->GetYaxis()->SetLabelOffset(999);
  histo_elefake->SetFillColor(kMagenta+1);

  TFile* f_WG = new TFile("ZnnG_WGJets.root");
  TH1F* histo_WG = (TH1F*)f_WG->Get(histname);
  histo_WG->SetStats(0);
  histo_WG->Scale(1.0/483496.0);
  histo_WG->Scale(12900.0);
  histo_WG->Scale(0.6637);
  cout<<"W+gamma: "<<(histo_WG->Integral()+histo_WG->GetBinContent(0)+histo_WG->GetBinContent(nBins+1))<<endl;
  histo_WG->SetTitle("");
  histo_WG->GetXaxis()->SetTitle("");
  histo_WG->GetXaxis()->SetTickLength(0);
  histo_WG->GetXaxis()->SetLabelOffset(999);
  histo_WG->GetYaxis()->SetTitle("");
  histo_WG->GetYaxis()->SetTickLength(0);
  histo_WG->GetYaxis()->SetLabelOffset(999);
  histo_WG->SetFillColor(kAzure+10);

//  for(int i = 1; i <= nBins; i++)
//  {
//    double bincontent = histo_WG->GetBinContent(i);
//    double scaleFactor = maxLikelihoodFactors[i-1];
//    bincontent = bincontent*scaleFactor;
//    histo_WG->SetBinContent(i,bincontent);
//  }

  TFile *f_GJets_40to100 = new TFile("ZnnG_GJets_HT-40To100.root");
  TFile *f_GJets_100to200 = new TFile("ZnnG_GJets_HT-100To200.root");
  TFile *f_GJets_200to400 = new TFile("ZnnG_GJets_HT-200To400.root");
  TFile *f_GJets_400to600 = new TFile("ZnnG_GJets_HT-400To600.root");
  TFile *f_GJets_600toInf = new TFile("ZnnG_GJets_HT-600ToInf.root");
  TH1F* histo_GJets_40to100 = (TH1F*)f_GJets_40to100->Get(histname);
  TH1F* histo_GJets_100to200 = (TH1F*)f_GJets_100to200->Get(histname);
  TH1F* histo_GJets_200to400 = (TH1F*)f_GJets_200to400->Get(histname);
  TH1F* histo_GJets_400to600 = (TH1F*)f_GJets_400to600->Get(histname);
  TH1F* histo_GJets_600toInf = (TH1F*)f_GJets_600toInf->Get(histname);
  //Lumi
  histo_GJets_40to100->Scale(12900.0);
  histo_GJets_100to200->Scale(12900.0);
  histo_GJets_200to400->Scale(12900.0);
  histo_GJets_400to600->Scale(12900.0);
  histo_GJets_600toInf->Scale(12900.0);
  //Events inspected
  histo_GJets_40to100->Scale(1.0/4468724.0);
  histo_GJets_100to200->Scale(1.0/5142782.0);
  histo_GJets_200to400->Scale(1.0/10281790.0);
  histo_GJets_400to600->Scale(1.0/2480186.0);
  histo_GJets_600toInf->Scale(1.0/2416755.0);
  //Cross section
  histo_GJets_40to100->Scale(20790);
  histo_GJets_100to200->Scale(9238);
  histo_GJets_200to400->Scale(2305);
  histo_GJets_400to600->Scale(274.4);
  histo_GJets_600toInf->Scale(93.46);
  cout<<"Gamma+Jets 40to100: "<<(histo_GJets_40to100->Integral()+histo_GJets_40to100->GetBinContent(0)+histo_GJets_40to100->GetBinContent(nBins+1))<<endl;
  cout<<"Gamma+Jets 100to200: "<<(histo_GJets_100to200->Integral()+histo_GJets_100to200->GetBinContent(0)+histo_GJets_100to200->GetBinContent(nBins+1))<<endl;
  cout<<"Gamma+Jets 200to400: "<<(histo_GJets_200to400->Integral()+histo_GJets_200to400->GetBinContent(0)+histo_GJets_200to400->GetBinContent(nBins+1))<<endl;
  cout<<"Gamma+Jets 400to600: "<<(histo_GJets_400to600->Integral()+histo_GJets_400to600->GetBinContent(0)+histo_GJets_400to600->GetBinContent(nBins+1))<<endl;
  cout<<"Gamma+Jets 600toInf: "<<(histo_GJets_600toInf->Integral()+histo_GJets_600toInf->GetBinContent(0)+histo_GJets_600toInf->GetBinContent(nBins+1))<<endl;
  histo_GJets_40to100->Add(histo_GJets_100to200);
  histo_GJets_40to100->Add(histo_GJets_200to400);
  histo_GJets_40to100->Add(histo_GJets_400to600);
  histo_GJets_40to100->Add(histo_GJets_600toInf);
  cout<<"Gamma+Jets combined: "<<(histo_GJets_40to100->Integral()+histo_GJets_40to100->GetBinContent(0)+histo_GJets_40to100->GetBinContent(nBins+1))<<endl;
  histo_GJets_40to100->SetTitle("");
  histo_GJets_40to100->GetXaxis()->SetTitle("");
  histo_GJets_40to100->GetXaxis()->SetTickLength(0);
  histo_GJets_40to100->GetXaxis()->SetLabelOffset(999);
  histo_GJets_40to100->GetYaxis()->SetTitle("");
  histo_GJets_40to100->GetYaxis()->SetTickLength(0);
  histo_GJets_40to100->GetYaxis()->SetLabelOffset(999);
  //histo_GJets_40to100->SetFillColor(TColor::GetColor("#FFFFCC"));
  histo_GJets_40to100->SetFillColor(43);

  TFile *f_ZllG = new TFile("ZnnG_ZLLGJets.root");
  TH1F* histo_ZllG = (TH1F*)f_ZllG->Get(histname);
  histo_ZllG->SetStats(0);
  histo_ZllG->Scale(12900.0);
  histo_ZllG->Scale(1.0/440320.0);
  histo_ZllG->Scale(0.143);
  cout<<"Z(ll)+Gamma: "<<(histo_ZllG->Integral()+histo_ZllG->GetBinContent(0)+histo_ZllG->GetBinContent(nBins+1))<<endl;
  histo_ZllG->SetTitle("");
  histo_ZllG->GetXaxis()->SetTitle("");
  histo_ZllG->GetXaxis()->SetTickLength(0);
  histo_ZllG->GetXaxis()->SetLabelOffset(999);
  histo_ZllG->GetYaxis()->SetTitle("");
  histo_ZllG->GetYaxis()->SetTickLength(0);
  histo_ZllG->GetYaxis()->SetLabelOffset(999);
  histo_ZllG->SetFillColor(kSpring-8);

  TFile *f_TTG = new TFile("ZnnG_TTGJets.root");
  TH1F* histo_TTG = (TH1F*)f_TTG->Get(histname);
  histo_TTG->SetStats(0);
  histo_TTG->Scale(12900.0);
  histo_TTG->Scale(1.0/4874091.0);
  histo_TTG->Scale(3.697);
  cout<<"tt+Gamma: "<<(histo_TTG->Integral()+histo_TTG->GetBinContent(0)+histo_TTG->GetBinContent(nBins+1))<<endl;
  histo_TTG->SetTitle("");
  histo_TTG->GetXaxis()->SetTitle("");
  histo_TTG->GetXaxis()->SetTickLength(0);
  histo_TTG->GetXaxis()->SetLabelOffset(999);
  histo_TTG->GetYaxis()->SetTitle("");
  histo_TTG->GetYaxis()->SetTickLength(0);
  histo_TTG->GetYaxis()->SetLabelOffset(999);
  histo_TTG->SetFillColor(kOrange);

  TFile *f_TG = new TFile("ZnnG_TGJets.root");
  TH1F* histo_TG = (TH1F*)f_TG->Get(histname);
  histo_TG->SetStats(0);
  histo_TG->Scale(12900.0);
  histo_TG->Scale(1.0/1449833.0);
  histo_TG->Scale(2.967);
  cout<<"t+Gamma: "<<(histo_TG->Integral()+histo_TG->GetBinContent(0)+histo_TG->GetBinContent(nBins+1))<<endl;
  histo_TG->SetTitle("");
  histo_TG->GetXaxis()->SetTitle("");
  histo_TG->GetXaxis()->SetTickLength(0);
  histo_TG->GetXaxis()->SetLabelOffset(999);
  histo_TG->GetYaxis()->SetTitle("");
  histo_TG->GetYaxis()->SetTickLength(0);
  histo_TG->GetYaxis()->SetLabelOffset(999);
  histo_TG->SetFillColor(kBlue-4);

  TFile* f_WWG = new TFile("ZnnG_WWG.root");
  TH1F* histo_WWG = (TH1F*)f_WWG->Get(histname);
  histo_WWG->SetStats(0);
  histo_WWG->Scale(12900.0);
  histo_WWG->Scale(1.0/999400.0);
  histo_WWG->Scale(0.2147);
  cout<<"WWG: "<<(histo_WWG->Integral()+histo_WWG->GetBinContent(0)+histo_WWG->GetBinContent(nBins+1))<<endl;
  histo_WWG->SetFillColor(kBlue-4);

// TFile* f_Diphoton = new TFile("ZnnG_DiPhotonJets.root");
// TH1F* histo_diphoton = (TH1F*)f_Diphoton->Get(histname);
// histo_diphoton->SetStats(0);
// histo_diphoton->Scale(12900.0);
// histo_diphoton->Scale(1.0/3648565.0);
// histo_diphoton->Scale(135.1);
// cout<<"Diphoton: "<<(histo_diphoton->Integral()+histo_diphoton->GetBinContent(0)+histo_diphoton->GetBinContent(nBins+1))<<endl;
// histo_diphoton->SetFillColor(kOrange+10);

  histo_WWG->Add(histo_TG);
  histo_WWG->Add(histo_jetfake);

  THStack *stackHisto = new THStack("stackHisto","Title");
  stackHisto->Add(histo_WWG);
  stackHisto->Add(histo_elefake);
  stackHisto->Add(histo_ZllG);
  stackHisto->Add(histo_GJets_40to100);
  stackHisto->Add(histo_TTG);
//  stackHisto->Add(histo_diphoton);
  stackHisto->Add(histo_WG);
  stackHisto->SetTitle("");
  stackHisto->Draw("HIST");

  c->Update();
  double ymax_background = c->GetUymax();

  //Accommodate both the data and background plots
  if(ymax_background > ymax_data)
    ymax_data = ymax_background;
  ymax_data /= yaxis_compress;
  histo_data->GetYaxis()->SetRangeUser(0.,ymax_data);
  histo_data->Draw();
  stackHisto->Draw("HIST SAME");
  histo_data->SetLineColor(kBlack);
  histo_data->SetMarkerColor(kBlack);
  histo_data->Draw("SAME");

  //Central location of leg defined to be location of leg in phoPt plot
  TLegend* leg = new TLegend(0.5+leg_xoffset,0.58075+leg_yoffset,0.885387+leg_xoffset,0.862969+leg_yoffset,"");
  leg->AddEntry(histo_data,"Data");
  leg->AddEntry(histo_WG,"W#gamma#rightarrowl#nu#gamma","F");
//  leg->AddEntry(histo_diphoton,"#gamma#gamma","F");
  leg->AddEntry(histo_TTG,"tt#gamma","F");
  leg->AddEntry(histo_GJets_40to100,"#gamma+jet","F");
  leg->AddEntry(histo_ZllG,"Z(ll)#gamma","F");
  leg->AddEntry(histo_elefake,"e#rightarrow#gamma MisID","F");
  leg->AddEntry(histo_WWG,"jet#rightarrow#gamma MisID,WW#gamma, t#gamma","F");
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

  c->SaveAs(TString("znng_"+plotname+".png"));
  c->SaveAs(TString("znng_"+plotname+".pdf"));
  delete(c);
}

void znng_plotter()
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

  histnames.push_back(TString("h_photon_SCEta_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #eta"));
  plotnames.push_back(TString("phoEta"));

  histnames.push_back(TString("h_photon_SCPhi_4"));
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

  histnames.push_back(TString("h_leptonPt_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("muPt"));

  histnames.push_back(TString("h_leptonEta_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon #eta"));
  plotnames.push_back(TString("muEta"));

  histnames.push_back(TString("h_leptonPhi_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon #phi"));
  plotnames.push_back(TString("muPhi"));

  histnames.push_back(TString("h_phoPT_over_photonicRecoil_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} / Recoil"));
  plotnames.push_back(TString("phoPtOverRecoil"));

  histnames.push_back(TString("h_leptonPt_over_pfMET_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon #it{p}_{T} / pfMET"));
  plotnames.push_back(TString("muPtOverpfMET"));

  histnames.push_back(TString("h_lepMET_MT_4"));
  yaxis_compressions.push_back(1.);
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Muon + pfMET m_{T} [GeV]"));
  plotnames.push_back(TString("muMETmT"));

  for(int i = 0; i < histnames.size(); i++)
  {
    plot(histnames[i],yaxis_compressions[i],leg_xoffsets[i],leg_yoffsets[i],xaxis_titles[i],plotnames[i]);
  }
}