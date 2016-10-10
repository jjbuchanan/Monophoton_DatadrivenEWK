////postAnalyzer.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom postAnalyzer analyze
//
//To run, assuming this is compiled to an executable named 'analyze':
//$ ./analyze /hdfs/store/user/jjbuch/LatestNtuples/ /afs/hep.wisc.edu/user/jjbuchanan/private/CMSSW_7_4_9/src/output.root -1 10000
//Runs over every event in the folder LatestNtuples, reporting progress every 10000 events
//and storing the resulting histograms in the file output.root.
//
//To plot, for example, single photon trigger efficiency as a function of photon pt:
//$ root -l
//root[0] TFile *f = new TFile("output.root");
//root[1] TGraphAsymmErrors *efficiency = new TGraphAsymmErrors((TH1F*)f->Get("Photon_Et_300_2"),(TH1F*)f->Get("Photon_Et_300_1"));
//root[2] efficiency->Draw("AP")
//root[3] efficiency->SetTitle("Single photon trigger efficiency")
//root[4] efficiency->GetXaxis()->SetTitle("Photon p_{T}")
//root[5] efficiency->Draw("AP")
//

#define postAnalyzer_cxx
#include "postAnalyzer_ZllG_mc_wg.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <set>
using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
  {
    std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
    return 1;
  }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
  {
    std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
    return 1;
  }
  postAnalyzer t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void postAnalyzer::Loop(Long64_t maxEvents, int reportEvery)
{
  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;
  float passed170;
  passed170 = 0.0;
  float totalInspectedSum_QCDscale_muR1_muF1 = 0.0;
  float sum170_QCDscale_muR1_muF1, sum170_QCDscale_muR1_muF2, sum170_QCDscale_muR1_muF0p5, sum170_QCDscale_muR2_muF1, sum170_QCDscale_muR2_muF2, sum170_QCDscale_muR2_muF0p5, sum170_QCDscale_muR0p5_muF1, sum170_QCDscale_muR0p5_muF2, sum170_QCDscale_muR0p5_muF0p5;
  sum170_QCDscale_muR1_muF1 = sum170_QCDscale_muR1_muF2 = sum170_QCDscale_muR1_muF0p5 = sum170_QCDscale_muR2_muF1 = sum170_QCDscale_muR2_muF2 = sum170_QCDscale_muR2_muF0p5 = sum170_QCDscale_muR0p5_muF1 = sum170_QCDscale_muR0p5_muF2 = sum170_QCDscale_muR0p5_muF0p5 = 0;
  float central_sum, central_sum_passing;
  central_sum = central_sum_passing = 0.0;

  int initialIndex = -1;//Should end up being 9.
  bool initialIndexNotSet = true;
  int nMCreplicas = 101;
  string initialID = "111"; //111 for basically everything: NNPDF30_lo_as_0130_nf_4 (LHAID 263400)
  std::vector<int> vecIndices;//The vector indices of all the MC replicas. Should just increase sequentially after the initial value.
  vecIndices.clear();
  std::vector<float> sum_passing;//Sum of weights of passing events for each MC replica.
  sum_passing.clear();

  std::vector<int> phoCand1;
  phoCand1.clear();

  std::vector<int> jetveto;
  jetveto.clear();

  bool debug=true;
  Long64_t nentries = fChain->GetEntries();
  std::cout<<"nentries:"<<nentries<<std::endl;
  //Look at up to maxEvents events, or all if maxEvents == -1.
  Long64_t nentriesToCheck = nentries;
  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;

  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  TStopwatch sw;
  sw.Start();
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
  {
    
    event_.clear();
    event_info.clear();
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //=1.0 for real data
    double event_weight=1.0;
    double EWK_corrected_weight=1.0;
    double NNLO_weight = 1.0;

    if(initialIndexNotSet)
    {
      //Match the initial vector index to the specified MC replica ID.
      for(int i = 0; i < lheWeightIDs->size(); i++)
      {
        if(lheWeightIDs->at(i) == initialID)
        {
          initialIndex = i;
          break;
        }
      }
      //DEBUG
      cout<<"initialIndex = "<<initialIndex<<endl;
      //Set up all the vector indices, and sum_passing while we're at it.
      for(int i = 0; i < nMCreplicas; i++)
      {
        vecIndices.push_back(initialIndex + i);
        sum_passing.push_back(0.0);
      }
      //Now vecIndices.size() and sum_passing.size() should both == nMCreplicas == 101
      //DEBUG
      for(int i = 0; i < vecIndices.size(); i++)
        cout<<"vector index: "<<vecIndices[i]<<", ID: "<<lheWeightIDs->at(vecIndices[i])<<endl;
      cout<<endl;
      //Don't initialize more than once.
      initialIndexNotSet = false;
    }

    totalInspectedSum_QCDscale_muR1_muF1 += genWeight_QCDscale_muR1_muF1;
    central_sum += lheNormalizedWeights->at(0);
    
    //Sequential cuts
    
    phoCand1   = getPhoCand(175,1.4442,1);
    
    if(phoCand1.size()>0)
      {
        jetveto = JetVetoDecision(phoCand1[0]);
      }

    //if(metFilters==0)
    //{
    //  nMETFiltersPassed++;
    //  if((HLTPho>>7&1 == 1)||(HLTPho>>8&1 == 1)|| (HLTPho>>9&1 == 1) || (HLTPho>>10&1 == 1) || (HLTPho>>11&1 == 1) || (HLTPho>>12&1 == 1) || (HLTPho>>22&1 == 1))
    //  {
    //    nHLTPassed++;
          if(phoCand1.size() >0)
          {
            // event_weight*=(1.013 - 0.0001168*phoEt->at(phoCand1[0]));
            if( TMath::Max( ( (*phoPFChWorstIso)[phoCand1[0]]  - rho*EAchargedworst((*phoSCEta)[phoCand1[0]]) ), 0.0) < 1.37 )
            {
              std::vector<int> elelist = electron_veto_looseID(phoCand1[0],10.0);
              std::vector<int> mulist = muon_veto_looseID(phoCand1[0],10.0);
              std::vector<int> leplist;
              leplist.clear();
              bool elePairSet = false;
              TLorentzVector e1, e2;
              bool muPairSet = false;
              TLorentzVector m1, m2;
              if(elelist.size() > 1)
              {
                for(int i=1; i<elelist.size(); ++i)
                {
                  if(eleCharge->at(0)*eleCharge->at(i) == -1)
                  {
                    e1.SetPtEtaPhiE(elePt->at(elelist[0]),eleEta->at(elelist[0]),elePhi->at(elelist[0]),eleEn->at(elelist[0]));
                    e2.SetPtEtaPhiE(elePt->at(elelist[i]),eleEta->at(elelist[i]),elePhi->at(elelist[i]),eleEn->at(elelist[i]));
                    elePairSet = true;
                    break;
                  }
                }
              }
              if(mulist.size() > 1)
              {
                for(int i=1; i<mulist.size(); ++i)
                {
                  if(muCharge->at(0)*muCharge->at(i) == -1)
                  {
                    m1.SetPtEtaPhiE(muPt->at(mulist[0]),muEta->at(mulist[0]),muPhi->at(mulist[0]),muEn->at(mulist[0]));
                    m2.SetPtEtaPhiE(muPt->at(mulist[i]),muEta->at(mulist[i]),muPhi->at(mulist[i]),muEn->at(mulist[i]));
                    muPairSet = true;
                    break;
                  }
                }
              }
              
              TLorentzVector l1, l2;
              if(elePairSet)
              {
                leplist = elelist;
                l1 = e1;
                l2 = e2;
                if(muPairSet)
                {
                  TLorentzVector ee = e1+e2;
                  Double_t dielectron_mass = ee.M();
                  TLorentzVector mm = m1+m2;
                  Double_t dimuon_mass = mm.M();
                  if(dimuon_mass > dielectron_mass)
                  {
                    leplist = mulist;
                    l1 = m1;
                    l2 = m2;
                    elePairSet = false; //Use elePairSet to discriminate whether an electron (true) or muon (false) pair was selected
                  }
                }
              }
              else if(muPairSet)
              {
                leplist = mulist;
                l1 = m1;
                l2 = m2;
              }

              TLorentzVector ll = l1+l2;
              Double_t dilepton_mass = ll.M();
              Double_t dilepton_pt = ll.Pt();
              TLorentzVector met_4vec;
              met_4vec.SetPtEtaPhiE(pfMET,0.,pfMETPhi,pfMET);
              TLorentzVector leptoMET_4vec = ll+met_4vec;
              Double_t leptoMET = leptoMET_4vec.Pt();
              Double_t leptoMET_phi = leptoMET_4vec.Phi();

              if(leptoMET > 170)
              {
                // Double_t EWK_percent_adjustment = ewkCorrection->GetBinContent(ewkCorrection->GetXaxis()->FindBin(phoEt->at(phoCand1[0])));
                // EWK_corrected_weight = event_weight*(1.0+.01*EWK_percent_adjustment);
                // NNLO_weight = event_weight*EWK_corrected_weight*NNLOCorrection(phoEt->at(phoCand1[0]));
                passed170+=NNLO_weight;
                fillHistos(nMCreplicas+9,NNLO_weight,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                //PDF
                central_sum_passing += NNLO_weight*lheNormalizedWeights->at(0);
                for(unsigned int i = 0; i < nMCreplicas; i++) //nMCreplicas is 101
                {
                  sum_passing[i] += NNLO_weight*lheNormalizedWeights->at(vecIndices[i]);
                  fillHistos(i,NNLO_weight*lheNormalizedWeights->at(vecIndices[i]),phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                }
                //Scale
                sum170_QCDscale_muR1_muF1 += NNLO_weight*genWeight_QCDscale_muR1_muF1;
                sum170_QCDscale_muR1_muF2 += NNLO_weight*genWeight_QCDscale_muR1_muF2;
                sum170_QCDscale_muR1_muF0p5 += NNLO_weight*genWeight_QCDscale_muR1_muF0p5;
                sum170_QCDscale_muR2_muF1 += NNLO_weight*genWeight_QCDscale_muR2_muF1;
                sum170_QCDscale_muR2_muF2 += NNLO_weight*genWeight_QCDscale_muR2_muF2;
                sum170_QCDscale_muR2_muF0p5 += NNLO_weight*genWeight_QCDscale_muR2_muF0p5;
                sum170_QCDscale_muR0p5_muF1 += NNLO_weight*genWeight_QCDscale_muR0p5_muF1;
                sum170_QCDscale_muR0p5_muF2 += NNLO_weight*genWeight_QCDscale_muR0p5_muF2;
                sum170_QCDscale_muR0p5_muF0p5 += NNLO_weight*genWeight_QCDscale_muR0p5_muF0p5;
                fillHistos(nMCreplicas,NNLO_weight*genWeight_QCDscale_muR1_muF1,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                fillHistos(nMCreplicas+1,NNLO_weight*genWeight_QCDscale_muR1_muF2,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                fillHistos(nMCreplicas+2,NNLO_weight*genWeight_QCDscale_muR1_muF0p5,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                fillHistos(nMCreplicas+3,NNLO_weight*genWeight_QCDscale_muR2_muF1,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                fillHistos(nMCreplicas+4,NNLO_weight*genWeight_QCDscale_muR2_muF2,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                fillHistos(nMCreplicas+5,NNLO_weight*genWeight_QCDscale_muR2_muF0p5,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                fillHistos(nMCreplicas+6,NNLO_weight*genWeight_QCDscale_muR0p5_muF1,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                fillHistos(nMCreplicas+7,NNLO_weight*genWeight_QCDscale_muR0p5_muF2,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
                fillHistos(nMCreplicas+8,NNLO_weight*genWeight_QCDscale_muR0p5_muF0p5,phoCand1[0],jetveto,leplist,dilepton_mass,dilepton_pt,leptoMET,leptoMET_phi,elePairSet);
              }
            }
          }
      //}
    //}

    
    
    tree->Fill();
    
    if (jentry%reportEvery == 0)
      {
        std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
      }
  }

  if((nentriesToCheck-1)%reportEvery != 0)
    std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
  sw.Stop();
  std::cout<<"All events checked."<<std::endl;
  //Report
  std::cout << "RealTime : " << sw.RealTime() / 60.0 << " minutes" << std::endl;
    std::cout << "CPUTime  : " << sw.CpuTime()  / 60.0 << " minutes" << std::endl;
  std::cout << std::endl;
  std::cout<<"Total number of events: "<<nTotal<<std::endl;
  std::cout << "Number of events inspected: " << nTotal << std::endl;
  std::cout<<std::endl;
  cout<<"Scale variations"<<endl;
  cout<<"passed170: "<<passed170<<endl;
  double central_estimate_of_acceptance = (double)passed170/nTotal;
  cout<<"Central estimate of acceptance: "<<central_estimate_of_acceptance<<endl;
  cout<<endl;
  cout<<"Excel format"<<endl;
  cout<<"Scale_variation: "<<"muR1_muF1"<<" "<<"muR1_muF2"<<" "<<"muR1_muF0p5"<<" "<<"muR2_muF1"<<" "<<"muR2_muF2"<<" "<<"muR2_muF0p5"<<" "<<"muR0p5_muF1"<<" "<<"muR0p5_muF2"<<" "<<"muR0p5_muF0p5"<<endl;
  cout<<"Total_inspected(unweighted): "<<totalInspectedSum_QCDscale_muR1_muF1<<" "<<totalInspectedSum_QCDscale_muR1_muF1<<" "<<totalInspectedSum_QCDscale_muR1_muF1<<" "<<totalInspectedSum_QCDscale_muR1_muF1<<" "<<totalInspectedSum_QCDscale_muR1_muF1<<" "<<totalInspectedSum_QCDscale_muR1_muF1<<" "<<totalInspectedSum_QCDscale_muR1_muF1<<" "<<totalInspectedSum_QCDscale_muR1_muF1<<" "<<totalInspectedSum_QCDscale_muR1_muF1<<endl;
  cout<<"Passing: "<<sum170_QCDscale_muR1_muF1<<" "<<sum170_QCDscale_muR1_muF2<<" "<<sum170_QCDscale_muR1_muF0p5<<" "<<sum170_QCDscale_muR2_muF1<<" "<<sum170_QCDscale_muR2_muF2<<" "<<sum170_QCDscale_muR2_muF0p5<<" "<<sum170_QCDscale_muR0p5_muF1<<" "<<sum170_QCDscale_muR0p5_muF2<<" "<<sum170_QCDscale_muR0p5_muF0p5<<endl;
  std::vector<float> scale_passing_physical = {sum170_QCDscale_muR1_muF1, sum170_QCDscale_muR1_muF2, sum170_QCDscale_muR1_muF0p5, sum170_QCDscale_muR2_muF1, sum170_QCDscale_muR2_muF2, sum170_QCDscale_muR0p5_muF1, sum170_QCDscale_muR0p5_muF0p5};
  std::sort(scale_passing_physical.begin(),scale_passing_physical.end());
  float central_scale_acc = sum170_QCDscale_muR1_muF1/totalInspectedSum_QCDscale_muR1_muF1;
  float max_acc = (scale_passing_physical.back())/totalInspectedSum_QCDscale_muR1_muF1;
  float min_acc = (scale_passing_physical[0])/totalInspectedSum_QCDscale_muR1_muF1;
  float scaleup_fractional = (max_acc-central_scale_acc)/central_scale_acc;
  float scaledown_fractional = (min_acc-central_scale_acc)/central_scale_acc;
  cout<<endl;
  cout<<"PDF variations"<<endl;
  cout<<"central_sum_passing: "<<central_sum_passing<<endl;
  cout<<endl;
  double sum_of_sums = 0.0;
  for(int i = 1; i < sum_passing.size(); i++)
  {
    cout<<"sum_passing["<<i<<"]: "<<sum_passing[i]<<endl;
    sum_of_sums += sum_passing[i];
  }
  double mean_Npassing = sum_of_sums/(nMCreplicas-1);
  cout<<"mean_Npassing: "<<mean_Npassing<<endl;
  double mean_acceptance = mean_Npassing/nTotal;
  cout<<"Mean acceptance: "<<mean_acceptance<<endl;
  double sum_of_squared_residuals = 0.0;
  vector<double> sum_passing_variations;
  sum_passing_variations.clear();
  for(int i = 1; i < sum_passing.size(); i++)
  {
    sum_of_squared_residuals += pow((sum_passing[i] - mean_Npassing),2.0);
    sum_passing_variations.push_back(sum_passing[i]);
  }
  double rms_error = sqrt(sum_of_squared_residuals/(nMCreplicas-2))/nTotal;
  cout<<"RMS uncertainty in mean acceptance: "<<rms_error<<endl;
  double rms_error_fractional = rms_error/mean_acceptance;
  cout<<"RMS uncertainty in mean acceptance (fractional): "<<rms_error_fractional<<endl;
  cout<<endl;
  std::sort(sum_passing_variations.begin(),sum_passing_variations.end());
  double CL_midpoint_acceptance = (sum_passing_variations[83]+sum_passing_variations[15])/(2*nTotal);
  cout<<"68% CL midpoint acceptance: "<<CL_midpoint_acceptance<<endl;
  double CL_acceptance_uncertainty = (sum_passing_variations[83]-sum_passing_variations[15])/(2*nTotal);
  cout<<"68% CL acceptance uncertainty: "<<CL_acceptance_uncertainty<<endl;
  double CL_acceptance_uncertainty_fractional = CL_acceptance_uncertainty/CL_midpoint_acceptance;
  cout<<"68% CL acceptance uncertainty (fractional): "<<CL_acceptance_uncertainty_fractional<<endl;
  cout<<endl;
  //Final summary
  cout<<"Acc*eff"<<" | Scale up variation"<<" | Scale down variation"<<" | PDF variation"<<endl;
  cout<<central_estimate_of_acceptance<<" "<<scaleup_fractional<<" "<<scaledown_fractional<<" "<<rms_error_fractional<<endl;
}

void postAnalyzer::BookHistos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("ADD","ADD");
  tree->Branch("event_","std::vector<unsigned int>",&event_);
  tree->Branch("event_info","std::vector<double>",&event_info);
  fileName->cd();
  
  Float_t PtBins[6]={175.,190.,250., 400., 700.0,1000.0};

  //Set up the histos to be filled with method fillHistos
  for(int i=0; i<125; i++)
  {
    char ptbins[100];
    sprintf(ptbins, "_%d", i);
    std::string histname(ptbins);
    h_photon_Et_range[i] = new TH1F(("Photon_Et_range"+histname).c_str(), "Photon_Et",5,PtBins);h_photon_Et_range[i]->Sumw2();
  }
}

//Fill the sequential histos at a particular spot in the sequence
void postAnalyzer::fillHistos(int histoNumber, double event_weight,int index,std::vector<int> jets,std::vector<int> leplist,Double_t dilepton_mass,Double_t dilepton_pt,Double_t leptoMET,Double_t leptoMET_phi,bool leptonsAreElectrons)
{
  Float_t uncorrectedPhoEt = ((*phoSCRawE)[index]/TMath::CosH((*phoSCEta)[index]));
  h_photon_Et_range[histoNumber]->Fill(uncorrectedPhoEt,event_weight);
}

void postAnalyzer::scaleHistos(int histoNumber, double scale_factor)
{}

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
double postAnalyzer::DeltaPhi(double phi1, double phi2)
{
  double pi = TMath::Pi();
  double dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}


//---------------------------------------------------
// get a photon candiate based on pt eta and isolation
//----------------------------------------------------

std::vector<int> postAnalyzer::getPhoCand(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons
  for(int p=0;p<nPho;p++)
  {
    Float_t uncorrectedPhoEt = ((*phoSCRawE)[p]/TMath::CosH((*phoSCEta)[p]));
    bool kinematic = uncorrectedPhoEt > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
    //Short_t IsoPass;
    //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
    bool photonId = (
         ((*phoHoverE)[p]                <  0.05   ) &&
         ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0102 ) &&
         ((*phohasPixelSeed)[p]              ==  0      ) &&
         (TMath::Max(((*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p])), 0.0) < 1.37 )  &&
         (TMath::Max(((*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p])), 0.0) < (1.06 + (0.014 * uncorrectedPhoEt) + (0.000019 * pow(uncorrectedPhoEt, 2.0))) )  &&
         (TMath::Max(((*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])), 0.0) < (0.28 + (0.0053 * uncorrectedPhoEt)) ) );
    
    //      bool noncoll = fabs((*phoseedTimeFull5x5)[p]) < 3. && (*phomipTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001;
    //      bool noncoll = fabs((*phoseedTimeFull5x5)[p]) < 3. && (*phomipTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001;

    if(photonId && kinematic )
    {
      tmpCand.push_back(p);
    }
  }

  return tmpCand;

}






std::vector<int> postAnalyzer::getPhoCand1(double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over photons
  for(int p=0;p<nPho;p++)
    {
      Float_t uncorrectedPhoEt = ((*phoSCRawE)[p]/TMath::CosH((*phoSCEta)[p]));
      bool kinematic = fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      bool photonId = (
           ((*phoHoverE)[p]                <  0.05   ) &&
           ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0101 ) &&
           //((*phohasPixelSeed)[p]              ==  0      ) &&
           ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.21 )  &&
           ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (0.65 + (0.014 * uncorrectedPhoEt) + (0.000019 * pow(uncorrectedPhoEt, 2.0))) )  &&
           ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (0.18 + (0.0053 * uncorrectedPhoEt)) ) );
      
      if(photonId && kinematic){
  tmpCand.push_back(p);
      }
    }

  return tmpCand;

}


std::vector<int> postAnalyzer::getQcdden(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons
  for(int p=0;p<nPho;p++)
    {
      bool upperBound=false;
      bool lowerBound =false;

      double  maxPFCharged= TMath::Min(5.0*(3.32) , 0.20*(*phoEt)[p]);
      double  maxPFPhoton = TMath::Min(5.0*(0.81 + (0.0053 * (*phoEt)[p])) , 0.20*(*phoEt)[p]);
      double  maxPFNeutral= TMath::Min(5.0*(1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))) , 0.20*(*phoEt)[p]);
      
      bool kinematic = (*phoEt)[p] > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      upperBound = (
           ((*phoHoverE)[p]                <  0.05   ) &&
           ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0106 ) &&
           ((*phohasPixelSeed)[p]              ==  0      ) &&
           ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < maxPFCharged )  &&
           ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < maxPFNeutral )  &&
           ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < maxPFPhoton ));

      lowerBound = (
        ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) > 3.32 )  ||
        ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) > (1.92 + (0.014* (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))  ||
        ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) > (0.81 + (0.0053 * (*phoEt)[p])) ));

      
      if(upperBound && lowerBound && kinematic){
  tmpCand.push_back(p);
      }
    }

  return tmpCand;

}

// Effective area to be needed in PF Iso for photon ID
// https://indico.cern.ch/event/455258/contribution/0/attachments/1173322/1695132/SP15_253rd.pdf -- slide-5
Double_t postAnalyzer::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0;

  return EffectiveArea;
}

Double_t postAnalyzer::EAneutral(Double_t eta){
  Float_t EffectiveArea = 0.;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0599;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0819;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0696;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0360;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0360;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0462;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0656;

  return EffectiveArea;
}

Double_t postAnalyzer::EAphoton(Double_t eta){
  Float_t EffectiveArea = 0.;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.1271;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1101;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0756;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.1175;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.1498;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.1857;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.2183;

  return EffectiveArea;
}

//Returns true if veto passed
//Veto failed if an electron is found that passes Loose Electron ID and elePtcut, and does not overlap the candidate photon within dR of 0.5
//Always true if no electrons with |SC eta| < 2.5, since ID always fails for |SC eta| > 2.

std::vector<int> postAnalyzer::electron_veto_looseID(int pho_index, float elePtCut)
{
  std::vector<int> ele_cands;
  ele_cands.clear();

  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false;
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue
  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
  {
    //Make sure these get reset for every electron
    pass_SigmaIEtaIEtaFull5x5 = false;
    pass_dEtaIn = false;
    pass_dPhiIn = false;
    pass_HoverE = false;
    pass_iso = false;
    pass_ooEmooP = false;
    pass_d0 = false;
    pass_dz = false;
    pass_missingHits = false;
    pass_convVeto = false;
    //Find EA for corrected relative iso.
    if(abs(eleSCEta->at(i)) <= 1.0)
      EA = 0.1752;
    else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
      EA = 0.1862;
    else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
      EA = 0.1411;
    else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
      EA = 0.1534;
    else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
      EA = 0.1903;
    else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
      EA = 0.2243;
    else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
      EA = 0.2687;
    EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));

    if(abs(eleSCEta->at(i)) <= 1.479)
    {
      pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
      pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
      pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
      pass_HoverE = eleHoverE->at(i) < 0.104;
      pass_iso = EAcorrIso < 0.0893;
      pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
      pass_d0 = abs(eleD0->at(i)) < 0.0261;
      pass_dz = abs(eleDz->at(i)) < 0.41;
      pass_missingHits = eleMissHits->at(i) <= 2;
      pass_convVeto = eleConvVeto->at(i) == 1;
    }
    else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
    {
      pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
      pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
      pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
      pass_HoverE = eleHoverE->at(i) < 0.0897;
      pass_iso = EAcorrIso < 0.121;
      pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
      pass_d0 = abs(eleD0->at(i)) < 0.118;
      pass_dz = abs(eleDz->at(i)) < 0.822;
      pass_missingHits = eleMissHits->at(i) <= 1;
      pass_convVeto = eleConvVeto->at(i) == 1;
    }
      //Electron passes Loose Electron ID cuts
    if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
    {
      //Electron passes pt cut
      if(elePt->at(i) > elePtCut)
      {
        //Electron does not overlap photon
        if(dR(eleSCEta->at(i),eleSCPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
        {
          ele_cands.push_back(i);
        }
      }
    }
  }
  return ele_cands;
}



//Veto failed if a muon is found that passes Loose Muon ID, Loose Muon Isolation, and muPtcut, and does not overlap the candidate photon within dR of 0.5
std::vector<int> postAnalyzer::muon_veto_looseID(int pho_index, float muPtCut)
{
  std::vector<int> mu_cands;
  mu_cands.clear();

  bool veto_passed = true; //pass veto if no good muon found
  bool pass_PFMuon = false;
  bool pass_globalMuon = false;
  bool pass_trackerMuon = false;
  bool pass_iso = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue
  Float_t zero = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
  {
    pass_PFMuon = muIsPFMuon->at(i);
    pass_globalMuon = muIsGlobalMuon->at(i);
    pass_trackerMuon = muIsTrackerMuon->at(i);
    muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
    tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
    pass_iso = tightIso_combinedRelative < 0.25;
    //Muon passes Loose Muon ID and PF-based combined relative, dBeta-corrected Loose Muon Isolation cuts
  //      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
    if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon))
    {
      //Muon passes pt cut
      if(muPt->at(i) > muPtCut)
      {
        //Muon does not overlap photon
        if(dR(muEta->at(i),muPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
        {
          mu_cands.push_back(i);
        }
      }
    }
  }
  return mu_cands;
}


std::vector<int> postAnalyzer::JetVetoDecision(int pho_index) {

  bool jetVeto=true;
  std::vector<int> jetindex;
  float value =0.0;

  for(int i = 0; i < nJet; i++)
    {

      if(0.0 < abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.5)
        value =-0.8;
      if(2.5 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.75)
        value =-0.95;
      if(2.75 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <3.0)
        value =-0.97;
      if(3.00 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <5.0)
        value =-0.99;



      //      std::cout<<"Jet size: "<<nJet<<std::endl;
      double deltar = 0.0 ;
      //      std::cout<<"Jet no:"<<i<<"coming here pujetid: "<<pfJet_pt[i]<<std::endl;
      //if(OverlapWithMuon(jetEta->at(i),jetPhi->at(i)))     continue;
      //      std::cout<<"Jet no:"<<i<<"coming here OverlapWithMuon: "<<pfJet_pt[i]<<std::endl;
      //      if(OverlapWithElectron(jetEta->at(i),jetPhi->at(i)))   continue;
      if(pho_index>=0){
        deltar= dR(jetEta->at(i),jetPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index));
        //std::cout<<"Delta R between photon and jet="<<dR<< "jet pt"<<pfJet_pt[i]<<std::endl;
      }
      if(deltar>0.4 && jetPt->at(i) >30.0 && jetPFLooseId->at(i)==1 && jetPUidFullDiscriminant->at(i)>value)
        {
          jetindex.push_back(i);
        }

      
    }


  //  std::cout<<"Jet size: "<< jetindex.size()<<std::endl;
  //if(jetindex.size()>1)jetVeto = false;
  return jetindex;

}



bool postAnalyzer::OverlapWithMuon(double eta, double phi){
  
  bool overlap = false;
  //  std::cout<<"No of muon:"<<Muon_n<<std::endl;
  bool pass_PFMuon = false;
  bool pass_globalMuon = false;
  bool pass_trackerMuon = false;
  bool pass_iso = false;

  Float_t zero1 = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
    {
      pass_PFMuon = muIsPFMuon->at(i);
      pass_globalMuon = muIsGlobalMuon->at(i);
      pass_trackerMuon = muIsTrackerMuon->at(i);
      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero1,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.25;
      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
        {
          if(muPt->at(i) > 10.)
            {
              if(dR(muEta->at(i),muPhi->at(i),eta,phi) < 0.5)
                {
                  overlap = true;
                  break;
                }
            }
        }
    }

  return overlap;


}



bool postAnalyzer::OverlapWithElectron(double eta, double phi){
  bool overlap = false;

  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false;
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;

  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
    {
      pass_SigmaIEtaIEtaFull5x5 = false;
      pass_dEtaIn = false;
      pass_dPhiIn = false;
      pass_HoverE = false;
      pass_iso = false;
      pass_ooEmooP = false;
      pass_d0 = false;
      pass_dz = false;
      pass_missingHits = false;
      pass_convVeto = false;

      if(abs(eleSCEta->at(i)) <= 1.0)
        EA = 0.1752;
      else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
        EA = 0.1862;
      else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
        EA = 0.1411;
      else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
        EA = 0.1534;
      else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
        EA = 0.1903;
      else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
        EA = 0.2243;
      else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
        EA = 0.2687;
      EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));

      if(abs(eleSCEta->at(i)) <= 1.479)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
          pass_HoverE = eleHoverE->at(i) < 0.104;
          pass_iso = EAcorrIso < 0.0893;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
          pass_d0 = abs(eleD0->at(i)) < 0.0261;
          pass_dz = abs(eleDz->at(i)) < 0.41;
          pass_missingHits = eleMissHits->at(i) <= 2;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }
      else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
          pass_HoverE = eleHoverE->at(i) < 0.0897;
          pass_iso = EAcorrIso < 0.121;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
          pass_d0 = abs(eleD0->at(i)) < 0.118;
          pass_dz = abs(eleDz->at(i)) < 0.822;
          pass_missingHits = eleMissHits->at(i) <= 1;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }
      //Electron passes Loose Electron ID cuts
      if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
        {
          //Electron passes pt cut
          if(elePt->at(i) > 10.)
            {
              //Electron does not overlap photon
              if(dR(eleSCEta->at(i),eleSCPhi->at(i),eta,phi) < 0.5)
                {
                  overlap = true;
                  break;
                }
            }
        }
    }





}


double postAnalyzer::dR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}





bool postAnalyzer::dPhiJetMET_veto(std::vector<int> jets)
{
  //pass veto if no jet is found within DeltaPhi(jet,MET) < 0.5
  bool passes = false;
  int njetsMax = jets.size();
  //Only look at first four jets
  if(njetsMax > 4)
    njetsMax = 4;
  int j = 0;
  for(; j < njetsMax; j++)
    {
      //fail veto if a jet is found with DeltaPhi(jet,MET) < 0.5
      if(DeltaPhi(jetPhi->at(jets[j]),pfMETPhi) < 0.5)
  break;
    }
  if(j==njetsMax)
    passes = true;

  return passes;
}


Double_t postAnalyzer::EAchargedworst(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.080721;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.063986;
  return EffectiveArea;
}
