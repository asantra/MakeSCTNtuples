#include <fstream>
#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <vector>

using namespace std;
void MakeTree(TString inFileName, TString outFileName){
    cout << "It is working" << endl;
    ifstream infile;   

    
    infile.open(inFileName);// file containing numbers in 3 columns 
    if(infile.fail()){ 
        cout << "error" << endl; 
        return; // no point continuing if the file didn't open...
    }
    
    Long64_t event_number;
    vector<Double_t> pT;
    vector<Double_t> trkEta;
    vector<Double_t> trkPhi;
    vector<Double_t> phiToWafer;
    vector<Int_t>    nStrip;
    vector<Int_t>    bec;
    vector<Int_t>    layer;
    vector<Int_t>    etaModule;
    vector<Int_t>    phiModule;

    TFile *f = new TFile(outFileName,"RECREATE");
    f->cd();
    TTree *tree = new TTree("tree","An example of ROOT tree with a few branches");
    tree->Branch("event_number", & event_number,     "event_number/L");
    tree->Branch("pT",           "vector<Double_t>",  & pT);
    tree->Branch("trkEta",       "vector<Double_t>",  & trkEta);
    tree->Branch("trkPhi",       "vector<Double_t>",  & trkPhi);
    tree->Branch("phiToWafer",   "vector<Double_t>",  & phiToWafer);
    tree->Branch("nStrip",       "vector<Int_t>",     & nStrip);
    tree->Branch("bec",          "vector<Int_t>",     & bec);
    tree->Branch("layer",        "vector<Int_t>",     & layer);
    tree->Branch("etaModule",    "vector<Int_t>",     & etaModule);
    tree->Branch("phiModule",    "vector<Int_t>",     & phiModule);
    
    pT.clear(); trkEta.clear(); trkPhi.clear(); phiToWafer.clear(); nStrip.clear();
    bec.clear(); layer.clear(); etaModule.clear(); phiModule.clear();
    vector<Long64_t> vecEventNumber; vecEventNumber.clear();
    
    int count = 0;
    
    while (!infile.eof()){
       
       string timeStamp, namePattern;
       Double_t ptHere, trkEtaHere, trkPhiHere, phiToWaferHere;
       Int_t nStripHere, becHere, layerHere, etaModuleHere, phiModuleHere;
       Long64_t event_numberHere;
       
       infile >> timeStamp >> namePattern >> event_numberHere >> pTHere >> trkEtaHere >> trkPhiHere >> phiToWaferHere >> nStripHere >> becHere >> layerHere >> etaModuleHere >> phiModuleHere;
       
       vecEventNumber.push_back(event_numberHere);
       
       size_t vecSize = vecEventNumber.size();
       
       if(vecSize = 1){
            pT.push_back(pTHere);
            trkEta.push_back(trkEtaHere);
            trkPhi.push_back(trkPhiHere);
            phiToWafer.push_back(phiToWaferHere);
            nStrip.push_back(nStripHere);
            bec.push_back(becHere);
            layer.push_back(layerHere);
            etaModule.push_back(etaModuleHere);
            phiModule.push_back(phiModuleHere); 
       }
       else if((vecSize > 1) && (vecEventNumber.at(vecSize) == vecEventNumber.at(vecSize-1))){
            pT.push_back(pTHere);
            trkEta.push_back(trkEtaHere);
            trkPhi.push_back(trkPhiHere);
            phiToWafer.push_back(phiToWaferHere);
            nStrip.push_back(nStripHere);
            bec.push_back(becHere);
            layer.push_back(layerHere);
            etaModule.push_back(etaModuleHere);
            phiModule.push_back(phiModuleHere);   
       }
       else if(vecEventNumber.at(vecSize) != vecEventNumber.at(vecSize-1)){
           
       }
       
      
       std::cout << timeStamp << " " << namePattern << " " << event_number << " " << pT << " " << trkEta << " " << trkPhi << " " << phiToWafer << " " << nStrip << " " << bec << " " << layer << " " << etaModule << " " << phiModule << std::endl;
       count++;
       
       if(count >= 200)break;
       tree->Fill();
    }
    f->cd();
    f->Write();
    delete tree;
    delete f;
  
} 