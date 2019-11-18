#include <fstream>
#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
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
    Double_t pT;
    Double_t trkEta;
    Double_t trkPhi;
    Double_t phiToWafer;
    Int_t    nStrip;
    Int_t    bec;
    Int_t    layer;
    Int_t    etaModule;
    Int_t    phiModule;
    Int_t    side;
    Double_t charge;

    TFile *f = new TFile(outFileName,"RECREATE");
    f->cd();
    TTree *tree = new TTree("tree","An example of ROOT tree with a few branches");
    tree->Branch("event_number", & event_number,   "event_number/L");
    tree->Branch("pT",           & pT,             "pT/D");
    tree->Branch("trkEta",       & trkEta,         "trkEta/D");
    tree->Branch("trkPhi",       & trkPhi,         "trkPhi/D");
    tree->Branch("phiToWafer",   & phiToWafer,     "phiToWafer/D");
    tree->Branch("nStrip",       & nStrip,         "nStrip/I");
    tree->Branch("bec",          & bec,            "bec/I");
    tree->Branch("layer",        & layer,          "layer/I");
    tree->Branch("etaModule",    & etaModule,      "etaModule/I");
    tree->Branch("phiModule",    & phiModule,      "phiModule/I");
    tree->Branch("side",         & side,           "side/I");
    tree->Branch("charge",       & charge,         "charge/D");
    
    
    
    int count = 0;
    while (!infile.eof()){
       
       string timeStamp, namePattern;
       
       infile >> timeStamp >> namePattern >> event_number >> pT >> trkEta >> trkPhi >> phiToWafer >> nStrip >> bec >> layer >> etaModule >> phiModule >> side >> charge;
       
       
       //std::cout << timeStamp << " " << namePattern << " " << event_number << " " << pT << " " << trkEta << " " << trkPhi << " " << phiToWafer << " " << nStrip << " " << bec << " " << layer << " " << etaModule << " " << phiModule << std::endl;
       count++;
       
       if(count >= 1000000)break;
       tree->Fill();
    }
    f->cd();
    f->Write();
    delete tree;
    delete f;
  
} 