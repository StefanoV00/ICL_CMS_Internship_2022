/// \author Stefano Veroni
/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
root Read1Event.C -q
 
*/
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric> 
#include <string>
#include <stdio.h>
#include <vector>

// WHEN RUNNING A ROOT SESSION FROM CONSOLE, ALL WILL BE
// AVAILABLE EVEN IF NOT INCLUDED HERE.
// THE INCLUSION IS JUST TO HAVE AUTHOMATIC COMPLETION
// BY VISUAL STUDIO CODE
// BEFORE RUNNING, THOSE .hxx NOT IN ROOT/INCLUDE MUST BE
// COMMENTED OUT
//#include "TApplication.h"
//#include <TArrow.h>
//#include <TBenchmark.h>
#include "TBox.h"
#include "TButton.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TColor.h"
#include "TDatabasePDG.h"
#include "TF1.h"
#include "TFile.h"
#include <TFrame.h>
#include "TGraphErrors.h"
#include "TH1.h"
#include <TH2.h>
//#include <THStack.h>
#include <TInterpreter.h>
#include "TKey.h"
//#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
//#include "TLorentzVector.h"
//#include "TLorentzRotation.h"
#include "TMath.h"
#include <TNtuple.h>
//#include "TPad.h"
//#include "TParticlePDG.h"
//#include <TPave.h>
//#include <TPaveText.h>
//#include "TPostScript.h"
//#include <TProfile.h>
//#include "TRandom.h"
//#include <TRandom3.h>
#include "TROOT.h"
//#include "TRotation.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"
//#include "TVector3.h"
//#include "TVirtualPad.h"
//#include "TVirtualPS.h"
//#include "Riostream.h"
//#include <RDataFrame.hxx>

using namespace std;

/////////////////////////////////////////////////////////////////////////
//-1. Functions' Definitions/////////////////////////////////////////////
//______________________________________________________________________
#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"




void Read1Event(bool save_csv = false) 
{
    ////////////////////////////////////////////////////////////////////
    //0. Set Variables   ///////////////////////////////////////////////
    //__________________________________________________________________
    string path    = "../myTTrees/";
//none if none
    string rootname= "step3ticl_pt50_eta21_run0_FlatTracksters";//root file
    string friendrootname= "treeFriends"; //root file with friends
    string dirname = "ticlTree"; //"none" if no directory in between
    string treename= "TSTree_SimFromCP";
    int i = 177;
    
    ////////////////////////////////////////////////////////////////////
    //1. Get TTree & Branches///////////////////////////////////////////
    //__________________________________________________________________
    //NOTE: FOR EACH ENTRY, LC_QUANTITY IS A VECTOR OF "DETECTIONS"
    //FOR EXAMPLE: lc_x AND lc_layer WILL HAVE N ENTRIES, ONE PER 
    //CLUSTER
    TFile* file = open_file(path, rootname);
    TTree* tree = get_tree (file, dirname, treename);
    // GET FRIEND: Normally, would just do
    //TTree* treef = tree->GetFriend((treename+"_friend").c_str());
    //But there was an error, could not save as friend TTree, so
    TFile* filef = open_file(path, friendrootname);
    const char* coords_name = (treename + "_Coords").c_str();
    TTree* treeC = get_tree(filef, "none", coords_name);
    const char* energy_name = (treename + "_Energy").c_str();
    TTree* treeE = get_tree(filef, "none", energy_name);
    const char* shape_name = (treename + "_Shape").c_str();
    TTree* treeS = get_tree(filef, "none", shape_name);
    
    
    // Get Entries Number
    int Nentries = tree->GetEntries();
    
    //POINTERS MUST BE INITIALISED
    vector <double>* x     = 0;
    vector <double>* y     = 0;
    vector <double>* z     = 0;
    vector <double>* eta   = 0;
    vector <double>* phi   = 0;
    vector <int>* which_g = 0;
    vector <int>* layerN  = 0;
    vector <int>* Nrechits= 0;
    tree->SetBranchAddress("lc_x",        &x);
    tree->SetBranchAddress("lc_y",        &y);
    tree->SetBranchAddress("lc_z",        &z);
    tree->SetBranchAddress("lc_eta",      &eta);
    tree->SetBranchAddress("lc_phi",      &phi);
    tree->SetBranchAddress("lc_layer",    &layerN);
    tree->SetBranchAddress("lc_nrechits", &Nrechits );
    tree->SetBranchAddress("lc_TSidx",    &which_g );
    //From Coords TTree
    vector <double>* xpr = 0;
    vector <double>* ypr = 0;
    treeC->SetBranchAddress("lc_x_pr",  &xpr);
    treeC->SetBranchAddress("lc_y_pr",  &ypr);
    //From Energy of TTree
    vector <double>* energy= 0;
    vector <double>* gCentreX_pr = 0;
    vector <double>* gCentreY_pr = 0;
    double centreTotX_pr = 0;
    double centreTotY_pr = 0;
    treeE->SetBranchAddress("lc_energyTrue",    &energy);
    treeE->SetBranchAddress("lc_gCentreX_pr",   &gCentreX_pr);
    treeE->SetBranchAddress("lc_gCentreY_pr",   &gCentreY_pr);
    treeE->SetBranchAddress("lc_centreTotX_pr", &centreTotX_pr);
    treeE->SetBranchAddress("lc_centreTotY_pr", &centreTotY_pr);
    cout<<"Addresses set"<<endl;
    
    ////////////////////////////////////////////////////////////////////
    //2. Prepare Hist //////////////////////////////////////////////////
    //__________________________________________________________________
    //ROOT::RDataFrame dataframe ("TSTree_SimFromCP", "ticlTree");
    //Prepare Vector of Histograms
    tree  -> GetEntry(i);
    treeC -> GetEntry(i);
    treeE -> GetEntry(i);

    double xprmin = xpr->at(argmin(*xpr))-0.02;
    double xprmax = xpr->at(argmax(*xpr))+0.02;
    double yprmin = ypr->at(argmin(*ypr))-0.02;
    double yprmax = ypr->at(argmax(*ypr))+0.02;
    cout<<"Extremi Taken"<<endl;
    TCanvas* c_projections = new TCanvas("call", "call");
        TH2D* h_projections = new TH2D("xy-projections", "xy-projections", 
                                        100, xprmin, xprmax,
                                        100, yprmin, yprmax); 
        h_projections->SetContour(99); //nicer, default is like 10
        h_projections->GetXaxis()->SetTitle("x projection");
        h_projections->GetYaxis()->SetTitle("y projection");
        h_projections->GetZaxis()->SetTitle("# E fraction");
    string id = "Event " + to_string(i);
    const char* id_text = id.c_str();
    TText* textpr = new TText(xprmax,  yprmax, id_text);
    textpr->SetTextAlign(31);
    textpr->SetTextSize(0.025);
    cout<<"Hist prep"<<endl;

    ////////////////////////////////////////////////////////////////////
    //3. Do stuff for chosen event /////////////////////////////////////
    //__________________________________________________________________
    //NOTE: RDataFrame.Define() could define new columns theta, etc...

    int N = x->size();
    vector <double> xpr0 = getAwhereB(*xpr, *which_g, 0);
    vector <double> ypr0 = getAwhereB(*ypr, *which_g, 0);
    vector <double> xpr1 = getAwhereB(*xpr, *which_g, 1);
    vector <double> ypr1 = getAwhereB(*ypr, *which_g, 1);
    vector <double> ene0 = getAwhereB(*energy, *which_g, 0);
    double E0 = myvecsum(ene0);
    vector <double> ene1 = getAwhereB(*energy, *which_g, 1);
    double E1 = myvecsum(ene1);
    double Etot = E0 + E1;

    // all hits, colour based on the photon they belong too
    TCanvas* c_2ph = new TCanvas("c_2ph", "c_2ph");
        c_2ph->DrawFrame(xprmin, yprmin, xprmax, yprmax);
        const int N0 = xpr0.size();
        double xpr0arr [N0]; double ypr0arr [N0];
        for (int i = 0; i<N0; i++)
            {xpr0arr[i] = xpr0[i]; ypr0arr[i] = ypr0[i];}
        auto gr0 = new TGraph(N0, xpr0arr, ypr0arr); 
        gr0->SetMarkerStyle(47);
        gr0->SetMarkerSize(0.5);
        gr0->SetMarkerColor(4); 
        gr0->SetTitle("Photons' Showers;x-projection;y-projection");

        const int N1 = xpr1.size();
        double xpr1arr [N1]; double ypr1arr [N1];
        for (int i = 0; i<N1; i++)
            {xpr1arr[i] = xpr1[i]; ypr1arr[i] = ypr1[i];}
        auto gr1 = new TGraph(N1, xpr1arr, ypr1arr);
        gr1->SetMarkerStyle(33);
        gr1->SetMarkerSize(0.6);
        gr1->SetMarkerColor(2);
        gr1->SetTitle("Photons' Showers;x-projection;y-projection");
    
    // all hits, colour based etc with efr > 0.005
    TCanvas* c_2phbig = new TCanvas("c_2phbig", "c_2phbig");
        c_2phbig->DrawFrame(xprmin, yprmin, xprmax, yprmax);
        vector <double> xpr0big = getAwhereB(xpr0, ene0, 
                            [E0](double ene){return ene/E0>0.005;});
        vector <double> ypr0big = getAwhereB(ypr0, ene0, 
                            [E0](double ene){return ene/E0>0.005;});
        const int N0big = xpr0big.size();
        double xpr0arrbig [N0big];
        double ypr0arrbig [N0big];
        for (int i = 0; i<N0big; i++)
            {xpr0arrbig[i] = xpr0big[i]; 
             ypr0arrbig[i] = ypr0big[i];}
        auto gr0big = new TGraph(N0big, xpr0arrbig, ypr0arrbig); 
        gr0big->SetMarkerStyle(47);
        gr0big->SetMarkerSize(0.5);
        gr0big->SetMarkerColor(4); 
        gr0big->SetTitle("Photons' Showers;x-projection;y-projection");

        vector <double> xpr1big = getAwhereB(xpr1, ene1, 
                            [E1](double ene){return ene/E1>0.005;});
        vector <double> ypr1big = getAwhereB(ypr1, ene1, 
                            [E1](double ene){return ene/E1>0.005;});
        const int N1big = xpr1big.size();
        double xpr1arrbig [N1big]; 
        double ypr1arrbig [N1big];
        for (int i = 0; i<N1big; i++)
            {xpr1arrbig[i] = xpr1big[i];
             ypr1arrbig[i] = ypr1big[i];}
        auto gr1big = new TGraph(N1big, xpr1arrbig, ypr1arrbig);
        gr1big->SetMarkerStyle(33);
        gr1big->SetMarkerSize(0.6);
        gr1big->SetMarkerColor(2);
        gr1big->SetTitle("Photons' Showers;x-projection;y-projection");
        const char* textup = "E-fraction > 0.005";
        TText* text = new TText(xprmin,  yprmax, textup);
        text->SetTextAlign(11);
        text->SetTextSize(0.025);

    cout<<"Graph prep"<<endl;
    for (int j = 0; j<N; j++)
    {
        //Find other quantities
        int layer_ij = layerN->at(j);
        double x_ij  = x->at(j);
        double y_ij  = y->at(j);
        double z_ij  = z->at(j);
        double xpr_ij  = xpr -> at(j);
        double ypr_ij  = ypr -> at(j);
        double efr_ij  = energy -> at(j) / Etot;

        h_projections->Fill(xpr_ij, ypr_ij, efr_ij);
    }
    
    ////////////////////////////////////////////////////////////////////
    //4. The end, ladies & gentlemen////////////////////////////////////
    //__________________________________________________________________
    c_projections->cd();
    gStyle->SetPalette(kRainBow); //nice colour option for colz
    gPad  ->SetLogz(1);
    h_projections-> SetStats(0);
    h_projections-> Draw("colz");
    textpr ->Draw();
    c_projections-> Update();

    c_2ph->cd();
    c_2ph->DrawFrame(xprmin, yprmin, xprmax, yprmax);
    gr0-> SetStats(0);
    gr0-> Draw("P"); //P = point only
    gr1-> SetStats(0);
    gr1-> Draw("P");
    textpr ->Draw();
    c_2ph-> Update();

    c_2phbig->cd();
    c_2phbig->DrawFrame(xprmin, yprmin, xprmax, yprmax);
    gr0big-> SetStats(0);
    gr0big-> Draw("P"); //P = point only
    gr1big-> SetStats(0);
    gr1big-> Draw("P");
    text  -> Draw();
    textpr ->Draw();
    c_2phbig-> Update();

    c_projections -> WaitPrimitive();
    
    return;
    file->Close();   
    filef->Close();
}