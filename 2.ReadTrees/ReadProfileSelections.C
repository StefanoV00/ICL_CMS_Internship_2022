/// \author Stefano Veroni
/// File to read information from a tree file
/// Also saves some in better known csv format.
/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
root ReadProfile.C -q

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <string>
#include <stdio.h>
#include <cmath>
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
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
//-1. Functions' Definitions/////////////////////////////////////////////
//______________________________________________________________________
#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"



void ReadProfileSelections() 
{
    ////////////////////////////////////////////////////////////////////
    //0. Set Variables   ///////////////////////////////////////////////
    //__________________________________________________________________
    string path    = "../myTTrees/";
//none if none
    string rootname= "treeSelections";//file
    string dirname = "CloseXY"; //"none" if no directory in between
    string treename= "TSTree_SimFromCP";


    ///////////////////////////////////////////////////////////////////
    // 1. GET TTREES AND SET ADDRESSES   //////////////////////////////
    ///////////////////////////////////////////////////////////////////
    TFile* file = open_file(path, rootname);
    TTree* tree = get_tree (file, dirname, treename);
    int Nentries = tree->GetEntries();

    const char* coords_name   = (treename + "_Coords").c_str();
    TTree* tfriendC = get_tree(file, dirname, coords_name);
    const char* deltas_name   = (treename + "_Deltas").c_str();
    TTree* tfriendD = get_tree(file, dirname, deltas_name);
    const char* energies_name = (treename + "_Energy").c_str();
    TTree* tfriendE = get_tree(file, dirname, energies_name);
    const char* shape_name    = (treename + "_Shape").c_str();
    TTree* tfriendS = get_tree(file, dirname, shape_name);
    

    
    //POINTERS MUST BE INITIALISED
    vector <double>* E_theory= 0;
    vector <double>* E_tot   = 0;
    vector <double>* energy  = 0;
    vector <int>* layer_gstart = 0;
    vector <int>* which_g      = 0;
    vector <int>*            totmult= 0;
    vector <vector<double>>* gmult  = 0;
    tree    ->SetBranchAddress("cp_energy",          &E_theory);
    tree    ->SetBranchAddress("ts_energy",          &E_tot);
    tree    ->SetBranchAddress("lc_TSidx",           &which_g);
    tfriendE->SetBranchAddress("lc_energyTrue",      &energy);
    tfriendS->SetBranchAddress("lc_layerFromgStart", &layer_gstart);
    tfriendS->SetBranchAddress("lc_totMult",         &totmult);
    tfriendS->SetBranchAddress("lc_gMult",           &gmult);
    
    ////////////////////////////////////////////////////////////////////
    //2. Prepare Hist //////////////////////////////////////////////////
    //__________________________________________________________________
    //ROOT::RDataFrame dataframe ("TSTree_SimFromCP", "ticlTree");
    //Prepare Vector of Histograms
    TCanvas* cEg = new TCanvas("EnergyProfile", "EnergyProfile");
    TH1F* hEg = new TH1F("EnergyProfile", "EnergyProfile", 
                                    30, 0, 30);
    hEg->GetXaxis()->SetTitle("Layer #");
    hEg->GetYaxis()->SetTitle("Energy Fraction");


    TH1F* hMult = new TH1F("hm", "hm", 30, 0, 30);

    TH1F* hgMult = new TH1F("hg", "hg", 30, 0, 30);

    
    
    ////////////////////////////////////////////////////////////////////
    //3. Loop: get event, find quantities, save data, fill hists ///////
    //__________________________________________________________________
    //NOTE: RDataFrame.Define() could define new columns theta, etc...
    double Eavg = 0;
    for (int i = 0; i<Nentries; i++)
    { 
        // Getting new entry assignes value of branch to branch address
        // overwriting at each iteration 
        tree     -> GetEntry(i);
        tfriendC -> GetEntry(i);
        tfriendE -> GetEntry(i);
        tfriendS -> GetEntry(i);
        int N = energy->size();

        for (int j = 0; j<N; j++)
        {
            //Find other quantities
            int gij = which_g->at(j);
            int layer_start_ij = layer_gstart->at(j);
            double E_ij_frac = energy->at(j) / E_tot->at(gij);
            
            hEg -> Fill(layer_start_ij, E_ij_frac);
        }  
        Eavg += E_theory->at(0) + E_theory->at(1); 

        for (int k = 0; k<totmult->size(); k++)
        {
            hMult->Fill(totmult->at(k));
            if (k<gmult->at(0).size()) hgMult->Fill(gmult->at(0)[k]);
            if (k<gmult->at(1).size()) hgMult->Fill(gmult->at(1)[k]);
        }
       
        /*vector <double> totvec = (double) *totmult;
        int* totarr = &totvec[0];
        vector <double> gvec;
        gvec.insert( gvec.end(), gmult->at(0).begin(), gmult->at(0).end() );
        gvec.insert( gvec.end(), gmult->at(1).begin(), gmult->at(1).end() );
        double* garr = &gvec[0];
        hMult->FillN(totvec.size(), totarr);*/
    }
    //hMult  -> Scale(1./hMult ->GetEntries());
    //hgMult -> Scale(1./hgMult->GetEntries());
    TH1* hcum =  hMult ->GetCumulative();
    TH1* hgcum = hgMult->GetCumulative();
    //TH1* hcum =  hMult;
    //TH1* hgcum = hgMult;   

    hEg->Scale((double)1./(2.*Nentries)); //normalise
    hEg->SetMinimum(0);
    Eavg /= (double) 2.*Nentries;

    cout<<"Checkpoint: Finished Loop"<<endl;

    
    ////////////////////////////////////////////////////////////////////
    //4. The end, ladies & gentlemen////////////////////////////////////
    //__________________________________________________________________
    cEg->cd();
    string E = "Energy: "+ to_string(lround(Eavg)) +" GeV";
    const char* E_text = E.c_str();
    TText* text = new TText(0, 1.1 * hEg->GetMaximum(), E_text);
    text->SetTextAlign(11);
    text->SetTextSize(0.03);
    hEg-> Draw("B");
    text->Draw();
    cEg-> Update();

    TCanvas* ccum = new TCanvas ();
    ccum -> cd();
    hcum->GetXaxis()->SetTitle("LCs per Layer");
    hcum->GetYaxis()->SetTitle("# Layers with #LCs <= x");
    hcum -> Draw();
    ccum -> Update();
    
    TCanvas* cgcum = new TCanvas ();
    cgcum -> cd();
    hgcum->GetXaxis()->SetTitle("gLCs per Layer");
    hgcum->GetYaxis()->SetTitle("# Layers with #gLCs <= x");
    hgcum -> Draw();
    cgcum -> Update();

    ccum -> WaitPrimitive();
    cEg -> WaitPrimitive();
    
    return;
}