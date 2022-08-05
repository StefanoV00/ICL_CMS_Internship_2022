/// \author Stefano Veroni
/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
root AnalysisCLUE.C -q

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cmath>
#include <vector>


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
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
//#include "TLorentzVector.h"
//#include "TLorentzRotation.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TPad.h"
#include "TROOT.h"
//#include "TRotation.h"
#include "TStyle.h"
#include "TSystem.h"
//#include "TText.h"
#include "TTree.h"

#include "../0.Helpers/MyCfunctions.h"
#include "../0.Helpers/MyRoot.h"
using namespace std;

/*
 COMPARE TRUTH AND RECONSTRUCTION
 
 Using methodology of a pretty comprehensive paper, available at
 https://twiki.cern.ch/twiki/pub/CALICE/CaliceAnalysisNotes/CAN-057.pdf:
 “The mixed event is considered to be successfully reconstructed if
 • it contains exactly two reconstructed EM showers; 
 • energies and X, Z barycenter coordinates agree within ±20% and ±5 mm,
   respectively, of the energies and the coordinates reconstructed 
   in the single shower events”.

 Here we work with energy only.
*/
void AnalysisCLUEhist ()
{
    ////////////////////////////////////////////////////////////////////
    //0. Set Stuff  ////////////////////////////////////////////////////
    //__________________________________________________________________
    string path     = "../myTTrees/";
//none if none
    string rootname = "step3ticl_pt50_eta21_run0_FlatTracksters";//root file
    
    string dirname1 = "ticlTree"; //"none" if no directory in between
    string treename1= "TSTree_SimFromCP";
    string rootname1f= "treeFriends"; //root file with associated
    
    string dirname2 = "ticlTree"; //"none" if no directory in between
    string treename2= "TSTree_CLUE3D3";

    int Nevents = 0;  //# events to consider, 0 if all
    float ievent = 177;//INFINITY;  
    
    ////////////////////////////////////////////////////////////////////
    //1. Get TTrees (original, associated, clue)////////////////////////
    //__________________________________________________________________
    
    // True info orginal tree
    TFile* file = open_file(path, rootname);
    TTree* ttrue = get_tree (file, dirname1, treename1);
    int Nentries = ttrue->GetEntries();
    
    // True info associated tree
    TFile* friends = open_file(path, rootname1f, "read");
    const char* coords_name   = (treename1 + "_Coords").c_str();
    TTree* ttrueC = get_tree(friends, "none", coords_name);
    const char* deltas_name   = (treename1 + "_Deltas").c_str();
    TTree* ttrueD = get_tree(friends, "none", deltas_name);
    const char* energies_name = (treename1 + "_Energy").c_str();
    TTree* ttrueE = get_tree(friends, "none", energies_name);
    const char* shape_name    = (treename1 + "_Shape").c_str();
    TTree* ttrueS = get_tree(friends, "none", shape_name);

    // Reconstructed tree - CLUE
    TTree* tclue = get_tree (file, dirname2, treename2);



    ////////////////////////////////////////////////////////////////////
    //2. Get Branches //////////////////////////////////////////////////
    //__________________________________________________________________
    
    vector <double>* energy =0;
    ttrueE->SetBranchAddress("lc_energyTrue",&energy);
    
    // To distinguish different cases
    double Dsz =0;ttrueD->SetBranchAddress("lc_deltaStartSeedZ" ,&Dsz);
    int    Dsl =0;ttrueD->SetBranchAddress("ts_deltaStartLayer" ,&Dsl);
    double Dsxy=0;ttrueD->SetBranchAddress("cp_deltaVtxXY",&Dsxy);

    // True and reconstructed Energy
    vector <double>* tsEnergy_true = 0;
    vector <double>* tsEnergy_clue = 0;
    ttrue->SetBranchAddress("ts_energy",   &tsEnergy_true);
    tclue->SetBranchAddress("ts_energy",   &tsEnergy_clue);


    ////////////////////////////////////////////////////////////////////
    //3. Prepare Plot(s)////////////////////////////////////////////////
    //__________________________________________________________________
    
    // 2D Hist (x->separation, y->deltastartlayer, z->%Reconstructed)
    //int Nbins = 10;
    double XYmax = 6;
    double Dlayermax = 4; 

    TH2D* hph   =new TH2D("Reconstruction Efficiency - Photon #",
                          "Reconstruction Efficiency - Photon #",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH2D* h20   =new TH2D("Reconstruction Efficiency - 20%% sigma_E",
                          "Reconstruction Efficiency - 20%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH2D* h10   =new TH2D("Reconstruction Efficiency - 10%% sigma_E",
                          "Reconstruction Efficiency - 10%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);

    TH2D* hnorm = new TH2D("Normalisation2D", "Normalisation2D",
                                        2*XYmax  , 0, XYmax,
                                        Dlayermax, 0, Dlayermax);
    



    ////////////////////////////////////////////////////////////////////
    //4. Loop, compare, fill ///////////////////////////////////////////
    //__________________________________________________________________
    
    for (int i = 0; i<Nentries; i++) //Loop & Fill
    {
        ttrue  -> GetEntry(i);
        ttrueC -> GetEntry(i); ttrueD -> GetEntry(i);
        ttrueE -> GetEntry(i); ttrueS -> GetEntry(i);
        tclue  -> GetEntry(i);
        
        int N = energy->size();

        hnorm -> Fill(Dsxy, Dsl);

        if (tsEnergy_clue->size() == tsEnergy_true->size())
        {
            hph -> Fill(Dsxy, Dsl);

            double Etrue0 = tsEnergy_true->at(0);
            double Etrue1 = tsEnergy_true->at(1);
            double Eclue0 = tsEnergy_clue->at(0);
            double Eclue1 = tsEnergy_clue->at(1);

            bool yes20_0 = abs( Etrue0 - Eclue0 ) < 0.2*Etrue0;
            bool yes20_1 = abs( Etrue1 - Eclue1 ) < 0.2*Etrue1;
            bool yes10_0 = abs( Etrue0 - Eclue0 ) < 0.1*Etrue0;
            bool yes10_1 = abs( Etrue1 - Eclue1 ) < 0.1*Etrue1;

            if      (yes10_0 && yes10_1)
                {h10 -> Fill(Dsxy, Dsl);
                 h20 -> Fill(Dsxy, Dsl);}
            
            else if (yes20_0 && yes20_1)
                {h20 -> Fill(Dsxy, Dsl);}
            
            else continue;
        }
    }
    cout<<"Loop done"<<endl;

    // Normalise
    TH2D hnorm_bit = (*hnorm)*0.000001;
    TH2D hphnorm = (*hph) + hnorm_bit;
    hphnorm = hphnorm / (*hnorm);
    hphnorm.SetContour(99); //nicer, default is like 10
    hphnorm.GetXaxis()->SetTitle("Conersion Vertex XY Distances");
    hphnorm.GetYaxis()->SetTitle("Delta Start Layer");
    hphnorm.GetZaxis()->SetTitle("Reconstruction Efficiency");

    TH2D h20norm = (*h20) + hnorm_bit;
    h20norm = h20norm / (*hnorm);
    h20norm.SetContour(99); //nicer, default is like 10
    h20norm.GetXaxis()->SetTitle("Conersion Vertex XY Distances");
    h20norm.GetYaxis()->SetTitle("Delta Start Layer");
    h20norm.GetZaxis()->SetTitle("Reconstruction Efficiency");
   
    TH2D h10norm = (*h10) + hnorm_bit;
    h10norm = h10norm / (*hnorm);
    h10norm.SetContour(99); //nicer, default is like 10
    h10norm.GetXaxis()->SetTitle("Conersion Vertex XY Distances");
    h10norm.GetYaxis()->SetTitle("Delta Start Layer");
    h10norm.GetZaxis()->SetTitle("Reconstruction Efficiency");



    ////////////////////////////////////////////////////////////////////
    //5. Draw //////////////////////////////////////////////////////////
    //__________________________________________________________________
    
    gStyle->SetPalette(kRainBow); //nice colour option for colz

    TCanvas* cph = new TCanvas("AllEvents-#Ph", "AllEvents-#Ph"); 
    cph -> Divide(2, 1);
    cph -> cd(1);
    hphnorm.SetStats(0);
    hphnorm.Draw("surf1");//("colz");
    cph -> cd(2);
    hphnorm.SetStats(0);
    hphnorm.Draw("colz");//("colz");
    cph-> Update();

    TCanvas* c20 = new TCanvas("AllEvents-20%", "AllEvents-20%"); 
    c20 -> Divide(2, 1);
    c20 -> cd(1);
    h20norm.SetStats(0);
    h20norm.Draw("surf1");//("colz");
    c20 -> cd(2);
    h20norm.SetStats(0);
    h20norm.Draw("colz");//("colz");
    c20-> Update();


    TCanvas* c10 = new TCanvas("AllEvents-10%", "AllEvents-10%");
    c10 -> Divide(2, 1);
    c10 -> cd(1);
    h10norm.SetStats(0);
    h10norm.Draw("surf1");
    c10 -> cd(2);
    h10norm.SetStats(0);
    h10norm.Draw("colz");
    c10-> Update();

    c20 -> WaitPrimitive();
    c10 -> WaitPrimitive();

    return;
}