/// \author Stefano Veroni
/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 3.Analysis
root AnalysisCLUE.C -q
cd ../

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cmath>
#include <vector>



#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include <TH2.h>
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"
using namespace std;

/*
 COMPARE TRUTH AND RECONSTRUCTION
 
 Get Accuracy values of CLUE3D reconstruction and print those as:
 - Event Accuracy: all photons in event are reconstructed 
   up to x accuracy
 - Photon Accuracy: count number of photons correctly reconstructed
   up to x accuracy. 
 - One-ph Accuracy: accuracy of single-photons events.
 - Two-ph Accuracy: (Event) accuracy of double-photons events.

 Plot hists of accuracy as function delta cp_vtx in xy plane and
 delta startlayer.
*/
void AnalysisCLUE ()
{
    ////////////////////////////////////////////////////////////////////
    //0. Set Stuff  ////////////////////////////////////////////////////
    //__________________________________________________________________
    string path     = "../myTTrees/";
    vector <string> paths = {
                             //path+"CloseByPhotons/"     ,
                             path+"CloseByPhotonsDR1d5/",
                             path+"CloseByPhotonsDR2d5/",
                             path+"CloseByPhotonsDR3/"  ,
                             path+"CloseByPhotonsDR3d5/",
                             path+"CloseByPhotonsDR4/"  ,
                             path+"CloseByPhotonsDR5/"  ,
                             path+"SinglePhotons/",
                             path+"SinglePhotons2/"
                             };

    string rootname= "step3ticl_pt50_eta21_run0_FlatTracksters";//root file
    string rootname1 = "step3ticl_";
    string rootname2 = "_eta";
    string rootname3 = "_run";
    string rootname4 = "_FlatTracksters";
    
    vector <string> pts  = {"En40to400", "En20to200"};
    vector <string> etas = {"21"};
    int             Nrun = 50;
    int istart   = 0;
    int ptstart  = 0;
    int etastart = 0;
    vector <int> runstart = {0, 0};
    vector <int> problematic_runs = {192, 279, 289};
    
    string dirname1 = "ticlTree"; //"none" if no directory in between
    string treename1= "TSTree_SimFromCP";
    
    string dirname2 = "ticlTree"; //"none" if no directory in between
    string treename2= "TSTree_CLUE3D3"; 
    
    ////////////////////////////////////////////////////////////////////
    //1. Prepare Plot(s)////////////////////////////////////////////////
    //__________________________________________________________________
    
    // 2D Hist (x->separation, y->deltastartlayer, z->%Reconstructed)
    //int Nbins = 10;
    double XYmax = 5.5;
    double XYmin = 1.5;
    int Nxy = (int) (XYmax - XYmin +1)* 2;
    double Dlayermax = 7; 

    //__________________________________________________________________
    TH1D* hph   =new TH1D("Event Reconstruction Efficiency - Photon #",
                          "Event Reconstruction Efficiency - Photon #",
                            Nxy  , XYmin, XYmax);
                            //Dlayermax, 0, Dlayermax);
    
    TH1D* h20   =new TH1D("Event Reconstruction Efficiency - 20%% sigma_E",
                          "Event Reconstruction Efficiency - 20%% sigma_E",
                            Nxy  , XYmin, XYmax);
                           // Dlayermax, 0, Dlayermax);
    
    TH1D* h10   =new TH1D("Event Reconstruction Efficiency - 10%% sigma_E",
                          "Event Reconstruction Efficiency - 10%% sigma_E",
                            Nxy  , XYmin, XYmax);
                            //Dlayermax, 0, Dlayermax);
    
    TH1D* h05   =new TH1D("Event Reconstruction Efficiency - 5%% sigma_E",
                          "Event Reconstruction Efficiency - 5%% sigma_E",
                            Nxy  , XYmin, XYmax);

    TH1D* hf20   =new TH1D("Event Reconstruction Efficiency - 20%% sigma_E",
                          "Event Reconstruction Efficiency - 20%% sigma_E",
                            Nxy  , XYmin, XYmax);
                           // Dlayermax, 0, Dlayermax);
    
    TH1D* hf10   =new TH1D("Event Reconstruction Efficiency - 10%% sigma_E",
                          "Event Reconstruction Efficiency - 10%% sigma_E",
                            Nxy  , XYmin, XYmax);
                            //Dlayermax, 0, Dlayermax);
    
    TH1D* hf05   =new TH1D("Event Reconstruction Efficiency - 5%% sigma_E",
                          "Event Reconstruction Efficiency - 5%% sigma_E",
                            Nxy  , XYmin, XYmax);
                            //Dlayermax, 0, Dlayermax);

    TH1D* hnorm = new TH1D("Event Normalisation2D", "Event Normalisation2D",
                            Nxy  , XYmin, XYmax);
                            //Dlayermax, 0, Dlayermax);
    

    //__________________________________________________________________
    TH1D* hdelta = new TH1D("DeltaReconstruction", "DeltaReconstruction",
                                        100, -100, 100);

    TH1D* hreldelta = new TH1D("RelDeltaReconstruct", "RelDeltaReconstruct",
                                        50, -1, 1);
    
    TH1D* hratio = new TH1D("RatioReconstruct", "RatioReconstruct",
                                        50, 0, 10);
    
    
////////////////////////////////////////////////////////////////////////
//////2. LOOP OVER TTREES //////////////////////////////////////////////
/////__________________________________________________________________
    //Guide to subscripts:
    // For accuracy
    // - ph : classification accuracy
    // - <N>: accuracy to N%
    // For counting what
    // - E : all info in event matches up to <N>
    // - P : photons treated individually
    // - _1: single photons events
    // - _2: double photons events

    int countphE = 0;  int countphP  = 0; 
    int countph_1= 0;  int countph_2 = 0;
   
    int count20E = 0;   int count20P  = 0;
    int count20_1= 0;   int count20_2 = 0;
    int countf20E = 0;  int countf20P  = 0;
    int countf20_1= 0;  int countf20_2 = 0;
   
    int count10E = 0;  int count10P  = 0;
    int count10_1= 0;  int count10_2 = 0;
    int countf10E = 0;  int countf10P  = 0;
    int countf10_1= 0;  int countf10_2 = 0;
   
    int count05E = 0;  int count05P  = 0;
    int count05_1= 0;  int count05_2 = 0;
    int countf05E = 0;  int countf05P  = 0;
    int countf05_1= 0;  int countf05_2 = 0;
   
    int countallE = 0; int countallP = 0; 
    int countall_1= 0; int countall_2= 0;
    
    for (int i_p = 0; i_p < paths.size(); i_p++)
    {
    string path = paths[i_p];
    for (int pt_i = ptstart; pt_i < pts.size(); pt_i++)
    {
    string rootname_pt = rootname1 + pts[pt_i];
    for (int e_i = etastart; e_i < etas.size(); e_i++)
    {
    string rootname_eta = rootname2 + etas[e_i];
    //Balance the single photons (last one)
    Nrun *= (i_p==paths.size()-1 && e_i==0 && pt_i==0) ? 5 : 1;
    runstart[pt_i] *= (i_p==paths.size()-1 && e_i==0) ? 5 : 1;
    for (int run = runstart[pt_i]; run < Nrun; run++)
    {
        if (contains(problematic_runs, run)) continue; 
        printf("Path: %s, ", path.c_str());
        printf("PT Index: %d, RUN N. is %2d --------\n", pt_i, run);
        string rootname_run= rootname3 + to_string(run) + rootname4;
        int new_i = int(run/10);

        string rootname_p_e_r = rootname_pt  +
                                rootname_eta +
                                rootname_run ; 
                     
        if (gSystem->AccessPathName((path+rootname_p_e_r+".root")
                                    .c_str()) == 0 )
        {
        TFile* file = new TFile((path+rootname_p_e_r+".root").c_str(),
                                "read");
        TDirectory* dir = (TDirectory*)file->Get(dirname1.c_str());
        if (dir)
        {

        // True info orginal tree
        TTree* ttrue = (TTree*) dir->Get(treename1.c_str());
        int Nentries = ttrue->GetEntries();
        // Reconstructed tree - CLUE
        TTree* tclue = (TTree*) dir->Get(treename2.c_str());


        ////////////////////////////////////////////////////////////
        //3. Get Branches //////////////////////////////////////////
        //__________________________________________________________
        // True and reconstructed Info
        vector <double>* emEnergy_true = 0;
        vector <double>* emEnergy_clue = 0;
        ttrue->SetBranchAddress("ts_emEnergy",   &emEnergy_true);
        tclue->SetBranchAddress("ts_emEnergy",   &emEnergy_clue);
        vector <double>* cp_vtxX = 0;
        vector <double>* cp_vtxY = 0;
        vector <double>* ts_firstlayer = 0;
        ttrue->SetBranchAddress("cp_vtxX",      &cp_vtxX);
        ttrue->SetBranchAddress("cp_vtxY",      &cp_vtxY);
        ttrue->SetBranchAddress("ts_firstLayer",&ts_firstlayer);
       


        ////////////////////////////////////////////////////////////
        //4. Loop, compare, fill ///////////////////////////////////
        //__________________________________________________________
        auto isecal   = [](int layN)     {return layN  <=28;};
        auto distinct = [](double tsmult){return tsmult<= 2;};
        for (int i = 0; i<Nentries; i++) //Loop & Fill
        {
            ttrue  -> GetEntry(i);
            tclue  -> GetEntry(i);
            int Nphs =  emEnergy_true->size();
            
            countallE  += 1;
            countallP  += Nphs;
            countall_1 += (Nphs == 1) ? 1 : 0;
            countall_2 += (Nphs == 2) ? 1 : 0;

            double Etrue0 = emEnergy_true->at(0);
            double Eclue0 = emEnergy_clue->at(0);
            double Etrue1 =(emEnergy_true->size()>1)
                           ? emEnergy_true->at(1)
                           : 0;
            double Eclue1 =(emEnergy_clue->size()>1)
                           ? emEnergy_clue->at(1)
                           : 0;
            double Etottrue = myvecsum(*emEnergy_true);
            double Etotclue = myvecsum(*emEnergy_clue);
            if (Etrue0 < Etrue1) //E0 must be bigger
                {double a = Etrue0*1.;Etrue0 = Etrue1*1.;Etrue1 = a*1.;}
            if (Eclue0 < Eclue1) //E0 must be bigger
                {double a = Eclue0*1.;Eclue0 = Eclue1*1.;Eclue1 = a*1.;} 
            double Eftrue0 = Etrue0/Etottrue;
            double Eftrue1 = Etrue1/Etottrue;
            double Efclue0 = Eclue0/Etotclue;
            double Efclue1 = Eclue1/Etotclue;

            double err0  = (Etrue0  - Eclue0 ) / Etrue0 ;
            double ferr0 = (Eftrue0 - Efclue0) / Efclue0;
            err0  = abs(err0);
            ferr0 = abs(ferr0);
            double err1  = (Nphs==2 && emEnergy_clue->size()==2) 
                            ? (Etrue1  - Eclue1 ) / Etrue1 : 0.;
            double ferr1 = (Nphs==2 && emEnergy_clue->size()==2) 
                            ? (Etrue1  - Eclue1 ) / Etrue1 : 0.;
            double errmax  = max(err0,  err1);
            double ferrmax = max(ferr0, ferr1);

            if (Nphs == emEnergy_clue->size())
            {
                // ACCURACY CALCULATIONS _______________________________
                countphE  += 1;
                countphP  += Nphs;
                countph_1 += (Nphs == 1) ? 1 : 0;
                countph_2 += (Nphs == 2) ? 1 : 0;

                //Now treat the absolute errors
                if (errmax < 0.05)
                {
                    count20E += 1;
                    count10E += 1;
                    count05E += 1;
                    count20P   += Nphs;
                    count20_1 += (Nphs == 1) ? 1 : 0;
                    count20_2 += (Nphs == 2) ? 1 : 0;
                    count10P  += Nphs;
                    count10_1 += (Nphs == 1) ? 1 : 0;
                    count10_2 += (Nphs == 2) ? 1 : 0;
                    count05P  += Nphs;
                    count05_1 += (Nphs == 1) ? 1 : 0;
                    count05_2 += (Nphs == 2) ? 1 : 0;
                }
                else if (errmax < 0.10)
                {
                    count20E  += 1;
                    count10E  += 1;
                    count20P  += Nphs;
                    count20_1 += (Nphs == 1) ? 1 : 0;
                    count20_2 += (Nphs == 2) ? 1 : 0;
                    count10P  += Nphs;
                    count10_1 += (Nphs == 1) ? 1 : 0;
                    count10_2 += (Nphs == 2) ? 1 : 0;
                    if (err0 < 0.05 || (err1 < 0.05 && Nphs == 2) )
                        {count05P += 1; }
                }
                else if (errmax < 0.20)
                {
                    count20E += 1;
                    count20P  += Nphs;
                    count20_1 += (Nphs == 1) ? 1 : 0;
                    count20_2 += (Nphs == 2) ? 1 : 0;
                    if (err0 < 0.05 || (err1 < 0.05 && Nphs == 2) )
                        {count10P += 1;
                         count05P += 1; }
                    else if (err0 < 0.10 || (err1 < 0.10 && Nphs == 2) )
                        {count10P += 1; }
                }
                else if (err0 < 0.05 || (err1 < 0.05 && Nphs == 2))
                    {count20P += 1;
                     count10P += 1;
                     count05P += 1; }
                else if (err0 < 0.10 || (err1 < 0.10 && Nphs == 2) )
                        {count20P += 1;
                         count10P += 1; }
                else if (err0 < 0.20 || (err1 < 0.20 && Nphs == 2) )
                        {count20P += 1;}

                
                // Now treat the fractional errors
                if (ferrmax < 0.05)
                {
                    countf20E += 1;
                    countf10E += 1;
                    countf05E += 1;
                    countf20P  += Nphs;
                    countf20_1 += (Nphs == 1) ? 1 : 0;
                    countf20_2 += (Nphs == 2) ? 1 : 0;
                    countf10P  += Nphs;
                    countf10_1 += (Nphs == 1) ? 1 : 0;
                    countf10_2 += (Nphs == 2) ? 1 : 0;
                    countf05P  += Nphs;
                    countf05_1 += (Nphs == 1) ? 1 : 0;
                    countf05_2 += (Nphs == 2) ? 1 : 0;
                }
                else if (ferrmax < 0.10)
                {
                    countf20E  += 1;
                    countf10E  += 1;
                    countf20P  += Nphs;
                    countf20_1 += (Nphs == 1) ? 1 : 0;
                    countf20_2 += (Nphs == 2) ? 1 : 0;
                    countf10P  += Nphs;
                    countf10_1 += (Nphs == 1) ? 1 : 0;
                    countf10_2 += (Nphs == 2) ? 1 : 0;
                    if (ferr0 < 0.05 || (ferr1 < 0.05 && Nphs == 2) )
                        {countf05P += 1; }
                }
                else if (ferrmax < 0.20)
                {
                    countf20E += 1;
                    countf20P  += Nphs;
                    countf20_1 += (Nphs == 1) ? 1 : 0;
                    countf20_2 += (Nphs == 2) ? 1 : 0;
                    if (ferr0 < 0.05 || (ferr1 < 0.05 && Nphs == 2) )
                        {countf10P += 1;
                         countf05P += 1; }
                    else if (ferr0 < 0.10 || (ferr1 < 0.10 && Nphs == 2))
                        {countf10P += 1; }
                }
                else if (ferr0 < 0.05 || (ferr1 < 0.05 && Nphs == 2))
                    {countf20P += 1;
                     countf10P += 1;
                     countf05P += 1; }
                else if (ferr0 < 0.10 || (ferr1 < 0.10 && Nphs == 2) )
                        {countf20P += 1;
                         countf10P += 1; }
                else if (ferr0 < 0.20 || (ferr1 < 0.20 && Nphs == 2) )
                        {countf20P += 1;}
            }
  
            // UPDATE HISTS_________________________________________
            double Dxy = 0; 
            int    Dl  = 0;
            if (Nphs == 2)
            {
              double dx = cp_vtxX->at(1) - cp_vtxX->at(0);
              double dy = cp_vtxY->at(1) - cp_vtxY->at(0);
              Dxy = sqrt( dx*dx + dy*dy );
              Dl  = abs(ts_firstlayer->at(1)-ts_firstlayer->at(0));

              hnorm->Fill(Dxy);//, Dl);
              if (emEnergy_clue->size()==2)
              { 
                hph -> Fill(Dxy);//, Dl);
                if (errmax < 0.05)
                {
                  h20->Fill(Dxy);//, Dl);
                  h10->Fill(Dxy);//, Dl);
                  h05->Fill(Dxy);//, Dl);
                }
                else if (errmax < 0.10)
                {
                  h20->Fill(Dxy);//, Dl);
                  h10->Fill(Dxy);//, Dl);
                }
                else if (errmax < 0.20)
                {
                  h20->Fill(Dxy);//, Dl); 
                }

                if (ferrmax < 0.05)
                {
                  hf20->Fill(Dxy);//, Dl);
                  hf10->Fill(Dxy);//, Dl);
                  hf05->Fill(Dxy);//, Dl);
                }
                else if (ferrmax < 0.10)
                {
                  hf20->Fill(Dxy);//, Dl);
                  hf10->Fill(Dxy);//, Dl);
                }
                else if (ferrmax < 0.20)
                {
                  hf20->Fill(Dxy);//, Dl); 
                }
              }   
                
              // CORRELATIONS CALCULATIONS ___________________________
              hdelta   ->Fill (Etrue0 - Eclue0);
              hreldelta->Fill ((Etrue0 - Eclue0)/Etrue0);
              hratio   ->Fill (Etrue0 / Eclue0);
              if (Nphs == 2)
              {
                hdelta   ->Fill (Etrue1 - Eclue1);
                hreldelta->Fill((Etrue1 - Eclue1)/Etrue1);
                hratio   ->Fill (Etrue1 / Eclue1);
              }
            }
        } //loop events
      dir->Close();
    }// if dir
    file -> Close();
    }// if file
    }// runs
    }// etas
    }// pts
    }// paths
    cout<<"Loop done"<<endl;

    printf("\nCLUE3D Accuracies\n");
    printf("Accuracy (single class)  : %4f \n",
          (double)countph_1/countall_1);
    printf("Accuracy (double class)  : %4f \n\n",
          (double)countph_2/countall_2);

    printf("CLUE3D Accuracies - Absolute Energy wise\n");
    printf("Event  Accuracy (0.20): %4f\n"  ,
          (double)count20E /countallE);
    printf("Photon Accuracy (0.20): %4f\n"  ,
          (double)count20P /countallP);
    printf("One-Ph Accuracy (0.20): %4f\n"  ,
          (double)count20_1/countall_1);
    printf("Two-Ph Accuracy (0.20): %4f\n\n",
          (double)count20_2/countall_2);

    printf("Event  Accuracy (0.10): %4f\n" ,
          (double)count10E /countallE);
    printf("Photon Accuracy (0.10): %4f\n" ,
          (double)count10P /countallP);
    printf("One-Ph Accuracy (0.10): %4f\n" ,
          (double)count10_1/countall_1);
    printf("Two-Ph Accuracy (0.10): %4f\n\n",
          (double)count10_2/countall_2);

    printf("Event  Accuracy (0.05): %4f \n" ,
          (double)count05E /countallE);
    printf("Photon Accuracy (0.05): %4f \n" ,
          (double)count05P /countallP);
    printf("One-Ph Accuracy (0.05): %4f \n" ,
          (double)count05_1/countall_1);
    printf("Two-Ph Accuracy (0.05): %4f \n\n",
          (double)count05_2/countall_2);


    printf("CLUE3D Accuracies - Fractional Energy wise\n");
    printf("Event  Accuracy (0.20): %4f\n"  ,
          (double)countf20E /countallE);
    printf("Photon Accuracy (0.20): %4f\n"  ,
          (double)countf20P /countallP);
    printf("One-Ph Accuracy (0.20): %4f\n"  ,
          (double)countf20_1/countall_1);
    printf("Two-Ph Accuracy (0.20): %4f\n\n",
          (double)countf20_2/countall_2);

    printf("Event  Accuracy (0.10): %4f\n" ,
          (double)countf10E /countallE);
    printf("Photon Accuracy (0.10): %4f\n" ,
          (double)countf10P /countallP);
    printf("One-Ph Accuracy (0.10): %4f\n" ,
          (double)countf10_1/countall_1);
    printf("Two-Ph Accuracy (0.10): %4f\n\n",
          (double)countf10_2/countall_2);

    printf("Event  Accuracy (0.05): %4f \n" ,
          (double)countf05E /countallE);
    printf("Photon Accuracy (0.05): %4f \n" ,
          (double)countf05P /countallP);
    printf("One-Ph Accuracy (0.05): %4f \n" ,
          (double)countf05_1/countall_1);
    printf("Two-Ph Accuracy (0.05): %4f\n\n",
          (double)countf05_2/countall_2);


    // Normalise
    //TH1D hnorm_bit = (*hnorm)*0.00001;
    
    TH1D hphnorm = (*hph) / (*hnorm);
    //hphnorm = hphnorm / (*hnorm);
    hphnorm.SetContour(50); //nicer, default is like 10
    hphnorm.GetXaxis()->SetTitle("Conversion Vtx XY Distances");
    //hphnorm.GetYaxis()->SetTitle("Delta Start Layer");
    //hphnorm.GetZaxis()->SetTitle("Reconstruction Efficiency (Class)");

    TH1D h20norm = (*h20 ) / (*hnorm);// + hnorm_bit;
    //h20norm = h20norm / (*hnorm);
    h20norm.SetContour(50); //nicer, default is like 10
    h20norm.GetXaxis()->SetTitle("Conersion Vtx XY Distances");
    //h20norm.GetYaxis()->SetTitle("Delta Start Layer");
    //h20norm.GetZaxis()->SetTitle("Reconstruction Efficiency (to 20%%)");
   
    TH1D h10norm = (*h10) / (*hnorm);// + hnorm_bit;
    //h10norm = h10norm / (*hnorm);
    h10norm.SetContour(50); //nicer, default is like 10
    h10norm.GetXaxis()->SetTitle("Conersion Vtx XY Distances");
    //h10norm.GetYaxis()->SetTitle("Delta Start Layer");
    //h10norm.GetZaxis()->SetTitle("Reconstruction Efficiency (to 10%%)");

    TH1D h05norm = (*h05) / (*hnorm);// + hnorm_bit;
    //h05norm = h05norm / (*hnorm);
    h05norm.SetContour(50); //nicer, default is like 10
    h05norm.GetXaxis()->SetTitle("Conersion Vtx XY Distances");
    //h05norm.GetYaxis()->SetTitle("Delta Start Layer");
    //h05norm.GetZaxis()->SetTitle("Reconstruction Efficiency (to 5%%)");

    TH1D hf20norm = (*hf20 ) / (*hnorm);// + hnorm_bit;
    //h20norm = h20norm / (*hnorm);
    hf20norm.SetContour(50); //nicer, default is like 10
    hf20norm.GetXaxis()->SetTitle("Conersion Vtx XY Distances");
    //h20norm.GetYaxis()->SetTitle("Delta Start Layer");
    //h20norm.GetZaxis()->SetTitle("Reconstruction Efficiency (to 20%%)");
   
    TH1D hf10norm = (*hf10) / (*hnorm);// + hnorm_bit;
    //h10norm = h10norm / (*hnorm);
    hf10norm.SetContour(50); //nicer, default is like 10
    hf10norm.GetXaxis()->SetTitle("Conersion Vtx XY Distances");
    //h10norm.GetYaxis()->SetTitle("Delta Start Layer");
    //h10norm.GetZaxis()->SetTitle("Reconstruction Efficiency (to 10%%)");

    TH1D hf05norm = (*hf05) / (*hnorm);// + hnorm_bit;
    //h05norm = h05norm / (*hnorm);
    hf05norm.SetContour(50); //nicer, default is like 10
    hf05norm.GetXaxis()->SetTitle("Conersion Vtx XY Distances");
    //h05norm.GetYaxis()->SetTitle("Delta Start Layer");
    //h05norm.GetZaxis()->SetTitle("Reconstruction Efficiency (to 5%%)");


    //gStyle->SetPalette(kRainBow); //nice colour option for colz
    TCanvas* c2d = new TCanvas("Double Events Accuracy", 
                               "Double Events Accuracy"); 
    c2d -> Divide(2, 2);
    c2d -> cd(1);
    hphnorm.SetStats(0);
    hphnorm.Draw("HIST");//("colz");
    c2d -> cd(2);
    h20norm.SetStats(0);
    h20norm.Draw("HIST");//("colz");
    c2d -> cd(3);
    h10norm.SetStats(0);
    h10norm.Draw("HIST");
    c2d -> cd(4);
    h05norm.SetStats(0);
    h05norm.Draw("HIST");
    c2d-> Update();

    TCanvas* cf2d = new TCanvas("Double Events Accuracy (fractional)", 
                               "Double Events Accuracy (fractional)"); 
    cf2d -> Divide(2, 2);
    cf2d -> cd(1);
    hphnorm.SetStats(0);
    hphnorm.Draw("HIST");//("colz");
    cf2d -> cd(2);
    hf20norm.SetStats(0);
    hf20norm.Draw("HIST");//("colz");
    cf2d -> cd(3);
    hf10norm.SetStats(0);
    hf10norm.Draw("HIST");
    cf2d -> cd(4);
    hf05norm.SetStats(0);
    hf05norm.Draw("HIST");
    cf2d-> Update();

    TCanvas* c0 = new TCanvas();
    c0 -> Divide(1, 3);
    c0 -> cd(1);
    hdelta->GetXaxis()->SetTitle("True - CLUE3D");
    hdelta->Draw();
    c0 -> cd(2);
    hreldelta->GetXaxis()->SetTitle("Relative True - CLUE3D");
    hreldelta->Draw();
    c0 -> cd(3);
    hratio->GetXaxis()->SetTitle("True - CLUE3D Ratio");
    hratio->Draw();
    c0 -> Update();

    c2d -> WaitPrimitive();
    cf2d-> WaitPrimitive();
    c0  -> WaitPrimitive();   
}






/*
//__________________________________________________________________
    TH1D* hphP   =new TH1D("Photons Reconstruction Efficiency - Photon #",
                          "Photons Reconstruction Efficiency - Photon #",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h20P   =new TH1D("Photons Reconstruction Efficiency - 20%% sigma_E",
                          "Photons Reconstruction Efficiency - 20%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h10P   =new TH1D("Photons Reconstruction Efficiency - 10%% sigma_E",
                          "Photons Reconstruction Efficiency - 10%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h05P   =new TH1D("Photons Reconstruction Efficiency - 5%% sigma_E",
                          "Photons Reconstruction Efficiency - 5%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);

    TH1D* hnormP = new TH1D("Photons Normalisation2D", "Photons Normalisation2D",
                                        2*XYmax  , 0, XYmax,
                                        Dlayermax, 0, Dlayermax);
    
    //__________________________________________________________________
    TH1D* hph_1   =new TH1D("Single Reconstruction Efficiency - Photon #",
                          "Single Reconstruction Efficiency - Photon #",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h20_1   =new TH1D("Single Reconstruction Efficiency - 20%% sigma_E",
                          "Single Reconstruction Efficiency - 20%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h10_1   =new TH1D("Single Reconstruction Efficiency - 10%% sigma_E",
                          "Single Reconstruction Efficiency - 10%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h05_1   =new TH1D("Single Reconstruction Efficiency - 5%% sigma_E",
                          "Single Reconstruction Efficiency - 5%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);

    TH1D* hnorm_1 = new TH1D("Single Normalisation2D", "Single Normalisation2D",
                                        2*XYmax  , 0, XYmax,
                                        Dlayermax, 0, Dlayermax);
    
    //__________________________________________________________________
    TH1D* hph_2   =new TH1D("Double Reconstruction Efficiency - Photon #",
                          "Double Reconstruction Efficiency - Photon #",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h20_2   =new TH1D("Double Reconstruction Efficiency - 20%% sigma_E",
                          "Double Reconstruction Efficiency - 20%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h10_2   =new TH1D("Double Reconstruction Efficiency - 10%% sigma_E",
                          "Double Reconstruction Efficiency - 10%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);
    
    TH1D* h05_2   =new TH1D("Double Reconstruction Efficiency - 5%% sigma_E",
                          "Double Reconstruction Efficiency - 5%% sigma_E",
                            2*XYmax  , 0, XYmax,
                            Dlayermax, 0, Dlayermax);

    TH1D* hnorm_2 = new TH1D("Double Normalisation2D", "Double Normalisation2D",
                                        2*XYmax  , 0, XYmax,
                                        Dlayermax, 0, Dlayermax);
    */