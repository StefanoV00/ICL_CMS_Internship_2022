/// \author Stefano Veroni
/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 3.Analysis
root AnalysisCLUEenergies.C -q
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
 Find distribution of fraction of energy of 2 most energetic tracksters 
 produced by CLUE3D algorithm when reconstructing double-photons events.
*/
void AnalysisCLUEenergies ()
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
                             };

    string rootname= "step3ticl_pt50_eta21_run0_FlatTracksters";//root file
    string rootname1 = "step3ticl_";
    string rootname2 = "_eta";
    string rootname3 = "_run";
    string rootname4 = "_FlatTracksters";
    
    vector <string> pts  = {"En20to200"};
    vector <string> etas = {"21"};
    int             Nrun = 100;
    int istart   = 0;
    int ptstart  = 0;
    int etastart = 0;
    vector <int> runstart = {0, 0};
    
    string dirname1 = "ticlTree"; //"none" if no directory in between
    string treename1= "TSTree_SimFromCP";
    
    string dirname2 = "ticlTree"; //"none" if no directory in between
    string treename2= "TSTree_CLUE3D3"; 
    
    ////////////////////////////////////////////////////////////////////
    //1. Prepare Plot //////////////////////////////////////////////////
    //__________________________________________________________________

    TH1D* henergy   =new TH1D("On Extra Tracksters Reconstruction",
                              "On Extra Tracksters Reconstruction",
                              20, 0, 1);
    
    TH1D* hreconstruct=new TH1D("N Reconstructed Photons",
                                "N Reconstructed Photons",
                                10, 0, 10);

       
////////////////////////////////////////////////////////////////////////
//////2. LOOP OVER TTREES //////////////////////////////////////////////
/////__________________________________________________________________
    int countall = 0;
    int count1   = 0;
    int count2   = 0;
    int countmore= 0;

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
    Nrun *= (i_p == paths.size()-1) ? 5 : 1;
    runstart[pt_i] *= (i_p == paths.size()-1) ? 5 : 1;
    for (int run = runstart[pt_i]; run < Nrun; run++)
    {
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
        vector <double>* energy_all= 0;
        vector <double>* tsmult_all= 0;
        vector <int>*    layer_all= 0;
        vector <int>*    which_g  = 0;
        ttrue->SetBranchAddress("lc_energy", &energy_all);
        ttrue->SetBranchAddress("lc_tsMult", &tsmult_all);
        ttrue->SetBranchAddress("lc_layer" , &layer_all);
        ttrue->SetBranchAddress("lc_TSidx" , &which_g);
        vector <double>* emEnergy_true = 0;
        vector <double>* emEnergy_clue = 0;
        ttrue->SetBranchAddress("ts_emEnergy",   &emEnergy_true);
        tclue->SetBranchAddress("ts_emEnergy",   &emEnergy_clue);


        ////////////////////////////////////////////////////////////
        //4. Loop, compare, fill ///////////////////////////////////
        //__________________________________________________________
        auto isecal   = [](int layN)     {return layN  <=28;};
        auto distinct = [](double tsmult){return tsmult<= 2;};
        for (int i = 0; i<Nentries; i++) //Loop & Fill
        {
            ttrue  -> GetEntry(i);
            tclue  -> GetEntry(i);
            int Nphs  =  emEnergy_true->size();
            int Nclue =  emEnergy_clue->size();
            
            countall  += 1;
            count1    += (Nclue== 1) ? 1 : 0;
            count2    += (Nclue== 2) ? 1 : 0;
            countmore += (Nclue > 2) ? 1 : 0;

            // Remove everything that is NOT E-Cal info
            vector<double> energyb=getAwhereB(*energy_all,*layer_all,
                                                            isecal);
            vector<double> tsmultb=getAwhereB(*tsmult_all,*layer_all,
                                                            isecal);
            vector<int>    layerb =getAwhereB(*layer_all ,*layer_all,
                                                               isecal);
            vector<int>    gb     =getAwhereB(*which_g   ,*layer_all,
                                                               isecal);
            // Get Energies
            for (int gj = 0; gj<emEnergy_true->size(); gj++)
              {emEnergy_true->at(gj) = 0;}
            for (int j = 0; j<gb.size(); j++)
              {emEnergy_true->at(gb[j]) += energyb[j]/tsmultb[j];}

            if (Nclue > 2)
            {
                std::sort (emEnergy_clue->begin(), emEnergy_clue->end()); 
                double Emax = emEnergy_clue->at(Nclue-1);
                double Epremax = emEnergy_clue->at(Nclue-2);
                double Etot  = myvecsum(*emEnergy_clue);
                double Efrac = (Emax + Epremax) / Etot; 

                henergy->Fill(Efrac);
            }  

            hreconstruct->Fill(Nclue); 

        } //loop events
    }// if dir
    file -> Close();
    }// if file
    }// runs
    }// etas
    }// pts
    }// paths
    cout<<"Loop done"<<endl;

    printf("\nCLUE3D Double-Events Reconstruction\n");
    printf(" - Double->Single  : %4f \n\n",(double)count1/countall);
    printf(" - Double->Double  : %4f \n\n",(double)count2/countall);
    printf(" - Double->More    : %4f \n\n",(double)countmore/countall);

    TCanvas* c0 = new TCanvas();
    double countalld = (double) countall;
    hreconstruct->Scale(1/countalld); 
    hreconstruct->GetXaxis()
                ->SetTitle("N Reconstructed Photons");
    hreconstruct->Draw();
    c0 -> Update();

    TCanvas* c1 = new TCanvas();
    henergy->GetXaxis()
           ->SetTitle("E-fraction of 2 Most Energetic Tracksters");
    henergy->Draw();
    c1 -> Update();

    c0 -> WaitPrimitive();
    c1 -> WaitPrimitive();

    return;  
}