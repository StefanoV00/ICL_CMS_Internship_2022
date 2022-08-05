/// \author Stefano Veroni
/// File to save events with energy between
/// 95 and 105 GeV in myTTrees/TreesEn95to100

/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 3.Analysis
root SaveSpecEn.C -q
cd ../

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <numeric>
#include <string>
#include <stdio.h>
#include <cmath>
#include <vector>


#include "TFile.h"
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"


using namespace std;

/////////////////////////////////////////////////////////////////////////
//-1. Functions' Definitions/////////////////////////////////////////////
//______________________________________________________________________
#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"




void SaveSpecEn (bool save_txt = true) 
{
////////////////////////////////////////////////////////////////////////
//////0. Set Variables   ///////////////////////////////////////////////
//////__________________________________________________________________
    string path    = "../myTTrees/";//none if none
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
    
    string rootname1 = "step3ticl_";
    string rootname2 = "_eta";
    string rootname3 = "_run";
    string rootname4 = "_FlatTracksters";
    
    vector <string> pts  = {"En20to200", "En40to400"};
    vector <string> etas = {"21"};
    int             Nrun = 100;

    int pathstart= 0;
    int istart   = 0;
    int ptstart  = 0;
    int etastart = 0;
    int runstart = 0;

    string dirname = "ticlTree";
    string treename= "TSTree_SimFromCP";
    

////////////////////////////////////////////////////////////////////////
/////1. PREPARE SAVING TTREES //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
    double Ed = 90. ;
    double Eu = 110.;
    //Get example tree
    string rootname= "step3ticl_pt50_eta21_run0_FlatTracksters";
    TFile* file = new TFile((path+rootname+".root").c_str(),
                                "read");
    TDirectory* dir = (TDirectory*)file->Get(dirname.c_str());
    TTree* tree0 = (TTree*) dir->Get(treename.c_str());

    
    string savepath1  = "../myTTrees/TreesEn"+
                       to_string((int)Ed)+"to"
                       +to_string((int)Eu)+"/Double/";
    string savepath2  = "../myTTrees/TreesEn"+
                       to_string((int)Ed)+"to"
                       +to_string((int)Eu)+"/Single/";
    string myrootname1 = "step3ticl_eta21_FlatTracksters";
    string myrootname2 = "step3ticl_eta21_FlatTracksters";
    TFile* twophs = new TFile((savepath1+myrootname1+".root").c_str(),
                                "recreate");
    twophs->cd();
    TTree* twotree = tree0->CloneTree(0);
    TFile* onephs = new TFile((savepath2+myrootname2+".root").c_str(),
                                "recreate");
    onephs->cd();
    TTree* onetree = tree0->CloneTree(0);
    dir -> Close();
    file -> Close();

    int onecounter = 0;
    int twocounter = 0;


////////////////////////////////////////////////////////////////////////
//////2. LOOP OVER TTREES //////////////////////////////////////////////
/////__________________________________________________________________
    for (int i_p = pathstart; i_p < paths.size(); i_p++)
    {
    string path = paths[i_p];
    for (int pt_i = ptstart; pt_i < pts.size(); pt_i++)
    {
    string rootname_pt = rootname1 + pts[pt_i];
    for (int e_i = etastart; e_i < etas.size(); e_i++)
    {
    string rootname_eta = rootname2 + etas[e_i];
    if (e_i == 0 && pt_i == 0) Nrun *= (i_p == paths.size()-1) ? 6 : 1;
    for (int run = runstart; run < Nrun; run++)
    {
        printf("Path Index: %d, ", i_p);
        printf("PT Index: %d, RUN N. is %2d --------\n", pt_i, run);
        string rootname_run= rootname3 + to_string(run) + rootname4;

        string rootname_p_e_r = rootname_pt  +
                                rootname_eta +
                                rootname_run ; 
                     
        if (gSystem->AccessPathName((path+rootname_p_e_r+".root")
                                    .c_str()) == 0 )
        {
        TFile* file = new TFile((path+rootname_p_e_r+".root").c_str(),
                                "read");
        TDirectory* dir = (TDirectory*)file->Get(dirname.c_str());
        if (dir)
        {
            ////////////////////////////////////////////////////////////
            //3. Get TTree & Branches///////////////////////////////////
            //////______________________________________________________
            TTree* tree = (TTree*) dir->Get(treename.c_str());
            vector <double>* totenergy = 0;
            tree -> SetBranchAddress("ts_energy", &totenergy);
            int Nentries = tree -> GetEntries();
            
            for (int i = 0; i<Nentries; i++)
            {
                tree->GetEntry(i);

                auto N  = totenergy->size();
                auto E0 = totenergy->at(0);
                auto E1 = (N == 2) ? totenergy->at(1) : 0.; 
                auto E  = E0 + E1;
                
                if (N==1)
                {
                    if (Ed < E && E < Eu)
                        {onecounter += 1;
                        cout<<E0<<" + "<<E1<<" -> "<<E<<endl;
                        tree->CopyAddresses(onetree);
                        onetree -> Fill();
                        }
                }
                else if (N==2)
                {
                    if (Ed < E && E < Eu)
                        {twocounter += 1;
                        cout<<E0<<" + "<<E1<<" -> "<<E<<endl;
                        tree->CopyAddresses(twotree);
                        twotree -> Fill();
                        }
                } 
            }//end of loop over entries
            file->Close();
            printf("Run %d saved \n", run);
        } //end of  (if dir)
        } //end of (if file)
    } // end of loop over runs
    runstart = 0;
    } // end of loop over etas
    etastart = 0;
    } // end of loop over pts
    ptstart = 0;
    }// end of paths 
    onetree->SetName(treename.c_str());
    twotree->SetName(treename.c_str());

    onephs-> Write("", TObject::kOverwrite);
    twophs-> Write("", TObject::kOverwrite);

    cout<<"------ FINISHED ------"<<endl;
    printf("TOT number of one events in (%f, %f): %d\n", 
    Ed, Eu, onecounter);
    printf("TOT number of two events in (%f, %f): %d\n", 
    Ed, Eu, twocounter);

}//end of main