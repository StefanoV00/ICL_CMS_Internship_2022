/// \author Stefano Veroni
/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 1.AddTTrees
root DoublesAddAAA.C -q
root DoublesAddDeltas.C -q
root DoublesAddEnergies.C -q
root DoublesAddShape.C -q
root DoublesAddTheSelections.C -q
cd ../

*/

#include <iostream>
#include <iomanip>
#include <fstream>
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
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "Math/Vector2D.h"

#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"
using namespace std;

/* 
 
 ADD TDIRECTORIES TO treeFriends.root FILE FOR:
 - entries with first seed DXY <= 3cm;
 - entries with first seed NOT DXY <= 3cm;
  - entries satisfying both of the above;
 - entries satisfying none of the above;
 - entries satisfying one of the above only (2);
 - entries with first detections in same layer;
 - entries with first detections NOT in same layer;
*/


void DoublesAddTheSelections ()
{
    string path    = "../myTTrees/TreesEn90to110/Double/"; //none if none
    string rootname= "step3ticl_eta21_FlatTracksters";//root file
    string friendrootname= "treeFriends"; //root file with friends
    string selrootname0= "treeSelections"; 
    string selrootname1= "treeSelectionsBad"; 
    string dirname = "none";
    string treename= "TSTree_SimFromCP";


    ///////////////////////////////////////////////////////////////////
    // 1. GET TTREES AND SET ADDRESSES   //////////////////////////////
    ///////////////////////////////////////////////////////////////////
    TFile* file = open_file(path, rootname);
    TTree* tree = get_tree (file, dirname, treename);
    int Nentries = tree->GetEntries();

    TFile* friends = open_file(path, friendrootname, "read");
    const char* coords_name   = (treename + "_Coords").c_str();
    TTree* tfriendC = get_tree(friends, "none", coords_name);
    const char* deltas_name   = (treename + "_Deltas").c_str();
    TTree* tfriendD = get_tree(friends, "none", deltas_name);
    const char* energies_name = (treename + "_Energy").c_str();
    TTree* tfriendE = get_tree(friends, "none", energies_name);
    const char* shape_name    = (treename + "_Shape").c_str();
    TTree* tfriendS = get_tree(friends, "none", shape_name);
    //const char* core_name    = (treename + "_ShapeCore").c_str();
    //TTree* tfriendSC = get_tree(friends, "none", core_name);

    double Dsz =0;tfriendD->SetBranchAddress("lc_deltaStartSeedZ" ,&Dsz);
    int    Dsl =0;tfriendD->SetBranchAddress("ts_deltaStartLayer" ,&Dsl);
    double Dsxy=0;tfriendD->SetBranchAddress("cp_deltaVtxXY",&Dsxy);


    ///////////////////////////////////////////////////////////////////
    // 2. PREPARE NEW tDIRECTORIES AND TTREES   ///////////////////////
    ///////////////////////////////////////////////////////////////////
    //TDirectory* dir_close_same= fr->mkdir("CloseXYSameStart",
    //                                      "CloseXYSameStart");
    TFile* fr0 = open_file(path, selrootname0, "recreate");
    TDirectory* dir_close_same= open_dir(fr0, "CloseXYSameStart");
    dir_close_same -> cd();
    TTree* close_same_0 = tree     -> CloneTree(0);
    TTree* close_same_c = tfriendC -> CloneTree(0);
    TTree* close_same_d = tfriendD -> CloneTree(0);
    TTree* close_same_e = tfriendE -> CloneTree(0);
    TTree* close_same_s = tfriendS -> CloneTree(0);
    //TTree* close_same_sc= tfriendSC-> CloneTree(0);

    //TDirectory* dir_close_diff= fr->mkdir("CloseXYDiffStart",
    //                                      "CloseXYDiffStart");
    TDirectory* dir_close_diff= open_dir(fr0, "CloseXYDiffStart");
    dir_close_diff -> cd();
    TTree* close_diff_0 = tree     -> CloneTree(0);
    TTree* close_diff_c = tfriendC -> CloneTree(0);
    TTree* close_diff_d = tfriendD -> CloneTree(0);
    TTree* close_diff_e = tfriendE -> CloneTree(0);
    TTree* close_diff_s = tfriendS -> CloneTree(0);
    //TTree* close_diff_sc= tfriendSC-> CloneTree(0);

    TDirectory* dir_close_samepm= open_dir(fr0, "CloseXYSameStart+-");
    dir_close_samepm -> cd();
    TTree* close_samepm_0 = tree     -> CloneTree(0);
    TTree* close_samepm_c = tfriendC -> CloneTree(0);
    TTree* close_samepm_d = tfriendD -> CloneTree(0);
    TTree* close_samepm_e = tfriendE -> CloneTree(0);
    TTree* close_samepm_s = tfriendS -> CloneTree(0);
    //TTree* close_samepm_sc= tfriendSC-> CloneTree(0);

    //TDirectory* dir_close_diff= fr->mkdir("CloseXYDiffStart",
    //                                      "CloseXYDiffStart");
    TDirectory* dir_close_diffpm= open_dir(fr0, "CloseXYDiffStart+-");
    dir_close_diffpm -> cd();
    TTree* close_diffpm_0 = tree     -> CloneTree(0);
    TTree* close_diffpm_c = tfriendC -> CloneTree(0);
    TTree* close_diffpm_d = tfriendD -> CloneTree(0);
    TTree* close_diffpm_e = tfriendE -> CloneTree(0);
    TTree* close_diffpm_s = tfriendS -> CloneTree(0);
    //TTree* close_diffpm_sc= tfriendSC-> CloneTree(0);

    //TDirectory* dir_far_same = fr->mkdir("FarXYSameStart",
    //                                     "FarXYSameStart");
    TDirectory* dir_far_same = open_dir(fr0, "FarXYSameStart");
    dir_far_same -> cd(); 
    TTree* far_same_0 = tree     -> CloneTree(0);
    TTree* far_same_c = tfriendC -> CloneTree(0);
    TTree* far_same_d = tfriendD -> CloneTree(0);
    TTree* far_same_e = tfriendE -> CloneTree(0);
    TTree* far_same_s = tfriendS -> CloneTree(0);
    //TTree* far_same_sc= tfriendSC-> CloneTree(0);


    //TDirectory* dir_close= fr->mkdir("CloseXY", "CloseXY");
    TDirectory* dir_close= open_dir(fr0, "CloseXY");
    dir_close -> cd();
    TTree* close_0 = tree     -> CloneTree(0);
    TTree* close_c = tfriendC -> CloneTree(0);
    TTree* close_d = tfriendD -> CloneTree(0);
    TTree* close_e = tfriendE -> CloneTree(0);
    TTree* close_s = tfriendS -> CloneTree(0);
    //TTree* close_sc= tfriendSC-> CloneTree(0);

    //TDirectory* dir_same = fr->mkdir("SameStart", "SameStart");
    TDirectory* dir_same = open_dir(fr0, "SameStart");
    dir_same -> cd();
    TTree* same_0 = tree     -> CloneTree(0);
    TTree* same_c = tfriendC -> CloneTree(0);
    TTree* same_d = tfriendD -> CloneTree(0);
    TTree* same_e = tfriendE -> CloneTree(0);
    TTree* same_s = tfriendS -> CloneTree(0);
    //TTree* same_sc= tfriendSC-> CloneTree(0);

    //TDirectory* dir_far = fr->mkdir("FarXY", "FarXY");
    TDirectory* dir_far = open_dir(fr0, "FarXY");
    dir_far -> cd();
    TTree* far_0 = tree     -> CloneTree(0);
    TTree* far_c = tfriendC -> CloneTree(0);
    TTree* far_d = tfriendD -> CloneTree(0);
    TTree* far_e = tfriendE -> CloneTree(0);
    TTree* far_s = tfriendS -> CloneTree(0);
    //TTree* far_sc= tfriendSC-> CloneTree(0);
   
    //TDirectory* dir_diff = fr->mkdir("DiffStart", "DiffStart");
    TDirectory* dir_diff = open_dir(fr0, "DiffStart");
    dir_diff -> cd();
    TTree* diff_0 = tree     -> CloneTree(0);
    TTree* diff_c = tfriendC -> CloneTree(0);
    TTree* diff_d = tfriendD -> CloneTree(0);
    TTree* diff_e = tfriendE -> CloneTree(0);
    TTree* diff_s = tfriendS -> CloneTree(0);
    //TTree* diff_sc= tfriendSC-> CloneTree(0);

    //TDirectory* dir_far_diff = fr->mkdir("FarXYDiffStart",
    //                                     "FarXYDiffStart");
    TDirectory* dir_far_diff = open_dir(fr0, "FarXYDiffStart");
    dir_far_diff -> cd(); 
    TTree* far_diff_0 = tree     -> CloneTree(0);
    TTree* far_diff_c = tfriendC -> CloneTree(0);
    TTree* far_diff_d = tfriendD -> CloneTree(0);
    TTree* far_diff_e = tfriendE -> CloneTree(0);
    TTree* far_diff_s = tfriendS -> CloneTree(0);
    //TTree* far_diff_sc= tfriendSC-> CloneTree(0);
   


    ///////////////////////////////////////////////////////////////////
    // 3. LOOP AND COPY ENTRIES    ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    bool closefirst = 0;
    bool samefirst  = 0;
    for (int i = 0; i< Nentries; i++)
    {
        cout<<"Entry: "<<i<<endl;
        tree     -> GetEntry(i);
        tfriendC -> GetEntry(i);
        tfriendD -> GetEntry(i);
        tfriendE -> GetEntry(i);
        tfriendS -> GetEntry(i);
        //tfriendSC-> GetEntry(i);

        bool close = Dsxy <= 3;
        bool same  = Dsl  == 0;
        bool samepm= abs(Dsl) <= 1;

        if (close)
        {
            close_0 -> Fill();
            close_c -> Fill();
            close_d -> Fill();
            close_e -> Fill();
            close_s -> Fill();
            //close_sc-> Fill();
            if (same) //close and same
            {
                same_0 -> Fill();
                same_c -> Fill();
                same_d -> Fill();
                same_e -> Fill();
                same_s -> Fill();
                //same_sc-> Fill();

                close_same_0 -> Fill();
                close_same_c -> Fill();
                close_same_d -> Fill();
                close_same_e -> Fill();
                close_same_s -> Fill();
                //close_same_sc-> Fill();
            }
            if (samepm) //close and same (plus or minus 1)
            {
                close_samepm_0 -> Fill();
                close_samepm_c -> Fill();
                close_samepm_d -> Fill();
                close_samepm_e -> Fill();
                close_samepm_s -> Fill();
                //close_samepm_sc-> Fill();
            }
            if (!same)     //close and with diff>=1
            {
            //continue;
                diff_0 -> Fill();
                diff_c -> Fill();
                diff_d -> Fill();
                diff_e -> Fill();
                diff_s -> Fill();
                //diff_sc-> Fill();

                close_diff_0 -> Fill();
                close_diff_c -> Fill();
                close_diff_d -> Fill();
                close_diff_e -> Fill();
                close_diff_s -> Fill();
                //close_diff_sc-> Fill();
            }
            if (!samepm)  // close and with diff>=2
            {

                close_diffpm_0 -> Fill();
                close_diffpm_c -> Fill();
                close_diffpm_d -> Fill();
                close_diffpm_e -> Fill();
                close_diffpm_s -> Fill();
                //close_diffpm_sc-> Fill();
            }
        }
        else //far
        {//continue;
            far_0 -> Fill();
            far_c -> Fill();
            far_d -> Fill();
            far_e -> Fill();
            far_s -> Fill();
            //far_sc-> Fill();
            if (same) //far and same
            {
                same_0 -> Fill();
                same_c -> Fill();
                same_d -> Fill();
                same_e -> Fill();
                same_s -> Fill();
                //same_sc-> Fill();

                far_same_0 -> Fill();
                far_same_c -> Fill();
                far_same_d -> Fill();
                far_same_e -> Fill();
                far_same_s -> Fill();
                //far_same_sc-> Fill();
            }
            else      //far and diff
            {
                diff_0 -> Fill();
                diff_c -> Fill();
                diff_d -> Fill();
                diff_e -> Fill();
                diff_s -> Fill();
                //diff_sc-> Fill();

                far_diff_0 -> Fill();
                far_diff_c -> Fill();
                far_diff_d -> Fill();
                far_diff_e -> Fill();
                far_diff_s -> Fill();
                //far_diff_sc-> Fill();
            }
        }
    }
    cout<<"Loop Finished"<<endl;
    fr0-> Write("", TObject::kOverwrite);
    //fr1-> Write("", TObject::kOverwrite);//->Write();
    cout<<"Files Written"<<endl;
} 