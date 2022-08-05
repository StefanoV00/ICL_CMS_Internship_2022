/// \author Stefano Veroni
/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 1.AddTTrees
root AddAAA.C -q
root AddDeltas.C -q
root AddEnergies.C -q
root AddEnergiesWrong.C -q
root AddShape.C -q
root AddShapeCore.C -q
root AddTheSelections.C -q
cd ../

*/
 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cmath>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include <TH2.h>
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"
using namespace std;


/*
 ADD DELTA INFORMATION: DISTANCES BETWEEN FIRST DETECTIONS:
 - Layer-wise
 - Coordinates-wise of lc first detection
 - Coordinates-wise of lc first seed detection
 - Coordinates-wise of sv vertex
 - Coordinates-wise of cp vertex (conversion)
*/
void AddDeltas () 
{
    string path    = "../myTTrees/TreesEn90to110/Double/"; //none if none
    string rootname= "step3ticl_eta21_FlatTracksters";//root file
    string friendrootname= "treeFriends"; //root file with friends
    string dirname = "none";
    string treename= "TSTree_SimFromCP";

    int icheck = 5; //Entry to check in loop to know where errors lie.


    //1. //////////////////////////////////////////////////////////////
    //Get stuff from the TTree we already have
    ///////////////////////////////////////////////////////////////////
    TFile *file = open_file(path, rootname);
    TTree *tree = get_tree (file, dirname, treename);
    //tree ->Scan("lc_TSidx:lc_nrechits:lc_layer:lc_mult:lc_tsMult:lc_seedEnergy");
    //return;
    int Nentries = tree->GetEntries();
    vector <double>* x    = 0;
    vector <double>* y    = 0;
    vector <double>* z    = 0;
    vector <double>* eta  = 0;
    vector <double>* phi   = 0;
    vector <double>* etaseed = 0;
    vector <double>* phiseed = 0;
    vector <double>* posX  = 0;
    vector <double>* posY  = 0;
    vector <double>* posZ  = 0;
    vector <double>* vtxX  = 0;
    vector <double>* vtxY  = 0;
    vector <double>* vtxZ  = 0;
    vector <double>* cpeta  = 0;
    vector <double>* cpphi  = 0;
    vector <int>* which_g = 0;
    vector <int>* layerN  = 0;
    vector <int>* layerfirst = 0;
    vector <int>* layerlast  = 0;
    tree->SetBranchAddress("lc_x",     &x);
    tree->SetBranchAddress("lc_y",     &y);
    tree->SetBranchAddress("lc_z",     &z);
    tree->SetBranchAddress("lc_eta",   &eta);
    tree->SetBranchAddress("lc_phi",   &phi);
    tree->SetBranchAddress("lc_seedEta", &etaseed);
    tree->SetBranchAddress("lc_seedPhi", &phiseed);
    tree->SetBranchAddress("lc_layer", &layerN);
    tree->SetBranchAddress("lc_TSidx", &which_g);
    tree->SetBranchAddress("sv_posX",  &posX);//SimVertex
    tree->SetBranchAddress("sv_posY",  &posY);
    tree->SetBranchAddress("sv_posZ",  &posZ);
    tree->SetBranchAddress("cp_vtxX",  &vtxX); //1st-conversion vertex
    tree->SetBranchAddress("cp_vtxY",  &vtxY);
    tree->SetBranchAddress("cp_vtxZ",  &vtxZ);
    tree->SetBranchAddress("cp_eta",   &cpeta);
    tree->SetBranchAddress("cp_phi",   &cpphi);
    tree->SetBranchAddress("cp_vtxY",  &vtxY);
    tree->SetBranchAddress("ts_firstLayer", &layerfirst);
    tree->SetBranchAddress("ts_lastLayer",  &layerlast);
    //tree->Scan("event:lc_TSidx:lc_x:lc_layer:sv_posX:cp_vtxX:cp_phi");


    //2. //////////////////////////////////////////////////////////////
    //Prepare new TTree, which will later become friend
    ///////////////////////////////////////////////////////////////////
    TFile* fr = open_file(path, friendrootname, "update");
    const char* coords_name = (treename + "_Coords").c_str();
    TTree* tfriendC = get_tree(fr, "none", coords_name);
    const char* friend_name = (treename + "_Deltas").c_str();
    TTree *tfriend = new TTree(friend_name, friend_name);
    tfriendC -> SetBranchStatus("*",0); //disable all existing branches
    tfriendC -> SetBranchStatus("lc_seedX",1);
    tfriendC -> SetBranchStatus("lc_seedY",1);
    tfriendC -> SetBranchStatus("lc_seedZ",1);
    vector<double>*xseed=0; tfriendC->SetBranchAddress("lc_seedX",&xseed);
    vector<double>*yseed=0; tfriendC->SetBranchAddress("lc_seedY",&yseed);
    vector<double>*zseed=0; tfriendC->SetBranchAddress("lc_seedZ",&zseed);
    //following use ts info
    int Dl;           tfriend -> Branch("ts_deltaStartLayer", &Dl); 
    //lc_Distances between first detections, use lc info
    double Dx;        tfriend -> Branch("lc_deltaStartX",  &Dx);       
    double Dy;        tfriend -> Branch("lc_deltaStartY",  &Dy);         
    double Dz;        tfriend -> Branch("lc_deltaStartZ",  &Dz);          
    double Dt;        tfriend -> Branch("lc_deltaStartXY", &Dt);
    double D ;        tfriend -> Branch("lc_deltaStart",   &D);          
    double Deta;      tfriend -> Branch("lc_deltaStartEta",&Deta);           
    double Dphi;      tfriend -> Branch("lc_deltaStartPhi",&Dphi);          
    double Dtheta;  //tfriend -> Branch("lc_deltaStartTheta", &Dtheta); 
    //seed_Distances between first seed detections, use true lc_info
    double Dxseed;    tfriend -> Branch("lc_deltaStartSeedX",  &Dxseed);
    double Dyseed;    tfriend -> Branch("lc_deltaStartSeedY",  &Dyseed);
    double Dzseed;    tfriend -> Branch("lc_deltaStartSeedZ",  &Dzseed);
    double Dtseed;    tfriend -> Branch("lc_deltaStartSeedXY", &Dtseed);
    double Dseed;     tfriend -> Branch("lc_deltaStartSeed",   &Dseed) ;
    double Detaseed;  tfriend -> Branch("lc_deltaStartSeedEta",&Detaseed);
    double Dphiseed;  tfriend -> Branch("lc_deltaStartSeedPhi",&Dphiseed);
    double Dthetaseed;    
    //sv_Distances between first detections, use true sv_info
    double Dxsv;      tfriend -> Branch("sv_deltaPosX",  &Dxsv);       
    double Dysv;      tfriend -> Branch("sv_deltaPosY",  &Dysv);          
    double Dzsv;      tfriend -> Branch("sv_deltaPosZ",  &Dzsv);            
    double Dtsv;      tfriend -> Branch("sv_deltaPosXY", &Dtsv);             
    double Dsv;       tfriend -> Branch("sv_deltaPoS",   &Dsv);
    //cp_Distances between first detections, use true cp_info
    double Dxcp;      tfriend -> Branch("cp_deltaVtxX",  &Dxcp);
    double Dycp;      tfriend -> Branch("cp_deltaVtxY",  &Dycp);
    double Dzcp;      tfriend -> Branch("cp_deltaVtxZ",  &Dzcp);
    double Dtcp;      tfriend -> Branch("cp_deltaVtxXY", &Dtcp);
    double Dcp;       tfriend -> Branch("cp_deltaVtx",   &Dcp);
    double Detacp;    tfriend -> Branch("cp_deltaVtxEta",&Detacp);
    double Dphicp;    tfriend -> Branch("cp_deltaVtxPhi",&Dphicp);
    double Dthetacp; 



    //3. //////////////////////////////////////////////////////////////
    // LOOP & Update Branches
    ///////////////////////////////////////////////////////////////////
    for (int i = 0; i<Nentries; i++)
    { 
        //icheck = 1;
        //Prepare for next cycle
        cout<<"Entry "<<i<<endl;
        tree   ->GetEntry(i);
        tfriendC->GetEntry(i);
        //tfriend->GetEntry(i);
        int N = x->size();

        // ts___________________________________________________________
        //Length and distance between first layers (from ts)
        int first0 = layerfirst->at(0); int last0 = layerlast->at(0);
        int first1 = layerfirst->at(1); int last1 = layerlast->at(1);
        Dl = abs (first1 - first0);
        //cout<<"Dl = "<<Dl<<endl;
        //if (i == 0) {cout<<"["<<layerfirst->at(0)<<","<<layerfirst->at(1)<<"]"<<endl;}

        // sv___________________________________________________________
        // Distance between first two SimVertices                                                                                                                                                                
        Dxsv = abs(posX -> at(0) - posX -> at(1));                                                                                                                                                         
        Dysv = abs(posY -> at(0) - posY -> at(1));                                                                                                                                                         
        Dzsv = abs(posZ -> at(0) - posZ -> at(1));                                                                                                                                                          
        Dtsv = sqrt(Dxsv * Dxsv + Dysv * Dysv);
        Dsv  = sqrt(Dtsv * Dtsv + Dzsv * Dzsv); 
        //cout<<"Dxsv = "<<Dxsv<<endl;
        if (i==icheck) {cout<<"found sv stuff"<<endl;}

        // cp___________________________________________________________
        // Distance between first two CPVertices                                                                                                                                                                
        Dxcp   = abs(vtxX -> at(0) - vtxX -> at(1));                                                                                                                                                         
        Dycp   = abs(vtxY -> at(0) - vtxY -> at(1));                                                                                                                                                         
        Dzcp   = abs(vtxZ -> at(0) - vtxZ -> at(1));                                                                                                                                                          
        Dtcp   = sqrt(Dxcp * Dxcp + Dycp * Dycp);
        Dcp    = sqrt(Dtcp * Dtcp + Dzcp * Dzcp); 
        Detacp = abs(cpeta -> at(0) - cpeta -> at(1)); 
        Dphicp = abs(cpphi -> at(0) - cpphi -> at(1)); 
        if (i==icheck) {cout<<"found cp stuff"<<endl;}
        
        // lc___________________________________________________________
        // Find first_q_i: list of q positions where photon number is i
        // and layer number is firsti, i.e the first one with detection
        vector <int> first_g_0 = getAwhereB(*which_g, *layerN, first0);
        vector <int> first_g_1 = getAwhereB(*which_g, *layerN, first1);
        //Prepare lc occurrence
        vector <double> x_0   = getAwhereB(*x,        *layerN, first0);
        vector <double> y_0   = getAwhereB(*y,        *layerN, first0);
        vector <double> z_0   = getAwhereB(*z,        *layerN, first0);
        vector <double> phi_0 = getAwhereB(*phi,      *layerN, first0);
        vector <double> eta_0 = getAwhereB(*eta,      *layerN, first0);
        vector <double> x_1   = getAwhereB(*x,        *layerN, first1);
        vector <double> y_1   = getAwhereB(*y,        *layerN, first1);
        vector <double> z_1   = getAwhereB(*z,        *layerN, first1);
        vector <double> phi_1 = getAwhereB(*phi,      *layerN, first1);
        vector <double> eta_1 = getAwhereB(*eta,      *layerN, first1);
        if (i==icheck) {cout<<"lc prepared"<<endl;}
        //Prepare lc seed occurrence
        vector <double> xs_0   = getAwhereB(*xseed,   *layerN, first0);
        vector <double> ys_0   = getAwhereB(*yseed,   *layerN, first0);
        vector <double> zs_0   = getAwhereB(*zseed,   *layerN, first0);
        vector <double> phis_0 = getAwhereB(*phiseed, *layerN, first0);
        vector <double> etas_0 = getAwhereB(*etaseed, *layerN, first0);
        vector <double> xs_1   = getAwhereB(*xseed,   *layerN, first1);
        vector <double> ys_1   = getAwhereB(*yseed,   *layerN, first1);
        vector <double> zs_1   = getAwhereB(*zseed,   *layerN, first1);
        vector <double> phis_1 = getAwhereB(*phiseed, *layerN, first1);
        vector <double> etas_1 = getAwhereB(*etaseed, *layerN, first1);
        if (i==icheck) {cout<<"Seed prepared"<<endl;}
        // Get what we are interested about - lc occurrence
        vector <double> first_x_0   = getAwhereB(x_0,   first_g_0, 0);
        vector <double> first_y_0   = getAwhereB(y_0,   first_g_0, 0);
        vector <double> first_z_0   = getAwhereB(z_0,   first_g_0, 0);
        vector <double> first_phi_0 = getAwhereB(phi_0, first_g_0, 0);
        vector <double> first_eta_0 = getAwhereB(eta_0, first_g_0, 0);
        vector <double> first_x_1   = getAwhereB(x_1,   first_g_1, 1);
        vector <double> first_y_1   = getAwhereB(y_1,   first_g_1, 1);
        vector <double> first_z_1   = getAwhereB(z_1,   first_g_1, 1);
        vector <double> first_phi_1 = getAwhereB(phi_1, first_g_1, 1);
        vector <double> first_eta_1 = getAwhereB(eta_1, first_g_1, 1);
        if (i==icheck) {cout<<"lc prepared 2"<<endl;}
        // Get what we are interested about - lc occurrence
        vector <double> first_xs_0   = getAwhereB(xs_0,   first_g_0, 0);
        vector <double> first_ys_0   = getAwhereB(ys_0,   first_g_0, 0);
        vector <double> first_zs_0   = getAwhereB(zs_0,   first_g_0, 0);
        vector <double> first_phis_0 = getAwhereB(phis_0, first_g_0, 0);
        vector <double> first_etas_0 = getAwhereB(etas_0, first_g_0, 0);
        vector <double> first_xs_1   = getAwhereB(xs_1,   first_g_1, 1);
        vector <double> first_ys_1   = getAwhereB(ys_1,   first_g_1, 1);
        vector <double> first_zs_1   = getAwhereB(zs_1,   first_g_1, 1);
        vector <double> first_phis_1 = getAwhereB(phis_1, first_g_1, 1);
        vector <double> first_etas_1 = getAwhereB(etas_1, first_g_1, 1);
        if (i==icheck) {cout<<"seed prepared 2"<<endl;}

        // Use them to get lc_Distances. 
        // NOTE: as there can be more than just one detection per first layer,
        // take smallest difference.
        //Occurrences
        Dx   = minAbsDiff(first_x_0, first_x_1);
        Dy   = minAbsDiff(first_y_0, first_y_1);
        Dz   = minAbsDiff(first_z_0, first_z_1);
        Dt   = sqrt(Dx*Dx + Dy*Dy);
        D    = sqrt(Dt*Dt + Dz*Dz);
        Deta = minAbsDiff(first_eta_0, first_eta_1);
        Dphi = minAbsDiff(first_phi_0, first_phi_1);
        //Seed
        Dxseed   = minAbsDiff(first_xs_0, first_xs_1);
        Dyseed   = minAbsDiff(first_ys_0, first_ys_1);
        Dzseed   = minAbsDiff(first_zs_0, first_zs_1);
        Dtseed   = sqrt(Dxseed*Dxseed + Dyseed*Dyseed);
        Dseed    = sqrt(Dtseed*Dtseed + Dzseed*Dzseed);
        Detaseed = minAbsDiff(first_etas_0, first_etas_1);
        Dphiseed = minAbsDiff(first_phis_0, first_phis_1);
        if (i==icheck) {cout<<"Found lc stuff"<<endl;}


        for (int j = 0; j<N; j++)
        {
            int layer_ij = layerN->at(j);
            int which_gj = which_g->at(j);
            int layer_ij_start = (which_gj== 0) ? 
                                layer_ij-first0 : layer_ij-first1;
            double x_ij  = x->at(j);
            double y_ij  = y->at(j);
            double z_ij  = z->at(j);
            double t_ij = 2 * atan (exp(- eta->at(j)));
            double x_pr_ij = x_ij / z_ij;
            double y_pr_ij = y_ij / z_ij;
            if (i==icheck && j==0) {cout <<"doubles work"<<endl;}

            //theta.push_back( t_ij );
            //x_pr.push_back ( x_pr_ij );
            //y_pr.push_back ( y_pr_ij );
            //if (i==icheck && j==0) {cout<<"push_back works"<<endl;}
        } 
        tfriend -> Fill();
        if (i==icheck) {cout<<"Fill works"<<endl;}
    }
    cout<<"Checkpoint: End of Loop"<<endl;
    tfriend -> Write();
    cout<<"Write works"<<endl;
    fr->Close();
    friend_name = (treename + "_Deltas").c_str();
    cout<<"Before Adding Friend, friend_name = "<<friend_name<<endl;
    //tree -> AddFriend(friend_name, "treeFriends.root");
    //cout<<"Friend added, done!"<<endl;
    file->Close();

    //tfriend -> Write("", TObject::kOverwrite);
}

//Garbage
    /*
    vector <double> *x_pr    = 0;
    vector <double> *y_pr    = 0;
    vector <double> *theta   = 0;
    double* Dy   = 0;
    double* Dz   = 0;
    double* Dt   = 0; //sqrt(Dx**2 + Dy**2)
    double* D    = 0;
    int* Dl   = 0;
    double* Deta = 0;
    double* Dphi = 0;
    double* Dtheta = 0;
    */


    //use DataFrame to only consider info at starting layer (filter) 
    //per which_g. Then do stuff
    /*
    ROOT::RDataFrame dtf = ROOT::RDataFrame(*tree);
    auto cut0 = Form("event == %i && lc_layer == %i && lc_TSidx == 0", (i+1, first0));
    auto cut1 = Form("event == %i && lc_layer == %i && lc_TSidx == 1", (i+1, first1));
    auto filtered0 = dtf.Filter(cut0);
    auto filtered1 = dtf.Filter(cut1);
    cout <<"hello there";
    auto d0 = filtered0.Display("lc_x"); d0 -> Print();
    auto d1 = filtered1.Display("lc_x"); d1 -> Print();
    auto filtered0x = filtered0.Take<vector<double>>("lc_x");
    auto filtered0y = filtered0.Take<vector<double>>("lc_y");
    */

   /*
   tfriend -> SetBranchAddress("lc_x_pr",         &x_pr);
    tfriend -> SetBranchAddress("lc_y_pr",         &y_pr);
    tfriend -> SetBranchAddress("lc_theta",        &theta);
    tfriend -> SetBranchAddress("ts_showerLength", &length);
    tfriend -> SetBranchAddress("ts_deltaStartLayer", &Dl);
    tfriend -> SetBranchAddress("ts_showerLength",    &length);
    tfriend -> SetBranchAddress("lc_x_pr",          &x_pr);
    tfriend -> SetBranchAddress("lc_y_pr",          &y_pr);
    tfriend -> SetBranchAddress("lc_theta",         &theta);
    tfriend -> SetBranchAddress("lc_deltaStartX",   &Dx);
    tfriend -> SetBranchAddress("lc_deltaStartY",   &Dy);
    tfriend -> SetBranchAddress("lc_deltaStartZ",   &Dz);
    tfriend -> SetBranchAddress("lc_deltaStartXY",  &Dt);
    tfriend -> SetBranchAddress("lc_deltaStart",    &D);
    tfriend -> SetBranchAddress("lc_deltaStartEta", &Deta);
    tfriend -> SetBranchAddress("lc_deltaStartPhi", &Dphi);
    //tfriend -> SetBranchAddress("lc_deltaStartTheta", &Dtheta);
    tfriend -> SetBranchAddress("sv_deltaPosX", &Dxsv);
    //cout << "Dxsv = "<<Dxsv <<",   Dxsv = "<<Dxsv<<endl; 
    tfriend -> SetBranchAddress("sv_deltaPosY", &Dysv);
    tfriend -> SetBranchAddress("sv_deltaPosZ", &Dzsv);
    tfriend -> SetBranchAddress("sv_deltaPosXY",&Dtsv);
    tfriend -> SetBranchAddress("sv_deltaPoS",  &Dsv);
    tfriend -> SetBranchAddress("cp_deltaVtxX",     &Dxcp);
    tfriend -> SetBranchAddress("cp_deltaVtxY",     &Dycp);
    tfriend -> SetBranchAddress("cp_deltaVtxZ",     &Dzcp);
    tfriend -> SetBranchAddress("cp_deltaVtxXY",    &Dtcp);
    tfriend -> SetBranchAddress("cp_deltaVtx",      &Dcp);
    tfriend -> SetBranchAddress("cp_deltaVtxEta",   &Dcp);
    tfriend -> SetBranchAddress("cp_deltaVtxPhi",   &Dcp);
    */