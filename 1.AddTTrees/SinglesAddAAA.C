/// \author Stefano Veroni
/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 1.AddTTrees
root SinglesAddAAA.C -q
root SinglesAddEnergies.C -q
root SinglesAddShape.C -q
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
#include "TFile.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"
using namespace std;

/*
 EXPLICIT SOME RELATIVELY STRAIGHTFORWARD INFO
 - POSITION RELATED (coordinates friend tree)
    - Theta, Rho, R
    - Projected x, y, rho
    - Seed positions: x, y, z, theta, rho, R, x_pr, y_pr, rho_pr
 - ENERGY RELATED (energy friend tree, placeholder)
    - Energy corrected with tsMult
    - Inverse of tsMult
*/
void SinglesAddAAA () 
{

    string path    = "../myTTrees/TreesEn90to110/Single/"; //none if none
    string rootname= "step3ticl_eta21_FlatTracksters";//root file
    string dirname = "none"; //"none" if no directory in between
    string treename= "TSTree_SimFromCP";
    int icheck = 5; //Entry to check in loop to know where errors lie.


    //1. //////////////////////////////////////////////////////////////
    //Get stuff from the TTree we already have
    ///////////////////////////////////////////////////////////////////
    cout<<"startmine"<<endl;
    TDirectory* myTTrees = new TDirectory();
    //TDirectory->cd();
    TFile *file = open_file(path, rootname);
    TTree *tree = get_tree (file, dirname, treename);
    //tree ->Scan("lc_TSidx:lc_x:lc_layer:lc_z:lc_seedEnergy");
    //return;
    int Nentries = tree->GetEntries();
    vector <double>* x    = 0;
    vector <double>* y    = 0;
    vector <double>* z    = 0;
    vector <double>* eta  = 0;
    vector <double>* phi = 0;
    vector <int>*    which_g = 0;
    vector <double>* energy = 0;
    vector <double>* energy_seed = 0;
    vector <double>* phi_seed = 0;
    vector <double>* eta_seed = 0;
    vector <double>* tsmult = 0;
    tree->SetBranchAddress("lc_x",        &x);
    tree->SetBranchAddress("lc_y",        &y);
    tree->SetBranchAddress("lc_z",        &z);
    tree->SetBranchAddress("lc_phi",      &phi);
    tree->SetBranchAddress("lc_eta",      &eta);
    tree->SetBranchAddress("lc_TSidx",    &which_g);
    tree->SetBranchAddress("lc_energy",   &energy);
    tree->SetBranchAddress("lc_seedEnergy",&energy_seed);
    tree->SetBranchAddress("lc_seedEta",   &eta_seed);
    tree->SetBranchAddress("lc_seedPhi",   &phi_seed);
    tree->SetBranchAddress("lc_tsMult",    &tsmult);
    //tree->Scan("event:lc_TSidx:lc_x:lc_layer:sv_posX:cp_vtxX:cp_phi");


    //2. //////////////////////////////////////////////////////////////
    //Prepare new TTree, which will later become friend
    ///////////////////////////////////////////////////////////////////
    cout<<"start friend"<<endl;
    TFile fr ((path+"treeFriends.root").c_str(), "recreate");
    const char* friend_name = (treename + "_Coords").c_str();
    TTree *tfriend = new TTree(friend_name, friend_name);
    const char* friendE_name = (treename + "_EnergyPlaceholder").c_str();
    TTree *tfriendE = new TTree(friendE_name, friendE_name);
    //"New" Coordinates
    vector <double> rho  ; tfriend -> Branch("lc_rho",  &rho);
    vector <double> R    ; tfriend -> Branch("lc_R",    &R);
    vector <double> theta; tfriend -> Branch("lc_theta",&theta); 
    //Projected Coordinates
    vector <double> x_pr  ; tfriend -> Branch("lc_x_pr",  &x_pr);
    vector <double> y_pr  ; tfriend -> Branch("lc_y_pr",  &y_pr);
    vector <double> rho_pr; tfriend -> Branch("lc_rho_pr",&y_pr);
    //Seed Coordinates
    vector <double> x_seed     ; tfriend -> Branch("lc_seedX",     &x_seed);
    vector <double> y_seed     ; tfriend -> Branch("lc_seedY",     &y_seed);
    vector <double> z_seed     ; tfriend -> Branch("lc_seedZ",     &z_seed);
    vector <double> rho_seed   ; tfriend -> Branch("lc_seedRho",   &rho_seed);
    vector <double> R_seed     ; tfriend -> Branch("lc_seedR",     &R_seed);
    vector <double> theta_seed ; tfriend -> Branch("lc_seedTheta", &theta_seed);
    vector <double> x_pr_seed  ; tfriend -> Branch("lc_seedX_pr",  &x_pr_seed);
    vector <double> y_pr_seed  ; tfriend -> Branch("lc_seedY_pr",  &y_pr_seed);
    vector <double> rho_pr_seed; tfriend -> Branch("lc_seedRho_pr",&rho_pr_seed);
    // Energy Stuff
    int Nlcs                   ; tfriendE-> Branch("lc_Nlcs",      &Nlcs);
    int Nshared                ; tfriendE-> Branch("lc_Nshared",   &Nshared);
    vector <double> true_energy; tfriendE-> Branch("lc_energyTrue",&true_energy);
    vector <double> tsfraction ; tfriendE-> Branch("lc_tsFraction",&tsfraction);
    vector<double>lcs_distinct ; tfriendE-> Branch("lc_distinctEs",&lcs_distinct);
    vector<vector<double>> lc_tsfracts;
                                 tfriendE-> Branch("lc_tsSameFractions",&lc_tsfracts);
    vector<int>lcs_sharingflags; tfriendE-> Branch("lc_shareFlags",&lcs_sharingflags);
    
    //3. //////////////////////////////////////////////////////////////
    // LOOP & Update Branches
    ///////////////////////////////////////////////////////////////////
    for (int i = 0; i<Nentries; i++)
    { 
        //icheck = i*1;
        //Prepare for next cycle
        tree->GetEntry(i);
        cout<<"Entry "<<i<<endl;
        int N = x->size();

        rho.clear()      ; R.clear()        ; theta.clear();
        x_pr.clear()     ; y_pr.clear()     ; rho_pr.clear();
        x_seed.clear()   ; y_seed.clear()   ; z_seed.clear();
        rho_seed.clear() ; R_seed.clear()   ; theta_seed.clear();
        x_pr_seed.clear(); y_pr_seed.clear(); rho_pr_seed.clear();

        true_energy.clear(); tsfraction.clear(); 
        Nlcs = 0;            Nshared = 0;
        lcs_distinct.clear();lc_tsfracts.clear(); lcs_sharingflags.clear();
        vector <int> kdoubles;
        

        if (i==icheck) {cout <<"Energy addresses set"<<endl;}

        for (int j = 0; j<N; j++)
        {
            //Coordinates (new, projected, seed)
            //_________________________________________________________
            double x_ij  = x->at(j);
            double y_ij  = y->at(j);
            double z_ij  = z->at(j);

            double rho_ij = sqrt(x_ij*x_ij + y_ij*y_ij);
            double R_ij   = sqrt(rho_ij*rho_ij + z_ij*z_ij);
            double t_ij   = 2 * atan (exp(- eta->at(j)));
            double x_pr_ij   = x_ij / z_ij;
            double y_pr_ij   = y_ij / z_ij;
            double rho_pr_ij = sqrt(x_pr_ij*x_pr_ij + y_pr_ij*y_pr_ij);

            double p_seed_ij   = phi_seed->at(j);
            double t_seed_ij   = 2 * atan (exp(- eta_seed->at(j)));
            double R_seed_ij   = z_ij / cos (t_seed_ij); //z_seed = lc_z
            double rho_seed_ij = z_ij * tan (t_seed_ij);
            double x_seed_ij   = rho_seed_ij * cos(p_seed_ij);
            double y_seed_ij   = rho_seed_ij * sin(p_seed_ij);
            double x_pr_seed_ij   = x_seed_ij   / z_ij;
            double y_pr_seed_ij   = y_seed_ij   / z_ij;
            double rho_pr_seed_ij = rho_seed_ij / z_ij;
            if (i==icheck && j==0) {cout <<"doubles work"<<endl;}

            rho  .push_back( rho_ij );
            R    .push_back( R_ij   );
            theta.push_back( t_ij   );
            x_pr  .push_back( x_pr_ij  );
            y_pr  .push_back( y_pr_ij  );
            rho_pr.push_back( rho_pr_ij);
            x_seed    .push_back( x_seed_ij  );
            y_seed    .push_back( y_seed_ij  );
            z_seed    .push_back( z_ij       );
            rho_seed  .push_back( rho_seed_ij);
            R_seed    .push_back( R_seed_ij  );
            theta_seed.push_back( t_seed_ij  );
            x_pr_seed  .push_back( x_pr_seed_ij  );
            y_pr_seed  .push_back( y_pr_seed_ij  );
            rho_pr_seed.push_back( rho_pr_seed_ij);
            if (i==icheck && j==0) {cout<<"push_back works"<<endl;}


            //Energy Stuff ///////////////////////////////////////////
            //________________________________________________________
            double Ej = energy->at(j);
            double Mj = tsmult->at(j);
            true_energy.push_back(Ej/Mj);
            tsfraction .push_back(1/Mj);

            if (contains(kdoubles, j)) {continue;} //nothing more to do
            Nlcs += 1;
            lcs_distinct.push_back(Ej);
            vector<double> tsfract_j (2);
            tsfract_j[which_g->at(j)] = 1/Mj;
            int flag = which_g->at(j);
            for (int k = j+1; k<N; k++)
            {
                if (x_ij==x->at(k) && y_ij==y->at(k) && z_ij==z->at(k))
                {
                    kdoubles.push_back(k);
                    double Mk = tsmult->at(k); 

                    Nshared += 1;
                    tsfract_j[which_g->at(j)== 0] = 1/Mk;
                    flag = 2;
                }
            }
            lc_tsfracts.push_back(tsfract_j);
            lcs_sharingflags.push_back(flag);
        } 

        tfriend -> Fill();
        tfriendE-> Fill();
        if (i==icheck) {cout<<"Fill works"<<endl;}
    }
    cout<<"Checkpoint: End of Loop"<<endl;
    tfriend -> Write();
    tfriendE-> Write();
    cout<<"Write works"<<endl;
    fr.Close();
    //friend_name = (treename + "_Coords").c_str();
    //cout<<"Before Adding Friend, friend_name = "<<friend_name<<endl;
    //tree -> AddFriend(friend_name, "treeFriends.root");
    //cout<<"Friend added, done!"<<endl;
    file->Close();
}

//Garbage
    /*
    TFile *fr = new TFile("treeFriends.root", "update");
    const char *friend_cycles = (treename +"_friend;*").c_str();
    fr -> Delete(friend_cycles);
    fr -> Close();
    //TFile fr ("treeFriends.root", "update");
    //const char* friend_name = (treename + "_friend").c_str();
    //TTree *tfriend = new TTree(friend_name, friend_name);

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