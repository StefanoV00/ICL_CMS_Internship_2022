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
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"
using namespace std;

/*
NEW TTREE W. INFO RELATED TO ENERGY DISTRIBUTIONS - CORRECTED WITH TSMULT
 - Energy deposited per layer per photon
 - Energy deposited per layer per photon since start of own shower

 - Layer # with maximum energy deposited by photon
 - Layer # since start of own shower with maximum energy by that photon

 - Barycentre positions: x, y, theta, phi, eta
   (the positions are 1001, 1001, 11, 11, 11 if no event in layer)

ALSO A SHAPE PLACEHOLDER
 - length of shower and stuff like that
*/
void SinglesAddEnergies () 
{
    string path    = "../myTTrees/TreesEn90to110/Single/"; //none if none
    string rootname= "step3ticl_eta21_FlatTracksters";//root file
    string friendrootname= "treeFriends"; //root file with friends
    string dirname = "none";
    string treename= "TSTree_SimFromCP";

    int icheck = 0; //Entry to check in loop to know where errors lie.


    ///////////////////////////////////////////////////////////////////
    //1. //////////////////////////////////////////////////////////////
    //Get stuff from the TTree we already have
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    TFile* file = open_file(path, rootname);
    TTree* tree = get_tree (file, dirname, treename);
    //tree->Print();
    //tree ->Scan("lc_TSidx:lc_nrechits:lc_layer:lc_mult:lc_tsMult:lc_seedEnergy");
    //return;
    int Nentries = tree->GetEntries();
    vector <double>* x    = 0;
    vector <double>* y    = 0;
    vector <double>* z    = 0;
    vector <double>* energyoverlapped = 0; //get the true one later
    vector <int>* nrechits = 0;
    vector <int>* layerN  = 0;
    vector <int>* layerfirst = 0;
    vector <int>* layerlast  = 0;
    tree->SetBranchAddress("lc_x",        &x);
    tree->SetBranchAddress("lc_y",        &y);
    tree->SetBranchAddress("lc_z",        &z);
    tree->SetBranchAddress("lc_energy", &energyoverlapped);
    tree->SetBranchAddress("lc_nrechits", &nrechits);
    tree->SetBranchAddress("lc_layer",    &layerN);
    tree->SetBranchAddress("ts_firstLayer", &layerfirst);
    tree->SetBranchAddress("ts_lastLayer",  &layerlast);
    //tree->Scan("event:lc_TSidx:lc_x:lc_layer:sv_posX:cp_vtxX:cp_phi");


    ///////////////////////////////////////////////////////////////////
    //2. //////////////////////////////////////////////////////////////
    //Prepare new TTree, which will later become friend
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    TFile* fr = open_file(path, friendrootname, "update");
    
    // Take some previous energy info (from the energy placeholder)
    const char* energies_name = (treename + "_EnergyPlaceholder").c_str();
    TTree* tfriende = get_tree(fr, "none", energies_name);
    vector <double>* energy = 0;
    tfriende -> SetBranchAddress("lc_energyTrue", &energy);
    
    // Create new ttrees: real energy one and shape placeholder
    TTree* tfriend =  tfriende -> CloneTree(0);
    const char* shapeholder_name = (treename + "_ShapePlaceholder").c_str();
    TTree* tfriendS= new TTree(shapeholder_name, shapeholder_name);
    //tfriend ->Scan("lc_energyTrue");

    // Take projection coordinates from Coords TTree
    const char* coords_name = (treename + "_Coords").c_str();
    TTree* tfriendC = get_tree(fr, "none", coords_name);
    vector <double>* xpr = 0; 
    tfriendC -> SetBranchAddress("lc_x_pr", &xpr);
    vector <double>* ypr = 0; 
    tfriendC -> SetBranchAddress("lc_y_pr", &ypr);
    
    // Shower info using ts info
    vector<int> length   ;tfriendS->Branch("ts_showerLength", &length);
    int whole_length = 0 ;tfriendS->Branch("ts_wholeLength",  &whole_length);
    
    // Effective Layer numbers, use lc info
    vector<int>layer_start ;tfriendS->Branch("lc_layerNFromStart",&layer_start);
    
    //Energy info, use lc
    vector <double> totEd;
    tfriend -> Branch("lc_EPerLayer",          &totEd);
    vector <double> totEd_start; 
    tfriend -> Branch("lc_EPerLayerFromStart", &totEd_start);  
    int LmaxETot = 0      ; 
    tfriend -> Branch("lc_layerMaxE",          &LmaxETot);
    int LmaxETot_start = 0; 
    tfriend -> Branch("lc_layerMaxEFromStart", &LmaxETot_start);
    double maxlcfraction = 0; //also overlapped
    tfriend -> Branch("lc_maxLCFraction",     &maxlcfraction);
    double maxLfraction = 0;
    tfriend -> Branch("lc_maxLayerFraction",  &maxLfraction);


    //Barycentre of entire shower
    double BCx = 0; tfriend->Branch("lc_BCx", &BCx);
    double BCy = 0; tfriend->Branch("lc_BCy", &BCy);
    double BCz = 0; tfriend->Branch("lc_BCz", &BCz);
    double BCl = 0; tfriend->Branch("lc_BCl", &BCl);
    double BCzstart = 0; tfriend->Branch("lc_BCzFromStart", &BCzstart);
    double BClstart = 0; tfriend->Branch("lc_BClFromStart", &BClstart);

    //Barycentre of projections
    double centreX_pr; tfriend->Branch("lc_BCx_pr", &centreX_pr);
    double centreY_pr; tfriend->Branch("lc_BCy_pr", &centreY_pr);
    
    //Barycentre of layer Information
    vector <double> x_totcentre  ;
    tfriend -> Branch("lc_centreX",     &x_totcentre);
    vector <double> y_totcentre  ;
    tfriend -> Branch("lc_centreY",     &y_totcentre);
    vector <double> eta_totcentre;
    tfriend -> Branch("lc_centreEta",   &eta_totcentre);
    vector <double> phi_totcentre;
    tfriend -> Branch("lc_centrePhi",   &phi_totcentre);
    
    //Barycentre of layer from start
    vector <double> x_totcentre_start  ;
    tfriend -> Branch("lc_centreXFromStart",    &x_totcentre_start);
    vector <double> y_totcentre_start  ;
    tfriend -> Branch("lc_centreYFromStart",    &y_totcentre_start);
    vector <double> eta_totcentre_start;
    tfriend -> Branch("lc_centreEtaFromStart",  &eta_totcentre_start);
    vector <double> phi_totcentre_start;
    tfriend -> Branch("lc_centrePhiFromStart",  &phi_totcentre_start);



    ///////////////////////////////////////////////////////////////////
    //3. //////////////////////////////////////////////////////////////
    // LOOP & Update Branches
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    for (int i = 0; i<Nentries; i++)
    { 
        //icheck = i*1;
        //Prepare for next cycle
        cout<<"Entry "<<i<<endl;
        tree     -> GetEntry(i);
        tfriendC -> GetEntry(i);
        tfriende -> GetEntry(i);
        layer_start .clear();
        int N = energy->size();
                
        // ts__________________________________________________________
        //Length and distance between first layers (from ts)
        int first_all = layerfirst->at(0);
        int last_all = layerlast->at(0);
        whole_length = last_all - first_all + 1;
        if (i==icheck) {cout <<"Lenght Info Set"<<endl
                            << "last_all   = "<<last_all<<endl;} 

        // lc : Barycentre of ALL shower & Projected
        //______________________________________________________________
        
        // Tot
        BCx = mymean(*x     , *energy);
        BCy = mymean(*y     , *energy);
        BCz = mymean(*z     , *energy);
        BCl = mymean(*layerN, *energy);
        BCzstart = BCz - z->at(argmin(*z));
        BClstart = BCl - first_all;
        centreX_pr = mymean(*xpr, *energy);
        centreY_pr = mymean(*ypr, *energy);
        if (i==icheck) {cout<<"Barycentre works"<<endl;}


        // lc: Energy-deposited vectors with number of entries
        // _____________________________________________________________
        double Etot = myvecsum(*energyoverlapped);
        double lcEmax = energyoverlapped->at(argmax(*energyoverlapped));
        maxlcfraction = lcEmax / Etot;
        ///*
        vector <double> Edi       (last_all);
        vector <double> Edi_start (whole_length);
        vector <vector<double>> Ed          (2, Edi); 
        vector <double>         totEd       (last_all); 
        vector <vector<double>> Ed_start    (2, Edi_start);
        vector <double>         totEd_start (whole_length);
        vector <vector<double>>* Ed_ptr          = &Ed; //0    
        vector <double>        * totEd_ptr       = &totEd;
        vector <vector<double>>* Ed_start_ptr    = &Ed_start;
        vector <double>        * totEd_start_ptr = &totEd_start; 
        tfriend -> SetBranchAddress("lc_EPerLayer",         &totEd_ptr);  
        tfriend -> SetBranchAddress("lc_EPerLayerFromStart",&totEd_start_ptr);
        if (i==icheck) {cout <<"Barycentre addresses set"<<endl;}


        // lc: Other Barycentre Information
        //______________________________________________________________
        // Initialise Barycentre vectors with number of entries
        ///*
        vector <double> x_totcentre     (last_all);
        vector <double> y_totcentre     (last_all);
        vector <double> z_totcentre     (last_all);
        vector <double> theta_totcentre (last_all);
        vector <double> eta_totcentre   (last_all);
        vector <double> phi_totcentre   (last_all);
        vector <vector<double>> x_centre      (2, x_totcentre);
        vector <vector<double>> y_centre      (2, y_totcentre);
        vector <vector<double>> z_centre      (2, z_totcentre);
        vector <vector<double>> theta_centre  (2, theta_totcentre);
        vector <vector<double>> eta_centre    (2, eta_totcentre);
        vector <vector<double>> phi_centre    (2, phi_totcentre); 
        vector <vector<double>>* xc_ptr     = &x_centre; 
        vector <vector<double>>* yc_ptr     = &y_centre; 
        vector <vector<double>>* zc_ptr     = &z_centre; 
        vector <vector<double>>* thetac_ptr = &theta_centre; 
        vector <vector<double>>* etac_ptr   = &eta_centre;
        vector <vector<double>>* phic_ptr   = &phi_centre;     
        vector <double>        * totxc_ptr  = &x_totcentre;
        vector <double>        * totyc_ptr  = &y_totcentre;
        vector <double>        * totzc_ptr  = &z_totcentre;
        vector <double>        * tottc_ptr  = &theta_totcentre;
        vector <double>        * totec_ptr  = &eta_totcentre;
        vector <double>        * totpc_ptr  = &phi_totcentre;
        tfriend -> SetBranchAddress("lc_centreX",     &totxc_ptr);
        tfriend -> SetBranchAddress("lc_centreY",     &totyc_ptr);
        tfriend -> SetBranchAddress("lc_centreEta",   &totec_ptr);
        tfriend -> SetBranchAddress("lc_centrePhi",   &totpc_ptr);  
        // Now same, but from start
        vector <double> x_totc_start     (whole_length);
        vector <double> y_totc_start     (whole_length);
        vector <double> z_totc_start     (whole_length);
        vector <double> theta_totc_start (whole_length);
        vector <double> eta_totc_start   (whole_length);
        vector <double> phi_totc_start   (whole_length);
        vector <double> footprint        (whole_length);
        vector <vector<double>> xc_start     (2, footprint);
        vector <vector<double>> yc_start     (2, footprint);
        vector <vector<double>> zc_start     (2, footprint);
        vector <vector<double>> thetac_start (2, footprint);
        vector <vector<double>> etac_start   (2, footprint);
        vector <vector<double>> phic_start   (2, footprint); 
        vector <vector<double>>* xcs_ptr     = &xc_start; 
        vector <vector<double>>* ycs_ptr     = &yc_start; 
        vector <vector<double>>* zcs_ptr     = &zc_start; 
        vector <vector<double>>* thetacs_ptr = &thetac_start; 
        vector <vector<double>>* etacs_ptr   = &etac_start;
        vector <vector<double>>* phics_ptr   = &phic_start;     
        vector <double>        * totxcs_ptr  = &x_totc_start;
        vector <double>        * totycs_ptr  = &y_totc_start;
        vector <double>        * totzcs_ptr  = &z_totc_start;
        vector <double>        * tottcs_ptr  = &theta_totc_start;
        vector <double>        * totecs_ptr  = &eta_totc_start;
        vector <double>        * totpcs_ptr  = &phi_totc_start;
        tfriend -> SetBranchAddress("lc_centreXFromStart",    &totxcs_ptr);
        tfriend -> SetBranchAddress("lc_centreYFromStart",    &totycs_ptr);
        tfriend -> SetBranchAddress("lc_centreEtaFromStart",  &totecs_ptr);
        tfriend -> SetBranchAddress("lc_centrePhiFromStart",  &totpcs_ptr); 
        //*/
        if (i==icheck) {cout <<"Energy addresses set"<<endl;}



        for (int j = 0; j<N; j++)
        {
            int layer_ij = layerN ->at(j);
            int layer_ij_start_all = layer_ij - first_all;
            layer_start .push_back(layer_ij_start_all);
            double E_ij  = energy->at(j);
            
            totEd_ptr->at(layer_ij-1)                  += E_ij;
            totEd_start_ptr->at(layer_ij_start_all)    += E_ij;
            if (i==icheck && j==0) {cout<<"Energy adding works"<<endl;}

            totxc_ptr -> at (layer_ij-1)      += x->at(j) * E_ij;
            totyc_ptr -> at (layer_ij-1)      += y->at(j) * E_ij;
            totzc_ptr -> at (layer_ij-1)      += z->at(j) * E_ij;
            totxcs_ptr -> at (layer_ij_start_all)   += x->at(j) * E_ij;
            totycs_ptr -> at (layer_ij_start_all)   += y->at(j) * E_ij;
            totzcs_ptr -> at (layer_ij_start_all)   += z->at(j) * E_ij; 
            if (i==icheck && j==0) {cout<<"Coords adding works"<<endl;} 
        } 
        double LEmax = totEd_start_ptr->at(argmax(*totEd_start_ptr));
        maxLfraction = LEmax / (Etot);


        // Finish Barycentres (normal ones)
        for (int lj=0; lj < last_all; lj++)
        {
            //Tot Ones
            double Etotij = totEd_ptr -> at(lj);
            if (Etotij == 0)
            {
                totxc_ptr -> at(lj) = 1001; //TMath::QuietNaN();
                totyc_ptr -> at(lj) = 1001; //TMath::QuietNaN();
                totzc_ptr -> at(lj) = 1001; //TMath::QuietNaN();
                totpc_ptr -> at(lj) =   11; //TMath::QuietNaN();
                tottc_ptr -> at(lj) =   11; //TMath::QuietNaN();
                totec_ptr -> at(lj) =   11; //TMath::QuietNaN();
            }
            else
            {
                totxc_ptr -> at(lj) /= Etotij;
                totyc_ptr -> at(lj) /= Etotij;
                totzc_ptr -> at(lj) /= Etotij;
 
                double xij = totxc_ptr -> at(lj);
                double yij = totyc_ptr -> at(lj);
                double zij = totzc_ptr -> at(lj);
                double rhoij = sqrt(xij*xij + yij*yij);
                double phiij = acos( xij / rhoij );
                double thetaij = atan (rhoij / zij);
                double etaij = - log ( tan(thetaij / 2) );
                totpc_ptr->at(lj) = phiij;
                tottc_ptr->at(lj) = thetaij;
                totec_ptr->at(lj) = etaij;
            }   
        }
        if (i==icheck) {cout<<"Barycentre Standard works"<<endl;} 
        

        // Finish Barycentres From Start
        for (int ljs=0; ljs < whole_length; ljs++)
        {           
            // Total Ones 
            double Etotsij = totEd_start_ptr -> at(ljs);
            if (Etotsij == 0)
            {
                totxcs_ptr -> at(ljs) = 1001; //TMath::QuietNaN();
                totycs_ptr -> at(ljs) = 1001; //TMath::QuietNaN();
                totzcs_ptr -> at(ljs) = 1001; //TMath::QuietNaN();
                totpcs_ptr -> at(ljs) =   11; //TMath::QuietNaN();
                tottcs_ptr -> at(ljs) =   11; //TMath::QuietNaN();
                totecs_ptr -> at(ljs) =   11; //TMath::QuietNaN();
            }
            else
            {
                totxcs_ptr -> at(ljs) /= Etotsij;
                totycs_ptr -> at(ljs) /= Etotsij;
                totzcs_ptr -> at(ljs) /= Etotsij;
                double xij = totxcs_ptr -> at(ljs);
                double yij = totycs_ptr -> at(ljs);
                double zij = totzcs_ptr -> at(ljs);
                double rhoij = sqrt(xij*xij + yij*yij);
                double phiij = acos( xij / rhoij );
                double thetaij = atan (rhoij / zij);
                double etaij = - log ( tan(thetaij / 2) );
                totpcs_ptr->at(ljs) = phiij;
                tottcs_ptr->at(ljs) = thetaij;
                totecs_ptr->at(ljs) = etaij;
            }    
        }
        if (i==icheck) {cout<<"Barycentre from Start works"<<endl;} 
        
        
        //Layer Max
        LmaxETot = argmax(totEd);
        LmaxETot_start = argmax(totEd_start);
        if (i==icheck) {cout<<"Lmaxs work"<<endl;} 
        
        tfriend -> Fill();
        tfriendS-> Fill();
        if (i==icheck) {cout<<"Fill works"<<endl;}
    }
    cout<<"Checkpoint: End of Loop"<<endl;
    const char* friend_name = (treename + "_Energy").c_str();
    tfriend -> SetName(friend_name);
    tfriend -> Write();//Write("", TObject::kOverwrite);
    tfriendS-> Write();
    cout<<"Write works"<<endl;
    const char* placeholder_namecycle = (treename + "_EnergyPlaceholder;*").c_str();
    fr->Delete(placeholder_namecycle);
    cout<<"Delete works"<<endl;
    fr->Close();
    //friend_name = (treename + "_Energies").c_str();
    //cout<<"Before Adding Friend, friend_name = "<<friend_name<<endl;
    //tree -> AddFriend(friend_name, "treeFriends.root");
    //cout<<"Friend added, done!"<<endl;
    file->Close();
    //const char* friend_name = (treename + "_friend").c_str();
    //cout<<"Before Adding Friend, friend_name = "<<friend_name<<endl;
    //tree -> AddFriend(friend_name, "treeFriends.root");
    //cout<<"Friend added, done!"<<endl;
    //file->Close();
}