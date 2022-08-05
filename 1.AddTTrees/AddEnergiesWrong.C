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
 ADD TO FRIEND INFO RELATED TO ENERGY DISTRIBUTIONS - NOT CORRECTED WITH TSMULT
 - Energy deposited per layer per photon
 - Energy deposited per layer in total
 - Energy deposited per layer per photon since start of own shower
 - Energy deposited per layer in total since start of very first shower

 - Layer # with maximum energy deposited by one photon
 - Layer # with maximum energy deposited in total
 - Layer # since start of own shower with maximum energy by that photon
 - Layer # since start of first shower with maximum total energy

 - Barycentre gpositions (per photon): x, y, theta, phi, eta
 - Barycentre positions total: x, y, theta, phi, eta
   (the positions are 1001, 1001, 11, 11, 11 if no event in layer)
*/
void AddEnergiesWrong () 
{
    string path    = "../myTTrees/TreesEn90to110/Double/"; //none if none
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
    vector <double>* eta  = 0;
    vector <double>* phi = 0;
    vector <double>* energy = 0;
    vector <int>* which_g = 0;
    vector <int>* layerN  = 0;
    vector <int>* layerfirst = 0;
    vector <int>* layerlast  = 0;
    vector <int>* nrechits = 0;
    vector <int>* mult = 0;  // 100 entries: 50 for whihc_g = 0, 50 for 1.
    vector <double>* tsmult = 0;
    tree->SetBranchAddress("lc_x",        &x);
    tree->SetBranchAddress("lc_y",        &y);
    tree->SetBranchAddress("lc_z",        &z);
    tree->SetBranchAddress("lc_phi",      &phi);
    tree->SetBranchAddress("lc_eta",      &eta);
    tree->SetBranchAddress("lc_energy",&energy);
    tree->SetBranchAddress("lc_layer", &layerN);
    tree->SetBranchAddress("lc_TSidx", &which_g);
    tree->SetBranchAddress("lc_nrechits", &nrechits);
    tree->SetBranchAddress("lc_mult",  &mult);
    tree->SetBranchAddress("lc_tsMult",  &tsmult);
    tree->SetBranchAddress("ts_firstLayer", &layerfirst);
    tree->SetBranchAddress("ts_lastLayer",  &layerlast);
    //tree->Scan("event:lc_TSidx:lc_x:lc_layer:sv_posX:cp_vtxX:cp_phi");


    ///////////////////////////////////////////////////////////////////
    //2. //////////////////////////////////////////////////////////////
    //Prepare new TTree, which will later become friend
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    TFile* fr = open_file(path, friendrootname, "update");
    const char* friend_name = (treename + "_EnergiesWrong").c_str();
    TTree* tfriend = new TTree(friend_name, friend_name);

    // Take projection coordinates from Coords TTree
    const char* coords_name = (treename + "_Coords").c_str();
    TTree* tfriendC = get_tree(fr, "none", coords_name);
    vector <double>* xpr = 0; 
    tfriendC -> SetBranchAddress("lc_x_pr", &xpr);
    vector <double>* ypr = 0; 
    tfriendC -> SetBranchAddress("lc_y_pr", &ypr);
    
    // Shower info using ts info
    vector<int> length(2); tfriend -> Branch("ts_showerLength", &length);
    int whole_length = 0 ; tfriend -> Branch("ts_wholeLength",  &whole_length);
    
    // Effective Layer numbers, use ln info
    vector<int>layerg_start;tfriend->Branch("lc_layerFromgStart",&layerg_start);
    vector<int>layer_start ;tfriend->Branch("lc_layerNFromStart",&layer_start);
    
    //Energy info, use lc
    vector <vector<double>> Ed (10); // Ed[i] energy deposited in layers by photon i
    tfriend -> Branch("lc_gEPerLayer",            &Ed);  
    vector <double> totEd;
    tfriend -> Branch("lc_totEPerLayer",          &totEd);
    vector <vector<double>> Ed_start; //_start prefix means "from start of shower"
    tfriend -> Branch("lc_gEPerLayerFromStart",   &Ed_start);  
    vector <double> totEd_start; 
    tfriend -> Branch("lc_totEPerLayerFromStart", &totEd_start);  
    vector <int> LmaxE (2); 
    tfriend -> Branch("lc_layerMaxE",             &LmaxE);
    int LmaxETot = 0      ; 
    tfriend -> Branch("lc_layerMaxTotE",          &LmaxETot);
    vector <int> LmaxE_start (2); 
    tfriend -> Branch("lc_layerMaxEFromStart",    &LmaxE_start); 
    int LmaxETot_start = 0; 
    tfriend -> Branch("lc_layerMaxTotEFromStart", &LmaxETot_start);

    //Barycentre Information, use lc
    vector <vector<double>> x_centre  ;
    tfriend -> Branch("lc_gCentreX",     &x_centre);
    vector <vector<double>> y_centre  ;
    tfriend -> Branch("lc_gCentreY",     &y_centre);
    vector <vector<double>> z_centre  ;
    tfriend -> Branch("lc_gCentreZ",     &z_centre);
    vector <vector<double>> theta_centre;
    tfriend -> Branch("lc_gCentreTheta", &theta_centre);
    vector <vector<double>> eta_centre;
    tfriend -> Branch("lc_gCentreEta",   &eta_centre);
    vector <vector<double>> phi_centre;
    tfriend -> Branch("lc_gCentrePhi",   &phi_centre);
    //Barycentre Tot Information
    vector <double> x_totcentre  ;
    tfriend -> Branch("lc_centreTotX",     &x_totcentre);
    vector <double> y_totcentre  ;
    tfriend -> Branch("lc_centreTotY",     &y_totcentre);
    vector <double> z_totcentre  ;
    tfriend -> Branch("lc_centreTotZ",     &z_totcentre);
    vector <double> theta_totcentre;
    tfriend -> Branch("lc_centreTotTheta", &theta_totcentre);
    vector <double> eta_totcentre;
    tfriend -> Branch("lc_centreTotEta",   &eta_totcentre);
    vector <double> phi_totcentre;
    tfriend -> Branch("lc_centreTotPhi",   &phi_totcentre);
    
    // Barycentre from start
    vector <vector<double>> x_centre_start  ;
    tfriend -> Branch("lc_gCentreXFromStart",    &x_centre_start);
    vector <vector<double>> y_centre_start  ;
    tfriend -> Branch("lc_gCentreYFromStart",    &y_centre_start);
    vector <vector<double>> z_centre_start  ;
    tfriend -> Branch("lc_gCentreZFromStart",    &z_centre_start);
    vector <vector<double>> theta_centre_start;
    tfriend -> Branch("lc_gCentreThetaFromStart",&theta_centre_start);
    vector <vector<double>> eta_centre_start;
    tfriend -> Branch("lc_gCentreEtaFromStart",  &eta_centre_start);
    vector <vector<double>> phi_centre_start;
    tfriend -> Branch("lc_gCentrePhiFromStart",  &phi_centre_start);
    //Barycentre Tot from start
    vector <double> x_totcentre_start  ;
    tfriend -> Branch("lc_centreTotXFromStart",    &x_totcentre_start);
    vector <double> y_totcentre_start  ;
    tfriend -> Branch("lc_centreTotYFromStart",    &y_totcentre_start);
    vector <double> z_totcentre_start  ;
    tfriend -> Branch("lc_centreTotZFromStart",    &z_totcentre_start);
    vector <double> theta_totcentre_start;
    tfriend -> Branch("lc_centreTotThetaFromStart",&theta_totcentre_start);
    vector <double> eta_totcentre_start;
    tfriend -> Branch("lc_centreTotEtaFromStart",  &eta_totcentre_start);
    vector <double> phi_totcentre_start;
    tfriend -> Branch("lc_centreTotPhiFromStart",  &phi_totcentre_start);

    //Barycentre of projections
    double centreX_pr; tfriend->Branch("lc_centreTotX_pr", &centreX_pr);
    double centreY_pr; tfriend->Branch("lc_centreTotY_pr", &centreY_pr);
    vector <double> gcentreX_pr(2);
    tfriend->Branch("lc_gCentreX_pr", &gcentreX_pr);
    vector <double> gcentreY_pr(2);
    tfriend->Branch("lc_gCentreY_pr", &gcentreY_pr);


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
        layer_start .clear();
        layerg_start.clear();
        int N = energy->size();
                
        // ts__________________________________________________________
        //Length and distance between first layers (from ts)
        int first0 = layerfirst->at(0); int last0 = layerlast->at(0);
        int first1 = layerfirst->at(1); int last1 = layerlast->at(1);
        int first_all = (first0 < first1) ? first0 : first1;
        int last_all  = (last0  > last1 ) ? last0  : last1;
        length[0] = last0 - first0 + 1;
        length[1] = last1 - first1 + 1;
        int length_max = (length[0] > length [1]) ? length[0] : length[1]; 
        whole_length = last_all - first_all + 1;
        if (i==icheck) {cout <<"Lenght Info Set"<<endl
                            << "last_all   = "<<last_all<<endl
                            << "length_max = "<<length_max<<endl;} 


        // lc__________________________________________________________
        // Initialise Energy-deposited vectors with number of entries
        ///*
        vector <double> Edi       (last_all);
        vector <double> Edi_start (length_max);
        vector <vector<double>> Ed          (2, Edi); 
        vector <double>         totEd       (last_all); 
        vector <vector<double>> Ed_start    (2, Edi_start);
        vector <double>         totEd_start (whole_length);
        vector <vector<double>>* Ed_ptr          = &Ed; //0    
        vector <double>        * totEd_ptr       = &totEd;
        vector <vector<double>>* Ed_start_ptr    = &Ed_start;
        vector <double>        * totEd_start_ptr = &totEd_start;
        tfriend -> SetBranchAddress("lc_gEPerLayer",           &Ed_ptr);  
        tfriend -> SetBranchAddress("lc_totEPerLayer",         &totEd_ptr);
        tfriend -> SetBranchAddress("lc_gEPerLayerFromStart",  &Ed_start_ptr);  
        tfriend -> SetBranchAddress("lc_totEPerLayerFromStart",&totEd_start_ptr);
        if (i==icheck) {cout <<"Barycentre addresses set"<<endl;}

        // lc__________________________________________________________
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
        tfriend -> SetBranchAddress("lc_gCentreX",     &xc_ptr);
        tfriend -> SetBranchAddress("lc_gCentreY",     &yc_ptr);
        tfriend -> SetBranchAddress("lc_gCentreZ",     &zc_ptr);
        tfriend -> SetBranchAddress("lc_gCentreTheta", &thetac_ptr);
        tfriend -> SetBranchAddress("lc_gCentreEta",   &etac_ptr);
        tfriend -> SetBranchAddress("lc_gCentrePhi",   &phic_ptr);
        tfriend -> SetBranchAddress("lc_centreTotX",     &totxc_ptr);
        tfriend -> SetBranchAddress("lc_centreTotY",     &totyc_ptr);
        tfriend -> SetBranchAddress("lc_centreTotZ",     &totzc_ptr);
        tfriend -> SetBranchAddress("lc_centreTotTheta", &tottc_ptr);
        tfriend -> SetBranchAddress("lc_centreTotEta",   &totec_ptr);
        tfriend -> SetBranchAddress("lc_centreTotPhi",   &totpc_ptr);  
        // Now same, but from start
        vector <double> x_totc_start     (whole_length);
        vector <double> y_totc_start     (whole_length);
        vector <double> z_totc_start     (whole_length);
        vector <double> theta_totc_start (whole_length);
        vector <double> eta_totc_start   (whole_length);
        vector <double> phi_totc_start   (whole_length);
        vector <double> footprint        (length_max);
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
        tfriend -> SetBranchAddress("lc_gCentreXFromStart",    &xcs_ptr);
        tfriend -> SetBranchAddress("lc_gCentreYFromStart",    &ycs_ptr);
        tfriend -> SetBranchAddress("lc_gCentreZFromStart",    &zcs_ptr);
        tfriend -> SetBranchAddress("lc_gCentreThetaFromStart",&thetacs_ptr);
        tfriend -> SetBranchAddress("lc_gCentreEtaFromStart",  &etacs_ptr);
        tfriend -> SetBranchAddress("lc_gCentrePhiFromStart",  &phics_ptr);
        tfriend -> SetBranchAddress("lc_centreTotXFromStart",    &totxcs_ptr);
        tfriend -> SetBranchAddress("lc_centreTotYFromStart",    &totycs_ptr);
        tfriend -> SetBranchAddress("lc_centreTotZFromStart",    &totzcs_ptr);
        tfriend -> SetBranchAddress("lc_centreTotThetaFromStart",&tottcs_ptr);
        tfriend -> SetBranchAddress("lc_centreTotEtaFromStart",  &totecs_ptr);
        tfriend -> SetBranchAddress("lc_centreTotPhiFromStart",  &totpcs_ptr); 
        //*/
        if (i==icheck) {cout <<"Energy addresses set"<<endl;}


        for (int j = 0; j<N; j++)
        {
            int layer_ij = layerN ->at(j);
            int which_gj = which_g->at(j);
            int layer_ij_gstart = (which_gj == 0) ? 
                                layer_ij-first0 : layer_ij-first1;
            int layer_ij_start_all = layer_ij - first_all;
            layerg_start.push_back(layer_ij_gstart);
            layer_start .push_back(layer_ij_start_all);
            double E_ij  = energy->at(j);
            
            Ed_ptr   ->at(which_gj)[layer_ij-1]        += E_ij;
            totEd_ptr->at(layer_ij-1)                  += E_ij;
            Ed_start_ptr->at(which_gj)[layer_ij_gstart]+= E_ij;//no -1  
            totEd_start_ptr->at(layer_ij_start_all)    += E_ij;
            if (i==icheck && j==0) {cout<<"Energy adding works"<<endl;}

            xc_ptr ->at(which_gj)[layer_ij-1] += x->at(j) * E_ij;
            yc_ptr ->at(which_gj)[layer_ij-1] += y->at(j) * E_ij;
            zc_ptr ->at(which_gj)[layer_ij-1] += z->at(j) * E_ij;
            totxc_ptr -> at (layer_ij-1)      += x->at(j) * E_ij;
            totyc_ptr -> at (layer_ij-1)      += y->at(j) * E_ij;
            totzc_ptr -> at (layer_ij-1)      += z->at(j) * E_ij;
            xcs_ptr ->at(which_gj)[layer_ij_gstart] += x->at(j) * E_ij;
            ycs_ptr ->at(which_gj)[layer_ij_gstart] += y->at(j) * E_ij;
            zcs_ptr ->at(which_gj)[layer_ij_gstart] += z->at(j) * E_ij;
            totxcs_ptr -> at (layer_ij_start_all)   += x->at(j) * E_ij;
            totycs_ptr -> at (layer_ij_start_all)   += y->at(j) * E_ij;
            totzcs_ptr -> at (layer_ij_start_all)   += z->at(j) * E_ij; 
            if (i==icheck && j==0) {cout<<"Coords adding works"<<endl;} 
        } 
        //print_vector(Ed_ptr->at(0));
        //print_vector(xc_ptr->at(0));
        //print_vector(xc_ptr->at(1));
        //print_vector(*totxcs_ptr);

        // Finish Barycentres (normal ones)
        for (int lj=0; lj < last_all; lj++)
        {
            //Individual Ones
            for (int g=0; g<2; g++)
            {
                double Eigj = Ed_ptr -> at(g)[lj];
                if (Eigj == 0) 
                {
                    xc_ptr -> at(g)[lj] = 1001; //TMath::QuietNaN();
                    yc_ptr -> at(g)[lj] = 1001; //TMath::QuietNaN();
                    zc_ptr -> at(g)[lj] = 1001; //TMath::QuietNaN();
                    phic_ptr  ->at(g)[lj] = 11;
                    thetac_ptr->at(g)[lj] = 11;
                    etac_ptr  ->at(g)[lj] = 11;
                }
                else 
                {
                    xc_ptr -> at(g)[lj] /= Eigj;
                    yc_ptr -> at(g)[lj] /= Eigj;
                    zc_ptr -> at(g)[lj] /= Eigj;

                    double xigj = xc_ptr -> at(g)[lj];
                    double yigj = yc_ptr -> at(g)[lj];
                    double zigj = zc_ptr -> at(g)[lj];
                    double rhoigj = sqrt(xigj*xigj + yigj*yigj);
                    double phiigj = acos( xigj / rhoigj );
                    double thetaigj = atan (rhoigj / zigj);
                    double etaigj = - log ( tan(thetaigj / 2) );
                    phic_ptr  ->at(g)[lj] = phiigj;
                    thetac_ptr->at(g)[lj] = thetaigj;
                    etac_ptr  ->at(g)[lj] = etaigj;
                }
            }
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
        if (i==icheck) {cout<<"Barycentre standard works"<<endl;} 
        //print_vector(xc_ptr->at(0));
        //print_vector(yc_ptr->at(0));

        // Finish Barycentres From Start
        for (int ljs=0; ljs < whole_length; ljs++)
        {   
            if (ljs < length_max) // Individual Ones
            {
                for (int g=0; g<2; g++)
                {
                    double Esigj = Ed_start_ptr -> at(g)[ljs];
                    if (Esigj == 0)
                    {
                        xcs_ptr -> at(g)[ljs] = 1001; //TMath::QuietNaN();
                        ycs_ptr -> at(g)[ljs] = 1001; //TMath::QuietNaN();
                        zcs_ptr -> at(g)[ljs] = 1001; //TMath::QuietNaN();
                        phics_ptr  ->at(g)[ljs] = 11;
                        thetacs_ptr->at(g)[ljs] = 11;
                        etacs_ptr  ->at(g)[ljs] = 11;
                    }
                    else 
                    {
                        xcs_ptr -> at(g)[ljs] /= Esigj;
                        ycs_ptr -> at(g)[ljs] /= Esigj;
                        zcs_ptr -> at(g)[ljs] /= Esigj;
                        double xigj = xcs_ptr -> at(g)[ljs];
                        double yigj = ycs_ptr -> at(g)[ljs];
                        double zigj = zcs_ptr -> at(g)[ljs];
                        double rhoigj = sqrt(xigj*xigj + yigj*yigj);
                        double phiigj = acos( xigj / rhoigj );
                        double thetaigj = atan (rhoigj / zigj);
                        double etaigj = - log ( tan(thetaigj / 2) );
                        phics_ptr->at(g)[ljs]   = phiigj;
                        thetacs_ptr->at(g)[ljs] = thetaigj;
                        etacs_ptr->at(g)[ljs]   = etaigj;
                    } 
                }
            }
           
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
        //print_vector(xc_ptr->at(0));
        //print_vector(xc_ptr->at(1));
        //print_vector(*totxcs_ptr);
        
        // Barycentre of projected events
        centreX_pr = mymean(*xpr, *energy);
        centreY_pr = mymean(*ypr, *energy);
        vector <double> xpr0 = getAwhereB(*xpr   , *which_g, 0);
        vector <double> ypr0 = getAwhereB(*ypr   , *which_g, 0);
        vector <double> ene0 = getAwhereB(*energy, *which_g, 0);
        gcentreX_pr[0] = mymean(xpr0, ene0);
        gcentreY_pr[0] = mymean(ypr0, ene0);
        vector <double> xpr1 = getAwhereB(*xpr   , *which_g, 1);
        vector <double> ypr1 = getAwhereB(*ypr   , *which_g, 1);
        vector <double> ene1 = getAwhereB(*energy, *which_g, 1);
        gcentreX_pr[1] = mymean(xpr1, ene1);
        gcentreY_pr[1] = mymean(ypr1, ene1);


        //Layer Max
        LmaxE[0] = argmax(Ed[0]) + 1;
        LmaxE[1] = argmax(Ed[1]) + 1;
        LmaxETot = argmax(totEd);
        LmaxE_start[0] = argmax(Ed_start[0]);
        LmaxE_start[1] = argmax(Ed_start[1]);
        LmaxETot_start = argmax(totEd_start);
        if (i==icheck) {cout<<"Lmaxs work"<<endl;} 
        
        tfriend -> Fill();
        if (i==icheck) {cout<<"Fill works"<<endl;}
    }
    cout<<"Checkpoint: End of Loop"<<endl;
    tfriend -> Write();
    cout<<"Write works"<<endl;
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