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
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "Math/Vector2D.h"

#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"
using XYVector = ROOT::Math::XYVector;
using XYZVector = ROOT::Math::XYZVector;
using namespace std;

/* 
 
 ADD TO FRIEND INFO RELATED TO SHAPE OF HITS
 - Layer # with maximum # of clusters by one photon
 - Layer # with maximum # of clusters in total
 - Layer # since start of own shower with maximum # of clusters 
    by that photon
 - Layer # since start of first shower with maximum # of clusters

 - max distance between lcs per layer
 - max distance between lcs of same photon
 - max distance between lcs of different photons
 - max distance between lcs of different photons minus initial distance
 - tag max distance ij: i = 0(1): for same(different) photons
                        j = 0(1): for at least one(no) lc with rechits<3
    This gives combinations: 00=0, 01=1, 10, 11. Plus, -1 if no event.
 - max distance from barycentre 
 - max distance from barycentre excluding small clusters (<3)

 - area of hits of individual photons
 - sum of those two areas
 - area defined by all entries
*/


void AddShape () 
{
    string path    = "../myTTrees/TreesEn90to110/Double/"; //none if none
    string rootname= "step3ticl_eta21_FlatTracksters";//root file
    string friendrootname= "treeFriends"; //root file with friends
    string dirname = "none";
    string treename= "TSTree_SimFromCP";

    int icheck = 5; //Entry to check in loop to know where errors lie.
    
    double lcoeff = 10.; // for normalised long 2nd moment
    double rcoeff =  4.; // for normalised tran 2nd moment 


    ///////////////////////////////////////////////////////////////////
    //1. //////////////////////////////////////////////////////////////
    //Get stuff from the TTree we already have
    ///////////////////////////////////////////////////////////////////
    TFile* file = open_file(path, rootname);
    TTree* tree = get_tree (file, dirname, treename);
    //tree->Print();
    //tree ->Scan("lc_TSidx:lc_nrechits:lc_layer:lc_mult:lc_tsMult:lc_seedEnergy");
    //return;
    int Nentries = tree->GetEntries();
    vector <double>* x = 0;
    vector <double>* y = 0;
    vector <double>* z = 0;
    vector <double>* cpx = 0;
    vector <double>* cpy = 0;
    vector <double>* cpz = 0;
    vector <double>* eneoverlap = 0;
    vector <int>*    which_g = 0;
    vector <int>*    layerN  = 0;
    vector <int>*    layerfirst = 0;
    vector <int>*    layerlast  = 0;
    vector <int>*    nrechits = 0;
    vector <int>*    mult = 0;//100 entries: 50 for whihc_g = 0, 50 for 1.
    tree->SetBranchAddress("lc_x",        &x);
    tree->SetBranchAddress("lc_y",        &y);
    tree->SetBranchAddress("lc_z",        &z);
    tree->SetBranchAddress("cp_vtxX",     &cpx);
    tree->SetBranchAddress("cp_vtxY",     &cpy);
    tree->SetBranchAddress("cp_vtxZ",     &cpz);
    tree->SetBranchAddress("lc_energy",   &eneoverlap);
    tree->SetBranchAddress("lc_layer",      &layerN);
    tree->SetBranchAddress("lc_TSidx",      &which_g);
    tree->SetBranchAddress("lc_nrechits",   &nrechits);
    tree->SetBranchAddress("lc_mult",       &mult);
    tree->SetBranchAddress("ts_firstLayer", &layerfirst);
    tree->SetBranchAddress("ts_lastLayer",  &layerlast);
    //tree->Scan("event:lc_TSidx:lc_x:lc_layer:sv_posX:cp_vtxX:cp_phi");
    //tree->Scan("lc_x:lc_y:lc_z:lc_energy:lc_tsMult");

    ///////////////////////////////////////////////////////////////////
    //2. //////////////////////////////////////////////////////////////
    //Prepare new TTree, which will later become friend
    ///////////////////////////////////////////////////////////////////
    TFile* fr = open_file(path, friendrootname, "update");
    const char* friend_namecycle = (treename + "_Shape;*").c_str();
    fr->Delete(friend_namecycle);

    // Create new shape ttree
    const char* placeholder_name = (treename + "_ShapePlaceholder").c_str();
    TTree* tfriendS = get_tree(fr, "none", placeholder_name);
    TTree* tfriend =  tfriendS -> CloneTree(0);

    // Take projection coordinates from Coords TTree
    const char* coords_name = (treename + "_Coords").c_str();
    TTree* tfriendC = get_tree(fr, "none", coords_name);
    vector <double>* xpr = 0; 
    tfriendC -> SetBranchAddress("lc_x_pr", &xpr);
    vector <double>* ypr = 0; 
    tfriendC -> SetBranchAddress("lc_y_pr", &ypr);

    // Take centres coordinates from Energy TTree.
    const char* energies_name = (treename + "_Energy").c_str();
    TTree* tfriendE = get_tree(fr, "none", energies_name);
    vector <double>* energy = 0; 
    vector <double>* tsfraction = 0;
    vector <double>* centrex = 0; //for each layer
    vector <double>* centrey = 0; //for each layer
    double BCx = 0;
    double BCy = 0;
    double BCz = 0;
    double BCl = 0;
    vector <double>* gBCx = 0; 
    vector <double>* gBCy = 0; 
    vector <double>* gBCz = 0; 
    vector <double>* gBCl = 0; 
    double centrex_pr = 0; 
    double centrey_pr = 0; 
    vector<double>* gcentrex_pr = 0; 
    vector<double>* gcentrey_pr = 0; 
    tfriendE -> SetBranchAddress("lc_energyTrue", &energy);
    tfriendE -> SetBranchAddress("lc_tsFraction", &tsfraction);
    tfriendE -> SetBranchAddress("lc_centreTotXFromStart", &centrex);
    tfriendE -> SetBranchAddress("lc_centreTotYFromStart", &centrey);
    tfriendE -> SetBranchAddress("lc_BCx",   &BCx);
    tfriendE -> SetBranchAddress("lc_BCy",   &BCy);
    tfriendE -> SetBranchAddress("lc_BCz",   &BCz);
    tfriendE -> SetBranchAddress("lc_BCl",   &BCl);
    tfriendE -> SetBranchAddress("lc_gBCx",  &gBCx);
    tfriendE -> SetBranchAddress("lc_gBCy",  &gBCy);
    tfriendE -> SetBranchAddress("lc_gBCz",  &gBCz);
    tfriendE -> SetBranchAddress("lc_gBCl",  &gBCl);
    tfriendE -> SetBranchAddress("lc_BCx_pr", &centrex_pr);
    tfriendE -> SetBranchAddress("lc_BCy_pr", &centrey_pr);
    tfriendE -> SetBranchAddress("lc_gBCx_pr",   &gcentrex_pr);
    tfriendE -> SetBranchAddress("lc_gBCy_pr",   &gcentrey_pr);
    //tfriend -> Print();



    ////////////////////////////////////////////////////////////////////
    // SET NEW BRANCHES ////////////////////////////////////////////////
    //__________________________________________________________________

    //Layers - Hits info, use lc
    vector <int> totmult;
    vector <vector<double>> gmult (2);
    vector <int> totmultbig;
    vector <vector<double>> gmultbig (2);
    tfriend->Branch("lc_totMult",    &totmult);
    tfriend->Branch("lc_gMult",      &gmult);
    tfriend->Branch("lc_totMultBig", &totmultbig);
    tfriend->Branch("lc_gMultBig",   &gmultbig);
    // FOLLOWING ACTUALLY WRONG, NOT WHAT I WANTED
    vector <int> totcummult;
    vector <vector<double>> gcummult (2);
    // tfriend->Branch("lc_totCumMult", &totcummult);  
    // tfriend->Branch("lc_gCumMult",   &gcummult); 

    vector <int> LmaxMult (2); 
    int LmaxMultTot;
    vector <int> LmaxMult_start (2); 
    int LmaxMultTot_start;
    tfriend->Branch("lc_layerMaxMult",             &LmaxMult);
    tfriend->Branch("lc_layerMaxTotMult",          &LmaxMultTot);
    tfriend->Branch("lc_layerMaxMultFromStart",    &LmaxMult_start);  
    tfriend->Branch("lc_layerMaxTotMultFromStart", &LmaxMultTot_start); 


    // Longitudinal & Transverse 2nd Moments 
    double lambda2 = 0;
    double erre2   = 0;
    double lnorm   = 0;
    double rnorm   = 0;
    vector <double> glambda2 (2);
    vector <double> gerre2 (2);
    vector <double> glnorm (2);
    vector <double> grnorm (2);
    tfriend->Branch("lc_2ndLongMoment", &lambda2);
    tfriend->Branch("lc_2ndTranMoment", &erre2);
    tfriend->Branch("lc_2ndLongMomentNorm", &lnorm);
    tfriend->Branch("lc_2ndTranMomentNorm", &rnorm); 
    tfriend->Branch("lc_g2ndLongMoment",&glambda2);
    tfriend->Branch("lc_g2ndTranMoment",&gerre2);
    tfriend->Branch("lc_g2ndLongMomentNorm",&glnorm);
    tfriend->Branch("lc_g2ndTranMomentNorm",&grnorm); 

    

    //Distances between hits - from start
    vector <double> maxdist_start; 
    tfriend->Branch("lc_maxDistFromStart",         &maxdist_start);
    vector <int> maxdist_tag_start; 
    tfriend->Branch("lc_maxDistTagFromStart",      &maxdist_tag_start);
    vector <double> maxdist_big_start; 
    tfriend->Branch("lc_maxDistBigFromStart",      &maxdist_big_start);
    vector <vector<double>> maxdist_same_start (2); 
    tfriend->Branch("lc_maxDistSameFromStart",     &maxdist_same_start);
    vector <vector<double>> maxdist_same_big_start (2); 
    tfriend->Branch("lc_maxDistSameBigFromStart",  &maxdist_same_big_start);
    vector <double> maxdist_diff_start; 
    tfriend->Branch("lc_maxDistDiffFromStart",     &maxdist_diff_start); 
    vector <double> maxdist_diff_net_start; 
    tfriend->Branch("lc_maxDistDiffNetFromStart",  &maxdist_diff_net_start);
    vector <double> maxdist_diff_big_start; 
    tfriend->Branch("lc_maxDistDiffBigFromStart",   &maxdist_diff_big_start); 
    vector <double> maxdist_diff_net_big_start; 
    tfriend->Branch("lc_maxDistDiffBigNetFromStart",&maxdist_diff_net_big_start); 
    vector <double> maxdist_centre_start; 
    tfriend->Branch("lc_maxDistCentreFromStart",   &maxdist_centre_start);  
    vector <double> maxdist_centre_big_start; 
    tfriend->Branch("lc_maxDistCentreBigFromStart",&maxdist_centre_big_start);  
    

    // Areas of layers' sections (as convex hulls of detections)
    vector <double> areatot_start;
    tfriend->Branch("lc_areaTotFromStart", &areatot_start);
    vector <vector<double>> areag_start (2);
    tfriend->Branch("lc_areagFromStart",   &areag_start);


    //Distances between hits - projections
    double maxdist_pr; 
    tfriend->Branch("lc_maxDist_pr",          &maxdist_pr);
    int maxdist_tag_pr; 
    tfriend->Branch("lc_maxDistTag_pr",       &maxdist_tag_pr);
    double maxdist_big_pr; 
    tfriend->Branch("lc_maxDistBig_pr",       &maxdist_big_pr);
    vector<double> maxdist_same_pr (2); 
    tfriend->Branch("lc_maxDistSame_pr",      &maxdist_same_pr);
    vector<double> maxdist_same_big_pr (2); 
    tfriend->Branch("lc_maxDistSameBig_pr",   &maxdist_same_big_pr);
    double maxdist_diff_pr; 
    tfriend->Branch("lc_maxDistDiff_pr",      &maxdist_diff_pr); 
    double maxdist_diff_net_pr; 
    tfriend->Branch("lc_maxDistDiffNet_pr",   &maxdist_diff_net_pr);
    double maxdist_diff_big_pr; 
    tfriend->Branch("lc_maxDistDiffBig_pr",   &maxdist_diff_big_pr); 
    double maxdist_diff_net_big_pr; 
    tfriend->Branch("lc_maxDistDiffBigNet_pr",&maxdist_diff_net_big_pr); 
    double maxdist_centre_pr; 
    tfriend->Branch("lc_maxDistCentre_pr",    &maxdist_centre_pr);  
    double maxdist_centre_big_pr; 
    tfriend->Branch("lc_maxDistCentreBig_pr", &maxdist_centre_big_pr);
    vector <double> maxdist_gcentre_pr (2); 
    tfriend->Branch("lc_maxDistGCentre_pr",   &maxdist_gcentre_pr);  
    vector <double> maxdist_gcentre_big_pr (2); 
    tfriend->Branch("lc_maxDistGCentreBig_pr",&maxdist_gcentre_big_pr);
    vector <double> maxdist_gVtx_pr (2); 
    tfriend->Branch("lc_maxDistGVtx_pr",      &maxdist_gVtx_pr);  
    vector <double> maxdist_gVtx_big_pr (2); 
    tfriend->Branch("lc_maxDistGVtxBig_pr",   &maxdist_gVtx_big_pr);  


    // Areas of projections (as convex hulls of all detections)
    double Atot_pr;
    vector<double> Ag_pr(2);
    double AtotBig_pr;
    vector<double> AgBig_pr(2);
    tfriend->Branch("lc_A_pr",    &Atot_pr);
    tfriend->Branch("lc_gA_pr",   &Ag_pr); 
    tfriend->Branch("lc_ABig_pr", &AtotBig_pr);
    tfriend->Branch("lc_gABig_pr",&AgBig_pr);
    double Astartot_pr;
    vector<double> Astarg_pr(2);
    double AstartotBig_pr;
    vector<double> AstargBig_pr(2);
    tfriend->Branch("lc_Astar_pr",    &Astartot_pr);
    tfriend->Branch("lc_gAstar_pr",   &Astarg_pr); 
    tfriend->Branch("lc_AstarBig_pr", &AstartotBig_pr);
    tfriend->Branch("lc_gAstarBig_pr",&AstargBig_pr);


    // Ratio Distance2-Areas of projections
    double Atodm2_ratio_pr;
    tfriend->Branch("lc_AtoDmax2Ratio_pr",    &Atodm2_ratio_pr);
    double Atodm2_ratioBig_pr;
    tfriend->Branch("lc_AtoDmax2BigsRatio_pr",&Atodm2_ratioBig_pr);
    double Atodc2_ratio_pr;
    tfriend->Branch("lc_AtoDcen2Ratio_pr",    &Atodc2_ratio_pr);
    double Atodc2_ratioBig_pr;
    tfriend->Branch("lc_AtoDcen2BigsRatio_pr",&Atodc2_ratioBig_pr);
    vector <double> gAtodm2_ratio_pr   (2);
    tfriend->Branch("lc_gAtoDmax2Ratio_pr",    &gAtodm2_ratio_pr);
    vector <double> gAtodm2_ratioBig_pr(2);
    tfriend->Branch("lc_gAtoDmax2BigsRatio_pr",&gAtodm2_ratioBig_pr);
    vector <double> gAtodc2_ratio_pr   (2);
    tfriend->Branch("lc_gAtoDcen2Ratio_pr",    &gAtodc2_ratio_pr);
    vector <double> gAtodc2_ratioBig_pr(2);
    tfriend->Branch("lc_gAtoDcen2BigsRatio_pr",&gAtodc2_ratioBig_pr);


    // Radii with 68, 90 & 95% energy
    double radius68_pr;
    tfriend->Branch("lc_RadiusTot68_pr", &radius68_pr);
    double radius90_pr;
    tfriend->Branch("lc_RadiusTot90_pr", &radius90_pr);
    double radius95_pr;
    tfriend->Branch("lc_RadiusTot95_pr", &radius95_pr);
    vector <double> radius68g_pr (2);
    tfriend->Branch("lc_RadiusG68_pr", &radius68g_pr);
    vector <double> radius90g_pr (2);
    tfriend->Branch("lc_RadiusG90_pr", &radius90g_pr);
    vector <double> radius95g_pr (2);
    tfriend->Branch("lc_RadiusG95_pr", &radius95g_pr);
    double ratio9068_pr;
    tfriend->Branch("lc_RatioTot9068_pr", &ratio9068_pr);
    double ratio9568_pr;
    tfriend->Branch("lc_RatioTot9568_pr", &ratio9568_pr);
    double ratio9590_pr;
    tfriend->Branch("lc_RatioTot9590_pr", &ratio9590_pr);
    vector <double> ratio9068g_pr (2);
    tfriend->Branch("lc_RatioG9068_pr", &ratio9068g_pr);
    vector <double> ratio9568g_pr (2);
    tfriend->Branch("lc_RatioG9568_pr", &ratio9568g_pr);
    vector <double> ratio9590g_pr (2);
    tfriend->Branch("lc_RatioG9590_pr", &ratio9590g_pr);   


    ///////////////////////////////////////////////////////////////////
    //3. //////////////////////////////////////////////////////////////
    // LOOP & Update Branches
    ///////////////////////////////////////////////////////////////////

    // First a lambda function needed for big hits
    auto bighits = [] (int hit, int th = 3) {return hit>=3;};

    // Now loop
    for (int i = 0; i<Nentries; i++)
    { 
        //Prepare for next cycle
        cout<<"Entry "<<i<<endl;
        tree     -> GetEntry(i);
        tfriendC -> GetEntry(i);
        tfriendE -> GetEntry(i);
        tfriendS -> GetEntry(i);
        int N = x->size();

        if (i!=0) //clear vectors
        {
            totmult.clear();
            gmult[0].clear();
            gmult[1].clear();
            totmultbig.clear();
            gmultbig[0].clear();
            gmultbig[1].clear();
            //totcummult.clear();
            //gcummult[0].clear();
            //gcummult[1].clear();

            maxdist_start.clear();  
            maxdist_big_start.clear(); 
            maxdist_same_start[0].clear();
            maxdist_same_start[1].clear(); 
            maxdist_same_big_start[0].clear();
            maxdist_same_big_start[1].clear(); 

            maxdist_diff_start.clear(); 
            maxdist_diff_net_start.clear(); 
            maxdist_diff_big_start.clear(); 
            maxdist_diff_net_big_start.clear(); 

            maxdist_centre_start.clear(); 
            maxdist_centre_big_start.clear(); 
            maxdist_tag_start.clear();
            
            areatot_start.clear();
            areag_start[0].clear();
            areag_start[1].clear();
        }

        // ts__________________________________________________________
        //Length and distance between first layers (from ts)
        int first0 = layerfirst->at(0); int last0 = layerlast->at(0);
        int first1 = layerfirst->at(1); int last1 = layerlast->at(1);
        int first_all = (first0 < first1) ? first0 : first1;
        int last_all  = (last0  > last1 ) ? last0  : last1;
        vector <int> length = {last0 - first0 + 1, last1 - first1 + 1};
        


        //cp___________________________________________________________
        XYVector cpvtx0 = XYVector(cpx->at(0), cpy->at(0));
        XYVector cpvtx1 = XYVector(cpx->at(1), cpy->at(1));
        //cout<<"cp done"<<endl;


        // lc : Moments _______________________________________________
        double Etot = 0;
        double E0tot = 0;
        double E1tot = 0;

        auto b = eneoverlap->begin(); auto e = eneoverlap->end();
        auto itmax = max_element(b, e);
        double Eomax = eneoverlap->at(distance(b, itmax));
        auto it2max = max_element(b, e,
                                [Eomax](double &a, double &b) {
                                    if (a == Eomax) return true;
                                    if (b == Eomax) return false;                                         
                                    return a < b;
                                });
        double Eo2max = eneoverlap->at(distance(b, it2max));
        //cout<<Eomax<<","<<Eo2max<<endl;
        
        vector <double> ene0  = getAwhereB(*energy, *which_g, 0);
        sort(ene0.begin(), ene0.end());
        double E0max    = ene0[ene0.size()-1];
        double E02max = ene0[ene0.size()-2];
        vector <double> ene1  = getAwhereB(*energy, *which_g, 1);
        sort(ene1.begin(), ene1.end());
        double E1max    = ene1[ene1.size()-1];
        double E12max = ene1[ene1.size()-2];

        double el2 = 0; double elmax = 0;
        vector <double> gel2(2); vector <double> gelmax(2);
        double er2 = 0; double ermax = 0;
        vector <double> ger2(2); vector <double> germax(2);

        //cout<<"Entry "<<i<<endl;
        vector <int> kdoubles;
        for (int j = 0; j<N; j++) // for each entry
        {
            double Ej  = energy->at(j);
            double Eoj = eneoverlap->at(j);
            double lj2 = (z->at(j) - BCz) * (z->at(j) - BCz) * Ej;
            double rj2 = ( (x->at(j) - BCx) * (x->at(j) - BCx)
                        +(y->at(j) - BCy) * (y->at(j) - BCy) ) * Ej;
           
            Etot    += Ej;
            lambda2 += lj2;
            erre2   += rj2;

            // NO NEED DOUBLE CHECK, CAUSE ADD FRACTION TWICE
            //if (contains(kdoubles, j)) {continue;}
            if ( Eoj == Eomax || Eoj == Eo2max ) 
            {
                //cout<<Ej<<endl;
                elmax += lcoeff*lcoeff*Ej;
                ermax += rcoeff*rcoeff*Ej;
            }
            else
            {
                el2 += lj2*Ej;
                er2 += rj2*Ej;
            } 

            int gj = which_g->at(j);
            if (gj == 0)
            {
                E0tot += Ej;
                glambda2[0] += lj2;
                gerre2  [0] += rj2;
                if ( Ej == E0max || Ej == E02max ) 
                    {
                        gelmax[0] += lcoeff*lcoeff*Ej;
                        germax[0] += rcoeff*rcoeff *Ej;
                    }
                else
                    {
                        gel2[0] += lj2*Ej;
                        ger2[0] += rj2*Ej;
                    }
            } 
            else
            {   
                E1tot += Ej;
                glambda2[1] += lj2;
                gerre2  [1] += rj2;
                if ( Ej == E1max || Ej == E12max ) 
                    {
                        gelmax[1] += lcoeff*lcoeff*Ej;
                        germax[1] += rcoeff*rcoeff *Ej;
                    }
                else
                    {
                        gel2[1] += lj2*Ej;
                        ger2[1] += rj2*Ej;
                    }
            }    
        }
        lambda2    /= Etot; erre2    /= Etot;
        glambda2[0]/=E0tot; gerre2[0]/=E0tot;
        glambda2[1]/=E0tot; gerre2[1]/=E0tot;

        lnorm = el2 / (el2 + elmax);
        rnorm = er2 / (er2 + ermax);
        glnorm[0] = gel2[0] / (gel2[0] + gelmax[0]);
        glnorm[1] = gel2[1] / (gel2[1] + gelmax[1]);
        grnorm[0] = ger2[0] / (ger2[0] + germax[0]);
        grnorm[1] = ger2[1] / (ger2[1] + germax[1]);



        //lc__________________________________________________________
        for (int ln = first_all; ln<=last_all; ln++) // for each layer
        {
            vector <double> xn = getAwhereB(*x, *layerN, ln);
            vector <double> yn = getAwhereB(*y, *layerN, ln);
            vector <int> hitsn = getAwhereB(*nrechits  , *layerN, ln);
            vector <int> gn    = getAwhereB(*which_g   , *layerN, ln);
            vector <double> fn = getAwhereB(*tsfraction, *layerN, ln);

            vector <double> xnbig = getAwhereB(xn, hitsn, bighits);
            vector <double> ynbig = getAwhereB(yn, hitsn, bighits);
            vector <int>    gnbig = getAwhereB(gn, hitsn, bighits);
            vector <double> fnbig = getAwhereB(fn, hitsn, bighits);
            

            // Mult ____________________________________________________
            double totcount    = myvecsum(fn);
            int totcount_i = round(totcount)+0.001;
            totmult   .push_back(totcount_i);
            double totcountbig = myvecsum(fnbig);
            int totcountbig_i = round(totcountbig)+0.01;
            totmultbig.push_back(totcountbig_i);
            //int cumcount = (ln == first_all)? 
            //               totcount:totcount + totcummult[ln-first_all-1];
            //totcummult.push_back(cumcount);
            if (ln>=first0 && ln<=last0)
            {
                int countg = count(gn.begin(), gn.end(), 0);
                gmult[0]   .push_back(countg);
                int countgbig = count(gnbig.begin(), gnbig.end(), 0);
                gmultbig[0].push_back(countgbig);
                //int cumcountg = (ln == first0)? 
                //                countg:countg + gcummult[0][ln-first0-1];
                //gcummult[0].push_back(cumcountg);
            }
            if (ln>=first1 && ln<=last1)
            {
                int countg = count(gn.begin(), gn.end(), 1);
                gmult[1]   .push_back(countg);
                int countgbig = count(gnbig.begin(), gnbig.end(), 1);
                gmultbig[1].push_back(countgbig);
                //int cumcountg = (ln == first1)? 
                //                countg:countg + gcummult[1][ln-first1-1];
                //gcummult[1].push_back(cumcountg);
            }       



            // Distances _______________________________________________
            vector <XYVector> posvec (xn.size());
            for (int j = 0; j<xn.size(); j++)
                {posvec[j] = XYVector(xn[j], yn[j]);}
            double cx = centrex->at(ln-first_all);
            double cy = centrey->at(ln-first_all); 
            XYVector centre = XYVector(cx, cy); 

            
            vector <double> soln=get_max_distances(posvec, gn, hitsn, 
                                                   cpvtx0, cpvtx1, 
                                                   centre);
            maxdist_start             .push_back(soln[0]);
            maxdist_big_start         .push_back(soln[1]);
            maxdist_same_start[0]     .push_back(soln[2]);
            maxdist_same_start[1]     .push_back(soln[3]);
            maxdist_same_big_start[0] .push_back(soln[4]);
            maxdist_same_big_start[1] .push_back(soln[5]);
            
            maxdist_diff_start        .push_back(soln[6]);
            maxdist_diff_big_start    .push_back(soln[7]);
            maxdist_diff_net_start    .push_back(soln[8]);
            maxdist_diff_net_big_start.push_back(soln[9]);
            
            maxdist_centre_start      .push_back(soln[10]);
            maxdist_centre_big_start  .push_back(soln[11]);
            maxdist_tag_start         .push_back(soln[12]);


            // areas __________________________________________________
            vector <double> xn0 = getAwhereB(xn, gn, 0);
            vector <double> yn0 = getAwhereB(yn, gn, 0);
            vector <double> xn1 = getAwhereB(xn, gn, 1);
            vector <double> yn1 = getAwhereB(yn, gn, 1);

            areatot_start.push_back(Shoelace(ConvexHull(xn, yn)));
            double A0_i = Shoelace(ConvexHull(xn0, yn0));
            areag_start[0].push_back(A0_i);
            double A1_i = Shoelace(ConvexHull(xn1, yn1));
            areag_start[1].push_back(A1_i);

        }
        //cout<<"Individual layers' stuff done"<<endl;

        // ____________________________________________________________
        // Distances in projected space _______________________________
        XYVector cpvtx0_pr = XYVector(cpx->at(0)/cpz->at(0), 
                                      cpy->at(0)/cpz->at(0));
        XYVector cpvtx1_pr = XYVector(cpx->at(1)/cpz->at(1), 
                                      cpy->at(1)/cpz->at(1));
        XYVector centre0_pr = XYVector(gcentrex_pr->at(0), 
                                       gcentrey_pr->at(0));
        XYVector centre1_pr = XYVector(gcentrex_pr->at(1), 
                                       gcentrey_pr->at(1));
        XYVector centre_pr  = XYVector(centrex_pr, centrey_pr);
        vector <double> soln = get_max_distances(*xpr, *ypr, 
                                                 *which_g,   *nrechits, 
                                                 cpvtx0_pr,  cpvtx1_pr,
                                                 centre0_pr, centre1_pr,
                                                 centre_pr);
        maxdist_pr             = soln[0];
        maxdist_big_pr         = soln[1];
        maxdist_same_pr        = {soln[2], soln[3]};
        maxdist_same_big_pr    = {soln[4], soln[5]};
        maxdist_diff_pr        = soln[6];
        maxdist_diff_big_pr    = soln[7];
        maxdist_diff_net_pr    = soln[8];
        maxdist_diff_net_big_pr= soln[9];
        maxdist_centre_pr      = soln[10];
        maxdist_centre_big_pr  = soln[11];
        maxdist_gcentre_pr     = {soln[12], soln[13]};
        maxdist_gcentre_big_pr = {soln[14], soln[15]};
        maxdist_tag_pr         = soln[16];
        maxdist_gVtx_pr        = {soln[17], soln[18]};
        maxdist_gVtx_big_pr    = {soln[19], soln[20]};
        //cout<<"get projected max done"<<endl;


        // ____________________________________________________________
        // Radii with % of energy (projected space) ____________________
        vector<double> radii =get_radii(*xpr, *ypr, *which_g, *energy,
                                    centre0_pr,centre1_pr,centre_pr);
        radius68_pr  = radii[0];
        radius90_pr  = radii[1];
        radius95_pr  = radii[2];
        radius68g_pr = {radii[3], radii[6]};
        radius90g_pr = {radii[4], radii[7]};
        radius95g_pr = {radii[5], radii[8]};

        ratio9068_pr  = radius90_pr / radius68_pr;
        ratio9568_pr  = radius95_pr / radius68_pr;
        ratio9590_pr  = radius95_pr / radius90_pr;
        ratio9068g_pr = {radius90g_pr[0] / radius68g_pr[0],
                         radius90g_pr[1] / radius68g_pr[1]};
        ratio9568g_pr = {radius95g_pr[0] / radius68g_pr[0],
                         radius95g_pr[1] / radius68g_pr[1]};
        ratio9590g_pr = {radius95g_pr[0] / radius90g_pr[0],
                         radius95g_pr[1] / radius90g_pr[1]};
        
        if (radius68g_pr[argmax(radius68g_pr)] > 0.01 ||
            radius90g_pr[argmax(radius90g_pr)] > 0.05 ||
            radius95g_pr[argmax(radius95g_pr)] > 0.05)
            {print_vector(radii);}

        if (ratio9068_pr     < 1 || ratio9590_pr     < 1 ||
            ratio9068g_pr[0] < 1 || ratio9068g_pr[1] < 1 ||
            ratio9590g_pr[0] < 1 || ratio9590g_pr[1] < 1 )
          {cout<<"A ratio is impossibly smaller than 1!"<<endl;return;}


        // ____________________________________________________________
        // Areas of projections _______________________________________
        vector <double> xpr0 = getAwhereB(*xpr     , *which_g, 0);
        vector <double> ypr0 = getAwhereB(*ypr     , *which_g, 0);
        vector <double> xpr1 = getAwhereB(*xpr     , *which_g, 1);
        vector <double> ypr1 = getAwhereB(*ypr     , *which_g, 1);

        Atot_pr  = Shoelace(ConvexHull(*xpr, *ypr));
        Ag_pr[0] = Shoelace(ConvexHull(xpr0, ypr0));
        Ag_pr[1] = Shoelace(ConvexHull(xpr1, ypr1));

        auto shull = StarHull(*xpr, *ypr,
                              centrex_pr, centrey_pr);
        Astartot_pr  = Shoelace(shull[0], shull[1]);
        shull = StarHull(xpr0, ypr0,
                         gcentrex_pr->at(0), gcentrey_pr->at(0));
        Astarg_pr[0] = Shoelace(shull[0], shull[1]);
        shull = StarHull(xpr1, ypr1,
                         gcentrex_pr->at(1), gcentrey_pr->at(1));
        Astarg_pr[1] = Shoelace(shull[0], shull[1]);
        //cout<<"areas of projections done"<<endl;

        //Big
        vector <double> xprbig =getAwhereB(*xpr    ,*nrechits, bighits);
        vector <double> yprbig =getAwhereB(*ypr    ,*nrechits, bighits);
        vector <int>    gbig   =getAwhereB(*which_g,*nrechits, bighits);
        AtotBig_pr = Shoelace(ConvexHull(xprbig, yprbig));
        vector <double> xpr0big = getAwhereB(xprbig, gbig, 0);
        vector <double> ypr0big = getAwhereB(yprbig, gbig, 0);
        AgBig_pr[0] = Shoelace(ConvexHull(xpr0big, ypr0big));
        vector <double> xpr1big = getAwhereB(xprbig, gbig, 1);
        vector <double> ypr1big = getAwhereB(yprbig, gbig, 1);
        AgBig_pr[1] = Shoelace(ConvexHull(xpr1big, ypr1big));
        shull = StarHull(xprbig, yprbig,
                         centrex_pr, centrey_pr);
        AstartotBig_pr  = Shoelace(shull[0], shull[1]);
        shull = StarHull(xpr0big, ypr0big,
                         gcentrex_pr->at(0), gcentrey_pr->at(0));
        AstargBig_pr[0] = Shoelace(shull[0], shull[1]);
        shull = StarHull(xpr1big, ypr1big,
                         gcentrex_pr->at(1), gcentrey_pr->at(1));
        AstargBig_pr[1] = Shoelace(shull[0], shull[1]);


        //______________________________________________________________
        // Area-Distances2 Ratios of projections _______________________
        Atodc2_ratio_pr    =Atot_pr/
                          (maxdist_centre_pr*maxdist_centre_pr);
        Atodc2_ratioBig_pr =AtotBig_pr/
                          (maxdist_centre_big_pr*maxdist_centre_big_pr);
        Atodm2_ratio_pr    =Atot_pr/
                          (maxdist_pr*maxdist_pr);
        Atodm2_ratioBig_pr =AtotBig_pr/
                          (maxdist_big_pr*maxdist_big_pr);

        gAtodc2_ratio_pr[0]   = Ag_pr[0]/
                             (maxdist_gcentre_pr[0]
                              *maxdist_gcentre_pr[0]);
        gAtodc2_ratio_pr[1]   = Ag_pr[1]/
                             (maxdist_gcentre_pr[1]
                              *maxdist_gcentre_pr[1]);
        gAtodc2_ratioBig_pr[0]= AgBig_pr[0]/
                             (maxdist_gcentre_big_pr[0]
                              *maxdist_gcentre_big_pr[0]);
        gAtodc2_ratioBig_pr[1]= AgBig_pr[1]/
                              (maxdist_gcentre_big_pr[1]
                              *maxdist_gcentre_big_pr[1]);
        gAtodm2_ratio_pr[0]   = Ag_pr[0]/
                              (maxdist_same_pr[0]
                              *maxdist_same_pr[0]);
        gAtodm2_ratio_pr[1]   = Ag_pr[1]/
                              (maxdist_same_pr[1]
                              *maxdist_same_pr[1]);
        gAtodm2_ratioBig_pr[0]= AgBig_pr[0]/
                              (maxdist_same_big_pr[0]
                              *maxdist_same_big_pr[0]);
        gAtodm2_ratioBig_pr[1]= AgBig_pr[1]/
                              (maxdist_same_big_pr[1]
                              *maxdist_same_big_pr[1]);


        // Layer number with maximum number of hits ____________________
        vector <int> mult_both = sum_portions(*mult, 2); //sum halves
        LmaxMult[0] = argmax(*mult, 0, 50) + 1;
        LmaxMult[1] = argmax(*mult, 50, 99) + 1 - 50;
        LmaxMultTot = argmax(mult_both) + 1;
        if (LmaxMultTot > 50) {LmaxMultTot -= 50;}
        LmaxMult_start[0] = LmaxMult[0] - first0;
        LmaxMult_start[1] = LmaxMult[1] - first1;
        LmaxMultTot_start = LmaxMultTot - first_all;
        //cout<<"lmax done"<<endl;

        tfriend -> Fill();
        //cout<<"Fill works"<<endl;
    }
    cout<<"Checkpoint: End of Loop"<<endl;
    const char* friend_name = (treename + "_Shape").c_str();
    tfriend -> SetName(friend_name);
    tfriend -> Write("", TObject::kOverwrite);
    cout<<"Write works"<<endl;
    const char* placeholder_namecycle = (treename + "_ShapePlaceholder;*").c_str();
    fr->Delete(placeholder_namecycle);
    cout<<"Delete works"<<endl;
    fr->Close();
    return;
    cout<<"Before Adding Friend, friend_name = "<<friend_name<<endl;
    tree -> AddFriend(friend_name, "treeFriends.root");
    cout<<"Friend added, done!"<<endl;
    file->Close();
    //const char* friend_name = (treename + "_friend").c_str();
    //cout<<"Before Adding Friend, friend_name = "<<friend_name<<endl;
    //tree -> AddFriend(friend_name, "treeFriends.root");
    //cout<<"Friend added, done!"<<endl;
    //file->Close();
}

 /*//Distances between hits
    vector <double> maxdist; 
    tfriend->Branch("lc_MaxDist",         &maxdist);
    vector <int> maxdist_tag; 
    tfriend->Branch("lc_MaxDistTag",      &maxdist_tag);
    vector <double> maxdist_big; 
    tfriend->Branch("lc_MaxDistBig",      &maxdist_big);
    vector <vector<double>> maxdist_same (2); 
    tfriend->Branch("lc_MaxDistSame",     &maxdist_same);
    vector <vector<double>> maxdist_same_big (2); 
    tfriend->Branch("lc_MaxDistSameBig",  &maxdist_same_big);
    vector <double> maxdist_diff; 
    tfriend->Branch("lc_MaxDistDiff",     &maxdist_diff); 
    vector <double> maxdist_diff_net; 
    tfriend->Branch("lc_MaxDistDiffNet",  &maxdist_diff_net);
    vector <double> maxdist_diff_big; 
    tfriend->Branch("lc_MaxDistDiffBig",  &maxdist_diff_big); 
    vector <double> maxdist_diff_net_big; 
    tfriend->Branch("lc_MaxDistDiffNetBig",&maxdist_diff_net_big); 
    vector <double> maxdist_centre; 
    tfriend->Branch("lc_MaxDistCentre",   &maxdist_centre);  
    vector <double> maxdist_centre_big; 
    tfriend->Branch("lc_MaxDistCentreBig",&maxdist_centre_big);  */