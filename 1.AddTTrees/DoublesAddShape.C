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
 - Layer # with maximum # of clusters in total
 - Layer # since start of first shower with maximum # of clusters

 - max distance between lcs per layer
 - tag max distance ij: i = 0(1): for same(different) photons
                        j = 0(1): for at least one(no) lc with rechits<3
    This gives combinations: 00=0, 01=1, 10, 11. Plus, -1 if no event.
 - max distance from barycentre 
 - max distance from barycentre excluding small clusters (<3)

 - area defined by all entries
*/


void DoublesAddShape () 
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
    tfriendE -> SetBranchAddress("lc_centreXFromStart", &centrex);
    tfriendE -> SetBranchAddress("lc_centreYFromStart", &centrey);
    tfriendE -> SetBranchAddress("lc_BCx",   &BCx);
    tfriendE -> SetBranchAddress("lc_BCy",   &BCy);
    tfriendE -> SetBranchAddress("lc_BCz",   &BCz);
    tfriendE -> SetBranchAddress("lc_BCl",   &BCl);
    tfriendE -> SetBranchAddress("lc_BCx_pr", &centrex_pr);
    tfriendE -> SetBranchAddress("lc_BCy_pr", &centrey_pr);
    //tfriend -> Print();



    ////////////////////////////////////////////////////////////////////
    // SET NEW BRANCHES ////////////////////////////////////////////////
    //__________________________________________________________________

    //Layers - Hits info, use lc
    vector <int> totmult;
    vector <int> totmultbig;
    tfriend->Branch("lc_TotMult",    &totmult);
    tfriend->Branch("lc_MultBig", &totmultbig);
    // FOLLOWING ACTUALLY WRONG, NOT WHAT I WANTED
    vector <int> totcummult;
    // tfriend->Branch("lc_totCumMult", &totcummult);  
    // tfriend->Branch("lc_gCumMult",   &gcummult); 
 
    int LmaxMultTot;
    int LmaxMultTot_start;
    tfriend->Branch("lc_layerMaxMult",          &LmaxMultTot); 
    tfriend->Branch("lc_layerMaxMultFromStart", &LmaxMultTot_start); 


    // Longitudinal & Transverse 2nd Moments 
    double lambda2 = 0;
    double erre2   = 0;
    double lnorm   = 0;
    double rnorm   = 0;
    tfriend->Branch("lc_2ndLongMoment", &lambda2);
    tfriend->Branch("lc_2ndTranMoment", &erre2);
    tfriend->Branch("lc_2ndLongMomentNorm", &lnorm);
    tfriend->Branch("lc_2ndTranMomentNorm", &rnorm);  

    

    //Distances between hits - from start
    vector <double> maxdist_start; 
    tfriend->Branch("lc_maxDistFromStart",         &maxdist_start);
    vector <int> maxdist_tag_start; 
    tfriend->Branch("lc_maxDistTagFromStart",      &maxdist_tag_start);
    vector <double> maxdist_big_start; 
    tfriend->Branch("lc_maxDistBigFromStart",      &maxdist_big_start);
    vector <double> maxdist_centre_start; 
    tfriend->Branch("lc_maxDistCentreFromStart",   &maxdist_centre_start);  
    vector <double> maxdist_centre_big_start; 
    tfriend->Branch("lc_maxDistCentreBigFromStart",&maxdist_centre_big_start);  
    

    // Areas of layers' sections (as convex hulls of detections)
    vector <double> areatot_start;
    tfriend->Branch("lc_areaFromStart", &areatot_start);
    vector <vector<double>> areag_start (2);


    //Distances between hits - projections
    double maxdist_pr; 
    tfriend->Branch("lc_maxDist_pr",          &maxdist_pr);
    int maxdist_tag_pr; 
    tfriend->Branch("lc_maxDistTag_pr",       &maxdist_tag_pr);
    double maxdist_big_pr; 
    tfriend->Branch("lc_maxDistBig_pr",       &maxdist_big_pr);
    double maxdist_centre_pr; 
    tfriend->Branch("lc_maxDistCentre_pr",    &maxdist_centre_pr);  
    double maxdist_centre_big_pr; 
    tfriend->Branch("lc_maxDistCentreBig_pr", &maxdist_centre_big_pr); 


    // Areas of projections (as convex hulls of all detections)
    double Atot_pr;
    double AtotBig_pr;
    tfriend->Branch("lc_A_pr",    &Atot_pr);
    tfriend->Branch("lc_ABig_pr", &AtotBig_pr);
    double Astartot_pr;
    double AstartotBig_pr;
    tfriend->Branch("lc_Astar_pr",    &Astartot_pr);
    tfriend->Branch("lc_AstarBig_pr", &AstartotBig_pr);


    // Ratio Distance2-Areas of projections
    double Atodm2_ratio_pr;
    tfriend->Branch("lc_AtoDmax2Ratio_pr",    &Atodm2_ratio_pr);
    double Atodm2_ratioBig_pr;
    tfriend->Branch("lc_AtoDmax2BigsRatio_pr",&Atodm2_ratioBig_pr);
    double Atodc2_ratio_pr;
    tfriend->Branch("lc_AtoDcen2Ratio_pr",    &Atodc2_ratio_pr);
    double Atodc2_ratioBig_pr;
    tfriend->Branch("lc_AtoDcen2BigsRatio_pr",&Atodc2_ratioBig_pr);


    // Radii with 68, 90 & 95% energy
    double radius68_pr;
    tfriend->Branch("lc_Radius68_pr", &radius68_pr);
    double radius90_pr;
    tfriend->Branch("lc_Radius90_pr", &radius90_pr);
    double radius95_pr;
    tfriend->Branch("lc_Radius95_pr", &radius95_pr);
    double ratio9068_pr;
    tfriend->Branch("lc_Ratio9068_pr", &ratio9068_pr);
    double ratio9568_pr;
    tfriend->Branch("lc_Ratio9568_pr", &ratio9568_pr);
    double ratio9590_pr;
    tfriend->Branch("lc_Ratio9590_pr", &ratio9590_pr);  


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
            totmultbig.clear();
            //totcummult.clear();
            //gcummult[0].clear();
            //gcummult[1].clear();

            maxdist_start.clear();  
            maxdist_big_start.clear(); 

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
        double er2 = 0; double ermax = 0;

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
            } 
            else
            {   
                E1tot += Ej;
            }    
        }
        lambda2    /= Etot; erre2    /= Etot;

        lnorm = el2 / (el2 + elmax);
        rnorm = er2 / (er2 + ermax);


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
                int countgbig = count(gnbig.begin(), gnbig.end(), 0);
                //int cumcountg = (ln == first0)? 
                //                countg:countg + gcummult[0][ln-first0-1];
                //gcummult[0].push_back(cumcountg);
            }
            if (ln>=first1 && ln<=last1)
            {
                int countg = count(gn.begin(), gn.end(), 1);
                int countgbig = count(gnbig.begin(), gnbig.end(), 1);
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

            maxdist_centre_start      .push_back(soln[10]);
            maxdist_centre_big_start  .push_back(soln[11]);
            maxdist_tag_start         .push_back(soln[12]);

            areatot_start.push_back(Shoelace(ConvexHull(xn, yn)));

        }

        // ____________________________________________________________
        // Distances in projected space _______________________________
        XYVector cpvtx0_pr = XYVector(cpx->at(0)/cpz->at(0), 
                                      cpy->at(0)/cpz->at(0));
        XYVector cpvtx1_pr = XYVector(cpx->at(1)/cpz->at(1), 
                                      cpy->at(1)/cpz->at(1));
        XYVector centre_pr  = XYVector(centrex_pr, centrey_pr);
        vector <double> soln = get_max_distances(*xpr, *ypr, 
                                                 *which_g,   *nrechits, 
                                                 cpvtx0_pr,  cpvtx1_pr,
                                                 centre_pr);
        maxdist_pr             = soln[0];
        maxdist_big_pr         = soln[1];

        maxdist_centre_pr      = soln[10];
        maxdist_centre_big_pr  = soln[11];

        maxdist_tag_pr         = soln[12];
        //cout<<"get projected max done"<<endl;

        // ____________________________________________________________
        // Radii with % of energy (projected space) ____________________
        vector<double> radii =get_radii(*xpr, *ypr, *which_g, *energy,
                                         centre_pr);
        radius68_pr  = radii[0];
        radius90_pr  = radii[1];
        radius95_pr  = radii[2];

        ratio9068_pr  = radius90_pr / radius68_pr;
        ratio9568_pr  = radius95_pr / radius68_pr;
        ratio9590_pr  = radius95_pr / radius90_pr;

        // ____________________________________________________________
        // Areas of projections _______________________________________

        Atot_pr  = Shoelace(ConvexHull(*xpr, *ypr));

        auto shull = StarHull(*xpr, *ypr,
                              centrex_pr, centrey_pr);
        Astartot_pr  = Shoelace(shull[0], shull[1]);

        //Big
        vector <double> xprbig =getAwhereB(*xpr    ,*nrechits, bighits);
        vector <double> yprbig =getAwhereB(*ypr    ,*nrechits, bighits);
        vector <int>    gbig   =getAwhereB(*which_g,*nrechits, bighits);
        AtotBig_pr = Shoelace(ConvexHull(xprbig, yprbig));

        shull = StarHull(xprbig, yprbig,
                         centrex_pr, centrey_pr);
        AstartotBig_pr  = Shoelace(shull[0], shull[1]);


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


        // Layer number with maximum number of hits ____________________
        vector <int> mult_both = sum_portions(*mult, 2); //sum halves
        LmaxMultTot = argmax(mult_both) + 1;
        if (LmaxMultTot > 50) {LmaxMultTot -= 50;}
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