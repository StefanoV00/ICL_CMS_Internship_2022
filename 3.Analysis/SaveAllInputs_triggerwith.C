/// \author Stefano Veroni
/// File to read following information from a root TTree to .txt.
///******************************************************
///  - General   : Etot, energy, layer
///  - Energy    : length, BCl, maxEnerLay, maxMultLay,
///                Edepo, mult (to be resized)
///  - #Photon   : maxEfr_lcs, maxEfr_lay, l2norm,r2norm, 
///                maxdist, area
///  - #Photon_pr: maxdist_pr, maxdistcentre_pr,
///                rad68_pr, rad90_pr, rad95_pr & ratios
///                acnvx_pr, astar_pr
///  - For both  : Nlcs
///  - True Inf  : Etot0, Etot1
// *****************************************************
/// Don't worry if for some reason it crashes: just change values
/// in start variables (lines 84-88) to reprise from last saved.

/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 3.Analysis
root SaveAllInputs_triggerwith.C -q
cd ../

*/

#include <iostream>
#include <iomanip>
#include <fstream>
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



void SaveAllInputs_triggerwith (bool save_txt = true) 
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
    
    //string rootname= "step3ticl_pt50_eta21_run0_FlatTracksters";//root file
    string rootname1 = "step3ticl_";
    string rootname2 = "_eta";
    string rootname3 = "_run";
    string rootname4 = "_FlatTracksters";
    
    vector <string> pts  = {"En20to200", "En40to400"};
    vector <string> etas = {"21"};
    int             Nrun = 600;

    int pathstart= 0;
    int istart   = 0;
    int ptstart  = 0;
    int etastart = 0;
    int runstart = 0;
    vector <int> runskip = {192, 279, 289, 320, 327, 483};
   
    int alsonorm = 0;
    string dirname = "ticlTree";//"none" if no directory in between
    string treename= "TSTree_SimFromCP";
    string cluename= "TSTree_CLUE3D3";
    
    string savepath = "../myData/";
    string ext      = ".txt";



////////////////////////////////////////////////////////////////////////
//////1. PREPARE SAVING FILES //////////////////////////////////////////
/////__________________________________________________________________
    // Prepare CSV Files where to save data
    // Use ";" as separator
    vector <ofstream> savefiles;savefiles.resize(6);
    vector <int>      saveNs = {100, 150, 200};
    saveNs.insert( saveNs.end(), saveNs.begin(), saveNs.end() );
    vector <string>   savenames = {savepath + treename + "100max" + ext,
                                   savepath + treename + "200max" + ext,
                                   savepath + treename + "200min" + ext,
                                   savepath+treename+"100max"+"_norm"+ext,
                                   savepath+treename+"200max"+"_norm"+ext,
                                   savepath+treename+"200min"+"_norm"+ext}; 
    if (save_txt && pathstart == 0 && istart == 0 && runstart == 0)
    {
        for (int si = 0;si<3*(1+alsonorm);si++)
        {
            remove (savenames[si].c_str());
            savefiles[si].open(savenames[si]);
            
            savefiles[si]<<"DeltaCPVtxXY;";
            savefiles[si]<<"DeltaLayerStart;";
            
            savefiles[si]<<"TotE;";
            for (int layn = 1;layn<=saveNs[si];layn++)//N biggestenergy LCs
            {savefiles[si]<<"E" +to_string(layn)+";";
             savefiles[si]<<"LN"+to_string(layn)+";";}
            
            for (int layn = 1;layn<=15;layn++)//First 15 lays from start 
            {savefiles[si]<<"EnerInLayer"  + to_string(layn)+";";
             savefiles[si]<<"NMultInLayerN"+ to_string(layn)+";";}
            savefiles[si]<<"Lenght(E-Cal);";
            savefiles[si]<<"BC_layerN;"    ;
            savefiles[si]<<"EnerMaxLayerN;";
            savefiles[si]<<"MultMaxLayerN;";

            savefiles[si]<<"MaxLCEfract;"   ;
            savefiles[si]<<"MaxLayEfract;"  ;
            savefiles[si]<<"2ndLMomentNorm;";
            savefiles[si]<<"2ndTMomentNorm;";
            for (int layn = 1;layn<=15;layn++)//First 15 lays from start 
            {savefiles[si]<<"AConvex"+ to_string(layn)+";";
             savefiles[si]<<"MaxDist"+ to_string(layn)+";";}

            savefiles[si]<<"ProjMaxDist;"      ;
            savefiles[si]<<"ProjMaxDistCentre;";
            savefiles[si]<<"ProjACnvx;"        ;
            savefiles[si]<<"ProjAstar;"        ;
            savefiles[si]<<"ProjRad68;"        ;
            savefiles[si]<<"ProjRad90;"        ;
            savefiles[si]<<"ProjRad95;"        ;
            savefiles[si]<<"ProjRatio6890;"    ;
            savefiles[si]<<"ProjRatio6895;"    ;
            savefiles[si]<<"ProjRatio9095;"    ;

            for (int lcs = 0;lcs<=14;lcs++)//First 11 NLCs values 
            {savefiles[si]<<"NLaysWNLCs"+to_string(lcs)+";";
             savefiles[si]<<"NLaysWNLCssmaller"+to_string(lcs)+";";}
            savefiles[si]<<"TotNLCs;";

            savefiles[si]<<"Score1ph;" ;
            savefiles[si]<<"Score2ph;" ;
            savefiles[si]<<"Etot0frac;";
            savefiles[si]<<"Etot1frac"  ;

            savefiles[si]<<endl;
        }
    }
    else if (save_txt)
        for (int si = 0; si<3*(1+alsonorm) ;si++)
            savefiles[si].open(savenames[si], ios::app);




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
    Nrun *= (i_p==paths.size()-1 && pt_i==0 && e_i==0) ? 6 : 1;
    for (int run = runstart; run < Nrun; run++)
    {
        if (contains(runskip, run)) continue;
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
            TTree* tclue = (TTree*)dir->Get(cluename.c_str());
            //tree->Scan("lc_energy:lc_tsMult");
            
            // Get Entries Number
            int Nentries = tree -> GetEntries();

            vector <double>* true_emEs = 0;
            tree->SetBranchAddress("ts_emEnergy", &true_emEs);
            vector <double>* cp_vtxx = 0;
            vector <double>* cp_vtxy = 0;
            tree->SetBranchAddress("cp_vtxX", &cp_vtxx);
            tree->SetBranchAddress("cp_vtxY", &cp_vtxy);
            vector <int>* glayerfirst = 0;
            vector <int>* glayerlast  = 0;
            tree->SetBranchAddress("ts_firstLayer", &glayerfirst);
            tree->SetBranchAddress("ts_lastLayer",  &glayerlast); 
            
            vector <double>* clue_emEs = 0;
            tclue->SetBranchAddress("ts_emEnergy", &clue_emEs);
            vector <double>* energy_all= 0;
            vector <double>* x_all = 0;
            vector <double>* y_all = 0;
            vector <double>* z_all = 0;
            vector <int>* layer_all= 0;
            tclue->SetBranchAddress("lc_energy", &energy_all);
            tclue->SetBranchAddress("lc_x",     &x_all);
            tclue->SetBranchAddress("lc_y",     &y_all);
            tclue->SetBranchAddress("lc_z",     &z_all);
            tclue->SetBranchAddress("lc_layer", &layer_all);
            vector <int>* layerfirst = 0;
            vector <int>* layerlast  = 0;
            tclue->SetBranchAddress("ts_firstLayer", &layerfirst);
            tclue->SetBranchAddress("ts_lastLayer",  &layerlast);
            
            

            ////////////////////////////////////////////////////////////
            //4. Loop: get event, find quantities, save data, fill hists
            //__________________________________________________________
            auto isecal   = [](int layN)     {return layN  <=28;};
            auto distinct = [](double tsmult){return tsmult<= 2;};
            for (int i = 0;i<Nentries;i++)
            { 
                tree  -> GetEntry(i);
                tclue -> GetEntry(i);
                if (i%250 == 0) cout<<i<<endl;

                ////////////////////////////////////////////////////////
                // 4.1 PREPARE DATA
                //______________________________________________________    
                
                //True info for statistics
                double dx = 0;       double dy = 0; 
                double DeltaVtx = 0; double DeltaLayer = 0; 
                int lfirst = 0;      int llast = 0;
                if (glayerfirst->size()==2)
                {
                    // Take deltaVTx
                    dx = cp_vtxx->at(0) - cp_vtxx->at(1);
                    dy = cp_vtxy->at(0) - cp_vtxy->at(1);
                    DeltaVtx = sqrt(dx*dx + dy*dy);
                    
                    // Remove differentiation between photons
                    DeltaLayer = abs(glayerfirst->at(0)-glayerfirst->at(1));
                    lfirst = (glayerfirst->at(0)<glayerfirst->at(1))?
                                glayerfirst->at(0):glayerfirst->at(1) ;
                    llast  = (glayerlast ->at(0)>glayerlast ->at(1))?
                                glayerlast ->at(0):glayerlast->at(1)  ;
                }
                else if (glayerfirst->size()==1)
                {
                    lfirst = glayerfirst->at(0);
                    llast  = glayerlast->at(0);
                }

                //Clue info for post-trigger analysis
                lfirst = layerfirst->at(0);
                llast  = layerlast ->at(0);

                // Remove everything that is NOT E-Cal info
                vector<double> energy=getAwhereB(*energy_all,*layer_all,
                                                               isecal);
                vector<double> x     =getAwhereB(*x_all     ,*layer_all,
                                                               isecal);
                vector<double> y     =getAwhereB(*y_all     ,*layer_all,
                                                               isecal);
                vector<double> z     =getAwhereB(*z_all     ,*layer_all,
                                                               isecal);
                vector<int>    layer =getAwhereB(*layer_all ,*layer_all,
                                                               isecal);
                llast = (llast<28) ? llast : 28;

                // Get True Energies (the true info)
                double Etot0 = true_emEs->at(0);
                double Etot1 = (true_emEs->size()>1)
                               ? true_emEs->at(1)
                               : 0;
                double Etot = Etot0 + Etot1;
                if (Etot0 < Etot1) //Etot0 must be bigger
                {double a = Etot0*1.;Etot0 = Etot1*1.;Etot1 = a*1.;}              
                
                // Sort everything REVERSE-energy-wise
                x     = sortAwithB(x    , energy, true);
                y     = sortAwithB(y    , energy, true);
                z     = sortAwithB(z    , energy, true);
                layer = sortAwithB(layer, energy, true);
                sort(energy.rbegin(), energy.rend());

                // Get tot #LCs
                int    Nlcs  = energy.size();
                

                // Get Projections and shift layer number to start
                vector <double> xpr, ypr;
                for (int j = 0; j<Nlcs; j++)
                {xpr.push_back( x[j]/z[j] );
                 ypr.push_back( y[j]/z[j] );
                layer[j] -= lfirst-1;} //First layer is 1

                // Check everything is alright
                if (Nlcs != energy.size() || Nlcs != layer.size() 
                    || Nlcs != x.size()   || Nlcs != y.size() 
                    || Nlcs != z.size())
                {cout<<Nlcs    <<","
                     <<x.size()<<","
                     <<y.size()<<","
                     <<z.size()<<endl;
                  return;}
                
                //******************************************************
                // Now I have: 
                //  - General : Etot, energy, layer
                //  - For both: Nlcs
                //  - True Inf: Etot0, Etot1
                // *****************************************************


                ////////////////////////////////////////////////////////
                // 4.2 TAKE USEFUL INFORMATION - Energy
                //______________________________________________________
                
                int  length = llast - lfirst +1;
                double BCl  = mymean(layer, energy);
                vector <double> Edepo;
                vector <int>    mult;
                for (int lay = 1; lay<=15; lay++)
                {
                    vector<double> energy_lay = getAwhereB(energy,layer, 
                                                                  lay);
                    double Edepoj = myvecsum(energy_lay);
                    Edepo.push_back(Edepoj);
                    mult .push_back(energy_lay.size());
                } 
                //Edepo.resize(15, 0);
                //mult .resize(15, 0);
                int maxEnerLay = argmax(Edepo)+1;
                int maxMultLay = argmax(mult)+1;
                //******************************************************
                // Now I have: 
                //  - General : Etot, energy, layer
                //  - Energy  : length, BCl, maxEnerLay, maxMultLay,
                //              Edepo, mult (to be resized)
                //  - For both: Nlcs
                //  - True Inf: Etot0, Etot1
                // *****************************************************



                ////////////////////////////////////////////////////////
                // 4.3 TAKE USEFUL INFORMATION - PHOTON # 
                //______________________________________________________
                double maxEfr_lcs = vecmax(energy) / Etot;
                double maxEfr_lay = vecmax(Edepo)  / Etot;
                
                // get normalised moments
                double BCx  = mymean(x, energy);
                double BCy  = mymean(y, energy);
                double BCz  = mymean(z, energy);
                double ml2 = 0; double mlmax = 0;
                double mr2 = 0; double mrmax = 0;
                for (int j = 0; j<2; j++)
                {
                    double Ej  = energy[j];
                    mlmax += 100.*Ej;
                    mrmax += 16. *Ej;
                }
                for (int j = 2; j<Nlcs; j++)
                {
                    double Ej  = energy[j];
                    double lj2 = (z[j] - BCz) * (z[j] - BCz) * Ej;
                    double rj2 = ( (x[j] - BCx) * (x[j] - BCx)
                                +(y[j] - BCy) * (y[j] - BCy) ) * Ej;
                    ml2 += lj2*Ej;
                    mr2 += rj2*Ej;
                }
                double l2norm = ml2 / (ml2 + mlmax);
                double r2norm = mr2 / (mr2 + mrmax);

                // get areas and distances
                vector <double> maxdist;
                vector <double> area;
                for (int lay = 1; lay <= 15; lay++)
                {
                    vector <double> xlay = getAwhereB(x, layer, lay);
                    vector <double> ylay = getAwhereB(y, layer, lay);
                    
                    double maxdist_lay = 0;
                    for (int j = 0; j<xlay.size(); j++)
                    {
                        double xj = xlay[j];double yj = ylay[j];
                        double delta_layk = 0;
                        for (int k = j+1; k<xlay.size(); k++)
                        {
                            double dx = xlay[k]-xj;
                            double dy = ylay[k]-yj;
                            delta_layk=sqrt( dx*dx + dy*dy);
                            if (delta_layk > maxdist_lay) 
                                maxdist_lay=delta_layk*1.;
                        }
                    }
                    maxdist.push_back(maxdist_lay);
                    area.push_back(Shoelace(ConvexHull(xlay, ylay)));
                }

                //******************************************************
                // Now I have: 
                //  - General   : Etot, energy, layer
                //  - Energy    : length, BCl, maxEnerLay, maxMultLay,
                //                Edepo, mult (to be resized)
                //  - #Photon   : maxEfr_lcs, maxEfr_lay, l2norm,r2norm, 
                //                maxdist, area
                //  - #Photon_pr:
                //  - For both  : Nlcs
                //  - True Inf: Etot0, Etot1
                // *****************************************************



                ////////////////////////////////////////////////////////
                // 4.4 TAKE USEFUL INFORMATION - PHOTON# (Projections)
                //______________________________________________________
                double centrex = mymean(xpr, energy);
                double centrey = mymean(ypr, energy);
                double maxdist_pr = 0;
                double maxdistcentre_pr = 0;
                double rad68_pr = 0;
                double rad90_pr = 0;
                double rad95_pr = 0;
                double ratio6890_pr = 0;
                double ratio6895_pr = 0;
                double ratio9095_pr = 0;
                vector <double> cdistances;
                vector <double> esum(Nlcs);
                
                for (int j = 0; j<Nlcs; j++)
                {
                    double xj = xpr[j];
                    double yj = ypr[j];
                    double delta_jk = 0;
                    for (int k = j+1; k<Nlcs; k++)
                    {
                        double dx = xpr[k]-xj;
                        double dy = ypr[k]-yj;
                        delta_jk = sqrt( dx*dx + dy*dy);
                        if (delta_jk > maxdist_pr) maxdist_pr = delta_jk;
                    }
                    double deltacentre_j =sqrt((xj-centrex)*(xj-centrex)
                                            +(yj-centrey)*(yj-centrey));
                    
                    cdistances.push_back(deltacentre_j);
                    if (deltacentre_j > maxdistcentre_pr)
                        maxdistcentre_pr=deltacentre_j;
                }
        
                //Get Radii 68, 90, 95
                vector <double> enes = sortAwithB(energy, cdistances);
                sort(cdistances.begin(), cdistances.end());
                enes[0] /= Etot;
                partial_sum(enes.begin(), enes.end(), esum.begin(), 
                    [Etot](double a, double b){b/=Etot;return a + b;});
                for (int j = 0; j<Nlcs; j++)
                    {if (esum[j] >= 0.68 || j == Nlcs-1) 
                        {rad68_pr = cdistances[j];
                        for (int jj = j*1; jj<Nlcs; jj++)
                            {if (esum[jj] >= 0.90 || jj == Nlcs-1) 
                                {rad90_pr = cdistances[jj];
                                for (int jjj = jj*1; jjj<Nlcs; jjj++)
                                    {if (esum[jjj]>=0.95 || jjj==Nlcs-1) 
                                        {rad95_pr = cdistances[jjj];
                                        break;}
                                    }
                                break;}
                            }
                        break;}
                    }
                ratio6890_pr = rad68_pr/rad90_pr;
                ratio6895_pr = rad68_pr/rad95_pr;
                ratio9095_pr = rad90_pr/rad95_pr;

                // Get Areas, projected
                vector<vector<double>>chull =ConvexHull(xpr, ypr);
                vector<vector<double>>shull =StarHull(xpr, ypr, 
                                                      centrex, centrey);
                double acnvx_pr = Shoelace(chull);
                double astar_pr = Shoelace(shull[0], shull[1]);

                //******************************************************
                // Now I have: 
                //  - General   : Etot, energy, layer
                //  - Energy    : length, BCl, maxEnerLay, maxMultLay,
                //                Edepo, mult (to be resized)
                //  - #Photon   : maxEfr_lcs, maxEfr_lay, l2norm,r2norm, 
                //                maxdist, area
                //  - #Photon_pr: maxdist_pr, maxdistcentre_pr,
                //                rad68_pr, rad90_pr, rad95_pr & ratios
                //                acnvx_pr, astar_pr
                //  - For both  : Nlcs
                //  - True Inf  : Etot0, Etot1
                // *****************************************************
                
                
                ////////////////////////////////////////////////////////
                // 4.5 TAKE USEFUL INFORMATION - BOTH
                //______________________________________________________
                
                vector <int> nlcs_distr (15);
                vector <int> nlcs_cumd (15);
                int lcount = count(mult.begin(), mult.end(), 0);
                nlcs_distr[0] = lcount;
                nlcs_cumd [0] = lcount;

                for (int nlcs = 1; nlcs<=14; nlcs++)
                    {
                        lcount = count(mult.begin(), mult.end(), nlcs);
                        nlcs_distr[nlcs] = lcount;
                        nlcs_cumd [nlcs] = nlcs_cumd [nlcs-1] + lcount;
                    }
                //******************************************************
                // Now I have: 
                //  - General   : Etot, energy, layer
                //  - Energy    : Edepo, mult, 
                //                length, BCl, maxEnerLay, maxMultLay
                //  - #Photon   : maxEfr_lcs, maxEfr_lay, l2norm,r2norm, 
                //                maxdist, area
                //  - #Photon_pr: maxdist_pr, maxdistcentre_pr,
                //                rad68_pr, rad90_pr, rad95_pr & ratios
                //                acnvx_pr, astar_pr
                //  - For both  : nlcs_distr, nlcs_cumd, Nlcs
                //  - True Inf  : Etot0, Etot1
                // *****************************************************


                //Save in csv format
                if (save_txt)
                {   
                    // The Original Info
                    int si = min(Etot/100, 2.);
                    double score1 = (Etot0 == 0. || Etot1==0.);
                    double score2 = 1. - score1;
                    if (Nlcs<saveNs[si])
                        {for (int n = 0; n<saveNs[si]-Nlcs; n++)
                            {energy.push_back(0.);
                             layer .push_back(0.); }}
                    if (Edepo.size()<15)
                        {for (int n = 0; n<15-Edepo.size(); n++)
                            {Edepo   .push_back(0.);
                             mult    .push_back(0.);}}
                    if (area.size()<15)
                        {for (int n = 0; n<15-area.size(); n++)
                             area    .push_back(0.);
                             maxdist .push_back(0.);}
                    if (savefiles[si].is_open())
                    {  
                        savefiles[si]<<DeltaVtx<<";";
                        savefiles[si]<<DeltaLayer<<";";
                        
                        savefiles[si]<<Etot<<";";
                        for (int layn = 1;layn<=saveNs[si];layn++)
                            {savefiles[si]<<energy[layn-1]<<";";
                             savefiles[si]<<layer [layn-1]<<";";}

                        for (int layn = 1;layn<=15;layn++) 
                            {savefiles[si]<<Edepo[layn-1]<<";";
                            savefiles[si]<<mult [layn-1]<<";";}
                            savefiles[si]<<length      <<";";
                            savefiles[si]<<BCl         <<";";
                            savefiles[si]<<maxEnerLay  <<";";
                            savefiles[si]<<maxMultLay  <<";";
                        
                        savefiles[si]<<maxEfr_lcs<<";"  ;
                        savefiles[si]<<maxEfr_lay<<";"  ;
                        savefiles[si]<<l2norm    <<";";
                        savefiles[si]<<r2norm    <<";";
                        for (int layn = 1;layn<=15;layn++)
                            {savefiles[si]<<area[layn-1]   <<";";
                            savefiles [si]<<maxdist[layn-1]<<";";}
                        
                        savefiles[si]<<maxdist_pr      <<";";
                        savefiles[si]<<maxdistcentre_pr<<";";
                        savefiles[si]<<acnvx_pr        <<";";
                        savefiles[si]<<astar_pr        <<";";
                        savefiles[si]<<rad68_pr        <<";";
                        savefiles[si]<<rad90_pr        <<";";
                        savefiles[si]<<rad95_pr        <<";";
                        savefiles[si]<<ratio6890_pr    <<";";
                        savefiles[si]<<ratio6895_pr    <<";";
                        savefiles[si]<<ratio9095_pr    <<";";

                        for (int lcs = 0;lcs<=14;lcs++) 
                        {savefiles[si]<<nlcs_distr[lcs]<<";";
                        savefiles[si]<<nlcs_cumd[lcs] <<";";}
                        savefiles[si]<<Nlcs<<";";

                        savefiles[si]<<score1     <<";";
                        savefiles[si]<<score2     <<";";
                        savefiles[si]<<Etot0/Etot <<";";
                        savefiles[si]<<Etot1/Etot;

                        savefiles[si]<<endl;
                    }
                    else 
                    {
                        cout<< "At event: " << i<< 
                            ", unable to open csv file" << endl;
                    }


                    // The one with "Normalised" Info
                    if (alsonorm)
                    {
                    si += 3;
                    if (Nlcs<saveNs[si-3])
                        {for (int n = 0; n<saveNs[si-3]-Nlcs; n++)
                            {energy.push_back(0.);
                            layer.push_back(0.); }}
                    if (savefiles[si].is_open())
                    {  
                        savefiles[si]<<DeltaVtx<<";";
                        savefiles[si]<<DeltaLayer<<";";
                        
                        savefiles[si]<<Etot/1000.<<";";
                        for (int layn = 1;layn<=saveNs[si-3];layn++)
                            {savefiles[si]<<energy[layn-1]/Etot<<";";
                             savefiles[si]<<layer [layn-1]/Etot<<";";}

                        for (int layn = 1;layn<=15;layn++) 
                            {savefiles[si]<<Edepo[layn-1]/Etot<<";";
                            savefiles[si]<<mult [layn-1]/1000.<<";";}
                            savefiles[si]<<length/28.    <<";";
                            savefiles[si]<<BCl/28.       <<";";
                            savefiles[si]<<maxEnerLay/28.<<";";
                            savefiles[si]<<maxMultLay/28.<<";";
                        
                        savefiles[si]<<maxEfr_lcs<<";"  ;
                        savefiles[si]<<maxEfr_lay<<";"  ;
                        savefiles[si]<<l2norm    <<";";
                        savefiles[si]<<r2norm    <<";";
                        for (int layn = 1;layn<=15;layn++)
                            {savefiles[si]<<area[layn-1]/1000.  <<";";
                            savefiles [si]<<maxdist[layn-1]/100.<<";";}
                        
                        savefiles[si]<<maxdist_pr      <<";";
                        savefiles[si]<<maxdistcentre_pr<<";";
                        savefiles[si]<<acnvx_pr        <<";";
                        savefiles[si]<<astar_pr        <<";";
                        savefiles[si]<<rad68_pr        <<";";
                        savefiles[si]<<rad90_pr        <<";";
                        savefiles[si]<<rad95_pr        <<";";
                        savefiles[si]<<ratio6890_pr    <<";";
                        savefiles[si]<<ratio6895_pr    <<";";
                        savefiles[si]<<ratio9095_pr    <<";";

                        for (int lcs = 0;lcs<=14;lcs++) 
                        {savefiles[si]<<nlcs_distr[lcs]/28.<<";";
                        savefiles[si]<<nlcs_cumd[lcs]/28. <<";";}
                        savefiles[si]<<Nlcs/1000.<<";";

                        savefiles[si]<<score1     <<";";
                        savefiles[si]<<score2     <<";";
                        savefiles[si]<<Etot0/Etot <<";";
                        savefiles[si]<<Etot1/Etot;

                        savefiles[si]<<endl;
                    }
                    else 
                    {
                        cout<< "At event: " << i<< 
                            ", unable to open csv file" << endl;
                    }
                    }
                    
                       
                }
            }//end of loop over entries
            file->Close();
            printf("Run %d saved \n", run);
        } //end of  (if dir)
        } //end of (if file)
    } // end of loop over runs
    } // end of loop over etas
    } // end of loop over pts
    runstart = 0;
    ptstart = 0;
    }// end of paths 
    cout<<"------ FINISHED ------"<<endl;

}//end of main