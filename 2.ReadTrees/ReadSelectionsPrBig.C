/// \author Stefano Veroni
// DO NOT CARE ABOUT "MUST HAVE CONSTANT VALUE" ERRORS

/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat
 
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
//#include "TApplication.h"
//#include <TArrow.h>
//#include <TBenchmark.h>
#include "TBox.h"
#include "TButton.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TColor.h"
#include "TDatabasePDG.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include <TFrame.h>
#include "TGraphErrors.h"
#include "TH1.h"
#include <TH2.h>
//#include <THStack.h>
#include <TInterpreter.h>
#include "TKey.h"
//#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
//#include "TLorentzVector.h"
//#include "TLorentzRotation.h"
#include "TMath.h"
#include <TNtuple.h>
//#include "TPad.h"
//#include "TParticlePDG.h"
//#include <TPave.h>
//#include <TPaveText.h>
//#include "TPostScript.h"
//#include <TProfile.h>
//#include "TRandom.h"
//#include <TRandom3.h>
#include "TROOT.h"
//#include "TRotation.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"
//#include "TVector3.h"
//#include "TVirtualPad.h"
//#include "TVirtualPS.h"
//#include "Riostream.h"
//#include <RDataFrame.hxx>

using namespace std;

/////////////////////////////////////////////////////////////////////////
//-1. Functions' Definitions/////////////////////////////////////////////
//______________________________________________________________________
#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"




void ReadSelectionsPrBig() 
// ALL QAUNTITIES INVOLVED ARE "BIG"-RELATED ALTHOUGH NOT EXPLICITLY
// DCLARED WITH THE USUAL "BIG" SUFFIX
{
    ////////////////////////////////////////////////////////////////////
    //0. Set Variables   ///////////////////////////////////////////////
    //__________________________________________________________________
    string path    = "../myTTrees/";
//none if none
    string rootname= "step3ticl_pt50_eta21_run0_FlatTracksters";//root file
    string friendrootname= "treeFriends"; //root file with friends
    string dirname = "ticlTree"; //"none" if no directory in between
    string treename= "TSTree_SimFromCP";

    string selname = "treeSelections";
    //string selvers = "Good";
    string seldir  = "CloseXYSameStart+-";
    string seltree = "TSTree_SimFromCP";
    string selrootname = selname;//+selvers;

    int Nevents = 10;  //# events to consider, 0 if all
    float ievent = INFINITY;  //event of which display info, INFINITE if all
    
    ////////////////////////////////////////////////////////////////////
    //1. Get TTree & Branches///////////////////////////////////////////
    //__________________________________________________________________
    //NOTE: FOR EACH ENTRY, LC_QUANTITY IS A VECTOR OF "DETECTIONS"
    //FOR EXAMPLE: lc_x AND lc_layer WILL HAVE N ENTRIES, ONE PER 
    //CLUSTER
    TFile* file  = open_file(path, selrootname);
    TTree* tree0 = get_tree (file, seldir, seltree);
    TTree* treec = get_tree (file, seldir, seltree+"_Coords");
    TTree* treed = get_tree (file, seldir, seltree+"_Deltas");
    TTree* treee = get_tree (file, seldir, seltree+"_Energies");
    TTree* trees = get_tree (file, seldir, seltree+"_Shape");
    //treee -> Scan("lc_centreTotX_pr:lc_gCentreX_pr");
    //tree0 -> Scan("which_g:energy");
    //treec -> Scan("lc_x_pr:lc_y_pr");
    
    // Get Entries Number
    int Nentries = tree0 -> GetEntries();
    Nevents = (Nentries>Nevents && Nevents>0) ? Nevents : Nentries;
    cout << "The TTree has total "<< Nentries <<" entries. "<< endl;
    cout << "We are working with "<< Nevents  <<" events."  << endl;
    
    //POINTERS MUST BE INITIALISED
    // FROM TREE0
    vector <double>* energy= 0;
    vector <int>* which_g = 0;
    vector <int>* layerN  = 0;
    vector <int>* nrechits= 0;
    tree0->SetBranchAddress("lc_energy",   &energy);
    tree0->SetBranchAddress("lc_layer",    &layerN);
    tree0->SetBranchAddress("lc_nrechits", &nrechits );
    tree0->SetBranchAddress("lc_TSidx",    &which_g );
    
    //FROM TREEC
    vector <double>* xpr = 0;treec->SetBranchAddress("lc_x_pr",&xpr);
    vector <double>* ypr = 0;treec->SetBranchAddress("lc_y_pr",&ypr);
    //FROM TREEE
    vector <double>* gCentreX_pr = 0;
    vector <double>* gCentreY_pr = 0;
    double centreTotX_pr = 0;
    double centreTotY_pr = 0;
    treee->SetBranchAddress("lc_gCentreX_pr",   &gCentreX_pr);
    treee->SetBranchAddress("lc_gCentreY_pr",   &gCentreY_pr);
    treee->SetBranchAddress("lc_centreTotX_pr", &centreTotX_pr);
    treee->SetBranchAddress("lc_centreTotY_pr", &centreTotY_pr);
    //FROM TREES
    vector <double>* maxDistGCen_pr = 0;
    double  maxDistCen_pr = 0;
    vector <double>* maxDistG_pr = 0;
    double maxDist_pr = 0;
    vector <double>* Ag_pr = 0;
    vector <double>* AgBig_pr = 0;
    double Atot_pr = 0;
    double AtotBig_pr = 0;
    trees->SetBranchAddress("lc_maxDistGCentreBig_pr", &maxDistGCen_pr);
    trees->SetBranchAddress("lc_maxDistCentreBig_pr",  &maxDistCen_pr);
    trees->SetBranchAddress("lc_maxDistSameBig_pr",    &maxDistG_pr);
    trees->SetBranchAddress("lc_maxDistBig_pr",        &maxDist_pr);
    trees->SetBranchAddress("lc_maxDistGCentreBig_pr", &maxDistGCen_pr);

    // First a lambda function needed for big hits
    auto bighits = [] (int hit) {return hit>=3;};

    ////////////////////////////////////////////////////////////////////
    //2. Loop: get event, find quantities, save data, fill hists ///////
    //__________________________________________________________________
    //NOTE: RDataFrame.Define() could define new columns theta, etc...
    for (int i = 0; i<Nevents; i++)
    { 
        tree0  ->GetEntry(i);
        treec -> GetEntry(i);
        treed -> GetEntry(i);
        treee -> GetEntry(i);
        trees -> GetEntry(i);
        
        int N = xpr->size();
        vector <double> xpr0 = getAwhereB(*xpr     , *which_g, 0);
        vector <double> ypr0 = getAwhereB(*ypr     , *which_g, 0);
        vector <int> hits0   = getAwhereB(*nrechits, *which_g, 0);
        xpr0 = getAwhereB(xpr0, hits0, bighits);
        ypr0 = getAwhereB(ypr0, hits0, bighits);
        vector <double> xpr1 = getAwhereB(*xpr     , *which_g, 1);
        vector <double> ypr1 = getAwhereB(*ypr     , *which_g, 1);
        vector <int> hits1   = getAwhereB(*nrechits, *which_g, 1);
        xpr1 = getAwhereB(xpr1, hits1, bighits);
        ypr1 = getAwhereB(ypr1, hits1, bighits);


        vector <double> ene0 = getAwhereB(*energy, *which_g, 0);
        double E0 = myvecsum(ene0);
        vector <double> ene1 = getAwhereB(*energy, *which_g, 1);
        double E1 = myvecsum(ene1);
        double Etot = E0 + E1;


        double xalld = centreTotX_pr - maxDistCen_pr;
        double xallu = centreTotX_pr + maxDistCen_pr;
        double yalld = centreTotY_pr - maxDistCen_pr;
        double yallu = centreTotY_pr + maxDistCen_pr;

        double x0d = gCentreX_pr->at(0) - maxDistGCen_pr->at(0);
        double x0u = gCentreX_pr->at(0) + maxDistGCen_pr->at(0);
        double y0d = gCentreY_pr->at(0) - maxDistGCen_pr->at(0);
        double y0u = gCentreY_pr->at(0) + maxDistGCen_pr->at(0);

        double x1d = gCentreX_pr->at(1) - maxDistGCen_pr->at(1);
        double x1u = gCentreX_pr->at(1) + maxDistGCen_pr->at(1);
        double y1d = gCentreY_pr->at(1) - maxDistGCen_pr->at(1);
        double y1u = gCentreY_pr->at(1) + maxDistGCen_pr->at(1);

        vector <double> xs = {xalld, x0d, x1d, xallu, x0u, x1u};
        double xprmin = vecmin(xs); double xprmax = vecmax(xs);
        vector <double> ys = {yalld, y0d, y1d, yallu, y0u, y1u};
        double yprmin = vecmin(ys); double yprmax = vecmax(ys);
        double extens = (xprmax-xprmin>yprmax-yprmin)?
                         xprmax-xprmin:yprmax-yprmin;
        xprmin -= extens*0.2;
        yprmin -= extens*0.2;
        xprmax  = xprmin + extens*1.2; 
        yprmax  = yprmin + extens*1.2;


        TCanvas* c_prtot = new TCanvas("cAll", "cAll");
        double csize=c_prtot->GetWindowWidth();
        c_prtot->SetWindowSize(csize, csize);
        TH2D* h_prtot    = new TH2D("Projections-All (Big)",
                                    "Projections-All (Big)", 
                                    100, xprmin, xprmax,
                                    100, yprmin, yprmax);
        h_prtot->SetContour(99); //nicer, default is like 10
        h_prtot->GetXaxis()->SetTitle("x projection");
        h_prtot->GetYaxis()->SetTitle("y projection");
        h_prtot->GetZaxis()->SetTitle("# Energy");
        string dist = "Max Distance: "+to_string(maxDist_pr);
        const char* dist_text = dist.c_str();
        TText* textall0 = new TText(xprmin,  yprmax, dist_text);
        textall0 -> SetTextAlign(11);
        textall0 -> SetTextSize(0.025);
        string distcen = "Max Distance from Centre: "+to_string(maxDistCen_pr);
        const char* distcen_text = distcen.c_str();
        TText* textall1 = new TText(xprmax,  yprmax, distcen_text);
        textall1 -> SetTextAlign(31);
        textall1 -> SetTextSize(0.025);

        string cen = "("+ to_string(centreTotX_pr)+", "
                        + to_string(centreTotY_pr) + ")";
        const char* cen_text = cen.c_str();
        TText* textall2 = new TText(xprmax, yprmax, cen_text);
        textall2->SetTextAlign(33);
        textall2->SetTextSize(0.025);

        TEllipse* elall = new TEllipse(centreTotX_pr, 
                                       centreTotY_pr,
                                       maxDistCen_pr,
                                       maxDistCen_pr);
        elall->SetFillStyle(0);



        TCanvas* c_pr0 = new TCanvas("cpr0", "cpr0");
        c_pr0->SetWindowSize(csize, csize);
        TH2D* h_pr0    = new TH2D("Projection-0th (Big)", 
                                  "Projection-0th (Big)", 
                                        100, xprmin, xprmax,
                                        100, yprmin, yprmax);
        h_pr0->SetContour(99); //nicer, default is like 10
        h_pr0->GetXaxis()->SetTitle("x projection");
        h_pr0->GetYaxis()->SetTitle("y projection");
        h_pr0->GetZaxis()->SetTitle("# Energy");
        dist = "Max Distance: "+to_string(maxDistG_pr->at(0));
        dist_text = dist.c_str();
        
        TText* text00 = new TText(xprmin, yprmax, dist_text);
        text00->SetTextAlign(11);
        text00->SetTextSize(0.025);
        distcen = "Max Distance from Centre: "+to_string(maxDistGCen_pr->at(0));
        distcen_text = distcen.c_str();
        TText* text01 = new TText(xprmax, yprmax, distcen_text);
        text01->SetTextAlign(31);
        text01->SetTextSize(0.025);

        cen = "("+ to_string(gCentreX_pr->at(0))+", "
                 + to_string(gCentreY_pr->at(0)) + ")";
        cen_text = cen.c_str();
        TText* text02 = new TText(xprmax, yprmax, cen_text);
        text02->SetTextAlign(33);
        text02->SetTextSize(0.025);

        TEllipse* el0 = new TEllipse(gCentreX_pr   ->at(0), 
                                     gCentreY_pr   ->at(0),
                                     maxDistGCen_pr->at(0),
                                     maxDistGCen_pr->at(0));
        el0->SetFillStyle(0);



        TCanvas* c_pr1 = new TCanvas("cpr1", "cpr1");
        c_pr1->SetWindowSize(csize, csize);
        TH2D* h_pr1   = new TH2D("Projection-1st (Big)", 
                                 "Projection-1st (Big)", 
                                        100, xprmin, xprmax,
                                        100, yprmin, yprmax);
        h_pr1->SetContour(99); //nicer, default is like 10
        h_pr1->GetXaxis()->SetTitle("x projection");
        h_pr1->GetYaxis()->SetTitle("y projection");
        h_pr1->GetZaxis()->SetTitle("# Energy");
        
        dist = "Max Distance: "+to_string(maxDistG_pr->at(1));
        dist_text = dist.c_str();
        TText* text10 = new TText(xprmin,  yprmax, dist_text);
        text10->SetTextAlign(11);
        text10->SetTextSize(0.025);
        
        distcen = "Max Distance from Centre: "+to_string(maxDistGCen_pr->at(1));
        distcen_text = distcen.c_str();
        TText* text11 = new TText(xprmax,  yprmax, distcen_text);
        text11->SetTextAlign(31);
        text11->SetTextSize(0.025);

        cen = "("+ to_string(gCentreX_pr->at(1))+", "
                 + to_string(gCentreY_pr->at(1)) + ")";
        cen_text = cen.c_str();
        TText* text12 = new TText(xprmax, yprmax, cen_text);
        text12->SetTextAlign(33);
        text12->SetTextSize(0.025);

        TEllipse* el1 = new TEllipse(gCentreX_pr   ->at(1), 
                                     gCentreY_pr   ->at(1),
                                     maxDistGCen_pr->at(1),
                                     maxDistGCen_pr->at(1));
        el1->SetFillStyle(0);



        TCanvas* c_2ph = new TCanvas("c_2ph", "c_2ph");
        c_2ph->SetWindowSize(csize, csize);
        c_2ph->DrawFrame(xprmin, yprmin, xprmax, yprmax);
        //auto pad = new TPad ();
        //pad -> SetTitle("Differentiating Photons (Big)");

        const int N0 = xpr0.size();
        double xpr0arr [N0]; double ypr0arr [N0];
        for (int i = 0; i<N0; i++)
            {xpr0arr[i] = xpr0[i]; ypr0arr[i] = ypr0[i];}
        auto gr0 = new TGraph(N0, xpr0arr, ypr0arr); 
        gr0->SetMarkerStyle(47);
        gr0->SetMarkerSize(0.5);
        gr0->SetMarkerColor(4); 
        gr0->SetTitle("Photons' Showers;x-projection;y-projection");
        //gr0->GetXaxis()->SetRange(xprmin, xprmax);
        //gr0->GetXaxis()->SetRange(yprmin, yprmax);

        const int N1 = xpr1.size();
        double xpr1arr [N1]; double ypr1arr [N1];
        for (int i = 0; i<N1; i++)
            {xpr1arr[i] = xpr1[i]; ypr1arr[i] = ypr1[i];}
        auto gr1 = new TGraph(N1, xpr1arr, ypr1arr);
        gr1->SetMarkerStyle(33);
        gr1->SetMarkerSize(0.6);
        gr1->SetMarkerColor(2);
        gr1->SetTitle("Photons' Showers;x-projection;y-projection");
        //gr0->SetTitle("Differentiating Photons");
        //gr0->GetXaxis()->SetTitle("x projection");
        //gr0->GetYaxis()->SetTitle("y projection");
        //gr1->GetXaxis()->SetTitle("x projection");
        //gr1->GetYaxis()->SetTitle("y projection");
        //gr1->GetXaxis()->SetRange(xprmin, xprmax);
        //gr1->GetXaxis()->SetRange(yprmin, yprmax);
        


        // FILL HISTS
        for (int j = 0; j<N; j++)
        {
            int g_ij = which_g->at(j);
            double Ei = (g_ij)? E1 : E0;
            double xpr_ij  = xpr->at(j);
            double ypr_ij  = ypr->at(j);
            double efr_ij  = energy->at(j)/Etot;
            double efrg_ij = energy->at(j)/Ei;

            // Fill Hists
            if (isinf(ievent) && nrechits->at(j)>=3) 
            {
                h_prtot->Fill(xpr_ij, ypr_ij, efr_ij);
                if (g_ij) 
                    {h_pr1 -> Fill(xpr_ij, ypr_ij, efrg_ij);}
                else  
                    {h_pr0 -> Fill(xpr_ij, ypr_ij, efrg_ij);}
            }
            else if (i==ievent && nrechits->at(j)>=3) 
            {
                h_prtot->Fill(xpr_ij, ypr_ij, efr_ij);
                if (g_ij) 
                    {h_pr1  ->Fill(xpr_ij, ypr_ij, efrg_ij);}
                else  
                    {h_pr0->Fill(xpr_ij, ypr_ij, efrg_ij);}
            } 
        }

        gStyle->SetPalette(kRainBow); //nice colour option for colz
        
        c_prtot->cd();
        h_prtot-> SetStats(0);
        h_prtot-> Draw("colz");
        textall0->Draw();
        textall1->Draw();
        textall2->Draw();
        elall  -> Draw();
        c_prtot-> Update();

        c_2ph->cd();
        c_2ph->DrawFrame(xprmin, yprmin, xprmax, yprmax);
        gr0-> SetStats(0);
        gr0-> Draw("P"); //P = point only
        gr1-> SetStats(0);
        gr1-> Draw("P");
        c_2ph-> Update();

        c_pr0  ->cd();
        h_pr0  -> SetStats(0);
        h_pr0  -> Draw("colz");
        text00 -> Draw();
        text01 -> Draw();
        text02 -> Draw();
        el0    -> Draw();
        c_pr0  -> Update();

        c_pr1   ->cd();
        h_pr1  -> SetStats(0);
        h_pr1  -> Draw("colz");//h_pr1  ->Draw("colz");
        text10 -> Draw();
        text11 -> Draw();
        text12 -> Draw();
        el1    -> Draw();
        c_pr1  -> Update();

        c_prtot -> WaitPrimitive();
        
        string x;
        while (1)
        {
            cout<<"\nType <go> to go to next event, <exit> to close:";
            cin >> x;
            if (x=="go") break;
            else if (x=="exit") return;
        }
        c_prtot->Close();
        c_pr0  ->Close();
        c_pr1  ->Close();
        c_2ph  ->Close();

        delete c_pr0, c_pr1, c_prtot;
        delete h_pr0, h_pr1, h_prtot;
        delete textall0, text00, text10;
        delete textall1, text01, text11;
        delete textall2, text02, text12;
        delete elall, el0, el1;
    

    }
    cout<<"Checkpoint: Finished Loop"<<endl;

    return;
    file->Close();   
}







/*double xprmin = xpr->at(argmin(*xpr));
        double xprmax = xpr->at(argmax(*xpr));
        double yprmin = ypr->at(argmin(*ypr));
        double yprmax = ypr->at(argmax(*ypr));
        double extens = (xprmax-xprmin>yprmax-yprmin)?
                         xprmax-xprmin:yprmax-yprmin;
        xprmin -= extens*0.1;
        yprmin -= extens*0.1;
        xprmax  = xprmin + extens*1.1; 
        yprmax  = yprmin + extens*1.1;*/