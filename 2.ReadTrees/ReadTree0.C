/// \author Stefano Veroni
/// File to read information from a tree file
/// Also saves some in better known csv format.

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
#include "TMath.h"
#include <TNtuple.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"
//#include <RDataFrame.hxx>

using namespace std;

/////////////////////////////////////////////////////////////////////////
//-1. Functions' Definitions/////////////////////////////////////////////
//______________________________________________________________________
#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"




void ReadTree0(bool save_csv = false) 
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
    int Nevents = 10;  //# events to consider, 0 if all
    float ievent = 177;//INFINITY;  //event of which display info, INFINITE if all
    
    ////////////////////////////////////////////////////////////////////
    //1. Get TTree & Branches///////////////////////////////////////////
    //__________________________________________________________________
    //NOTE: FOR EACH ENTRY, LC_QUANTITY IS A VECTOR OF "DETECTIONS"
    //FOR EXAMPLE: lc_x AND lc_layer WILL HAVE N ENTRIES, ONE PER 
    //CLUSTER
    TFile* file = open_file(path, rootname);
    TTree* tree = get_tree (file, dirname, treename);
    tree->Scan("lc_energy:lc_tsMult:lc_nrechits");

    // GET FRIEND: Normally, would just do
    //TTree* treef = tree->GetFriend((treename+"_friend").c_str());
    //But there was an error, could not save as friend TTree, so
    TFile* filef = open_file(path, friendrootname);
    const char* energy_name = (treename + "_Energy").c_str();
    TTree* treef = get_tree(filef, "none", energy_name);
    //treef->Print();
    //tree ->Scan("event:lc_TSidx:lc_x:lc_layer:sv_posX:cp_vtxX:cp_phi");
    //treef->Scan("lc_startDeltaX:lc_x_pr");
    //cout<<"\nNOTE return is here after the Scan, at beginning!!!!"<<endl;
    
    // Get Entries Number
    int Nentries = tree->GetEntries();
    Nevents = (Nentries>Nevents && Nevents>0) ? Nevents : Nentries;
    cout << "The TTree has total "<< Nentries <<" entries. "<< endl;
    cout << "We are working with "<< Nevents  <<" events."  << endl;
    
    //POINTERS MUST BE INITIALISED
    vector <double>* x     = 0;
    vector <double>* y     = 0;
    vector <double>* z     = 0;
    vector <double>* eta   = 0;
    vector <double>* phi   = 0;
    vector <double>* energy= 0;
    vector <int>* which_g = 0;
    vector <int>* layerN  = 0;
    vector <int>* Nrechits= 0;
    vector <int>* tsMult  = 0;
    tree->SetBranchAddress("lc_x",        &x);
    tree->SetBranchAddress("lc_y",        &y);
    tree->SetBranchAddress("lc_z",        &z);
    tree->SetBranchAddress("lc_eta",      &eta);
    tree->SetBranchAddress("lc_phi",      &phi);
    tree->SetBranchAddress("lc_energy",   &energy);
    tree->SetBranchAddress("lc_layer",    &layerN);
    tree->SetBranchAddress("lc_nrechits", &Nrechits );
    tree->SetBranchAddress("lc_TSidx",    &which_g );
    tree->SetBranchAddress("lc_tsMult",   &tsMult );
    //From Friend of TTree
    vector <double>* x_pr = 0;
    vector <double>* y_pr = 0;
    vector <double>* theta = 0;
    treef->SetBranchAddress("lc_x_pr",  &x_pr);
    treef->SetBranchAddress("lc_y_pr",  &y_pr);
    treef->SetBranchAddress("lc_theta", &theta);

    
    ////////////////////////////////////////////////////////////////////
    //2. Prepare Hist & Saving File/////////////////////////////////////
    //__________________________________________________________________
    //ROOT::RDataFrame dataframe ("TSTree_SimFromCP", "ticlTree");
    //Prepare Vector of Histograms
    TCanvas* c_projections = new TCanvas("c0", "c0");
    TH2D* h_projections = new TH2D("xy-projections", "xy-projections", 
                                    200, -0.6, 0.6,
                                    200, -0.6, 0.6);
    h_projections->SetContour(99); //nicer, default is like 10
    h_projections->GetXaxis()->SetTitle("x projection");
    h_projections->GetYaxis()->SetTitle("y projection");
    h_projections->GetZaxis()->SetTitle("# Events");

    vector <TH2D*> hist_vec; 
    vector <TCanvas*> canvas_vec; 
    for (int i=1; i<=50; i++) 
    {
        //TCanvas* c = new TCanvas(Form("c%d", i), Form("c%d", i));
        TH2D* h = new TH2D(Form("xyHist%d", i), Form("xy-Layer%d", i),
                            200, -100, 100,
                            200, -100, 100); 
        h->SetContour(99); //nicer, default is like 10
        h->GetXaxis()->SetTitle("x [cm]");
        h->GetYaxis()->SetTitle("y [cm]");
        h->GetZaxis()->SetTitle("Events");
        hist_vec.push_back(h); 
        //canvas_vec.push_back(c); //GIVES ERRORS
    }
    //Prepare CSV File where to save data
    ofstream save_file;
    if (save_csv)
    {
    string name = treename;
    string ext  = ".csv";
    string fullname = name+ext; 
    remove (fullname.c_str());
    save_file.open (fullname);
    save_file << "EventN; WhichPhoton;" //Entry number & Photon Flag
            << "x; y; z; x_pr; y_pr;"
            << "eta; theta; phi; energy;"
            << "layerN; Nrechits; tsMult"
            <<endl;
    }
    
    
    ////////////////////////////////////////////////////////////////////
    //3. Loop: get event, find quantities, save data, fill hists ///////
    //__________________________________________________________________
    //NOTE: RDataFrame.Define() could define new columns theta, etc...
    for (int i = 0; i<Nevents; i++)
    { 
        // Getting new entry assignes value of branch to branch address
        // overwriting at each iteration 
        tree  ->GetEntry(i);
        treef -> GetEntry(i);
        int N = x->size();

        for (int j = 0; j<N; j++)
        {
            //Find other quantities
            int layer_ij = layerN->at(j);
            double x_ij  = x->at(j);
            double y_ij  = y->at(j);
            double z_ij  = z->at(j);
            double x_pr_ij  = x_pr -> at(j);
            double y_pr_ij  = y_pr -> at(j);

            // Fill Hists
            if (isinf(ievent)) 
            {
                h_projections->Fill(x_pr_ij, y_pr_ij);
                hist_vec[layer_ij-1]->Fill(x_ij, y_ij);
            }
            else if (i==ievent) 
            {
                h_projections->Fill(x_pr_ij, y_pr_ij);
                hist_vec[layer_ij-1]->Fill(x_ij, y_ij);
            }

            //Save in csv format
            //Also, might (or not) want to try save_file<<tree->Scan()
            if (save_csv)
            {
                if (save_file.is_open())
                {  
                save_file<<            i << ";" //evnt N
                        << which_g->at(j)<< ";" //photon identifier
                        <<          x_ij << ";" 
                        <<          y_ij << ";"
                        <<          z_ij << ";"
                        <<       x_pr_ij << ";"
                        <<       y_pr_ij << ";"
                        <<   theta->at(j)<< ";"
                        <<     eta->at(j)<< ";"
                        <<     phi->at(j)<< ";"
                        <<  energy->at(j)<< ";"
                        <<      layer_ij << ";"
                        <<Nrechits->at(j)<< ";"
                        <<  tsMult->at(j)<< endl;
                }
                else 
                {
                    cout<< "At event: " << i<<", step "<<j << 
                        ", unable to open csv file" << endl;
                }
            }   
        }
    }
    save_file.close();
    cout<<"Checkpoint: Finished Loop"<<endl;
    
    ////////////////////////////////////////////////////////////////////
    //4. The end, ladies & gentlemen////////////////////////////////////
    //__________________________________________________________________
    gStyle->SetPalette(kRainBow); //nice colour option for colz
    h_projections->Draw("colz");
    //hist_vec[14].Draw();
    if (isinf(ievent)) //plot all events together
    {
        c_projections -> cd();
        h_projections ->Draw("colz");
        for (int i = 1; i < 20; i+=3) 
        {
            //canvas_vec[i].cd();
            TCanvas* c = new TCanvas(Form("c%d", i), Form("c%d", i));
            hist_vec[i-1]->SetStats(0);
            hist_vec[i-1]->Draw("colz");
        }
    }
    else //plot just ith event, specified by ievent
    {
        c_projections -> cd();
        h_projections ->Draw("colz");//Draw("SURF3");
        for (int i = 1; i < 20; i+=3) 
        {
            //canvas_vec[i].cd();
            TCanvas* c = new TCanvas(Form("c%d", i), Form("c%d", i));
            hist_vec[i-1]->SetStats(0);
            hist_vec[i-1]->Draw("colz");
        }
    }
    return;
    file->Close();   
    filef->Close();
}