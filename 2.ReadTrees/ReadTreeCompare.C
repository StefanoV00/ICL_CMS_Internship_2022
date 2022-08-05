/// \author Stefano Veroni
/// File to read information from a tree file

/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 2.ReadTrees
root ReadTreeCompare.C -q
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


#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include <TFrame.h>
#include "TGraphErrors.h"
#include "TH1.h"
#include <TH2.h>
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
//-1. Functions' Definitions/////////////////////////////////////////////
//______________________________________________________________________
#include "../0.Helpers/MyRoot.h"
#include "../0.Helpers/MyCfunctions.h"



void ReadTreeCompare() 
{
    ////////////////////////////////////////////////////////////////////
    //0. Set Variables   ///////////////////////////////////////////////
    //__________________________________________________________________
    string pathd    = "../myTTrees/TreesEn90to110/Double/";
    string paths    = "../myTTrees/TreesEn90to110/Single/";//none if none
    string rootname= "step3ticl_eta21_FlatTracksters";//root file
    string frootname = "treeFriends";
    //string rootname= "treeSelections";//file
    //string dirname = "CloseXY"; //"none" if no directory in between
    string dirname = "none";
    string treename= "TSTree_SimFromCP";
    
    string friendtype  = "_Shape";
    const char* title  = "2nd Normalised Transverse Moment";
    const char* totname= "lc_2ndTranMomentNorm";
    const char* gname  = "lc_2ndTranMomentNorm";
    int rebinint = 1;

    auto getxmax = [](double xmaxt, double xmaxg)
                   {return max(xmaxt, xmaxg);};
    auto getxmin = [](double xmint, double xming)
                   {return min(xmint, xming);};
    auto getymax = [](double ymaxt, double ymaxg)
                   {return 2*max(ymaxt, ymaxg);};
    
    
    const char* cumtitle  = "%Layers with #lcs<x (from start)";
    bool drawcum = 0;
    bool normcum = 1;


    ///////////////////////////////////////////////////////////////////
    // 1. GET TTREES AND SET ADDRESSES   //////////////////////////////
    /////////////////////////////////////////////////////////////////// 
    TTree* myfriendd; int dNentries;
    TTree* myfriends; int sNentries;
    if (friendtype != "none")
    {
        TFile* filed = new TFile((pathd+frootname+".root").c_str());
        //open_file(pathd, frootname);
        cout<<"doubles file obtained"<<endl;
        cout<<treename<<endl;
        cout<<friendtype<<endl;
        myfriendd = (TTree*) filed
                                ->Get((treename + friendtype).c_str());
        dNentries = myfriendd->GetEntries();
        cout<<"doubles ttree obtained"<<endl;

        TFile* files = new TFile((paths+frootname+".root").c_str());
        cout<<"singles file obtained"<<endl;
        myfriends = (TTree*) files
                                ->Get((treename + friendtype).c_str());
        sNentries = myfriends->GetEntries();
        cout<<"singles ttree obtained"<<endl;
    }
    else
    {
        TFile* filed = new TFile((pathd+rootname+".root").c_str());
        //open_file(pathd, frootname);
        cout<<"doubles file obtained"<<endl;
        cout<<treename<<endl;
        cout<<friendtype<<endl;
        myfriendd = (TTree*) filed
                                ->Get((treename).c_str());
        dNentries = myfriendd->GetEntries();
        cout<<"doubles ttree obtained"<<endl;

        TFile* files = new TFile((paths+rootname+".root").c_str());
        cout<<"singles file obtained"<<endl;
        myfriends = (TTree*) files
                                ->Get((treename).c_str());
        sNentries = myfriends->GetEntries();
        cout<<"singles ttree obtained"<<endl;
    }
    

////////////////////////////////////////////////////////////////////////
//////2. Prepare Hist //////////////////////////////////////////////////
//////__________________________________________________________________
    
    // Prepare separate ones
    TCanvas* c0 = new TCanvas("c0", "c0"); 
    myfriendd->Draw(totname);
    auto htot = (TH1F*)gPad->GetPrimitive("htemp"); // 1D

    TCanvas* cg = new TCanvas("c1", "c1"); 
    myfriends->Draw(gname);
    auto hg = (TH1F*)gPad->GetPrimitive("htemp"); // 1D

    // Get separate ones' parameters - xaxis
    htot -> Scale(1./dNentries);
    hg   -> Scale(1./sNentries);
    htot ->Rebin(rebinint);
    hg   ->Rebin(rebinint);
    double xmint = htot->GetXaxis()->GetXmin();
    double xmaxt = htot->GetXaxis()->GetXmax();
    int    Ntot  = htot->GetNbinsX();
    double widtht= (xmaxt - xmint) / (double) Ntot;
    double xming = hg  ->GetXaxis()->GetXmin();
    double xmaxg = hg  ->GetXaxis()->GetXmax(); 
    int    Ng    = hg  ->GetNbinsX();
    double widthg= (xmaxg - xming) / (double) Ng;
    // Get ps for general plotting -xaxis
    auto xmin  = getxmin(xmint, xming);
    auto xmax  = getxmax(xmaxt, xmaxg);
    int  nx    = (int)((xmax-xmin)/min(widtht, widthg));
    // Get separate ones' parameters - yaxis 
    auto   ymaxt = htot->GetMaximum();
    auto   ymaxg = hg  ->GetMaximum();
    // Get ps for general plotting - yaxis
    double totmultiplier = (double) Ntot / (double)nx;
    double gmultiplier   = (double) Ng / (double)nx;
    ymaxt *= 1.1 * totmultiplier *rebinint;
    ymaxg *= 1.1 * gmultiplier   *rebinint;
    auto ymax  = getymax(ymaxt, ymaxg);
    
    // Prepare separate ones with new parameters
    //TCanvas* c0_bis = new TCanvas("c0_bis", "c0_bis"); 
    string histtotname (totname);
    histtotname += (string)">>htot"+(string)"("
                   +to_string(nx)  +","
                   +to_string(xmin)+","
                   +to_string(xmax)+")";
    myfriendd->Draw(histtotname.c_str());
    htot = (TH1F*)gPad->GetPrimitive("htot"); // 1D

    //TCanvas* c1_bis = new TCanvas("c1_bis", "c1_bis");
    string histgname (gname);
    histgname   += (string)">>hg"  +(string)"("
                   +to_string(nx)  +","
                   +to_string(xmin)+","
                   +to_string(xmax)+")";
    myfriends->Draw(histgname.c_str());
    hg = (TH1F*)gPad->GetPrimitive("hg"); // 1D

    // Do combined one
    TCanvas* creal = new TCanvas("creal", "creal");
    htot -> Scale(1./dNentries);
    htot -> Rebin(rebinint);
    htot -> SetLineColor(2);
    htot -> SetStats(0);
    hg   -> Rebin(rebinint);
    hg   -> SetStats(0);
    hg   -> Scale(1./sNentries);
    htot->SetMaximum(ymax*1.1); //htot->SetAxisRange(xmin, xmax);
    hg  ->SetMaximum(ymax*1.1); //hg  ->SetAxisRange(xmin, xmax);

    creal->DrawFrame(xmin, 0, xmax, ymax*1.1);
    creal->cd();
    hg   -> Draw("HIST");
    htot -> Draw("HIST SAME");
    auto legend = new TLegend(0.7, 0.8, 0.88, 0.88);
    legend->AddEntry(htot,"Doubles");
    legend->AddEntry(hg  ,"Singles");
    legend->Draw();
    hg   ->SetXTitle(title);
    htot ->SetXTitle(title);
    hg   ->SetTitle(title);
    htot ->SetTitle(title);
    creal->Update();



    if (drawcum)
    {
        TCanvas* ccum= new TCanvas(); 
        auto htotcum = htot->GetCumulative();
        auto hgcum   = hg  ->GetCumulative();

        htotcum ->Rebin(rebinint);
        htotcum -> SetLineColor(2);
        htotcum -> SetStats(0);

        hgcum   ->Rebin(rebinint);
        hgcum   -> SetStats(0);

        //double xmint = htot->FindFirstBinAbove(); 
        //double xming = hg  ->FindLastBinAbove();
        double xmint = htotcum->GetXaxis()->GetXmin();
        double xming = hgcum  ->GetXaxis()->GetXmin();
        auto xmaxt = htotcum->GetXaxis()->GetXmax();
        auto xmaxg = hgcum  ->GetXaxis()->GetXmax(); 
        auto xmin = (xmint<xming)? xmint : xming;
        auto xmax = (xmaxt>xmaxg)? xmaxt : xmaxg;
        htotcum->SetAxisRange(xmin, xmax);
        hgcum  ->SetAxisRange(xmin, xmax);
        if (normcum)
        {
            double ymaxt = htotcum->GetBinContent(htotcum->FindLastBinAbove());
            double ymaxg = hgcum  ->GetBinContent(hgcum  ->FindLastBinAbove());
            auto ymax = (ymaxt>ymaxg) ? ymaxt : ymaxg;
            htotcum -> Scale(1./ymaxt);
            hgcum   -> Scale(1./ymaxg);
        }
        else
        {
            double ymaxt = htotcum->GetBinContent(htotcum->FindLastBinAbove());
            double ymaxg = hgcum  ->GetBinContent(hgcum  ->FindLastBinAbove());
            auto ymax = (ymaxt>ymaxg) ? ymaxt : ymaxg;
            htotcum ->SetMaximum(ymax*1.1);
            hgcum   ->SetMaximum(ymax*1.1);     
        }

        //ccum->DrawFrame(xmin, 0, xmax, ymax*1.2);
        hgcum   -> Draw("HIST");
        hgcum   ->SetTitle(cumtitle);
        hgcum   ->SetXTitle("x");
        htotcum ->SetTitle(cumtitle);
        htotcum ->SetXTitle("x");
        htotcum -> Draw("HIST SAME");
        ccum->Update();

        ccum   ->WaitPrimitive();
    }

    c0   ->WaitPrimitive();
    cg   ->WaitPrimitive();
    creal->WaitPrimitive();
    return;
}


/*TCanvas* creal2 = new TCanvas();
    creal2 ->cd();
    htot -> SetLineColor(2);
    TH1* htotnorm = htot->DrawNormalized("HIST");
    //htot -> Draw("HIST");
    creal->Update();
    creal2->cd();
    TH1* hgnorm = hg-> DrawNormalized("HIST SAME");
    creal2->Update();*/

/*if (ymaxt<ymaxg)
    {
        htot -> Draw("HIST");
        hg   -> Draw("HIST SAME");
    }
    else
    {
        hg   -> Draw("HIST");
        htot -> Draw("HIST SAME");
    }
*/