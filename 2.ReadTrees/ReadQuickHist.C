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




void ReadQuickHist() 
{
    string path    = "../myTTrees/TressEn95to105/Double/";
//none if none
    string rootname= "step3ticl_eta21_FlatTracksters";//root file
    string friendrootname= "treeFriends"; //root file with friends
    string dirname = "ticlTree"; //"none" if no directory in between
    string treename= "TSTree_SimFromCP";

    TFile* file = open_file(path, rootname);
    TTree* tree = get_tree (file, dirname, treename);
    TFile* filef = open_file(path, friendrootname);
    const char* coords_name = (treename + "_Coords").c_str();
    TTree* treeC = get_tree(filef, "none", coords_name);
    const char* energy_name = (treename + "_Energy").c_str();
    TTree* treeE = get_tree(filef, "none", energy_name);
    const char* shape_name = (treename + "_Shape").c_str();
    TTree* treeS = get_tree(filef, "none", shape_name);

    treeE->Scan("lc_tsFraction");
    return;

    TCanvas* c0 = new TCanvas();
    treeS->Draw("lc_areaTotFromStart");

    TCanvas* c1 = new TCanvas();
    treeS->Draw("lc_areagFromStart");
    c0 -> WaitPrimitive();
    return;
}