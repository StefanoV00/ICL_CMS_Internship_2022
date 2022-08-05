/// \author Stefano Veroni

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cmath>
#include <vector>
#include "MyCfunctions.h"

// WHEN RUNNING A ROOT SESSION FROM CONSOLE, ALL WILL BE
// AVAILABLE EVEN IF NOT INCLUDED HERE.
// THE INCLUSION IS JUST TO HAVE AUTHOMATIC COMPLETION
// BY VISUAL STUDIO CODE
// BEFORE RUNNING, THOSE .hxx NOT IN ROOT/INCLUDE MUST BE
// COMMENTED OUT
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "Math/Vector2D.h"

#include "MyRoot.h"
using XYVector = ROOT::Math::XYVector;
using XYZVector = ROOT::Math::XYZVector;
using namespace std;


void MyTestRootStuff ()
{
    /*
    string path    = "../myTTrees/";
//none if none
    string rootname= "step3ticl_pt50_eta21_run0_FlatTracksters";//file
    string friendrootname= "treeFriends"; //root file with friends
    string dirname = "ticlTree"; //"none" if no directory in between
    string treename= "TSTree_SimFromCP";

    TFile* fr = open_file(path, friendrootname, "update");
    TDirectory* dir = (TDirectory*)fr->Get("CloseXYSameStart");
    bool a = dir == nullptr;
    bool b = dir == 0;
    cout <<dir<<endl;
    cout<<a<<endl;
    cout<<b<<endl;
    
    vector <double> xs = {-2.7, 0  , -2.1, 6.1, 3.3, 3.4, 1.5, 3.4, 1.1};
    vector <double> ys = { 3.2, 0.1, -1.2, 2.3,-2.1, 4.1,-6.2, 2  , 2.9};
    vector <double>* xpr = &xs;
    vector <double>* ypr = &ys;
    vector <int> gs =    {   1,   1,    0,   1,   0,   0,   1,   1,   0};
    vector <int> hs =    {   6,  10,    9,   5,   2,   3,   2,   4,   4};
    vector <int>* which_g = &gs;
    vector <int>* nrechits = &hs;
    auto cpvtx0_pr = XYVector(1, 0);
    auto cpvtx1_pr = XYVector(4.1, -1);
    auto centre0_pr = XYVector(1.2, 1.7);
    auto centre1_pr = XYVector(2  , 1.1);
    auto centre_pr  = 0.4 * centre1_pr + 0.6 * centre0_pr;

    vector <double> soln = get_max_distances(*xpr, *ypr, 
                                                 *which_g,   *nrechits, 
                                                 cpvtx0_pr,  cpvtx1_pr,
                                                 centre0_pr, centre1_pr,
                                                 centre_pr);
    //Cheked vs Geogebra 
    
    // Max distances between points
    cout <<soln[0]<<" vs 10.4738"<<endl;
    cout <<soln[1]<<" vs 8.91572"<<endl;
    cout <<soln[2]<<" vs 7.63806"<<endl;
    cout <<soln[3]<<" vs 10.2944"<<endl;
    cout <<soln[4]<<" vs 7.63806"<<endl;
    cout <<soln[5]<<" vs 8.84596"<<endl;
    cout <<soln[6]<<" vs 10.4738"<<endl;
    cout <<soln[7]<<" vs 8.91572"<<endl;

    // Centres
    cout <<soln[10]<<" vs 7.6600"<<endl;
    cout <<soln[11]<<" vs 4.6564"<<endl;
    cout <<soln[12]<<" vs 4.3932"<<endl;
    cout <<soln[13]<<" vs 7.3171"<<endl;
    cout <<soln[14]<<" vs 4.3932"<<endl;
    cout <<soln[15]<<" vs 5.1478"<<endl;

    // Vtxs
    cout <<soln[17]<<" vs 4.7508"<<endl;
    cout <<soln[18]<<" vs 7.9924"<<endl;
    cout <<soln[19]<<" vs 4.7508"<<endl;
    cout <<soln[20]<<" vs 7.9924"<<endl;
    */
    
    vector <double> rhos = {2,1,3,4,5,6,7,9,8,
                            11,10,12,13,14,15,16,17,20,18,19};
    vector <double> phis = {2,3,1,2,1,1,0.4,0.3,1,
                            2,2,3,2.3,1.2,2.3,1.5,1.8,0.8,0.9,0.3};
    vector <double> xs;
    vector <double> ys;
    for (int i = 0; i<rhos.size(); i++)
    {
        xs.push_back(rhos[i]*cos(phis[i]));
        ys.push_back(rhos[i]*sin(phis[i]));
    }
                                     //         //        //
    vector <double> ene = {12,  31,  29,  15,   4,   2,   2,   3,   2,
                                     ///  //   /// //-///      //
                           15,  17,  28,  10,  10,  10,   4,   2,   2,   1,   1};
    vector <int> gs =    {  1,   1,   1,   1,   1,   1,   1,   1,   1,
                            0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};
    auto centre  = XYVector( 0, 0);
    auto centre0 = XYVector( 1, 1);
    auto centre1 = XYVector( 0, 0);
    cout<<"Prepared get radii"<<endl;
    print_vector(get_radii(xs, ys, gs, ene, centre0, centre1, centre));

    return;
}