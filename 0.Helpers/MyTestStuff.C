/// \author Stefano Veroni

/*
cd C:\root_v6.26.02\root-on-vscode 
C:\root_v6.26.02\bin\thisroot.bat 
cd 0.Helpers
root MyTestStuff.C -q
cd ../

*/

#include "MyCfunctions.h"
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
#include <TInterpreter.h>
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include <TNtuple.h>
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
#include "MyRoot.h"


void MyTestStuff ()
{

    /*
    TEST NAN
    double nan = NAN;
    cout<<nan;

    TEST getIndexes
    vector <double> myvec (100);
    myvec[10] = 10.2;
    myvec[50] = 10.2;
    myvec[99] = 10.2;
    vector <int> sol = getIndexes(myvec, 10.2);
    cout<<sol.size()<<endl;
    print_vector(sol); 

    TEST getAManywhereB
    int a = 5;
    int b = 100;
    vector <double> prototype (b);
    vector< vector<double> > myVec(a, prototype);
    vector <float> myU (b);
    float c = 10.2;
    myU[10] = c;
    myU[50] = c;
    myU[99] = c;
    for (int i = 0; i < myVec.size(); i++)
    {
        for (int j = 0; j < myVec[i].size(); j++)
        {
            myVec[i][j] = i * 100 + j;
        }
    }
    vector <vector<double>> sol = getAwhereB(myVec, myU, c);
    cout<<sol.size()<<endl;
    print_vector(sol[1]);

    TEST minAbsDiff
    vector <int> v1 = {1,2,3};
    vector <float> v2 = {19,17,21,8};
    cout<<minAbsDiff(v1, v2)<<endl;

    TEST argmax
    vector <double> myvec = {1, 6, 4, 8, 5, 9};
    cout<<argmax(myvec, 0, 0)<<endl;

    TEST sum_portions 
    vector <double> myvec = {1, 6, 4, 8, 5, 9};
    auto sol = sum_portions(myvec, 3);
    print_vector(sol); 

    int i = 4;
    if (i==NAN) {cout<<"yes"<<endl;}
     
    
    TEST sortAwithB 
    vector <string> a = {"b", "c", "a", "e", "d"};
    vector <double> b = {2.4, 3.5, 1, 67, 65};
    vector <string> anew = sortAwithB(a, b);
    print_vector(anew); 

    //TEST CONVEXHULL
    vector <vector<double>> myvec = {{16,3}, {12,17}, {0,6}, {-4,-6}, 
    {16,6}, {16,-7}, {16,-3}, {17,-4}, {5,19}, {19,-8}, 
    {3,16}, {12,13}, {3,-4}, {17,5}, {-3,15}, {-3,-9}, 
    {0,11}, {-9,-3}, {-4,-2}, {12,10}};

    const int N = myvec.size();
    vector<Point> points (N);
        for (int i = 0; i<N; i++)
        {
            points[i] = {myvec[i][0], myvec[i][1]};
        }
    //cout<<points[0].x<<endl;

    auto hull = ConvexHull(myvec);
    print_vector(hull);
    cout<<Shoelace(hull)<<endl;

    double* a = nullptr;
    cout<<a<<endl; 

    //TEST MEAN
    vector <int> myvec = {1,2,3,4};
    vector <double> myw = {0.5,0.5,0.55,0.20};

    cout<<myvecsum(myw)<<endl;

    cout<<mymean(myvec, myw)<<endl;
    
    double s = 0;
    double ws = 0;
    for (int i = 0; i<myvec.size(); i++)
    {
        s += myvec[i]*myw[i];
        ws += myw[i];
    }
    cout<<s/ws<<endl;

    //TEST ARGMIN E MAX
    vector <double> myvec = {.2,.1,.2,0,.3,.4,.5, .8, .6,.5};
    cout<<argmin(myvec)<<endl;
    cout<<argmax(myvec)<<endl;

    //TEST getAwhereB & LAMBDA function
    auto mylambda = [](string s) {return s=="yes";};
    vector<double> v1 = {1.2,    2.4,  3.1,   6.2,   7.8,    7.9, 7.8};
    vector<string> v2 = {"yes", "yes", "no", "yes", "yes", "no", "yes"};
    print_vector(getAwhereB(v1, v2, mylambda));

    //TEST CONTAINS
    vector <double> myvec = {1, 2, 5, 3.5, 6};
    cout<<contains(myvec, 4)<<endl; 

    // TEST PARTIAL SUM
    vector <double> enes = {1,2,3,4,5,6.1};
    double E = myvecsum(enes);
    vector <double> esum (enes.size());
    enes[0] /= E;
    partial_sum(enes.begin(), enes.end(), esum.begin(), [E](double a, double &b){b/=E; return a + b;});
    print_vector(esum);
    print_vector(enes); 

    // TEST FUNCTION MEAN
    vector <double> myvec = {1.2, 4.5, 6.5, 2.3};
    vector <double> weight = {1.3, 6.1, 0.1, 2.3};
    cout<<mymean(myvec, [](double x){return x*x;}, weight)
        <<endl; 

    // TEST RESIZE and inverse sortAwithB
    vector <double> myvec = {1.2, 4.5, 6.5, 2.3};
    vector <double> order = {1.2, 4.5, 3.1, 1};
    auto newvec = sortAwithB(myvec, order, false);
    print_vector(newvec);
    newvec.resize(2);
    print_vector(newvec); 

    // TEST MULT CUT
    vector <int> mult = {0, 0, 0, 1, 2};
    int first = 4;
    mult.erase(mult.begin(), mult.begin() + first - 1);
    print_vector(mult);

    int a = 2;
    cout<<a/3<<endl;

    vector <double> myvec = {1.3, 4.5, 6.7, 3.2};
    cout<<*(max_element(myvec.begin(), myvec.end()))<<endl; */

    if (2<3)
    {
        cout<<"a\n";
    }
    else if (2<1.5)
    {
        cout<<"b\n";
    }
    else if (2<5)
    {
        cout<<"c\n";
    }

    cout<<"Test Ended"<<endl;
    return;  
}