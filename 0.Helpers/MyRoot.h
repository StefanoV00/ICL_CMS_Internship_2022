/// \author Stefano Veroni

#pragma once

#include "MyCfunctions.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
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
#include "TVector3.h"
#include "Math/Vector3D.h"
#include "Math/Vector2D.h"
//#include "Riostream.h"
//#include <ROOT/RDataFrame.hxx>


using XYVector = ROOT::Math::XYVector;
using polarVector = ROOT::Math::Polar2DVector;
using XYZVector = ROOT::Math::XYZVector;
using namespace std;


inline
TFile* open_file (std::string path, 
                  std::string rootname, 
                  const char* option = "read")
{
    std::string path_file = (path != "none") ? path + rootname : rootname;
    path_file += ".root";
    TFile* file = new TFile(path_file.c_str(), option);
    if (file->IsZombie()) 
    {
        std::cout << "ZOMBIE: Error opening TFile" << std::endl;
    }
    return file;
}

inline
TDirectory* open_dir (TFile* file, 
                 std::string dirname) 
{
    const char* mydirname = dirname.c_str();
    TDirectory *mydir = (TDirectory*)file->Get(mydirname);  
    if (mydir == nullptr) 
        {mydir = file->mkdir(mydirname, mydirname);}
    return mydir;
}

inline
TTree* get_tree (TFile* file, 
                 std::string dirname, 
                 std::string treename) {
    if (dirname != "none") {
        TDirectory *mydir = (TDirectory*)file->Get(dirname.c_str());
        TTree *tree = (TTree*)mydir->Get(treename.c_str());  
        return tree;  
    }
    else {
        TTree *tree = (TTree*)file->Get(treename.c_str());
        return tree;
    }   
}

inline
TTree* get_myfriend (TFile* file, 
                     std::string dirname, 
                     std::string treename) {
    const char* friend_name = (treename + "_friend").c_str();
    try 
    {
    if (dirname != "none") {
        TDirectory* mydir = (TDirectory*)file->Get(dirname.c_str());
        TTree* treef = (TTree*)mydir->Get(friend_name);  
        return treef;  
    }
    else {
        TTree* treef = (TTree*)file->Get(friend_name);
        return treef;
    }  
    }
    catch (...) {
        std::cout<<"\nError opening Friend, returning nullptr."<<std::endl;
        TTree* treef = nullptr;
        return treef;
    } 
}

inline
int remove_TObject (TFile * file, 
                    std::string dirname, 
                    std::string objname)
/** @brief NOTE: objname = "name;cycle".
 * **/
{
    if (dirname != "none") {
        TDirectory *mydir = (TDirectory*)file->Get(dirname.c_str());
        mydir->Delete(objname.c_str());   
        return 1; 
    }
    else {
        file->Delete(objname.c_str());
        return 1;
    }   
}



////////////////////////////////////////////////////////////////////////
// Get Max Distances in set of lcs (for AddShape.C)
////////////////////////////////////////////////////////////////////////
vector <double> get_max_distances (vector <XYVector> posvec, 
                                    vector <int> which_g, 
                                    vector <int> nrechits,
                                    XYVector vtx0, XYVector vtx1,
                                    XYVector centre)
{
    XYVector deltavtx = vtx1 - vtx0;
    int Nentries = posvec.size();
    double maxdist              = .0;
    double maxdist_big          = .0; 
    double maxdist_same0        = .0;
    double maxdist_same1        = .0;
    double maxdist_same0_big    = .0;
    double maxdist_same1_big    = .0;

    double maxdist_diff         = .0;
    double maxdist_diff_big     = .0;
    double maxdist_diff_net     = .0;
    double maxdist_diff_net_big = .0;
    
    double maxdist_centre       = .0;
    double maxdist_centre_big   = .0;
    double maxdist_tag          =-1.;

    double maxdist_vtx0         = .0;
    double maxdist_vtx1         = .0;
    double maxdist_vtx0_big     = .0;
    double maxdist_vtx1_big     = .0;

    if (Nentries == 0) 
    {   vector <double> sol= {
        maxdist              ,
        maxdist_big          , 
        maxdist_same0        ,
        maxdist_same1        ,
        maxdist_same0_big    ,
        maxdist_same1_big    ,
        maxdist_diff         ,
        maxdist_diff_big     ,
        maxdist_diff_net     ,
        maxdist_diff_net_big ,
        maxdist_centre       ,
        maxdist_centre_big   ,
        maxdist_tag          , //NOTE: NO EVENTS MEANS TAG = -1
        maxdist_vtx0         ,
        maxdist_vtx1         ,
        maxdist_vtx0_big     ,
        maxdist_vtx1_big     }; 
        return sol;}

    for (int i=0; i<Nentries; i++)
    {
        XYVector pos  = posvec[i];
        int gflag     = which_g[i];
        int nhits     = nrechits[i];

        for (int j=i+1; j<Nentries; j++)
        {
            XYVector deltapos_ij = pos - posvec[j];
            double delta_ij = sqrt(deltapos_ij.mag2());
            if (delta_ij>maxdist) 
                {maxdist = delta_ij;}
            if (nhits>=3 && nrechits[j]>=3 && maxdist_big < delta_ij) 
                {maxdist_big = delta_ij;}

            if (gflag!=which_g[j])
            {
                if (delta_ij>maxdist_diff)
                    {maxdist_diff = delta_ij;
                    int ori = get_sign(deltapos_ij.Dot(deltavtx));
                    XYVector deltapos_ij_net = deltapos_ij - ori*deltavtx;
                    double  delta_ij_net = sqrt(deltapos_ij_net.mag2());
                    maxdist_diff_net = delta_ij_net;}
                if (nhits>=3 && nrechits[j]>=3 && maxdist_diff_big<delta_ij)
                    {maxdist_diff_big = delta_ij;
                    int ori = get_sign(deltapos_ij.Dot(deltavtx));
                    XYVector deltapos_ij_net = deltapos_ij - ori*deltavtx;
                    double  delta_ij_net = sqrt(deltapos_ij_net.mag2());
                    maxdist_diff_net_big = delta_ij_net;}
            }

            else if (gflag==0 && which_g[j]==0)
            { 
                if (delta_ij>maxdist_same0)
                    {maxdist_same0 = delta_ij;}
                if (nhits>=3 && nrechits[j]>=3 && maxdist_same0_big<delta_ij)
                    {maxdist_same0_big = delta_ij;}
            }
            
            else if (gflag==1 && which_g[j]==1)
            {
                if (delta_ij>maxdist_same1)
                    {maxdist_same1 = delta_ij;}
                if (nhits>=3 && nrechits[j]>=3 && maxdist_same1_big<delta_ij)
                    {maxdist_same1_big = delta_ij;}
            }
            
        }
        
        XYVector deltapos_ij_centre = pos - centre;
        double delta_ij_centre = sqrt(deltapos_ij_centre.mag2());
        if (delta_ij_centre > maxdist_centre)
            {maxdist_centre = delta_ij_centre;}
        if (nhits >= 3 && maxdist_centre_big < delta_ij_centre)
            {maxdist_centre_big = delta_ij_centre;}
        
        
        if (gflag == 0)
        {
            XYVector deltapos_ij_vtx = pos - vtx0;
            double delta_ij_vtx = sqrt(deltapos_ij_vtx.mag2());
            if (delta_ij_vtx > maxdist_vtx0)
                {maxdist_vtx0 = delta_ij_vtx;}
            if (nhits >= 3 && maxdist_vtx0_big < delta_ij_vtx)
                {maxdist_vtx0_big = delta_ij_vtx;}
        }
        else
        {
            XYVector deltapos_ij_vtx = pos - vtx1;
            double delta_ij_vtx = sqrt(deltapos_ij_vtx.mag2());
            if (delta_ij_vtx > maxdist_vtx1)
                {maxdist_vtx1 = delta_ij_vtx;}
            if (nhits >= 3 && maxdist_vtx1_big < delta_ij_vtx)
                {maxdist_vtx1_big = delta_ij_vtx;}
        }   
    }

    std::string s_tag = "";
    s_tag += (maxdist == maxdist_diff) ? "1" : "0";
    s_tag += (maxdist == maxdist_big)  ? "1" : "0";
    maxdist_tag = atof(s_tag.c_str());
    
    vector <double> sol= {
        maxdist              ,
        maxdist_big          , 
        maxdist_same0        ,
        maxdist_same1        ,
        maxdist_same0_big    ,
        maxdist_same1_big    ,
        
        maxdist_diff         ,
        maxdist_diff_big     ,
        maxdist_diff_net     ,
        maxdist_diff_net_big ,
        
        maxdist_centre       ,
        maxdist_centre_big   ,
        maxdist_tag          , //NOTE: NO EVENTS MEANS TAG = -1
        
        maxdist_vtx0         ,
        maxdist_vtx1         ,
        maxdist_vtx0_big     ,
        maxdist_vtx1_big     }; 
    return sol;
} 

template <typename T>
vector <double> get_max_distances ( vector <T> xvec, vector <T> yvec, 
                                    vector <int> which_g, 
                                    vector <int> nrechits,
                                    XYVector vtx0, XYVector vtx1,
                                    XYVector centre)
{
    vector <XYVector> posvec (xvec.size());
    for (int j = 0; j<xvec.size(); j++)
        {posvec[j] = XYVector(xvec[j], yvec[j]);}
    return get_max_distances(posvec, which_g, nrechits,
                              vtx0, vtx1, centre);
}

vector <double> get_max_distances (vector <XYVector> posvec, 
                                    vector <int> which_g, 
                                    vector <int> nrechits,
                                    XYVector vtx0, XYVector vtx1,
                                    XYVector centre0, XYVector centre1,
                                    XYVector centre)
{
    XYVector deltavtx = vtx1 - vtx0;
    int Nentries = posvec.size();
    double maxdist              = .0; //0
    double maxdist_big          = .0; 
    double maxdist_same0        = .0;
    double maxdist_same1        = .0;
    double maxdist_same0_big    = .0;
    double maxdist_same1_big    = .0; //5

    double maxdist_diff         = .0;
    double maxdist_diff_big     = .0;
    double maxdist_diff_net     = .0;
    double maxdist_diff_net_big = .0; //9
    
    double maxdist_centre       = .0;
    double maxdist_centre_big   = .0;
    double maxdist_centre0      = .0;
    double maxdist_centre1      = .0;
    double maxdist_centre0_big  = .0;
    double maxdist_centre1_big  = .0;
    double maxdist_tag          =-1.; //16

    double maxdist_vtx0         = .0;
    double maxdist_vtx1         = .0;
    double maxdist_vtx0_big     = .0;
    double maxdist_vtx1_big     = .0; //20

    if (Nentries == 0) 
    {   vector <double> sol= {
        maxdist              ,
        maxdist_big          , 
        maxdist_same0        ,
        maxdist_same1        ,
        maxdist_same0_big    ,
        maxdist_same1_big    ,

        maxdist_diff         ,
        maxdist_diff_big     ,
        maxdist_diff_net     ,
        maxdist_diff_net_big ,

        maxdist_centre       ,
        maxdist_centre_big   ,
        maxdist_centre0      ,
        maxdist_centre1      ,
        maxdist_centre0_big  ,
        maxdist_centre1_big  ,
        maxdist_tag          , //NOTE: NO EVENTS MEANS TAG = -1
        
        maxdist_vtx0         ,
        maxdist_vtx1         ,
        maxdist_vtx0_big     ,
        maxdist_vtx1_big     }; 
        return sol;}

    vector <double> best_centre0 (3);
    vector <double> best_centre1 (3);
    vector <double> best_centre  (3);
    for (int i=0; i<Nentries; i++)
    {
        XYVector pos  = posvec[i];
        int gflag     = which_g[i];
        int nhits     = nrechits[i];

        for (int j=i+1; j<Nentries; j++)
        {
            XYVector deltapos_ij = pos - posvec[j];
            double delta_ij = sqrt(deltapos_ij.mag2());
            if (delta_ij>maxdist) 
                {maxdist = delta_ij;}
            if (nhits>=3 && nrechits[j]>=3 && maxdist_big < delta_ij) 
                {maxdist_big = delta_ij;}

            if (gflag!=which_g[j])
            {
                if (delta_ij>maxdist_diff)
                    {maxdist_diff = delta_ij;
                    int ori = get_sign(deltapos_ij.Dot(deltavtx));
                    XYVector deltapos_ij_net = deltapos_ij - ori*deltavtx;
                    double  delta_ij_net = sqrt(deltapos_ij_net.mag2());
                    maxdist_diff_net = delta_ij_net;}
                if (nhits>=3 && nrechits[j]>=3 && maxdist_diff_big<delta_ij)
                    {maxdist_diff_big = delta_ij;
                    int ori = get_sign(deltapos_ij.Dot(deltavtx));
                    XYVector deltapos_ij_net = deltapos_ij - ori*deltavtx;
                    double  delta_ij_net = sqrt(deltapos_ij_net.mag2());
                    maxdist_diff_net_big = delta_ij_net;}
            }

            else if (gflag==0 && which_g[j]==0)
            { 
                if (delta_ij>maxdist_same0)
                    {maxdist_same0 = delta_ij;}
                if (nhits>=3 && nrechits[j]>=3 && maxdist_same0_big<delta_ij)
                    {maxdist_same0_big = delta_ij;}
            }
            
            else if (gflag==1 && which_g[j]==1)
            {
                if (delta_ij>maxdist_same1)
                    {maxdist_same1 = delta_ij;}
                if (nhits>=3 && nrechits[j]>=3 && maxdist_same1_big<delta_ij)
                    {maxdist_same1_big = delta_ij;}
            }
            
        }

        XYVector deltapos_i_centre = pos - centre;
        double delta_i_centre = sqrt(deltapos_i_centre.mag2());
        if (delta_i_centre > maxdist_centre)
            {maxdist_centre = delta_i_centre;
             best_centre[0] = i;
             best_centre[1] = pos.x();
             best_centre[2] = pos.y();}
        if (nhits >= 3 && maxdist_centre_big < delta_i_centre)
            {maxdist_centre_big = delta_i_centre;}
        
        
        if (gflag == 0)
        {
            XYVector deltapos_i_vtx0 = pos - vtx0;
            double delta_i_vtx0 = sqrt(deltapos_i_vtx0.mag2());
            if (delta_i_vtx0 > maxdist_vtx0)
                {maxdist_vtx0 = delta_i_vtx0;}
            if (nhits>=3 && maxdist_vtx0_big < delta_i_vtx0)
                {maxdist_vtx0_big = delta_i_vtx0;}
                
            XYVector deltapos_i_centre0 = pos - centre0;
            double delta_i_centre0 = sqrt(deltapos_i_centre0.mag2());
            if (delta_i_centre0 > maxdist_centre0)
                {maxdist_centre0 = delta_i_centre0;
                best_centre0[0] = i;
                best_centre0[1] = pos.x();
                best_centre0[2] = pos.y();}
            if (nhits>=3 && maxdist_centre0_big < delta_i_centre0)
                {maxdist_centre0_big = delta_i_centre0;}
        }
        else
        {
            XYVector deltapos_i_vtx1 = pos - vtx1;
            double delta_i_vtx1 = sqrt(deltapos_i_vtx1.mag2());
            if (delta_i_vtx1 > maxdist_vtx1)
                {maxdist_vtx1 = delta_i_vtx1;}
            if (nhits>=3 && maxdist_vtx1_big < delta_i_vtx1)
                {maxdist_vtx1_big = delta_i_vtx1;}
                
            
            XYVector deltapos_i_centre1 = pos - centre1;
            double delta_i_centre1 = sqrt(deltapos_i_centre1.mag2());
            if (delta_i_centre1 > maxdist_centre1)
                {maxdist_centre1 = delta_i_centre1;
                best_centre1[0] = i;
                best_centre1[1] = pos.x();
                best_centre1[2] = pos.y();}
            if (nhits>=3 && maxdist_centre1_big < delta_i_centre1)
                {maxdist_centre1_big = delta_i_centre1;}
        }   
    }
    //cout<<centre0.x()<<","<<centre0.y()<<endl;
    //print_vector(best_centre0);
    //cout<<centre1.x()<<","<<centre1.y()<<endl;
    //print_vector(best_centre1);
    //cout<<centre .x()<<","<<centre .y()<<endl;
    //print_vector(best_centre );

    std::string s_tag = "";
    s_tag += (maxdist == maxdist_diff) ? "1" : "0";
    s_tag += (maxdist == maxdist_big)  ? "1" : "0";
    maxdist_tag = atof(s_tag.c_str());
    
    vector <double> sol= {
        maxdist              ,
        maxdist_big          , 
        maxdist_same0        ,
        maxdist_same1        ,
        maxdist_same0_big    ,
        maxdist_same1_big    ,
        
        maxdist_diff         ,
        maxdist_diff_big     ,
        maxdist_diff_net     ,
        maxdist_diff_net_big ,
        
        maxdist_centre       ,
        maxdist_centre_big   ,
        maxdist_centre0      ,
        maxdist_centre1      ,
        maxdist_centre0_big  ,
        maxdist_centre1_big  ,
        maxdist_tag          ,
        
        maxdist_vtx0         ,
        maxdist_vtx1         ,
        maxdist_vtx0_big     ,
        maxdist_vtx1_big     }; 
    return sol;
} 

template <typename T>
vector <double> get_max_distances ( vector <T> xvec, vector <T> yvec, 
                                    vector <int> which_g, 
                                    vector <int> nrechits,
                                    XYVector vtx0, XYVector vtx1,
                                    XYVector centre0, XYVector centre1,
                                    XYVector centre)
{
    vector <XYVector> posvec (xvec.size());
    for (int j = 0; j<xvec.size(); j++)
        {posvec[j] = XYVector(xvec[j], yvec[j]);}
    return get_max_distances(posvec, which_g, nrechits, vtx0, vtx1, 
                              centre0, centre1, centre);
}



vector <double> get_radii (vector <XYVector> posvec,
                           vector <int> which_g,
                           vector <double> true_energy,
                           XYVector centre)
{
    vector <double> radii (3); //60, 90, 95 for tot

    int N = posvec.size();
    double E = myvecsum(true_energy);
    vector <double> esum (N);
    vector <double> distances (N);
    for (int i = 0; i<N; i++)
        {distances[i] = sqrt( (posvec[i]-centre).mag2() );}

    vector <double> enes = sortAwithB(true_energy, distances);
    sort(distances.begin(), distances.end());

    // enes is made fraction, and esum is cumulative fraction
    enes[0] /= E;
    partial_sum(enes.begin(), enes.end(), esum.begin(), 
               [E](double a, double b){b/=E; return a + b;});
    for (int i = 0; i<N; i++)
    {
        if (esum[i] >= 0.68) 
        {
            radii[0] = distances[i];
            for (int j = i; j<N; j++)
            {
                if (esum[j] >= 0.90)
                {
                    radii[1] = distances[j];
                    for (int k = j; k<N; k++)
                    {
                        if (esum[k] >= 0.95)
                        {
                            radii[2] = distances[k];
                            break; //k
                        }
                    }
                    break; //j
                }
            }
            break; //i
        }
    }
   
    return radii;
}

vector <double> get_radii (vector <XYVector> posvec,
                           vector <int> which_g,
                           vector <double> true_energy,
                           XYVector centre0, XYVector centre1,
                           XYVector centre)
{
    vector <double> radii (9); //60, 90, 95 for tot, photon0, photon1

    int N = posvec.size();
    double E = myvecsum(true_energy);
    vector <double> esum (N);
    vector <double> distances (N);
    for (int i = 0; i<N; i++)
        {distances[i] = sqrt( (posvec[i]-centre).mag2() );}

    vector <double> enes = sortAwithB(true_energy, distances);
    sort(distances.begin(), distances.end());

    // enes is made fraction, and esum is cumulative fraction
    enes[0] /= E;
    partial_sum(enes.begin(), enes.end(), esum.begin(), 
               [E](double a, double b){b/=E; return a + b;});
    for (int i = 0; i<N; i++)
    {
        if (esum[i] >= 0.68) 
        {
            radii[0] = distances[i];
            for (int j = i; j<N; j++)
            {
                if (esum[j] >= 0.90)
                {
                    radii[1] = distances[j];
                    for (int k = j; k<N; k++)
                    {
                        if (esum[k] >= 0.95)
                        {
                            radii[2] = distances[k];
                            break; //k
                        }
                    }
                    break; //j
                }
            }
            break; //i
        }
    }


    for (int gi = 0; gi<2; gi++)
    {
        auto posveci  = getAwhereB(posvec, which_g, gi);
        auto enei    = getAwhereB(true_energy, which_g, gi);
        int Ni       = enei.size();
        double Ei    = myvecsum(enei);

        vector <double> esumi (Ni);
        vector <double> distancesi (Ni);
        XYVector centrei = (gi == 0)? centre0 : centre1;
        for (int i = 0; i<Ni; i++)
            {distancesi[i] = sqrt( (posveci[i]-centrei).mag2() );}
        
        vector <double> enesi = sortAwithB(enei, distancesi);
        sort(distancesi.begin(), distancesi.end());

        // enes is made fraction, and esum is cumulative fraction
        enesi[0] /= Ei;
        partial_sum(enesi.begin(), enesi.end(), esumi.begin(), 
                    [Ei](double a, double b){b/=Ei; return a + b;});
        
        for (int i = 0; i<Ni; i++)
        {
            if (esumi[i] >= 0.68) 
            {
                radii[(gi+1)*3 + 0] = distancesi[i];
                for (int j = i; j<Ni; j++)
                {
                    if (esumi[j] >= 0.90)
                    {
                        radii[(gi+1)*3 + 1] = distancesi[j];
                        for (int k = j; k<Ni; k++)
                        {
                            if (esumi[k] >= 0.95)
                            {
                                radii[(gi+1)*3 +2] = distancesi[k];
                                break; //k
                            }
                        }
                        break; //j
                    }
                }
                break; //i
            }
        }
    }
   
    return radii;
}

template <typename T>
vector <double> get_radii (vector <T> xvec, vector<T> yvec,
                           vector <int> which_g,
                           vector <double> true_energy,
                           XYVector centre)
{
    vector <XYVector> posvec (xvec.size());
    for (int j = 0; j<xvec.size(); j++)
        {posvec[j] = XYVector(xvec[j], yvec[j]);}
    return get_radii(posvec, which_g, true_energy, centre);
}

template <typename T>
vector <double> get_radii (vector <T> xvec, vector<T> yvec,
                           vector <int> which_g,
                           vector <double> true_energy,
                           XYVector centre0, XYVector centre1,
                           XYVector centre)
{
    vector <XYVector> posvec (xvec.size());
    for (int j = 0; j<xvec.size(); j++)
        {posvec[j] = XYVector(xvec[j], yvec[j]);}
    return get_radii(posvec, which_g, true_energy, 
                     centre0, centre1, centre);
}
/*
    int c68 = 0;
    int c90 = 0;
    for (int i = 0; i<N; i++)
    {
        if (enes[i] > 0.95) 
        {
            radii(2) = distances[i];
            break;
        }
        else if (enes[i] > 0.90 && c90 == 0) 
        {
            radii(1) = distances[i];
            c90 = 1;
            continue;
        }
        else if (enes[i] > 0.68 && c68 == 0) 
        {
            radii[0] = distances[i];
            c68 = 1;
            continue;
        }
    }
    */