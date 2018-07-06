#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TGraph.h"
using namespace Eigen;
using namespace std;

void compute_piecewise_function(double A, vector<double>& b){
  b = {0,0,0,0,0};
  if(A<1000){
    b[0] = 1. - A/1000.;
    b[1] = A/1000.;
  }
  else if(A<2000){
    b[1] = 1. - (A-1000.)/1000.;
    b[2] = (A-1000)/1000.;
  }
  else if(A<3000){
    b[2] = 1 - (A-2000.)/1000.;
    b[3] = (A-2000)/1000.;
  }
  else if(A<4000){
    b[3] = 1 - (A-3000.)/1000.;
    b[4] = (A-3000.)/1000.;
  }
}

class WfPoint {
public:
	double t;
	double A;
	double G;
};

class Calpha {
public:
	vector<double> _t;
	vector<int> _N;
	vector<double> _c; // (A+f)/G

  void Compute(vector<WfPoint>, vector<double>); // A, G and piecewise function {y0,y1,...,yn}
	double GetValue(double); // t <=> c
};

/* c = avg of {(A+f)/G} in the same time-tick*/
void Calpha::Compute(vector<WfPoint> recv, vector<double> yv){
  _t.clear();
  _N.clear();
  _c.clear();
  vector<double>::iterator it;
  for(int i=0; i<recv.size(); i++){
    WfPoint rec = recv[i];
    double time = rec.t;
    /* find a float time in vector */
    it = std::find_if(_t.begin(), _t.end(), [time](double b) {return abs(time-b)<1E-5;} );
    if(it == _t.end()){
      _t.push_back(time);
      _N.push_back(1);
      vector<double> b(5);
      compute_piecewise_function(rec.A, b); // FIXME: reuse the function for same A
      double D = b[0] * yv[0] + b[1] * yv[1] + b[2] * yv[2] + b[3] * yv[3] + b[4] * yv[4];
      _c.push_back( (rec.A + D)/rec.G ); 
    }
    else{
      int index = it - _t.begin();
      _N[index] += 1;
      vector<double> b(5);
      compute_piecewise_function(rec.A, b); // FIXME: reuse the function for same A
      double D = b[0] * yv[0] + b[1] * yv[1] + b[2] * yv[2] + b[3] * yv[3] + b[4] * yv[4];
      _c[index] += (rec.A + D)/rec.G; 
    }
  }
  for(int i=0; i<_t.size(); i++){
    _c[i] /= _N[i];
  }
}

double Calpha::GetValue(double time){
  vector<double>::iterator it;
   /* find a float time in vector */
   it = std::find_if(_t.begin(), _t.end(), [time](double b) {return abs(time-b)<1E-5;} );
   if(it == _t.end()){
    cout << "Error: out of range in time series!" << endl;
    return 0;
   }
   int index = it - _t.begin();
   // cout << "c= " << _c[index] << " ; t= " << time << endl; 
   return _c[index];
}


int main()
{
  vector<WfPoint> recv;
  WfPoint therec;
  /* save all data points */
  double thistime, thisA, thisG;
  TFile* ifile = TFile::Open("protodune-nl-MC.root");
  TTree* T = (TTree*)ifile->Get("T");
  for(int i=0; i<T->GetEntries(); i++){
  	T->GetEntry(i);
  	double Ameas = T->GetLeaf("Ameas")->GetValue();
  	double time = T->GetLeaf("time")->GetValue();
  	double Vdac = T->GetLeaf("Vdac")->GetValue();
  	double Gamp = T->GetLeaf("Gamp")->GetValue();
  	cout << Ameas << " " << time << " " << Vdac << " " << Gamp << endl;
  	if(Ameas>500 && Ameas<4000 && (Vdac*183.<100.) && (Vdac*183.*Gamp*1E-3<1.4)){
  		therec.t = time;
  		therec.A = Ameas;
  		therec.G = Vdac*Gamp;
  		recv.push_back(therec);
  	}
  }
  ifile->Close();
  cout << "size of recv: " << recv.size() << endl;

  vector<double> yv = {0,0,0,0,0};
  // vector<double> yv = {-10,-10,-10,-10,-10};
  // vector<double> yv = {-22, -16, -9.5, -3.2, 3.1};
  // cout << "   true pars: -22, -16, -9.5, -3.2, 3.1" << endl;
  cout << "   initial guess: " << yv[0] << " " << yv[1] << " " << yv[2] << " " << yv[3] << " " << yv[4] << endl;

  
  /* ... iterate for n times ... */
  int nloop = 10;
  Calpha theC;  
  for(int iloop=0; iloop<nloop; iloop++){
    /* re-compute the coef c */
    theC.Compute(recv, yv);
    cout << "... nloop = " << iloop << endl;
    int ndim = recv.size();
    VectorXd M(ndim);
    SparseMatrix<double> R(ndim,yv.size());
    for(int i=0; i<ndim; i++){
      therec = recv[i];
      double c = theC.GetValue(therec.t);
      M(i) = therec.A - c * therec.G;

      vector<double> b(5);
      compute_piecewise_function(therec.A, b); 
      for(int j=0; j<yv.size(); j++){
        R.insert(i,j) = -b[j];
      }  
    }

    SparseMatrix<double> RTR(5,5);
    VectorXd x(5), RTM(5);
    RTR = R.transpose() * R;
    RTM = R.transpose() * M;
    /* ... solve the linear eqn: RTR*x = RTM ... */
    BiCGSTAB<SparseMatrix<double> > solver;
    solver.compute(RTR);
    x = solver.solve(RTM);

    // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error()      << std::endl;
    cout << "   pars = " << x.transpose() << endl;
    VectorXd chi(ndim);
    chi = M - R * x;
    cout << "   Chi2 = " << chi.transpose() * chi << endl << endl;

    /* save values for next iteration */  
    yv = {x(0), x(1), x(2), x(3), x(4)};

  }

  theC.Compute(recv,yv);
  double x[5] = {0, 1000, 2000, 3000, 4000};
  TGraph* gNL = new TGraph(5, x, &(yv[0]));
  gNL->SetName("NLcorrection");

  TGraph* gResp = new TGraph();
  gResp->SetName("Response");
  // sort the time vector of theC
  vector<double> v = theC._t;
  sort(v.begin(), v.end(), [](double a, double b) {return a<b;} );
  for(int i=0; i<v.size(); i++){
  	cout << v[i] << " " << theC.GetValue(v[i]) << endl;
  	gResp->SetPoint(gResp->GetN(), v[i], theC.GetValue(v[i]) );
  }


  TFile* ofile = new TFile("protodune-nl-MC-fitres.root","recreate");
  gNL->Write();
  gResp->Write();
  ofile->Close();

  return 0;
}
