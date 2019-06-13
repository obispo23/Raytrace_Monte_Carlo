#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<random>
#include<ctime>
//#include"vec_bvh.h"

//Returns Energy yield, arguements: Primary energy, Incident Angle.
double source_f_energy(double A, double B, double& C, string spec1, string spec2){
  //---------------------------species 1 to 2-----------------------------
  if(spec1 == "Xe" && spec2 == "W"){
    C = 0;
    B = B * 180.0 / M_PI; 
    C = 31.38+0.01812*A-569.7/A+0.002789*B*B+0.00000009318*A*B*B*B-0.1*B-0.0000002533*B*A*A;
    return C;
  }
  if(spec1 == "Ar" && spec2 == "W"){
    double M1 = 39.948; // projectile atomic mass (Ar)
    double M2 = 183.84; // target atomic mass (W)  
    double gamma = 4*M1*M2/(pow((M1+M2),2.0)); //# energy transfer factor
    double Us = 8.9; //# sublimation energy [eV]
    double setrange = A/2.0;
    double value[1000];
    double sum[1000];
    double r = unifrnd(generator);
    double Energy;
    double total = 0;
    for(int i = 0; i < 1000; i++){
      Energy = setrange*double(i)/1000.0;
      sum[i] = total+(Energy)/pow((Energy+Us),3.0) * (1.0-sqrt((Energy+Us)/(gamma*A)));
      value[i] = (Energy)/pow((Energy+Us),3.0) * (1.0-sqrt((Energy+Us)/(gamma*A)));
      total = sum[i];
    }
    double hold = 0;
    double start = 1.0;
    for(int i  = 0; i < 1000; i++){
      double diff = sum[i]/total - r;
      if(abs(diff) < start){
	start = abs(diff);
	hold = i;
	Energy = setrange*double(hold)/1000.0;
      }
    }

    C = Energy;
    if(C < eps || isnan(C) )
      C = 0;
    return C;

  }
  //--------------------------species 2 to 2------------------------------  
  if(spec1 == "W" && spec2 == "W"){
    C = 0;
    B = B * 180.0 / M_PI; 
    if(A >= 70)
      C = 8.33518957053645 + 0.0707587988443273*A + 0.000119442008784669*pow(B,3) + -3.32685474492003e-5*pow(B,4)/(10.584377161381 + A) - 4.38325075068488e-7*pow(B,4) - 2.6108244601802e-6*B*A*A;
    if(A < 70)
      C = 0.21214123247341*A + 0.00997059596574171*A*B + 1.92676159301919e-6*A*A*B*B + sin(A + 5.02255050678406e-9*A*pow(B,4)*sin(1.92676159301919e-6*A*A*B*B) - 5.02255050678406e-9*A*pow(B,4)) - 0.000199809153144667*B*A*A;
    return C;
  }
    
}
//Returns Electron yield, arguments: Primary Energy, Indicent Angle.
double source_f_yield(double A, double B, double& C, string spec1, string spec2){
  //---------------------------species 1 to 2-----------------------------
  if(spec1 == "Xe" && spec2 == "W"){
    C = 0;
    B = B * 180.0 / M_PI; 
    C = 0.07246+0.002236*A+0.00004875*B*B+0.0000004775*A*B*B-0.005515*B-0.0000004702*A*A-0.00000000006214*A*B*B*B*B;
    return C;
  }
  if(spec1 == "Ar" && spec2 == "W"){
    double E = A;
    double alpha = B;
    double M1 = 39.948; // projectile atomic mass (Ar)
    double M2 = 183.84;// target atomic mass (W)
    double Z1 = 18.0; //# projectile atmonic number (Ar)
    double Z2 = 74.0; //# target atomic number (W)
    double K = 8.478*Z1*Z2/(pow((pow(Z1,(2.0/3.0))+pow(Z2,(2.0/3.0))),0.5))*M1/(M1+M2);
    double Q = 1.1;
    double Us = 8.9; //# sublimation energy [eV]
    double alpha_star = (0.08+0.164*pow((M2/M1),0.4)+0.0145*pow((M2/M1),1.29));
    double Eth = (1.9+3.8*(M1/M2)+0.134*pow((M2/M1),1.24))*Us;
    
    //def Y_Yamamura(E,alpha):
    //    def Y_Yamamura_energy(E):
    double epsilon = 0.03255/(Z1*Z2*pow((pow(Z1,(2.0/3.0))+pow(Z2,(2.0/3.0))),0.5))*M2/(M1+M2)*E;
    double k = 0.079*pow((M1+M2),1.5)/(pow(M1,1.5)*pow(M2,0.5))*pow(Z1,(2.0/3.0))*pow(Z2,0.5)/(pow((pow(Z1,(2.0/3.0))+pow(Z2,(2.0/3.0))),0.75));
    double sn = 3.441*sqrt(epsilon)*log(epsilon+2.718)/(1.0+6.355*sqrt(epsilon)+epsilon*(-1.708+6.882*sqrt(epsilon)));
    double se = k*sqrt(epsilon);
    double Y_Yamamura_energy = 0.42*alpha_star*Q*K*sn/(Us*(1.0+0.35*Us*se))*pow((1.0-sqrt(Eth/E)),2.8);
    //    def Y_Yamamura_angular(alpha):
    double fs = 1.28; //# Sigmund f
    double gamma = 4*M1*M2/(pow((M1+M2),2.0)); //# energy transfer factor
    double h = 0.834;
    double E_th = 1.5*Us/gamma*pow((1+1.38*pow((M1/M2),h)),2.0);
    double a_0 = 0.529; //     # Bohr radius [A]
    double a_L = 0.8853/sqrt(pow(Z1,(2.0/3.0))+pow(Z2,(2.0/3.0)))*a_0; //# Linhard screening radius [A]
    double R0 = 3.1652; //# average lattice constant [A]
    double phi = pow((a_L/R0),1.5)*sqrt(Z1*Z2/( sqrt(pow(Z1,(2.0/3.0))+pow(Z2,(2.0/3.0)) )*E));
    double theta_opt = M_PI/2-286.0/180.0*M_PI*pow(phi,0.45);
    double eta = 1-sqrt(E_th/E);
    double f = fs*(1+2.5*((1-eta)/eta));
    double SUM = f*cos(theta_opt);
    double Y_Yamamura_angular =  pow(cos(alpha),-f)*exp(f*cos(theta_opt)*(1.0-(pow(cos(alpha),-1.0))));
    C =  Y_Yamamura_energy*Y_Yamamura_angular;
    if(C < eps || isnan(C) )
      C = 0;
    //cout << C << " " << alpha << " " << E << endl;
    return Y_Yamamura_energy*Y_Yamamura_angular;
  }
  //--------------------------species 2 to 2------------------------------  
  if(spec1 == "W" && spec2 == "W"){
    C = 0;
    B = B * 180.0 / M_PI;                                                                                                                     
    if(A >= 70)
      C = 0.0019022408974498*A + 0.0315762558647962*cos(5.63199894492506 + 0.0593689474611619*B) + 0.0301302565807804*B/(A + 0.0593689474611619*B*cos(0.350081032392476 + B)) - 0.0824151246397562 - 0.000630089851589718*A*cos(5.39549024714073 + 0.0632916176908535*B);
    if(A < 70)
      C = 1.22455927394761e-5*A*A + 2.58328306854474e-7*A*B*B + 3.75197223368283e-9*A*A*B*B - sin(5.93895283013407e-6*A*B) - A*sin(6.07517493158528e-11*A*B*B*B) - 0.000214533716629772*A;
    if(C < eps)
      C = 0;
    return C;
  }
}

//Energy yield, arguements: Primary energy, Incident angle.
double dyieldf(double A, double B, double& C){
 
}
//Electron yield, arguments: Primary Energy, Indicent Angle. 
double dyieldd(double A, double B, double& C){
}
//---------------------------------------------------------------  

