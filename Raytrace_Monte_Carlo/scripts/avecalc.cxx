#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>
#include"vec_bvh.h"

using namespace std;

int main(int argc, char **argv)
{
  ofstream outfile;
  ifstream infile;
  int n =0;
  string temp1, temp2;
  double numt;
  double Nstep = 1E6; // total number of beams
  int par = 10;
  int tn = 0;
  char buffer[32];
  string input, first, last;
  string filename = "logfile";
  if(argc > 1){
    for(int i = 1; i < argc; i++){
      input = argv[i];
      grabber(input, first, last, '=');
      if(first == "filename")
        filename = last;
      if(first == "Nstep")
	Nstep = stod(last);   
      if(first == "par")
	par = stoi(last);

    }
  }
  double number[par];
  snprintf(buffer, sizeof(char)*32, filename.c_str());
  infile.open(buffer);
  if(infile.good() == false){
    cout << buffer << " not found " << endl;
    return 1;
  }

  //for logfile use only modify for other stuff
  /*  for(int i =0; i < par; i++){
    n=0;
    while (n < 12)
      {
	if( n == 1)
	  infile >> tn; 
	  if(n == 7)
	  {
	    infile >> number[i];
	    numt = numt + number[i];	    
	    getline(infile,temp1);
	  }	    
	else
	  getline(infile,temp1);
	
	n++;
      }
      }*/
  n = 0;
  numt = 0;
  while(!infile.eof()){
    infile >> temp1 >> temp2;
    if(temp2 == "yield"){
      number[n] = stod(temp1);
      numt = numt + number[n];
      n++;
    }
    if(temp2 == "triangles")
      tn = stoi(temp1);
    getline(infile,temp1);
  }

  double std = 0.0;
  for(int i = 0; i < par; i++){
    //    cout << number[i]/(Nstep/par) - numt/Nstep<< endl;
    std = std + pow(number[i]/(Nstep/par) - numt/Nstep,2);
  }
  std = sqrt(std/(double(par-1)));
  //  cout << numt << " secondary electrons emitted" << endl;
  //cout << numt/Nstep << " secondary electron yield" << endl;
  //cout << std << " standard deviation" << endl;
  //  cout << "SEY and STD" << endl;
  infile.close();
  double distsee = 0, distseeave = 0;
  infile.open("distribution.dat");
  if(infile.good() == false){
    //cout << "distribution.dat" << " not found " << endl;
    //return 1;
  }
  else{
      for(int i=0; i < int(par*Nstep); i++){
	infile >> temp1 >> temp1 >> distsee; 
	distseeave = distseeave + distsee;
	getline(infile,temp1);
      }
      distseeave = distseeave * 1/(double(par)*Nstep);
      infile.close();
  }
  cout << "    " << numt/double(par) << "    " << std << "    " << distseeave << endl;
  
  infile.open("in_out.dat");
  if(infile.good() == false){
    cout << "in_out.dat" << " not found " << endl;
    return 1;
  }

  struct in_out{
    int in=0; 
    int out=0;
  };
  int idhold = 0;
  int inhold = 0;
  int outhold = 0;
  in_out *data;
  data = new in_out[tn];
  for(int i = 0; i < tn; i++){
    data[i].in = 0;
  data[i].out = 0;
  }
  //  while(!infile.eof()){
  for(int i = 0; i < tn*par; i++){
    infile >> idhold >> inhold >> outhold;
    getline(infile, temp1);
    data[idhold].in += inhold;
    data[idhold].out += outhold;
  }

  infile.close();
  outfile.open("in_out.dat");

  for(int i = 0; i < tn; i++)
    outfile << i << "\t" << data[i].in << "\t" << data[i].out << endl;
  outfile.close();
  return 0;
}
