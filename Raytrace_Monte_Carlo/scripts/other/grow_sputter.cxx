#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>

using namespace std;

int main(int argc, char **argv){
  string filename, temp, temp2, temp3;
  char buffer[32];
  filename = "remeshlog";
  ofstream outfile;
  ifstream infile;
  double sputter;
  double deposit;

  snprintf(buffer, sizeof(char)*32, filename.c_str());
  infile.open(buffer);
  while(!infile.eof()){
    infile >> temp;
    if(temp == "sputter:")
      infile >> sputter;
    if(temp == "deposit:")
      infile >> deposit;
  }
  infile.close();

  int count;
  double yield=0;
  filename = "logfile";
  snprintf(buffer, sizeof(char)*32, filename.c_str());
  infile.open(buffer);
  while(!infile.eof()){
    infile >> temp >> temp2;
    getline(infile,temp3);
    if(temp2 == "yield"){
      count++;
      yield += stod(temp);
    }
  }
  infile.close()
  yield = yield/double(count);
  cout << (sputter-deposit)/sputter << "\t " << deposit/sputter << "\t" << yield << endl;
  return 0; 
}

