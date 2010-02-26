#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>


using namespace std;


int main() {


     double quadrature, phase;

     ifstream datafile("data.txt");
     // double[] tmp(2);
     
     if (!datafile) {
          cout << "Cannot open file" << endl;
          return 1;
     }

     vector<vector<double> > data;
     data.reserve(1000000);

     while (!datafile.eof()) {
          datafile >> quadrature >> phase;
          data.push_back(vector<double>(2));
          data.back().at(0) = quadrature;
          data.back().at(1) = phase;
     }

     for(vector<double>::iterator i=data.begin(); i != data.end(); ++i) {
          cout << "Quadrature: " << *i.at(0) <<  "  Phase:  " << *i.at(1) << endl;
     }


     cout << data.size() << endl;
     cout << data.max_size() << endl;
     cout << data.capacity() << endl;

     cout << "ok" << endl;

     datafile.close();

     return 0;
}
