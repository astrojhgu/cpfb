#include <memory>
#include <iostream>
#include <cpfb/array2d.hpp>
#include <cpfb/batch_fir.hpp>
#include <cpfb/ospfb.hpp>
#include <cpfb/windowed_fir.hpp>
#include <cpfb/oscillator.hpp>
using namespace cpfb;


int main(int argc,char* argv[]){
    if(argc!=4){
        std::cout<<"Usage: "<<argv[0]<<" <nch> <tap> <k>";
        return -1;
    }

    int nch=std::stoi(argv[1]);
    int tap=std::stoi(argv[2]);
    double k=std::stod(argv[3]);
    auto c=coeff<double>(nch, tap, k);
    for(auto& x: c){
        std::cout<<x<<std::endl;
    }
}