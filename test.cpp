#include <memory>
#include <iostream>
#include <cpfb/array2d.hpp>
#include <cpfb/batch_fir.hpp>
#include <cpfb/cspfb.hpp>
#include <cpfb/windowed_fir.hpp>
#include <cpfb/oscillator.hpp>
using namespace cpfb;



int main(){
    Array2D<double> a(3,4,std::vector<double>{1,2,3,4,5,6,7,8,9,10,11,12});
    std::cout<<a<<std::endl;
    Array2D<double> b(3,4);
    a.transpose(b);
    std::cout<<b<<std::endl;
}
