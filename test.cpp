#include <memory>
#include <iostream>
#include <array2d.hpp>
#include <batch_fir.hpp>
#include <cspfb.hpp>
#include <windowed_fir.hpp>
using namespace cpfb;



int main(){
    Array2D<double> a(3,4,std::vector<double>{1,2,3,4,5,6,7,8,9,10,11,12});
    std::cout<<a<<std::endl;
    Array2D<double> b(3,4);
    a.transpose(b);
    std::cout<<b<<std::endl;
}
