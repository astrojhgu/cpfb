#include <memory>
#include <iostream>
#include <array2d.hpp>
#include <batch_fir.hpp>
#include <cspfb.hpp>
#include <windowed_fir.hpp>
using namespace cpfb;



int main(){
    auto nch=8;
    auto tap=12;
    auto nsignal=128;
    auto coeff1=coeff<double>(nch, tap, 0.2);

    
    CsPFB<double> pfb(nch, coeff1);

    
    std::vector<std::complex<double>> signal(nch*(tap-1)*nsignal);

    auto dphi=2*PI<double>()*-1.0/nch+0.01;

    std::complex<double> phi;

    for(auto& x:signal){
        x=std::exp(std::complex<double>(0,1)*phi);
        phi+=dphi;
    }

    
    for(int i=0;i<nsignal-1;++i){
        auto x=pfb.analyze(signal.begin()+i*nch*(tap-1), signal.begin()+(i+1)*nch*(tap-1));
        //std::cout<<x.transform([](const std::complex<double>& x)->double{return std::arg(x);})<<std::endl;
        std::cout<<x.transform((double(*)(const std::complex<double>&))std::arg)<<std::endl;
    }
    
    auto x=pfb.analyze(signal.begin()+(nsignal-1)*nch*(tap-1), signal.begin()+nsignal*nch*(tap-1));;
    
    std::cout<<x.transform([](const std::complex<double>& x)->double{return std::arg(x);})<<std::endl;
    std::cout<<x<<std::endl;

    

    std::vector<std::complex<double>> signal1(nch*(tap-1));
    pfb.analyze(signal1);
    //std::cout<<pfb.fir<<std::endl;
    
}
