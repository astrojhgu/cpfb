#include <memory>
#include <iostream>
#include <array2d.hpp>
#include <batch_fir.hpp>
#include <cspfb.hpp>
#include <windowed_fir.hpp>
#include <oscillator.hpp>
using namespace cpfb;


std::vector<double> ampl_response(double dphi_dpt, size_t nch, size_t tap, double k, double nsignal=4){
    COscillator<double> osc(0.0, dphi_dpt);
    CsPFB<double> pfb(nch, coeff<double>(nch, tap, k));
    std::vector<std::complex<double>> signal(nch*(tap-1)*nsignal);
    osc.fill(signal);
    Array2D<std::complex<double>> x(nch, (tap-1));
    for(int i=0;i<nsignal;++i){
        x=pfb.analyze(signal.begin()+i*nch*(tap-1), signal.begin()+(i+1)*nch*(tap-1));
    }
    auto spec=x.transform([](const std::complex<double>& x)->double{return std::abs(x);});
    std::vector<double> result(spec.ncols());
    for(size_t i=0;i<spec.nrows();++i){
        for(size_t j=0;j<spec.ncols();++j){
            result[j]+=spec(i,j);
        }
    }
    return result;
}



int main(){
    size_t nch=8;
    
    for(auto dphi_dpt=-PI<double>();dphi_dpt<PI<double>();dphi_dpt+=0.01){
        auto spec=ampl_response(dphi_dpt, 16, 32, 0.5, 2);
        for(auto& x: spec){
            std::cout<<x<<" ";
        }
        std::cout<<std::endl;
    }
}
