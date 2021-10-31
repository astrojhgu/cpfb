#include <memory>
#include <iostream>
#include <cpfb/array2d.hpp>
#include <cpfb/batch_fir.hpp>
#include <cpfb/cspfb.hpp>
#include <cpfb/windowed_fir.hpp>
#include <cpfb/oscillator.hpp>
using namespace cpfb;

template <typename T>
std::vector<T> ampl_response(T dphi_dpt, size_t nch, size_t tap, T k, size_t nsignal=4){
    COscillator<T> osc(0.0, dphi_dpt);
    CsPFB<T> pfb(nch, coeff<T>(nch, tap, k));
    std::vector<std::complex<T>> signal(nch*(tap-1)*nsignal);
    osc.fill(signal);
    Array2D<std::complex<T>> x(nch, (tap-1));
    for(int i=0;i<nsignal;++i){
        x=pfb.analyze(signal.begin()+i*nch*(tap-1), signal.begin()+(i+1)*nch*(tap-1));
    }
    auto spec=x.transform([](const std::complex<T>& x)->T{return std::abs(x);});
    std::vector<T> result(spec.ncols());
    for(size_t i=0;i<spec.nrows();++i){
        for(size_t j=0;j<spec.ncols();++j){
            result[j]+=spec(i,j);
        }
    }
    return result;
}



int main(){
    size_t nch=8;
    using Tfloat=double;

    
    for(auto dphi_dpt=-PI<Tfloat>();dphi_dpt<PI<Tfloat>();dphi_dpt+=0.001){
        auto spec=ampl_response(dphi_dpt, nch, 64, (Tfloat)0.3, 32);
        for(auto& x: spec){
            std::cout<<x<<" ";
        }
        std::cout<<std::endl;
    }
    //auto spec=ampl_response(PI<double>()/nch*2, 4, 3, (Tfloat)0.5, 2);
}
