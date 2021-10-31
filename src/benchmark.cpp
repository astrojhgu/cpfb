#include <memory>
#include <iostream>
#include <cpfb/array2d.hpp>
#include <cpfb/batch_fir.hpp>
#include <cpfb/cspfb.hpp>
#include <cpfb/windowed_fir.hpp>
#include <cpfb/oscillator.hpp>
using namespace cpfb;
using Tfloat=double;
template <typename T>
void bench(std::vector<std::complex<double>> signal, size_t nch, size_t tap, T k){
    CsPFB<T> pfb(nch, coeff<T>(nch, tap, k));
    auto batch=tap-1;
    auto nsignal=signal.size()/nch/batch;
    assert(nch*batch*nsignal==signal.size());
    for(size_t i=0;i<nsignal;++i){
        if (i%100==0){
            std::cerr<<i<<" / "<<nsignal<<std::endl;
        }
        
        pfb.analyze_insitu(std::span<std::complex<double>>(signal.begin()+i*nch*batch, nch*batch));
    }
    
    
}



int main(){
    constexpr size_t nch=1024;
    constexpr size_t tap=9;
    constexpr size_t batch=tap-1;
    constexpr size_t nsignal=2000;
    std::cerr<<"signal len="<<nsignal*batch*nch/1e6<<" Mpts"<<std::endl;
    double dphi_dpt=PI<double>()/8;
    COscillator<double> osc(0.0, dphi_dpt);
    std::vector<std::complex<double>> signal(nch*batch*nsignal);
    osc.fill(signal);
    for(size_t i=0;i<10;++i){
        bench(signal, nch, tap, (Tfloat)0.3);
    }
    
    //auto spec=ampl_response(PI<double>()/nch*2, 4, 3, (Tfloat)0.5, 2);
}
