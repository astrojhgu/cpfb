#include <memory>
#include <iostream>
#include <fstream>
#include <cpfb/array2d.hpp>
#include <cpfb/batch_fir.hpp>
#include <cpfb/ospfb.hpp>
#include <cpfb/windowed_fir.hpp>
#include <cpfb/oscillator.hpp>
using namespace cpfb;

template <typename T>
std::vector<T> ampl_response(T dphi_dpt, size_t nch, size_t tap, T k, size_t nsignal=4){
    COscillator<T> osc(0.0, dphi_dpt);
    OsPFB<T> pfb(nch, coeff<T>(nch, tap, k));
    std::vector<std::complex<T>> signal(nch*(tap-1)*nsignal);
    //COscillator<T> osc1(0.0, PI<T>()/nch);
    osc.fill(signal);
    OsPFBOutput<double, false> channelized(nch, tap-1);
    for(int i=0;i<nsignal;++i){
        //x=pfb.analyze(signal.begin()+i*nch*(tap-1), signal.begin()+(i+1)*nch*(tap-1));
        channelized=pfb.analyze(std::span<std::complex<double>>(signal.begin()+i*nch*(tap-1), nch*(tap-1)));
    }
    auto spec_even=channelized.even.transform([](const std::complex<T>& x)->T{return std::norm(x);});
    auto spec_odd=channelized.odd.transform([](const std::complex<T>& x)->T{return std::norm(x);});

    auto phase_even=channelized.even.transform([](const std::complex<T>& x)->T{return std::arg(x);});
    auto phase_odd=channelized.odd.transform([](const std::complex<T>& x)->T{return std::arg(x);});

    std::vector<T> result(spec_even.ncols()*2);
    for(size_t i=0;i<spec_even.nrows();++i){
        for(size_t j=0;j<spec_even.ncols();++j){
            result[2*j]+=spec_even(i,j);
            result[2*j+1]+=spec_even(i,j);
        }
    }
    return result;
}



int main(){
    size_t nch=32;
    using Tfloat=double;
    std::vector<double> result;
    for(auto dphi_dpt=-PI<Tfloat>();dphi_dpt<PI<Tfloat>();dphi_dpt+=0.001){
        auto spec=ampl_response(dphi_dpt, nch,15, (Tfloat)0.5, 2);
        //std::cout<<spec.size()<<std::endl;
        std::copy(spec.begin(), spec.end(), std::back_inserter(result));
    }
    std::ofstream output("chirp_ospfb_output.bin");
    output.write((char*)result.data(), result.size()*sizeof(double));
    //auto spec=ampl_response(PI<double>()/nch*2, 4, 3, (Tfloat)0.5, 2);
}
