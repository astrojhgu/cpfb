#include <memory>
#include <iostream>
#include <fstream>
#include <cpfb/array2d.hpp>
#include <cpfb/batch_fir.hpp>
#include <cpfb/ospfb.hpp>
#include <cpfb/windowed_fir.hpp>
#include <cpfb/oscillator.hpp>
#include <cpfb/csp_channelizer.hpp>
#include <cpfb/chained_pfb.hpp>
using namespace cpfb;

template <typename T>
void phase_response(T dphi_dpt, size_t nch, size_t nch_fine, size_t tap, size_t tap_fine, T k, T k_fine, size_t nsignal=4){
    COscillator<T> osc(0.0, dphi_dpt);
    OsPFB<T> pfb(nch, coeff<T>(nch, tap, k));
    std::vector<std::complex<T>> signal(nch*(tap-1)*nsignal);
    //COscillator<T> osc1(0.0, PI<T>()/nch);
    osc.fill(signal);
    auto channelized=std::make_pair(
        Array2D<std::complex<T>, false>(nch, (tap-1), nullptr)
        ,Array2D<std::complex<T>, false>(nch, (tap-1), nullptr));

    CSPChannelizer<T> csp(nch, nch_fine*2, coeff<T>(nch_fine*2, tap_fine, k_fine), {4,5,6});
    
    for(int i=0;i<nsignal;++i){
        //x=pfb.analyze(signal.begin()+i*nch*(tap-1), signal.begin()+(i+1)*nch*(tap-1));
        channelized=pfb.analyze(std::span<std::complex<double>>(signal.begin()+i*nch*(tap-1), nch*(tap-1)));
        auto& even=channelized.first;
        auto& odd=channelized.second;
        csp.feed_station_data(even, odd);
        while(csp.fetch()){
            for(int j=0;j<csp.output_buffer.ncols();++j){
                for(auto k:std::vector<int>{1,2,3,4,5,6,7,8}){
                    std::cerr<<std::abs(csp.output_buffer(k,j))<<" ";
                }
                std::cerr<<std::endl;
            }
            
        }
    }
}



int main(){
    size_t nch=8;
    using Tfloat=double;
    phase_response(PI<double>()/nch*5+0.1, nch, 4, 15, 4, (Tfloat)0.5, (Tfloat)0.5, 120);
}
