#ifndef OSCILLATOR_HPP
#define OSCILLATOR_HPP
#include <concepts>
#include <complex>
#include <iterator>

namespace cpfb{
    template <std::floating_point T>
    struct COscillator{
        T phase{};
        T dphi_dpt;
        T ampl;
        COscillator(T phi0, T dphi_dpt1, T ampl1=(T)1)
        :phase(phi0), dphi_dpt(dphi_dpt1), ampl(ampl1)
        {}

        COscillator()=delete;
        COscillator(const COscillator&)=default;
        COscillator& operator=(const COscillator& rhs)=default;

        std::complex<T> get(){
            auto result=std::complex<T>(ampl*std::cos(phase), ampl*std::sin(phase));
            phase+=dphi_dpt;
            return result;
        }

        template <std::forward_iterator U>
        void fill(U b, U e){
            for(auto i=b;i!=e;++i){
                *i=get();
            }
        }

        template <std::ranges::range R>
        void fill(R& x){
            fill(x.begin(), x.end());
        }
    };

    template <std::floating_point T>
    struct HalfChShifter{
        size_t nch;
        std::vector<std::complex<T>> factor;
        size_t idx;

        HalfChShifter(size_t nch1, int direction=-1)
        :nch(nch1), factor(nch1*2), idx(0){
            COscillator<T> osc(0.0, direction*PI<T>()/nch);
            osc.fill(factor);
        }
        template <std::ranges::range R>
        void shift(R& x){
            for(auto& x1:x){
                x1*=factor[(idx++)%(2*nch)];
                if (idx>=2*nch){
                    idx=0;
                }
            }
        }
    };
}


#endif
