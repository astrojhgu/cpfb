#ifndef OSCILLATOR_HPP
#define OSCILLATOR_HPP
#include <concepts>
#include <complex>
#include <iterator>
#include "utils.hpp"
#include "array2d.hpp"

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
        COscillator(COscillator&&)=default;
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
    struct ROscillator{
        T phase{};
        T dphi_dpt;
        T ampl;
        ROscillator(T phi0, T dphi_dpt1, T ampl1=(T)1)
        :phase(phi0), dphi_dpt(dphi_dpt1), ampl(ampl1)
        {}

        ROscillator()=delete;
        ROscillator(const ROscillator&)=default;
        ROscillator& operator=(const ROscillator& rhs)=default;

        T get(){
            auto result=ampl*std::cos(phase);
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

        HalfChShifter(HalfChShifter&&)=default;
        HalfChShifter(const HalfChShifter&)=default;
        HalfChShifter<T>& operator=(const HalfChShifter&)=default;

        HalfChShifter(size_t nch1, int direction=-1)
        :nch(nch1), factor(nch1*2), idx(0){
            COscillator<T> osc(0.0, direction*PI<T>()/nch);
            osc.fill(factor);
        }
        template <std::ranges::range R>
        void shift(R& x){
            for(auto& x1:x){
                x1*=factor[idx++];
                if (idx>=2*nch){
                    idx=0;
                }
            }
        }

        template <bool owned>
        void shift(Array2D<std::complex<T>,owned>& x){
            for(size_t j=0;j<x.ncols();++j){
                for(size_t i=0;i<x.nrows();++i){            
                    x(i,j)*=factor[idx];
                }
                ++idx;
                if(idx>=2*nch){
                    idx=0;
                }
            }
        }
    };
}


#endif
