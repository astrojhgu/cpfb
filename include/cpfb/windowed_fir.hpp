#ifndef WINDOWED_FIR
#define WINDOWED_FIR
#include <cmath>
#include <vector>
#include <concepts>
#include <numeric>
#include "utils.hpp"

namespace cpfb{
    template <typename T>
    T blackman_window(size_t i, size_t n){
        constexpr T a0=0.3635819;
        constexpr T a1=0.4891775;
        constexpr T a2=0.1365995;
        constexpr T a3=0.0106411;
        auto x=(T)i/(T)n*PI<T>();
        return a0-a1*std::cos(2*x)+a2*std::cos(4*x)-a3*std::cos(6*x);
    }

    template <typename T>
    T raised_cosine(size_t i, size_t n){
        auto x=std::sin(PI<T>()*i/n);
        return x*x;
    }

    template <typename T>
    T blackman_nuttall(size_t i, size_t n){
        T a0=0.3635819;
        T a1=0.4891775;
        T a2=0.1365995;
        T a3=0.0106411;
        return a0-a1*cos(2*PI<T>()*i/n)+a2*cos(4*PI<T>()*i/n)-a3*cos(6*PI<T>()*i/n);
    }

    template <std::ranges::range R, typename F>
    void apply_window(R& workspace, F&& f){
        auto iter=workspace.begin();
        auto end=workspace.end();
        size_t n=end-iter;
        size_t i=0;
        for(;iter!=end;++i,++iter){
            *iter*=f(i, n);
        }
    }

template <typename T>
    std::vector<T> lp_coeff(size_t nch, size_t l, T k){
        std::vector<T> result(l*nch);
        for(size_t i=0;i<l*nch;++i){
            if((T)i<=l*k || (T)i>=(T)(l*nch)-l*k){
                result[i]=1.0;
            }else{
                result[i]=0.0;
            }
        }
        return result;
    }

    template <typename T>
    void symmetrize(std::vector<T>& workspace){
        auto n1=workspace.size();
        for(size_t n=0;n<=(n1/2-2);++n){
            workspace[n1-n-1]=workspace[n+1];
        }
    }

    template <typename T>
    std::vector<T> to_time_domain(std::vector<T>& input){
        std::vector<std::complex<T>> a(input.size());
        for(int i=0;i<input.size();++i){
            a[i]=std::complex<T>(input[i]);
        }
        auto plan=FftwTraits<T>::plan_dft_1d(input.size(), FFTW_BACKWARD, FFTW_ESTIMATE);
        FftwTraits<T>::execute_dft(plan, a.data(), a.data());
        std::vector<T> b;
        std::transform(a.begin(), a.end(), std::back_inserter(b), [](auto x){return x.real();});
        auto s=std::accumulate(b.begin(), b.end(), (T)0);
        std::transform(b.begin(), b.end(), b.begin(), [=](T x){return x/s;});
        return fftshift(std::span(b.begin(), b.size()));
    }

    template <typename T>
    std::vector<T> coeff(size_t nch, size_t tap_per_ch, T k=1.0){
        auto a=lp_coeff(nch, tap_per_ch, k);
        symmetrize(a);
        auto b=to_time_domain(a);
        //apply_blackman_window(b);
        //apply_window(b, (T(*)(size_t,size_t))(blackman_window));
        apply_window(b, (T(*)(size_t,size_t))(blackman_nuttall));
        return b;
    }    
}


#endif
