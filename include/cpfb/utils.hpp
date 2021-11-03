#ifndef UTILS_HPP
#define UTILS_HPP
#include <ranges>
#include "array2d.hpp"

namespace cpfb{
    template <typename T>
    constexpr T PI()
    {
        return 3.141592653589793238462643383279502884L;
    }

    template <typename T>
    std::vector<T> fftshift(std::span<T> x){
        std::vector<T> y(x.size());
        auto n2=x.size()/2;
        assert(x.size()%2==0);
        for(size_t i=0;i<n2;++i){
            y[i]=x[i+n2];
            y[i+n2]=x[i];
        }
        return y;
    }

    template <typename T>
    void fftshift_insitu(std::span<T> x){
        auto n2=x.size()/2;
        assert(x.size()%2==0);
        for(size_t i=0;i<n2;++i){
            std::swap(x[i],x[i+n2]);
        }
    }

    template <typename T, bool owned>
    void fftshift_insitu(Array2D<T, owned>& x){
        size_t ncols=x.ncols();
        size_t n2=ncols/2;
        assert(ncols%2==0);
        for(int i=0;i<x.nrows();++i){
            for(size_t j=0;j<n2;++j){
                std::swap(x(i, j), x(i, j+n2));
            }
        }
    }
}

#endif