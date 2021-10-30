#ifndef UTILS_HPP
#define UTILS_HPP

namespace cpfb{
    template <typename T>
    constexpr T PI()
    {
        return 3.141592653589793238462643383279502884L;
    }

    template <typename T>
    std::vector<T> fftshift(const std::vector<T>& x){
        std::vector<T> y(x.size());
        auto n2=x.size()/2;
        assert(x.size()%2==0);

#ifdef FINE_PAR
#pragma omp parallel for
#endif
        for(size_t i=0;i<n2;++i){
            y[i]=x[i+n2];
            y[i+n2]=x[i];
        }
        return y;
    }
}


#endif