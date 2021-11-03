#ifndef OSPFB_HPP
#define OSPFB_HPP
#include "cspfb.hpp"
#include "oscillator.hpp"
#include <concepts>

namespace cpfb{
    template <std::floating_point T>
    std::complex<T> to_complex(const T& x){
        return std::complex<T>(x);
    }

    template <std::floating_point T>
    std::complex<T> to_complex(const std::complex<T>& x){
        return x;
    }

    template <std::floating_point T>
    struct OsPFB{
        CsPFB<T> even_pfb;
        CsPFB<T> odd_pfb;
        std::vector<std::complex<T>> even_buffer;
        std::vector<std::complex<T>> odd_buffer;
        HalfChShifter<T> shifter;

        OsPFB()=delete;
        OsPFB(const OsPFB<T>&)=delete;
        OsPFB(OsPFB&& rhs)
        :even_pfb(std::move(rhs.even_pfb))
        ,odd_pfb(std::move(rhs.odd_pfb))
        , even_buffer(std::move(rhs.even_buffer))
        , odd_buffer(std::move(rhs.odd_buffer))
        , shifter(std::move(rhs.shifter))
        {}
        OsPFB<T>& operator=(const OsPFB<T>& rhs)=delete;



        template <std::ranges::range H>
        OsPFB(size_t nch1, const H& h0)
        :even_pfb(nch1, h0), odd_pfb(nch1, h0),
        even_buffer(even_pfb.size_per_shoot()),
        odd_buffer(even_pfb.size_per_shoot()),
        shifter(nch1)
        {
        }

        template <typename U>
        std::pair<Array2D<std::complex<T>, false>,
        Array2D<std::complex<T>, false>> analyze(std::span<U> data){
            assert(data.size()==even_buffer.size());
            //std::copy(data.begin(), data.end(), even_buffer.begin());
            std::transform(data.begin(), data.end(), even_buffer.begin(), (std::complex<T>(*)(const U&))to_complex);
            std::copy(even_buffer.begin(), even_buffer.end(), odd_buffer.begin());
            shifter.shift(odd_buffer);
            auto result_even=even_pfb.analyze_insitu(even_buffer);
            auto result_odd=odd_pfb.analyze_insitu(odd_buffer);
            return std::make_pair(std::move(result_even), std::move(result_odd));
        }
    };
}

#endif
