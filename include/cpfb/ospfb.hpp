#ifndef OSPFB_HPP
#define OSPFB_HPP
#include "cspfb.hpp"
#include "oscillator.hpp"

namespace cpfb{
    
    template <std::floating_point T>
    struct OsPFB{
        CsPFB<T> even_pfb;
        CsPFB<T> odd_pfb;
        std::vector<std::complex<T>> even_buffer;
        std::vector<std::complex<T>> odd_buffer;
        HalfChShifter<T> shifter;

        OsPFB()=delete;
        OsPFB(const OsPFB<T>&)=delete;
        OsPFB<T>& operator=(const OsPFB<T>& rhs)=delete;



        template <std::ranges::range H>
        OsPFB(size_t nch1, const H& h0)
        :even_pfb(nch1, h0), odd_pfb(nch1, h0),
        even_buffer(even_pfb.size_per_shoot()),
        odd_buffer(even_pfb.size_per_shoot()),
        shifter(nch1)
        {
        }

        std::pair<Array2D<std::complex<T>, false>,
        Array2D<std::complex<T>, false>> analyze(std::span<std::complex<T>> data){
            assert(data.size()==even_buffer.size());
            std::copy(data.begin(), data.end(), even_buffer.begin());
            std::copy(data.begin(), data.end(), odd_buffer.begin());
            shifter.shift(odd_buffer);
            auto result_even=even_pfb.analyze_insitu(even_buffer);
            auto result_odd=odd_pfb.analyze_insitu(odd_buffer);
            return std::make_pair(std::move(result_even), std::move(result_odd));
        }
    };
}

#endif
