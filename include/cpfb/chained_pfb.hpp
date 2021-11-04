#ifndef CHAINED_PFB_HPP
#define CHAINED_PFB_HPP
#include <cpfb/array2d.hpp>
#include <cpfb/batch_fir.hpp>
#include <cpfb/ospfb.hpp>
#include <cpfb/windowed_fir.hpp>
#include <cpfb/oscillator.hpp>
#include <cpfb/csp_channelizer.hpp>



namespace cpfb{

    template <std::floating_point T>
    struct ChainedPFB{
        std::vector<size_t> ch_selected;
        std::vector<OsPFB<T>> pfb;
        CSPChannelizer<T> csp;
        std::vector<std::span<std::complex<T>>> input_ptrs;
        size_t input_size;


        template <std::ranges::range R>
        ChainedPFB(size_t nants, size_t nch_coarse, const R& h0_coarse, size_t nch_fine_per_coarse_kept, const R& h0_fine, const std::vector<size_t>& ch_selected1)
        :ch_selected(ch_selected1)
        ,csp(nch_coarse, nch_fine_per_coarse_kept*2, h0_fine, ch_selected)
        , input_ptrs(nants)
        , input_size(0)
        {
            for(size_t i=0;i<nants;++i){
                pfb.emplace_back(nch_coarse, h0_fine);
            }
        }

        ChainedPFB(ChainedPFB&&)=default;
        ChainedPFB(const ChainedPFB&)=delete;
        ChainedPFB& operator=(const ChainedPFB&)=delete;

        Array2D<std::complex<T>>& output_buffer(){
            return csp.output_buffer;
        }

        size_t chunk_size()const{
            return pfb[0].size_per_shoot();
        }

        void feed(std::vector<std::span<std::complex<T>>>&& input_ptrs){
            assert(input_ptrs[0].size()%chunk_size()==0);
            this->input_ptrs=input_ptrs;
            input_size=0;
        }

        bool fetch(){
            if(csp.fetch()){
                return true;
            }
            if(input_size==input_ptrs[0].size()){
                return false;
            }else{
                auto channelized=pfb[0].analyze(input_ptrs[0]);
                auto& even=channelized.first;
                auto& odd=channelized.second;
                for(size_t n=1;n<pfb.size();++n){
                    auto channelized1=pfb[n].analyze(input_ptrs[n]);
                    auto& even1=channelized1.first;
                    auto& odd1=channelized1.second;
                    for(size_t i=0;i<even.size();++i){
                        even.data[i]+=even1.data[i];
                        odd.data[i]+=odd1.data[i];
                    }
                }

                csp.feed_station_data(even, odd);
                csp.fetch();
                return true;
            }
        }

    };
}

#endif
