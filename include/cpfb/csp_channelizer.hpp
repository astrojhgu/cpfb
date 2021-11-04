#ifndef CSP_CHANNELIZER
#define CSP_CHANNELIZER

#include "ospfb.hpp"
#include "cspfb.hpp"
#include "oscillator.hpp"

namespace cpfb{
    template <typename T>
    struct RowBuffer{
        const size_t ncols;
        const size_t output_nrows;
        size_t input_idx;
        size_t output_idx;
        Array2D<T, false> source;
        Array2D<T, true> buffer;

        RowBuffer(size_t nrows1, size_t ncols1)
        :ncols(ncols1),output_nrows(nrows1), input_idx(0), output_idx(0), source(),buffer(nrows1, ncols1){
        }

        template <bool owned>
        void feed(Array2D<T, owned>& food){
            buffer.reshape(output_nrows, ncols);
            input_idx=0;
            source=food;
        }

        bool fetch(){
            buffer.reshape(output_nrows, ncols);
            size_t nrows_needed=output_nrows-output_idx;
            size_t nrows_available=source.nrows()-input_idx;
            size_t nrows_to_fetch=std::min(nrows_needed, nrows_available);
            //std::cerr<<nrows_to_fetch<<std::endl;
            std::copy(source.data.get()+input_idx*ncols, source.data.get()+(input_idx+nrows_to_fetch)*ncols, buffer.data.get()+output_idx*ncols);
            input_idx+=nrows_to_fetch;
            output_idx+=nrows_to_fetch;
            
            if(output_idx==output_nrows){
                output_idx=0;
                return true;
            }else{
                return false;
            }
        }
    };

    template <typename T>
    struct CSPChannelizer{
        size_t nch_coarse;//half of total ch of ospfb, including pos and neg channels,
        size_t nch_fine_per_coarse;//both kept and dropped are included
        std::vector<size_t> ch_coarse_selected;
        std::vector<CsPFB<T>> secondary_pfbs;
        size_t secondary_batch;
        HalfChShifter<T> even_half_ch_shifter;
        HalfChShifter<T> odd_half_ch_shifter;
        RowBuffer<std::complex<T>> even_buffer;
        RowBuffer<std::complex<T>> odd_buffer;
        Array2D<std::complex<T>> output_buffer;
        

        template <std::ranges::range H>
        CSPChannelizer(size_t nch_c, size_t all_fine_ch_per_coarse, const H h0, const std::vector<size_t>& ch_coarse_selected1)
        //all_fine_ch_per_coarse is the final kept channels, so should be x2 to build fft
        :nch_coarse(nch_c),nch_fine_per_coarse(all_fine_ch_per_coarse)
        , ch_coarse_selected(ch_coarse_selected1)
        , secondary_pfbs(), 
        secondary_batch((h0.end()-h0.begin())-all_fine_ch_per_coarse),
        even_half_ch_shifter(all_fine_ch_per_coarse), 
        odd_half_ch_shifter(all_fine_ch_per_coarse), 
        even_buffer(secondary_batch, nch_c),
        odd_buffer(secondary_batch, nch_c)
        , output_buffer(ch_coarse_selected.size()*all_fine_ch_per_coarse/2, secondary_batch/all_fine_ch_per_coarse)
        //output_buffer(ch_coarse_selected)
        {
            assert(all_fine_ch_per_coarse%4==0);
            for(auto& c:ch_coarse_selected1){
                secondary_pfbs.emplace_back(all_fine_ch_per_coarse, h0);
            }
            assert(secondary_pfbs[0].size_per_shoot()==even_buffer.output_nrows);
        }

        template <bool owned>
        void feed_station_data(OsPFBOutput<T, owned>& data){
            auto& even=data.even;
            auto& odd=data.odd;
            assert(even.ncols()==odd.ncols());
            assert(even.nrows()==odd.nrows());
            assert(even.ncols()==nch_coarse);
            even_buffer.feed(even);
            odd_buffer.feed(odd);
        }

        bool fetch(){
            auto be=even_buffer.fetch();
            auto bo=odd_buffer.fetch();
            assert(be==bo);
            if(!be){
                return false;
            }

            
            even_buffer.buffer.transpose_self();
            odd_buffer.buffer.transpose_self();
            assert(even_buffer.buffer.nrows()==nch_coarse);
            assert(even_buffer.buffer.ncols()==secondary_batch);
            
            even_half_ch_shifter.shift(even_buffer.buffer);
            odd_half_ch_shifter.shift(odd_buffer.buffer);
            output_buffer.reshape(ch_coarse_selected.size()*nch_fine_per_coarse/2, secondary_batch/nch_fine_per_coarse);
            for(size_t i=0;i<ch_coarse_selected.size();++i){
                
                auto ch=ch_coarse_selected[i];
                Array2D<std::complex<T>, false> y;
                if(ch%2==0){
                    y=secondary_pfbs[i].analyze_insitu(std::span(even_buffer.buffer.data.get()+ch/2*secondary_batch, secondary_batch));

                }else{
                    y=secondary_pfbs[i].analyze_insitu(std::span(odd_buffer.buffer.data.get()+ch/2*secondary_batch, secondary_batch));
                }
                fftshift_insitu(y);
                y.transpose_self();
                
                Array2D<std::complex<T>, false> y1(nch_fine_per_coarse/2, secondary_batch/nch_fine_per_coarse, 
                output_buffer.data.get()+i*secondary_batch/2
                );
                std::copy(y.data.get()+secondary_batch/4, y.data.get()+secondary_batch/4+y1.size(), y1.data.get());
            }
            output_buffer.transpose_self();
            assert(output_buffer.ncols()==nch_fine_per_coarse*nch_coarse/2);
            return true;
        }
    };
}

#endif
