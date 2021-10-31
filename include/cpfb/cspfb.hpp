#ifndef CSPFB_HPP
#define CSPFB_HPP
#include <complex>
#include <concepts>
#include "fft_wrapper.hpp"
#include "batch_fir.hpp"

//#define INPLACE_ALGO

namespace cpfb{
    template <std::floating_point T>
    struct CsPFB{
        size_t nch;
        size_t batch;
        BatchFIR<T, std::complex<T>> fir;
        typename FftwTraits<T>::UniquePtrPlan plan;
#ifndef INPLACE_ALGO
        Array2D<std::complex<T>> buffer;
#endif

        typename FftwTraits<T>::UniquePtrPlan init_plan(size_t nch, size_t batch){
            int n[]={(int)nch};
            int howmany=(int)batch;
            int idist=(int)nch;
            int odist=(int)nch;
            int istride=1;
            int ostride=1;
            int* inembed=n;
            int* onembed=n;
            int rank=1;
            return FftwTraits<T>::plan_many_dft(rank, n, howmany,
            inembed, 
            istride, idist,
            onembed, 
            ostride, odist, 
            FFTW_FORWARD, FFTW_ESTIMATE);
        }

        template <std::ranges::range H>
        Array2D<T> calc_coeffs(size_t nch, const H& h0){
            auto tap=(h0.end()-h0.begin())/nch;
            assert(tap*nch==(h0.end()-h0.begin()));
            Array2D<T> c(tap, nch, h0);
            c.transpose_self();
            c.reverse_row_self();
            return c;
        }

        CsPFB()=delete;
        CsPFB(const CsPFB<T>&)=delete;
        CsPFB<T>& operator=(const CsPFB<T>& rhs)=delete;



        template <std::ranges::range H>
        CsPFB(size_t nch1, const H& h0)
        :nch(nch1), batch((h0.end()-h0.begin())/nch1-1),
        fir(calc_coeffs(nch1, h0)), plan(init_plan(nch, batch))
#ifndef INPLACE_ALGO
        ,buffer(nch1, batch)
#endif
        {
        }

        size_t size_per_shoot()const{
            return nch*batch;
        }


        template <std::forward_iterator X>
        Array2D<std::complex<T>> analyze(X begin, X end){
            Array2D<std::complex<T>> x(batch, nch ,begin, end);
#ifndef INPLACE_ALGO
            x.transpose(buffer);
            //x.reshape(x.ncols(), x.rows());
            x.swap(buffer);
#else
            x.transpose_self();
#endif
            //x.swap(buffer);
            analyze_transposed(x);
            return x;
        }


        template <std::ranges::range X>
        Array2D<std::complex<T>> analyze(const X& x1){
            return analyze(x1.begin(), x1.end());
        }

        void analyze_transposed(Array2D<std::complex<T>>& x){
/*
 *elements of row-major 2D array x should be in order :
 0  8 16 24
 1  9 17 25
 2 10 18 26
 3 11 19 27
 4 12 20 28
 5 13 21 29
 6 14 22 30
 7 15 23 31
 */
            fir.filter(x);
            //x.transpose_self();
            //std::cout<<x<<std::endl;
#ifndef INPLACE_ALGO
            x.transpose(buffer);
            //x.reshape(x.ncols(), x.nrows());
            x.swap(buffer);
#else
            x.transpose_self();
#endif
            FftwTraits<T>::execute_dft(plan, x.data.get(), x.data.get());
        }
    };
}

#endif

