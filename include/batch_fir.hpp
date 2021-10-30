#ifndef BATCH_FIR_HPP
#define BATCH_FIR_HPP
#include "array2d.hpp"

namespace cpfb{
    template <typename C, typename D>
    struct BatchFIR{
        Array2D<C> coeff_rev;
        Array2D<D> state;

        BatchFIR(Array2D<C>&& c)
        :coeff_rev(std::move(c)), state(coeff_rev.nrows(), coeff_rev.ncols()-1)
        {
            coeff_rev.reverse_col_self();
            state.transform_self([](auto){return D();});
        }


        BatchFIR(BatchFIR&&)=default;
        BatchFIR(const BatchFIR&)=delete;
        BatchFIR& operator=(const BatchFIR&)=delete;
        BatchFIR& operator=(BatchFIR&&)=default;

        void filter(Array2D<D>& x){
            assert(x.nrows()==coeff_rev.nrows());
            assert(x.ncols()==state.ncols());
            auto coeff_rev_nrows=coeff_rev.nrows();
            auto state_ncols=state.ncols();
            auto x_ncols=x.ncols();
            auto coeff_rev_ncols=coeff_rev.ncols();
            for(size_t i=0;i<coeff_rev_nrows;++i){
                auto d1=&state(i,0);
                auto d2=&x(i,0);
                auto c1=&coeff_rev(i,0);
                for(size_t j=0;j<state_ncols;++j){
                    size_t k=0;
                    D x1=D();
                    for(;k+j<state_ncols;++k){
                        x1+=c1[k]*d1[k+j];
                    }
                    for(;k<coeff_rev_ncols;++k){
                        assert(k+j-state_ncols<x_ncols);
                        x1+=c1[k]*d2[k+j-state_ncols];
                    }
                    d1[j]=x1;
                }
            }
            state.swap(x);
        }
    };

    template <typename C, typename D>
    std::ostream& operator<<(std::ostream& os, const BatchFIR<C, D>& bf){
        os<<"BatchFIR {"<<std::endl;
        os<<"coeffs reversed:"<<std::endl;
        os<<bf.coeff_rev<<std::endl;
        os<<"state:"<<std::endl;
        os<<bf.state<<std::endl;
        os<<"}";
        return os;
    }
}

#endif
