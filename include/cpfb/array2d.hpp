#ifndef ARRAY2D_HPP
#define ARRAY2D_HPP
#include <memory>
#include <iostream>
#include <vector>
#include <type_traits>
#include <cassert>
#include <iterator>
#include <ranges>
#include <type_traits>
#include <functional>

namespace cpfb{
    template <typename T, bool owned=true>
    struct Array2D{
        size_t _nrows;
        size_t _ncols;
        std::unique_ptr<T[], std::function<void(T*)>> data;

        Array2D()
        :_nrows(0), _ncols(0), data(nullptr,  [](T* p){if(owned){delete p;}}){
            static_assert(!owned);
        }

        Array2D(size_t nrows1, size_t ncols1)
        :_nrows(nrows1), _ncols(ncols1), data(new T[nrows1*ncols1], [](T* p){if(owned){delete p;}}){
            static_assert(owned);
        }

        Array2D(size_t nrows1, size_t ncols1, T* p)
        :_nrows(nrows1), _ncols(ncols1), data(p, [](T* p){if(owned){delete p;}}){
            static_assert(!owned);
        }

        template <std::ranges::range U>
        Array2D(size_t nrows1, size_t ncols1, const U& rhs)
        :Array2D(nrows1, ncols1, rhs.begin(), rhs.end()){
            static_assert(owned);
        }

        template <std::forward_iterator U>
        Array2D(size_t nrows1, size_t ncols1, U begin, U end)
        :Array2D(nrows1, ncols1){
            std::copy(begin, end, data.get());
            assert(end-begin==nrows1*ncols1);
            static_assert(owned);
        }

        Array2D(const Array2D<T,owned>&)=delete;
        Array2D& operator= (const Array2D<T,owned>&) = delete;

        Array2D(Array2D<T,owned>&& rhs)
        :_nrows(rhs.nrows()), _ncols(rhs.ncols()), data(std::move(rhs.data))
        {
            rhs._nrows=0;
            rhs._ncols=0;
        }

        template <bool owned1>
        Array2D(Array2D<T, owned1>& rhs)
        :Array2D(rhs._nrows, rhs._ncols, rhs.data.get())
        {
        }


        Array2D& operator=(Array2D&& rhs){
            data=std::move(rhs.data);
            _nrows=rhs._nrows;
            rhs._nrows=0;
            _ncols=rhs._ncols;
            rhs._ncols=0;
            return *this;
        }

        Array2D& operator=(Array2D& rhs){
            static_assert(!owned);
            _nrows=rhs._nrows;
            _ncols=rhs._ncols;
            data.reset(rhs.data.get());
            return *this;
        }

        Array2D<T, true> clone()const{
            Array2D<T,true> result(_nrows, _ncols);
            std::copy(data.get(), data.get()+size(), result.data.get());
            return result;
        }

        void swap(Array2D<T,owned>& rhs){
            if(!owned){
                assert(0);
            }
            assert(rhs.size()==size());
            data.swap(rhs.data);
            std::swap(_nrows, rhs._nrows);
            std::swap(_ncols, rhs._ncols);
        }

        size_t nrows()const{
            return _nrows;
        }

        size_t ncols()const{
            return _ncols;
        }

        size_t size()const{
            return _nrows*_ncols;
        }

        T& get(size_t i, size_t j){
            return data[i*_ncols+j];
        }

        const T& get(size_t i, size_t j) const{
            return data[i*_ncols+j];
        }

        T& operator()(size_t i, size_t j){
            return data[i*_ncols+j];
        }

        const T& operator()(size_t i, size_t j) const{
            return data[i*_ncols+j];
        }

        Array2D& reverse_col_self(){
            for(size_t i=0;i<_nrows;++i){
                for(size_t j=0;j<_ncols/2;++j){
                    std::swap(get(i,j), get(i, _ncols-1-j));
                }
            }
            return *this;
        }

        Array2D& reverse_row_self(){
            for(size_t i=0;i<_nrows/2;++i){
                for(size_t j=0;j<_ncols;++j){
                    std::swap(get(_nrows-1-i,j), get(i, j));
                }
            }
            return *this;
        }


        Array2D<T, true> transpose()const{
            Array2D<T,true> result(ncols(), nrows());
            for(size_t i=0;i<_nrows;++i){
                for(size_t j=0;j<_ncols;++j){
                    result(j,i)=this->get(i,j);
                }
            }
            return result;
        }

        template <bool owned1>
        void transpose(Array2D<T,owned1>& result)const{
            assert(result.size()==size());
            result._ncols=_nrows;
            result._nrows=_ncols;
            for(size_t i=0;i<_nrows;++i){
                for(size_t j=0;j<_ncols;++j){
                    result.get(j,i)=get(i,j);
                }
            }
        }

        void reshape(size_t nrows1, size_t ncols1){
            assert(nrows1*ncols1==size());
            _nrows=nrows1;
            _ncols=ncols1;
        }

        void transpose_self()
        {
            auto first=data.get();
            auto last=first+_nrows*_ncols;
            const int mn1 = (last - first - 1);
            
            std::vector<bool> visited(_nrows*_ncols);
            auto cycle = first;
            while (++cycle != last) {
                if (visited[cycle - first])
                    continue;
                size_t a = cycle - first;
                do  {
                    a = a == mn1 ? mn1 : (_nrows * a) % mn1;
                    std::swap(*(first + a), *cycle);
                    visited[a] = true;
                } while ((first + a) != cycle);
            }
            std::swap(_nrows, _ncols);
        }

        template <typename F>
        Array2D& transform_self(F&& f){
            for(size_t i=0;i<size();++i){
                data[i]=f(data[i]);
            }
            return *this;
        }

        template <typename F>
        auto transform(F&& f)const ->Array2D<typename std::invoke_result<F, T>::type, true> {
            typedef decltype(f(get(0,0))) Tres;
            Array2D<Tres, true> result(_nrows, _ncols);
            for(int i=0;i<size();++i){
                result.data[i]=f(data[i]);
            }
            return result;
        }
    };

    template <typename T, bool owned>
    std::ostream& operator<<(std::ostream& os, const Array2D<T, owned>& x){
        os<<"Array ("<<x.nrows()<<" x "<<x.ncols()<<"): ["<<std::endl;
        for(size_t i=0;i<x.nrows();++i){
            for(size_t j=0;j<x.ncols();++j){
                os<<x(i,j)<<" ";
            }
            os<<std::endl;
        }
        os<<"]";
        return os;
    }
}

#endif
