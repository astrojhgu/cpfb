#ifndef FFT_WRAPPER
#define FFT_WRAPPER
#include <fftw3.h>
#include <memory>
#include <complex>
namespace cpfb{
    template <std::floating_point T>
    struct FftwTraits{
    };

    template <>
    struct FftwTraits<double>{
        using UniquePtrPlan=std::unique_ptr<std::remove_pointer<fftw_plan>::type, 
        void(*)(fftw_plan)>;

        static UniquePtrPlan plan_dft_1d(int n, int sign, unsigned int flags){
            return UniquePtrPlan(fftw_plan_dft_1d(n, nullptr, nullptr, sign, flags),
            fftw_destroy_plan
            );
        }

        static UniquePtrPlan plan_many_dft(int rank, const int* n, int howmany, 
                                            int* inembed, int istride, int idist,
                                            const int* onembed,
                                            int ostride, int odist, int sign, unsigned int flags)
        {
            return UniquePtrPlan(fftw_plan_many_dft(
                rank, n, howmany, nullptr, inembed, istride, idist, 
                nullptr, onembed, ostride, odist, sign, flags
            ),
            fftw_destroy_plan
            );

        }

        static void execute_dft(UniquePtrPlan& plan, std::complex<double>* input, std::complex<double>* output){
            fftw_execute_dft(plan.get(), (fftw_complex*)input, (fftw_complex*)output);
        }
        

    };

    template <>
    struct FftwTraits<float>{
        using UniquePtrPlan=std::unique_ptr<std::remove_pointer<fftwf_plan>::type, 
        void(*)(fftwf_plan)>;

        static UniquePtrPlan plan_dft_1d(int n, int sign, unsigned flags){
            return UniquePtrPlan(fftwf_plan_dft_1d(n, nullptr, nullptr, sign, flags),
            fftwf_destroy_plan
            );
        }

        static UniquePtrPlan plan_many_dft(int rank, const int* n, int howmany, 
                                            int* inembed, int istride, int idist,
                                            const int* onembed,
                                            int ostride, int odist, int sign, unsigned int flags)
        {
            return UniquePtrPlan(fftwf_plan_many_dft(
                rank, n, howmany, nullptr, inembed, istride, idist, 
                nullptr, onembed, ostride, odist, sign, flags
            ),
            fftwf_destroy_plan
            );
        }

        static void execute_dft(UniquePtrPlan& plan, std::complex<float>* input, std::complex<float>* output){
            fftwf_execute_dft(plan.get(), (fftwf_complex*)input, (fftwf_complex*)output);
        }
    };


}

#endif
