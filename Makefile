targets=test validate_cspfb_ampl_response

all:$(targets)

INC=-I include `pkg-config --cflags fftw3`
LIBS=`pkg-config --libs fftw3` `pkg-config --libs fftw3f`

test: test.cpp include/array2d.hpp include/batch_fir.hpp include/fft_wrapper.hpp include/cspfb.hpp include/windowed_fir.hpp
	g++ -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

validate_cspfb_ampl_response: src/validate_cspfb_ampl_response.cpp include/array2d.hpp include/batch_fir.hpp include/fft_wrapper.hpp include/cspfb.hpp include/windowed_fir.hpp include/oscillator.hpp
	g++ -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20



clean:
	rm -f $(targets)
