targets=test validate_cspfb_ampl_response benchmark

headers=include/cpfb/array2d.hpp include/cpfb/batch_fir.hpp include/cpfb/fft_wrapper.hpp include/cpfb/cspfb.hpp include/cpfb/windowed_fir.hpp include/cpfb/oscillator.hpp
all:$(targets)


INC=-I include `pkg-config --cflags fftw3`
LIBS=`pkg-config --libs fftw3` `pkg-config --libs fftw3f`

test: test.cpp $(headers)
	g++ -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

validate_cspfb_ampl_response: src/validate_cspfb_ampl_response.cpp $(headers)
	g++ -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

benchmark: src/benchmark.cpp $(headers)
	g++ -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

clean:
	rm -f $(targets)
