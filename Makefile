targets=test benchmark validate_ospfb_ampl_response validate_cspfb_ampl_response coeffs chirp_ospfb chirp_cspfb

headers=include/cpfb/array2d.hpp include/cpfb/batch_fir.hpp include/cpfb/fft_wrapper.hpp include/cpfb/cspfb.hpp include/cpfb/windowed_fir.hpp include/cpfb/oscillator.hpp
all:$(targets)
CXX=g++

INC=-I include `pkg-config --cflags fftw3`
LIBS=`pkg-config --libs fftw3` `pkg-config --libs fftw3f`

test: test.cpp $(headers)
	$(CXX) -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

validate_cspfb_ampl_response: src/validate_cspfb_ampl_response.cpp $(headers)
	$(CXX) -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

validate_ospfb_ampl_response: src/validate_ospfb_ampl_response.cpp $(headers)
	$(CXX) -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

coeffs: src/coeffs.cpp $(headers)
	$(CXX) -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

chirp_ospfb: src/chirp_ospfb.cpp $(headers)
	$(CXX) -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

chirp_cspfb: src/chirp_cspfb.cpp $(headers)
	$(CXX) -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

benchmark: src/benchmark.cpp $(headers)
	$(CXX) -o $@ $< -O3 -g $(INC) $(LIBS) --std=c++20

clean:
	rm -f $(targets)
