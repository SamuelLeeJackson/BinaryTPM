thermal: thermal.cpp thermal_GPU.cu
	nvcc thermal.cpp thermal_GPU.cu -o thermal -Xcompiler -fopenmp -O3

visibility: visibility.cpp visibility_GPU.cu
	nvcc visibility.cpp visibility_GPU.cu -o visibility -Xcompiler -fopenmp -O3

observe: wise_observe.cpp
	g++ wise_observe.cpp -o observe -std=c++1y -O3 -fopenmp

fitting: fitting.cpp
	g++ fitting.cpp -o fitting -std=c++1y -O3 -fopenmp