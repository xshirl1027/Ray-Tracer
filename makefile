CC=g++
CXXFLAGS = -g -O2 -fopenmp -std=gnu++11

all: RayTracer RayTracer2

RayTracer: RayTracer.o 
	g++ $(CXXFLAGS) svdDynamic.c RayTracer.c utils.c -lm -o RayTracer

RayTracer2: RayTracer2.o 
	g++ $(CXXFLAGS) svdDynamic.c RayTracer2.c utils.c -lm -o RayTracer2
