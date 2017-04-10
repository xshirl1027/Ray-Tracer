#!/bin/sh
g++ -O4 -g -O0 svdDynamic.c RayTracer.c utils.c -lm -o RayTracerA3

g++ -O4 -g -O0 svdDynamic.c RayTracer2.c utils.c -lm -o RayTracerA4

g++ -O4 -g -O0 svdDynamic.c RayTracerScene2.c utils.c -lm -o RayTracerScene
