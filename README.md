# tensor-noise

Generating tensor noise with isotropic autocovariance matrix.
The description of the problem is available in

https://gapolinario.github.io/physics/2020/04/24/generating-random-gaussian-correlated-tensors.html

Create a folder called `data`. Compile C code with

```
gcc tensor-noise.c -O3 -lm -lgsl -lgslcblas -o tensor-noise.x
```
