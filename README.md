# tight-binding
Density of states using tight-binding model, programmed in C with OpenMP parallel implementation.

With spin-orbit interaction.

Dependencies
------------

* OpenMP
* GSL (GNU Scientific Library)
* BLAS
* LAPACK
* CBLAS

Download /clone
---------------
$ git clone https://github.com/graguirre/tight-binding.git

Make and run
------------
* Make

```
$ make tight-binding
```

* Calculate DoS 

```
$ cat input/pt-pt.xyz | ./tight-binding -d
```

* Calculate DoS with spin-orbit (parameter lambda=0.1) 

```
$ cat input/pt-pt.xyz | ./tight-binding -d -l 0.1
```

* Plot DOS, using gnuplot

```
gnuplot> plot '<(cat input/cadenaPt.xyz | ./tight-binding -d)' u 1:2 w l
```

* Plot and compare DOS

```
$ cat input/cadenaPt.xyz | ./tight-binding -g > dos-green.dat
gnuplot> plot '<(cat input/cadenaPt.xyz | ./tight-binding -d)' u 1:2 w l, 'dos-green.dat' u 1:2 w l, '<(octave --silent extra/dos.m)' u 1:2 w l
```


* Plot Hamiltonian matrix, using gnuplot

```
gnuplot> set view map
gnuplot> plot '<(cat input/pt-pt.xyz | ./tight-binding -h)' matrix with image
```
