# tight-binding
Density of states using tight-binding model, programmed in C with OpenMP parallel implementation.

Pending spin-orbit interaction.

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
$ make tight-binding

Calculate DoS
$ echo platino.inp | ./tight-binding -d

