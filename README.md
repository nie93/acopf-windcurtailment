# Wind Curtailment Remedial Action Scheme (ACOPF)

## Development Environment

* Windows
* Python 2.7, 3.5.2
* numpy (==1.16.1+mkl)
* scipy (==0.19.0, 1.2.1)
* pypmu (==0.1.1)

*Notes: Experienced data type problem in `numpy-1.12.0` when calling `numpy.imag()` functions. Converting sparse matrix to array consumes more time significantly in the computation.*

<<<<<<< HEAD
## BeagleBone Black (armv7h) Development

### Using `ipopt` and Python interface `cyipopt`

[Install Scripts Reference](https://github.com/matthias-k/cyipopt)

#### Package Requirements

```sh
$ sudo apt install cython
$ sudo apt-get -y install liblapack-dev libblas-dev
$ sudo apt-get -y install gfortran
$ sudo apt install coinor-libipopt1v5 coinor-libipopt-dev
$ sudo apt install python-numpy python-six python-future
```

#### Build and Install `cyipopt`

```sh
$ git clone https://github.com/matthias-k/cyipopt.git
$ cd cyipopt
$ python setup.py build
```

Check that everything linked correctly with `ldd` (depends on the CPU architecture):

```sh
$ ldd build/lib.linux-x86_64-2.7/cyipopt.so
```

```sh
$ ldd build/lib.linux-armv7l-2.7/cyipopt.so
```

Then, install

```sh
$ sudo python setup.py install
```
=======
## Other Branches

* [`dev-older-scipy`](https://github.com/nie93/acopf-windcurtailment/tree/dev-older-scipy)
* [`dev-ipopt`](https://github.com/nie93/acopf-windcurtailment/tree/dev-ipopt)
* [`dev-latest-scipy`](https://github.com/nie93/acopf-windcurtailment/tree/dev-latest-scipy)
>>>>>>> dev-older-scipy
