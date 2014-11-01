DosExtraction
=============

This program allows to extract the Density of States (DoS) in polymer semiconductors, assessed by means of Capacitance-Voltage measurements on Metal-Insulator-Semiconductor capacitors.

Compile
=======

In order to generate the executables, first open the *CMakeLists.txt* file (in the top-level folder) and, if necessary, edit it to your needs.

Then create a build directory and move into it:

```
$ mkdir build
$ cd build
```

Now you're ready to configure your system:

```
$ cmake ..
```

or, should you want the compiler to produce also debug symbols:

```
$ cmake -DCMAKE_BUILD_TYPE=Debug ..
```

The further option:

```
$ -DCMAKE_INSTALL_PREFIX=my_dir
```

will allow you to customize the installation directory, which is by default */usr/local*.

Finally:

```
$ make
```

will build the *simulate\_dos* and *fit\_dos* executables and the *dosextraction* shared library under the *bin/* and *lib/*
directories (or the ones specified in *CMakeLists.txt*) respectively.

Install
=======

If you wish to install:

- the executable, into *${CMAKE_INSTALL_PREFIX}/bin/*;
- the shared library, into *${CMAKE_INSTALL_PREFIX}/lib/*;
- the header files, into *${CMAKE_INSTALL_PREFIX}/include/dosextraction/*;

just type:

```
# make install
```

while:

```
# make uninstall
```

will remove them.

Documentation
=============

If [Doxygen](http://www.doxygen.org) (version 3.8.6 or above) and [GraphViz](http://www.graphviz.org)
are found, the following command will generate the documentation under the *doc/* folder (or the one specified
in *CMakeLists.txt*):

```
$ make doc
```

Please **read it** for further information.
