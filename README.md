DosExtraction
=============

This program allows to extract the Density of States (DoS) in polymer semiconductors, assessed by means of Capacitance-Voltage measurements on Metal-Insulator-Semiconductor capacitors.

Compile
=======

In order to generate the executable, first open the *CMakeLists.txt* file (in the top-level folder) and, if necessary, edit it to your needs.

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

Finally:

```
$ make
```

will build the *simulate\_dos* executable and the *dosextraction* shared library under the *bin/* and *lib/*
directories (or the ones specified in *CMakeLists.txt*) respectively.

If you wish to install:

- the executable, into */usr/local/bin/*;
- the shared library, into */usr/local/lib/*;
- the header files, into */usr/local/include/dosextraction/*;

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