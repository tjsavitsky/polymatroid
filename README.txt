kpolyext.py is a reimplementation in Python of the algorithm for enumerating
integer polymatroids found in:

Thomas J. Savitsky. Enumeration of 2-polymatroids on up to seven elements.
SIAM J. Discrete Math., 28(4) (2014), 1641--1650.

It requires the python-igraph package for python.

The following commands will generate a catalog of matroids.

    $ echo '0,0' > n0.matroids.txt
    $ python3 kpolyext.py -c 1 n0.matroids.txt n1.matroids.txt
    $ python3 kpolyext.py -c 1 n1.matroids.txt n2.matroids.txt
    $ python3 kpolyext.py -c 1 n2.matroids.txt n3.matroids.txt
    $ python3 kpolyext.py -c 1 n3.matroids.txt n4.matroids.txt
    $ python3 kpolyext.py -c 1 n4.matroids.txt n5.matroids.txt
    $ python3 kpolyext.py -c 1 n5.matroids.txt n6.matroids.txt
    $ python3 kpolyext.py -c 1 n6.matroids.txt n7.matroids.txt
    $ python3 kpolyext.py -c 1 n7.matroids.txt n8.matroids.txt
    $ python3 kpolyext.py -c 1 n8.matroids.txt n9.matroids.txt

The following commands will generate a catalog of 2-polymatroids.

    $ echo '0,0' > n0.2pm.txt
    $ python3 kpolyext.py -c 2 n0.2pm.txt n1.2pm.txt 
    $ python3 kpolyext.py -c 2 n1.2pm.txt n2.2pm.txt
    $ python3 kpolyext.py -c 2 n2.2pm.txt n3.2pm.txt
    $ python3 kpolyext.py -c 2 n3.2pm.txt n4.2pm.txt
    $ python3 kpolyext.py -c 2 n4.2pm.txt n5.2pm.txt
    $ python3 kpolyext.py -c 2 n5.2pm.txt n6.2pm.txt
    $ python3 kpolyext.py -c 2 n6.2pm.txt n7.2pm.txt

GNU Parallel can be used to speed up the process.  For example,
the last command can be parallelized as follows:
    $ cat n6.2pm.txt | parallel --pipe --progress -N1 \
        'python3 kpolyext.py -c 2' > n7.2pm.txt

The last command may take several days to execute on a typical
desktop computer.  One may wish to use the shuf command instead
of cat when generating catalogs.

Small catalogs of k-polymatroids for k>2 can also be generated.

Polymatroids are represented as space-separated lists of flats and
their ranks.  Flats are given as bitmasks in hexadecimal without
the leading '0x', and the rank is set off from the bitmask by a comma.
For example: 0,0 1,2 2,2 3,3 is the 2-polymatroid consisting of
two lines placed freely in the plane.

This python implementation is significantly slower than the C version
described in the paper.  Hopefully, its improved readability makes
up for its lack of speed.

Two perl scripts are included.  One converts the strings
representing a polymatroid to a more human-readable format.
The other tallies polymatroids by their ranks.

A cython version of kpolyext.py is available in the cython directory.

--Thomas J. Savitsky
savitsky@gwmail.gwu.edu
July 14, 2020
