kpolyext.py is a reimplementation in Python3 of the algorithm for enumerating
integer polymatroids found in:

Thomas J. Savitsky. Enumeration of 2-polymatroids on up to seven elements.
SIAM J. Discrete Math., 28(4) (2014), 1641--1650.

It requires the ipython package for python.

The following commands will generate a catalog of 2-polymatroids.

    $ echo -n '0,0' > n0.2pm.txt
    $ python3 kpolyext.py -c 2 n0.2pm.txt n1.2pm.txt 
    $ python3 kpolyext.py -c 2 n1.2pm.txt n2.2pm.txt
    $ python3 kpolyext.py -c 2 n2.2pm.txt n3.2pm.txt
    $ python3 kpolyext.py -c 2 n3.2pm.txt n4.2pm.txt
    $ python3 kpolyext.py -c 2 n4.2pm.txt n5.2pm.txt
    $ python3 kpolyext.py -c 2 n5.2pm.txt n6.2pm.txt
    $ python3 kpolyext.py -c 2 n6.2pm.txt n7.2pm.txt

Polymatroids are represented as space-separated lists of flats and
their ranks.  Flats are given as bitmasks in hexadecimal without
the leading '0x', and the rank is set off from the bitmask by a comma.
For example: 0,0 1,2 2,2 3,3 is the 2-polymatroid consisting of
two lines placed freely in the plane.

The following commands will generate a catalog of matroids.

    $ echo -n '0,0' > n0.matroids.txt
    $ python3 kpolyext.py -c 1 n0.matroids.txt n1.matroids.txt
    $ python3 kpolyext.py -c 1 n1.matroids.txt n2.matroids.txt
    $ python3 kpolyext.py -c 1 n2.matroids.txt n3.matroids.txt
    $ python3 kpolyext.py -c 1 n3.matroids.txt n4.matroids.txt
    $ python3 kpolyext.py -c 1 n4.matroids.txt n5.matroids.txt
    $ python3 kpolyext.py -c 1 n5.matroids.txt n6.matroids.txt
    $ python3 kpolyext.py -c 1 n6.matroids.txt n7.matroids.txt
    $ python3 kpolyext.py -c 1 n7.matroids.txt n8.matroids.txt
    $ python3 kpolyext.py -c 1 n8.matroids.txt n9.matroids.txt

GNU Parallel can be used to speed up the process.  For example,
the last command can be parallelized as follows:
    $ cat n8.matroids.txt | parallel -N1 --pipe --progress \
        'python3 kpolyext.py -c 1' > n9.matroids.txt

Two perl scripts are present.  One converts the bitstrings
representing a polymatroid to a more human-readable format.
The other tallies polymatroids by their ranks.

--Thomas J. Savitsky
savitsky@gwmail.gwu.edu
July 10, 2020
