PYTHON3=/usr/bin/python3

kpolylib: setup.py kpolymatroid.pyx
	${PYTHON3} setup.py build_ext --inplace

clean:
	rm -rf build
