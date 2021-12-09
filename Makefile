all:
	cd src; make
clean:
	cd src; make clean; 
	cd src/expand; make clean
distclean:
	rm calliope expand *.nc *.dat *.tmp restart/*; make clean
print:
	cd src; make print
