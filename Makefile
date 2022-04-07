all:
	cd src; make
clean:
	cd src; make clean; 
	cd src/expand; make clean
expand:
	cd src/expand; make expand
distclean:
	rm calliope expand *.nc *.dat *.tmp restart/*; make clean
print:
	cd src; make print
