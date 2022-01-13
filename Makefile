FC=gfortran
CFLAGS=-c -g

parallel_tempering: main.f90
	$(FC) main.f90 -o parallel_tempering

run:
	make
	./parallel_tempering

clean:
	rm parallel_tempering

