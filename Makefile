#executable name
exec = main
main = main.o

#test if ifort is present
wres = $(shell which flang > /dev/null; echo $$?)
wres2 = $(shell which ifort > /dev/null; echo $$?)
ifeq "$(wres)" "0"
	fc = flang
	switchOPT = -O3 -march=znver2 -mavx -flto -funroll-loops -fremap-arrays
else ifeq "$(wres2)" "0"
	fc = ifort
	switchOPT = -O3 -ipo -ip -unroll -xHost -g
        switchOPT = -O3 -ipo -ipo-jobs10 -unroll -axAVX -xSSE4.1 -g -traceback
	switchDBG = -O0 -check all -warn all
	switchDBG += -fpe0 -u -traceback -warn nounused -g
	switchDBG += -init=snan,zero,arrays
	switchOMP = $(switchOPT) -qopenmp
	nowarn = -nowarn
else
	fc = gfortran
	switchOPT = -ffree-line-length-none -O3 -g -fallow-argument-mismatch
	switchDBG = -fbacktrace -g
	switchDBG += -ffpe-trap=zero,overflow,invalid
	switchDBG += -fbounds-check -ffree-line-length-none -O3
	nowarn = -w
endif


#NPROCS = $(shell sysctl hw.ncpu  | grep -o '[0-9]\+')
NPROCS = 4
MAKEFLAGS = -j$(NPROCS)

#default switch
switch = $(switchOPT)

#objects
objs = opkda2.o
objs += opkda1.o
objs += opkdmain.o
objs += kemimo_commons.o
objs += kemimo_sticking.o
objs += kemimo_reactionarray.o
objs += kemimo_gas_rates.o
objs += kemimo_dust_rates.o
objs += kemimo_rates.o
objs += kemimo_flux.o
objs += kemimo_ode.o
objs += kemimo.o

#lib = -I/usr/include/python2.7 -lpython2.7


#default target
all:	main

main:	$(objs) $(main)
	$(fc) $(objs) $(main) -o $(exec) $(switch) $(lib)
########################
## For parallel execution, explicit dependencies:
opkda2.o: opkda2.f
	$(fc) $(switch) $(nowarn) -c $< -o $@

objs2 = opkda2.o
opkda1.o: opkda1.f $(opjs2)
	$(fc) $(switch) $(nowarn) -c $< -o $@

objs2 += opkda1.o
opkdmain.o: opkdmain.f $(objs2)
	$(fc) $(switch) $(nowarn) -c $< -o $@

objs2 += opkdmain.o
kemimo_commons.o: kemimo_commons.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@

objs2 += kemimo_commons.o

kemimo_sticking.o: kemimo_sticking.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@

objs2 += kemimo_sticking.o
kemimo_reactionarray.o: kemimo_reactionarray.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@

objs2 += kemimo_reactionarray.o

kemimo_gas_rates.o: kemimo_gas_rates.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@
kemimo_dust_rates.o: kemimo_dust_rates.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@

objs2 += kemimo_gas_rates.o
objs2 += kemimo_dust_rates.o

kemimo_rates.o: kemimo_rates.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@

objs2 += kemimo_rates.o
kemimo_flux.o: kemimo_flux.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@

objs2 += kemimo_flux.o
kemimo_ode.o: kemimo_ode.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@

objs2 += kemimo_ode.o
kemimo.o: kemimo.f90 $(objs2)
	$(fc) $(switch) -c $< -o $@

main.o: main.f90 $(objs)
	$(fc) $(switch) -c $< -o $@

########################
#full debug target
debug: switch = $(switchDBG)
debug: all

#openmp target
omp: switch = $(switchOMP)
omp: all

#openmp target alias
openmp: omp

#shared library target
sharedlib: switch += -fPIC
sharedlib: $(objs)
	$(fc) $(objs) -o libkemimo.so $(switch) -shared $(lib)

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec) libkemimo.so

.PHONY: clean

#rule for f90
%.o:%.f90
	$(fc) $(switch) -c $^ -o $@

#rule for f
%.o:%.f
	$(fc) $(switch) $(nowarn) -c $^ -o $@

#rule for c
%.o:%.c
	$(cc) $(cswitch) $(lib) -c $^ -o $@
