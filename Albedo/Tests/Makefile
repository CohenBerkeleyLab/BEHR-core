FC = gfortran
MEX = /Applications/MATLAB_R2014b.app/bin/mexsh
OBJFILES = ../user_brdf_kernels.o brdf_kernel_test.o

mex: mex_brdf_kernels.F user_brdf_kernels.f
	$(MEX) $^

mexdry: mex_brdf_kernels.F user_brdf_kernels.f
	$(MEX) -n $^

brdf_test: $(OBJFILES)
	$(FC) -o brdf_kernel_test.exe $^

%.o: %.f
	$(FC) -c -o $@ $<

.PHONY: clean mex
clean:
	rm -rf *.exe mex_brdf_kernels.mexmaci64 $(OBJFILES)
