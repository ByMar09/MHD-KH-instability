OBJECTS =  parameters.o scalar.o threevectors.o fourvectors.o funciones.o imex_scheme.o  mp5.o  mp7.o  mp9.o \
	01_openfiles.o 01_open_backup.o 02_init_interface.o \
	03_init_var_primitive.o 03_restart_primitive.o\
	04_boundary.o 05_reconstruct.o 06_conserved_var.o 07_clean_var.o \
	08_swap1.o 09_sumint.o 10_varint.o 11_elecvarint.o \
	12_varprimitive.o 12_varprimitive_nr.o  13_intercon_var.o\
	14_flowinit.o 15_lax_flow.o 16_hll_flow_mcl.o 16_hll_flow_mpx.o \
	16_hllc_flow_mcl.o 16_hllc_flow_mpx.o 16_hlle_flow_mcl.o 16_hlle_flow_mpx.o \
	17_fluxandsource.o 18_sumfin.o 19_varfin.o \
	20_electricfield.o 21_swap2.o 22_writer.o 23_fin_interface.o  \
	24_boundary_conserved.o 24_boundary_flux_conserved.o \
	25_boundary_electric.o 26_source_boundary.o \
	main_cueva_code.o

MODULES = parameters.mod scalar.mod threevectors.mod fourvectors.mod \
	  funciones.mod

.PHONY: clean		

cueva.out: $(MODULES) $(OBJECTS)
	gfortran  -fopenmp -mcmodel=large -O5 -ftree-vectorize -funroll-all-loops -fno-unsafe-math-optimizations -frounding-math $(OBJECTS) -o cueva.out 

# export PATH=/opt/slurm/sbin:/opt/slurm/bin:$PATH
# export OMP_NUM_THREADS= 4
# -fbounds-check


%.mod: %.f95
	gfortran  -fopenmp -mcmodel=large -c  -O5 -ftree-vectorize -funroll-all-loops -fno-unsafe-math-optimizations -frounding-math $<

%.o: %.f95
	gfortran  -fopenmp -mcmodel=large -c  -O5 -ftree-vectorize -funroll-all-loops -fno-unsafe-math-optimizations -frounding-math $<


clean:
	rm -f $(OBJECTS) $(MODULES) cueva.out cueva.log

veryclean:
	rm -f $(OBJECTS) $(MODULES) cueva.out cueva.log *dat data_global/*dat data_local/*dat data_print/*dat data_restart/*dat
