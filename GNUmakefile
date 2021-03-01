include arch.mk

default: all

#all: clean_Destdir Common Wannier90 Gw Mqsgw_dmft Ucal Fullgw_dmft Analytical_continuation

ifdef USE_HDF5
all: com_script com_wannier90 com_comgw com_comlowh com_comdc com_comcoulomb com_comwann com_ctqmc
else
all: com_script com_wannier90 com_comgw com_comlowh com_comdc com_comcoulomb com_comwann com_ctqmc com_risb
endif


com_script:
	cd bin/flapwmbpt_input && $(MAKE) all && cd ../../
com_wannier90:  	
	cd wannier90_2.1 && $(MAKE) all && cd ../
com_comgw:
	cd gw && $(MAKE) && cd ../
com_comlowh:
	cd ComLowH && $(MAKE) && cd ../
com_comdc:
	cd ComDC && $(MAKE) && cd ../
com_comcoulomb:
	cd ComCoulomb && $(MAKE) && cd ../
com_comwann:
	cd ComWann && $(MAKE) && cd ../  
com_ctqmc:
	cd ComCTQMC && $(MAKE) && cp ./bin/EVALSIM ../bin && cp ./bin/CTQMC ../bin && cd ../ 
com_risb:
	cd ComRISB && $(MAKE) && cd ../


#clean: clean_Common clean_Wannier90 clean_Gw clean_Mqsgw_dmft clean_Ucal clean_Fullgw_dmft clean_Analytical_continuation clean_Destdir

ifdef USE_HDF5
clean: clean_script clean_wannier90 clean_comgw clean_comlowh clean_comdc clean_comcoulomb clean_comwann clean_ctqmc clean_Destdir
else
clean: clean_script clean_wannier90 clean_comgw clean_comlowh clean_comdc clean_comcoulomb clean_comwann clean_ctqmc clean_comrisb clean_Destdir
endif


clean_script:
	cd bin && cd flapwmbpt_input && $(MAKE) clean && cd ../ && cd ../
clean_ctqmc:
	cd ComCTQMC && $(MAKE) clean && rm ../bin/CTQMC && rm ../bin/EVALSIM && cd ../
clean_comlowh:
	cd ComLowH && $(MAKE) clean && cd ../
clean_comgw:
	cd gw && $(MAKE) clean && cd ../
clean_comdc:
	cd ComDC && $(MAKE) clean && cd ../
clean_comcoulomb:
	cd ComCoulomb && $(MAKE) clean && cd ../
clean_comwann:
	cd ComWann && $(MAKE) clean && cd ../
clean_Common:  	
	cd common && $(MAKE) clean && cd ../
clean_wannier90:  	
	cd wannier90_2.1 && $(MAKE) clean && cd ../
clean_Gw:  	
	cd gw && $(MAKE) clean && cd ../
clean_Ucal:  	
	cd ucal && $(MAKE) clean && cd ../
clean_Mqsgw_dmft:  	
	cd mqsgw_dmft && $(MAKE) clean && cd ../
clean_Fullgw_dmft:  	
	cd fullgw_dmft && $(MAKE) clean && cd ../
clean_Analytical_continuation:  	
	cd analytical_continuation && $(MAKE) clean && cd ../
clean_comrisb:
	cd ComRISB && $(MAKE) clean && cd ..
	cd bin && rm -rf CyGutz CyGutzB check_band_gaps.py complot_bands.py \
		comrisb.py create_wh_rl_init.py exe_spci_analysis exe_spci_j2_mott \
		exe_spci_j2_mott_analysis exe_spci_mott exe_spci_s2_mott \
		exe_spci_sjz_mott exe_spci_sz_mott gs_ed.py gs_idmrg.py gs_ml.py \
		gs_rspci_mott_onfly.py gs_syten.py gw_gutz gwannden.py gwannier.py \
		gwannier_plot_bands.py init_ga.py init_magnetism.py init_mott.py \
		init_subval_tsolver.py misc modify_ginit.py modify_gparam.py \
		modify_u_matrix_for_impurity.py plot_band_tf.py plot_dos_tf.py \
		post_process run_ga.py save_ldag stepin_wien_gutz.py switch_gparam.py \
		&& cd ..
clean_Destdir:
	rm ./bin/rspflapw.exe ./bin/ComLowH ./bin/ComDC ./bin/ComWann ./bin/ComCoulomb
