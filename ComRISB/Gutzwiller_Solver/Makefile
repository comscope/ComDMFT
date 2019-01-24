SRC_lib = external/primme 
SRC_exe = \
		SRC_cmplx/CYGUTZ  \
		SRC_cmplx/CYGUTZ_BANDS \
		SRC_cmplx/GSPCI_ANALYSIS  \
		SRC_cmplx/GSPCI_S2_MOTT  \
		SRC_cmplx/GSPCI_MOTT  \
		SRC_cmplx/GSPCI_SZ_MOTT  \
		SRC_cmplx/GSPCI_JZ_MOTT  \
		SRC_cmplx/GSPCI_J2JPJN_MOTT  \
		SRC_cmplx/GSPCI_J2_MOTT_ANALYSIS  

all:
	mkdir -p lib include
	for dir in ${SRC_lib} ${SRC_exe} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir && ${MAKE} ; ) ; done

clean:
	for dir in ${SRC_lib} ${SRC_exe} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir && ${MAKE} clean ; ) ; done

install :
	for dir in ${SRC_exe} - ; do \
    if [ "$$dir" = "-" ] ; then break ; fi ; \
    (cd $$dir && ${MAKE} install ; ) ; done

