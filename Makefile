# Generate intermediate data files

.PHONY : sav
sav : bx_z249.sav by_z249.sav bz_z249.sav n_z249.sav te_z249.sav \
	bx_z302.sav bz_z302.sav n_z302.sav te_z302.sav \
	bx_z357.sav by_z357.sav bz_z357.sav n_z357.sav te_z357.sav \
	bx_z416.sav by_z416.sav bz_z416.sav n_z416.sav te_z416.sa

idl_bin=/Applications/exelis/idl85/bin/idl

bx_z249.sav:
  cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot3a','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bx_z249.sav

by_z249.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot3a','y','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav by_z249.sav

bz_z249.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot3a','z','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bz_z249.sav

n_z249.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'3p1','z','n',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav n_z249.sav

te_z249.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'3p1','z','te',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav te_z249.sav

bx_z302.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot3a','x','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav bx_z302.sav

#by_z302.sav:
#	cd ../output/intermediate/; \
#	${idl_bin} -e "pro00710,'bdot3a','y','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
#	mv 11*.sav by_z302.sav

bz_z302.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot3a','z','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav bz_z302.sav

n_z302.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'3p1','z','n',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav n_z302.sav

te_z302.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'3p1','z','te',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav te_z302.sav

bx_z357.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot10','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bx_z357.sav

by_z357.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot10','y','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav by_z357.sav

bz_z357.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot10','z','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bz_z357.sav

n_z357.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'3p2','z','n',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav n_z357.sav

te_z357.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'3p2','z','te',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav te_z357.sav

bx_z416.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot10','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bz_z416.sav

by_z416.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot10','y','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav by_z416.sav

bz_z416.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'bdot10','z','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav bz_z416.sav

n_z416.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'3p2','z','te',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav n_z416.sav

te_z416.sav:
	cd ../output/intermediate/; \
	${idl_bin} -e "pro00710,'3p2','z','te',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav te_z416.sav

.PHONY : clean

clean :
	cd ../output/intermediate/; \
	rm *.tar; \
  rm *.txt
