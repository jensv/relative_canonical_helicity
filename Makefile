# Generate intermediate data files

.PHONY : vts
vts : b_z249.vts n_z249.vts te_z249.vts \
	b_z302.vts n_z302.vts te_z302.vts \
	b_z357.vts n_z357.vts te_z357.vts \
	b_z416.vts n_z416.vts te_z416.vts

b_z249.vts: bx_z249.sav by_z249.sav bz_z249.sav


b_z302.vts:

b_z357.vts:

b_z416.vts:

n_z249.vts:

n_z302.vts:

n_z357.vts:

n_z416.vts:

te_z249.vts:

te_z302.vts:

te_z357.vts:

te_z416.vts:


intermediate : 
	cd ../output/intermediate/

	idl -e "pro00710,'bdot3a','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* bx_z249.sav
	idl -e "pro00710,'bdot3a','y','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* by_z249.sav
	idl -e "pro00710,'bdot3a','z','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* bz_z249.sav
	idl -e "pro00710,'3p1','z','n',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* n_z249.sav
	idl -e "pro00710,'3p1','z','te',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* te_z249.sav


	idl -e "pro00710,'bdot3a','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* bx_z302.sav	
	idl -e "pro00710,'bdot3a','y','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* by_z302.sav
	idl -e "pro00710,'bdot3a','z','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* bz_z302.sav
	idl -e "pro00710,'3p1','z','n',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* n_z302.sav
	idl -e "pro00710,'3p1','z','te',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* te_z302.sav

	idl -e "pro00710,'bdot10','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* bx_z357.sav	
	idl -e "pro00710,'bdot10','y','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* by_z357.sav
	idl -e "pro00710,'bdot10','z','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* bz_z357.sav
	idl -e "pro00710,'3p2','z','n',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* n_z357.sav
	idl -e "pro00710,'3p2','z','te',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* te_z357.sav


	idl -e "pro00710,'bdot10','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* bx_z416.sav
	idl -e "pro00710,'bdot10','y','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* by_z416.sav
	idl -e "pro00710,'bdot10','z','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* bz_z416.sav
	idl -e "pro00710,'3p2','z','n',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* n_z416.sav
	idl -e "pro00710,'3p2','z','te',indgen(21)*1./20.,1,shotset='001',current_rise=0"
	mv 11* te_z416.sav

