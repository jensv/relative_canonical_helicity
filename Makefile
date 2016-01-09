# Generate intermediate data files

.PHONY : vts
vts : b_z249_000000.vts t_z249_000000.vts n_z249_00000.vts \
	b_z302_000000.vts t_z302_000000.vts n_z302_000000.vts \
	b_z357_000000.vts t_z357_000000.vts n_z257_000000.vts \
	b_z416_000000.vts t_z416_000000.vts n_z416_000000.vts

idl_bin=/Applications/exelis/idl85/bin/idl
python_script=idl_plane_to_structured_vtk.py
time_points=21

b_z249_000000.vts : bx_z249.sav by_z249.sav bz_z249.sav
	ipython ${python_script} vector $^ 0.249 ${time_points} b_z249 B 1 1 None 25 25
t_z249_000000.vts : t_z249.sav
	ipython ${python_script} scalar $^ 0.249 ${time_points} t_z249 t 1 1 None 25 25
n_z249_000000.vts : n_z249.sav
	ipython ${python_script} scalar $^ 0.249 ${time_points} n_z249 n 1 1 None 25 25

b_z302_000000.vts : bx_z249.sav by_z249.sav bz_z249.sav
	ipython ${python_script} vector $^ 0.302 ${time_points} b_z302 B 1 1 None 25 25
t_z302_000000.vts : bx_z249.sav
	ipython ${python_script} scalar $^ 0.302 ${time_points} t_z302 t 1 1 None 25 25
n_z302_000000.vts : bx_z249.sav
	ipython ${python_script} scalar $^ 0.302 ${time_points} n_z302 n 1 1 None 25 25

b_z357_000000.vts : bx_z249.sav by_z249.sav bz_z249.sav
	ipython ${python_script} vector $^ 0.249 ${time_points} b_z249 B 1 1 None 25 25
b_z357_000000.vts : bx_z249.sav
	ipython ${python_script} scalar $^ 0.249 ${time_points} t_z249 t 1 1 None 25 25
b_z357_000000.vts : bx_z249.sav
	ipython ${python_script} scalar $^ 0.249 ${time_points} n_z249 n 1 1 None 25 25

b_z416_000000.vts : bx_z249.sav by_z249.sav bz_z249.sav
	ipython ${python_script} vector $^ 0.249 ${time_points} b_z249 B 1 1 None 25 25
b_z416_000000.vts : bx_z249.sav by_z249.sav bz_z249.sav
	ipython ${python_script} scalar $^ 0.249 ${time_points} t_z249 t 1 1 None 25 25
b_z416_000000.vts : bx_z249.sav by_z249.sav bz_z249.sav
	ipython ${python_script} scalar $^ 0.249 ${time_points} n_z249 n 1 1 None 25 25

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
