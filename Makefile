# Generate intermediate data files

.PHONY : vts
vts : b_z249_000000.vts p_z249_000000.vts \
	p_z302_000000.vts \
	b_z357_000000.vts p_z357_000000.vts \
	b_z416_000000.vts p_z416_000000.vts

idl_bin=/Applications/exelis/idl85/bin/idl
python_script=../../writing_to_vtk/source/idl_plane_to_structured_vtk.py
time_points=21
kx=1
ky=1
x_points=30
y_points=30
smooth_factor=None

b_z249_000000.vts : bx_z249.sav by_z249.sav bz_z249.sav
	ipython ${python_script} vector $^ 0.249 ${time_points} b_z249 B ${kx} ${ky} ${smooth_factor} ${x_points} ${y_points}

p_z249_000000.vts : n_z249.sav te_z249.sav
	ipython ${python_script} scalar $^ 0.249 ${time_points} p_z249 n te ${kx} ${ky} ${smooth_factor} ${x_points} ${y_points}

#b_z302_000000.vts : bx_z302.sav by_z302.sav bz_z302.sav
#	ipython ${python_script} vector $^ 0.302 ${time_points} b_z302 B ${kx} ${ky} None ${x_points} ${y_points}

p_z302_000000.vts : n_z302.sav te_z302.sav
	ipython ${python_script} scalar $^ 0.302 ${time_points} p_z302 n te ${kx} ${ky} ${smooth_factor} ${x_points} ${y_points}

b_z357_000000.vts : bx_z357.sav by_z357.sav bz_z357.sav
	ipython ${python_script} vector $^ 0.357 ${time_points} b_z357 B ${kx} ${ky} ${smooth_factor} ${x_points} ${y_points}

p_z357_000000.vts : n_z357.sav te_z357.sav
	ipython ${python_script} scalar $^ 0.357 ${time_points} p_z357 n te ${kx} ${ky} ${smooth_factor} ${x_points} ${y_points}

b_z416_000000.vts : bx_z416.sav by_z416.sav bz_z416.sav
	ipython ${python_script} vector $^ 0.416 ${time_points} b_z416 B ${kx} ${ky} ${smooth_factor} ${x_points} ${y_points}

p_z416_000000.vts : n_z416.sav te_z416.sav
	ipython ${python_script} scalar $^ 0.416 ${time_points} p_z416 n te ${kx} ${ky} ${smooth_factor} ${x_points} ${y_points}

bx_z249.sav:
	${idl_bin} -e "pro00710,'bdot3a','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bx_z249.sav

by_z249.sav:
	${idl_bin} -e "pro00710,'bdot3a','y','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav by_z249.sav

bz_z249.sav:
	${idl_bin} -e "pro00710,'bdot3a','z','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bz_z249.sav

n_z249.sav:
	${idl_bin} -e "pro00710,'3p1','z','n',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav n_z249.sav

te_z249.sav:
	${idl_bin} -e "pro00710,'3p1','z','te',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav te_z249.sav

bx_z302.sav:
	${idl_bin} -e "pro00710,'bdot3a','x','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav bx_z302.sav

#by_z302.sav:
#	${idl_bin} -e "pro00710,'bdot3a','y','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
#	mv 11*.sav by_z302.sav

bz_z302.sav:
	${idl_bin} -e "pro00710,'bdot3a','z','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav bz_z302.sav

n_z302.sav:
	${idl_bin} -e "pro00710,'3p1','z','n',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav n_z302.sav

te_z302.sav:
	${idl_bin} -e "pro00710,'3p1','z','te',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav te_z302.sav

bx_z357.sav:
	${idl_bin} -e "pro00710,'bdot10','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bx_z357.sav

by_z357.sav:
	${idl_bin} -e "pro00710,'bdot10','y','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav by_z357.sav

bz_z357.sav:
	${idl_bin} -e "pro00710,'bdot10','z','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bz_z357.sav

n_z357.sav:
	${idl_bin} -e "pro00710,'3p2','z','n',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav n_z357.sav

te_z357.sav:
	${idl_bin} -e "pro00710,'3p2','z','te',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav te_z357.sav

bx_z416.sav:
	${idl_bin} -e "pro00710,'bdot10','x','b',indgen(21)*1./20.,1,shotset='001',current_rise=0"; \
	mv 11*.sav bx_z416.sav

by_z416.sav:
	${idl_bin} -e "pro00710,'bdot10','y','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav by_z416.sav

bz_z416.sav:
	${idl_bin} -e "pro00710,'bdot10','z','b',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav bz_z416.sav

n_z416.sav:
	${idl_bin} -e "pro00710,'3p2','z','n',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav n_z416.sav

te_z416.sav:
	${idl_bin} -e "pro00710,'3p2','z','te',indgen(21)*1./20.,1,shotset='003',current_rise=0"; \
	mv 11*.sav te_z416.sav

.PHONY : clean

clean :
	cd ../output/intermediate/; \
	rm *.tar; \
  rm *.txt
