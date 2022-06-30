require(myRtools)
require(myRplots)
require(Rflexpart)
require(ncdf4)


#	SETTINGS

icon.dir = "/project/ivme/MCH-1/icon-art-BRM/icon_scripts_output_code/icon_art_output/150418/structured/"
icon.fn = "ICON-ART-OEM_1504CDO_wsnow0_DOM01_00140000.nc"		#	file name gives time after release
icon.var = "testtr12"

flex.dir = "/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_nnb_interpolation/150418"
flex.fn = "grid_conc_20180415110000.nc"
tidx =3	#	use tidx to get output for specific time 1: 1 hour after start, etc

xlim = c(7.,9.)
ylim = c(46.75, 47.75)
zlim = c(1E1, 1E4)
mu.tracer = 28.97	#	in the xml file for the icon run it says 6.4E-2. I assume this means kg/mol. g/mol would not make sense. No such molecule exists. 


levels= c(5,10,15,20)		#	level to be plotted, counted from the bottom


#	constants
mu.air = 28.97
R.air = 287

#	END OF SETTINGS


#	data is in 'conc' as a 5D array. Dimensions 4 and 5 only contain one entry. 
flex = read.grid.output(file.path(flex.dir, flex.fn), release="BEROMUENSTER 12", tidx=tidx)
#	drop not used dimension of conc
flex$conc = array(flex$conc, dim=dim(flex$conc)[1:3])

#	release location 
release = list(x=flex$hdr$rel.lng1, y=flex$hdr$rel.lat1, pch=4)

#	default FLEXPART plot
plot(flex, xlim=xlim, ylim=ylim, map.db="worldHires", lev=1)
for(level in levels){
#	plot using the same function as for ICON data; flexpart levels start with 1 at surface!
fill.2dplot(flex, para="conc", xlim=xlim, ylim=ylim, zlim=zlim, lev=level, log.scale=TRUE, 
	borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE, 
	nlevels=30, pnts=release, key.title="C (ng/m3)")

#	load icon output variable 'icon.var' as a 'reg.grid' object. 'reg.grid' is a class defined in 
#	'myRtools'. The 3D field is stored as icon$zz.  
icon = read.reg.grid.from.ncdf(file.path(icon.dir, icon.fn), sn=icon.var)
#	flip vertical order
icon$zz = icon$zz[,,dim(icon$zz)[3]:1]

#	convert units. ICON units are actually mole fractions (volume mixing ratios). That's not the same as in FLEXPART (mass mixing ratio). Conversion requires molar weight of the tracer in ICON and that of air. 
icon$zz = icon$zz * mu.tracer/mu.air * 1E12 # mol/mol -> g/g -> ng/kg


#	default plot function for class 'reg.grid' translates to fill.2dplot. By default field 'zz' is plotted. If 3D array then level 1 is displayed. Since there is not flip of levels in 'read.reg.grid.from.ncdf', we need to specify level=80 to plot surface level.
fill.2dplot(icon, xlim=xlim, ylim=ylim, zlim=zlim, lev=level, log.scale=TRUE, borders=TRUE, 
	type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE, nlevels=30,key.title="C (ng/m3)")


#	single plot for both models with one color scale
png(paste0("/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_nnb_interpolation/150418/images/cdo_nnb_150418_14h_level_",level,".png"), height=12, width=24, units="cm", res=800)
layout(matrix(1:3, ncol=3), widths=c(1,1,lcm(2)))
par(mar=c(4,4,2,1)+.1)
# my_rmse = rmse(icon$zz, flex$conc, na.rm = FALSE)
col.scale= fill.2dplot(flex, para="conc", xlim=xlim, ylim=ylim, zlim=zlim,lev= level, log.scale=TRUE,
    borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE, 
    nlevels=30, pnts=release, plot.key=FALSE, main="FLEXPART")
col.scale= fill.2dplot(icon, xlim=xlim, ylim=ylim, zlim=zlim, lev=level, log.scale=TRUE,
    borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE, 
    nlevels=30, pnts=release, plot.key=FALSE, main="ICON")
plot.colorpalette(col.scale)
dev.off()
 }
#	We should try to compare the total mass in both simulations. 
###################################################################
#	1) Get air density from ICON output
icon.rho = read.reg.grid.from.ncdf(file.path(icon.dir, icon.fn), sn="rho")
icon.rho$zz = icon.rho$zz[,,dim(icon.rho$zz)[3]:1]

#	2) crop ICON and FLEXPART to same grid
icon.rho = crop.grid(icon.rho, xrng=xlim, yrng=ylim)
icon = crop.grid(icon, xrng=xlim, yrng=ylim)
flex = crop.grid(flex, para="conc", xrng=xlim, yrng=ylim)
#	3) get level heights 
icon.hhl = read.reg.grid.from.ncdf(file.path(icon.dir, icon.fn), sn="z_ifc")
#	flip levels
icon.hhl$zz = icon.hhl$zz[,,dim(icon.hhl$zz)[3]:1]
icon.hhl = crop.grid(icon.hhl, xrng=xlim, yrng=ylim)
#	level thickness
dz = aperm(apply(icon.hhl$zz, 1:2, diff), c(2,3,1))

#	surface are of grid
A = get.surface.area(grid=icon, units="m^2")
#	volume of grid cells
V = rep(A$zz, dim(dz)[3]) * dz
#	total air mass of grid cells
mass.air = V * icon.rho$zz

#	total mass 
mass.icon = sum(mass.air * icon$zz) * 1E-12	# ng -> kg
mass.flex = sum(mass.air[,,1:dim(flex$conc)[3]] * flex$conc) * 1E-12

cat("Total mass ICON:", mass.icon, "(kg)\n")
cat("Total mass FLEXPART:", mass.flex, "(kg)\n")


