require(myRtools)
require(myRplots)
require(Rflexpart)
require(ncdf4)


#	SETTINGS

flex.dir = "/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_nnb_interpolation/20180525"
flex2.dir = "/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_rbf_interpolation/20180525"
flex3.dir = "/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_byc_interpolation/20180525"
flex.fn = "grid_conc_20180525110000.nc"
tidx =6	#	use tidx to get output for specific time 1: 1 hour after start, etc

xlim = c(7.75,11)
ylim = c(47., 48.75)
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
flex2 = read.grid.output(file.path(flex2.dir, flex.fn), release="BEROMUENSTER 12", tidx=tidx)
#       drop not used dimension of conc
flex2$conc = array(flex2$conc, dim=dim(flex2$conc)[1:3])
flex3 = read.grid.output(file.path(flex3.dir, flex.fn), release="BEROMUENSTER 12", tidx=tidx)
#       drop not used dimension of conc
flex3$conc = array(flex3$conc, dim=dim(flex3$conc)[1:3])



#	default FLEXPART plot
plot(flex, xlim=xlim, ylim=ylim, map.db="worldHires", lev=1)
for(level in levels){
#	plot using the same function as for ICON data; flexpart levels start with 1 at surface!
fill.2dplot(flex, para="conc", xlim=xlim, ylim=ylim, zlim=zlim, lev=level, log.scale=TRUE, 
	borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE, 
	nlevels=30, pnts=release, key.title="C (ng/m3)")
fill.2dplot(flex2, para="conc", xlim=xlim, ylim=ylim, zlim=zlim, lev=level, log.scale=TRUE,
        borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE,
        nlevels=30, pnts=release, key.title="C (ng/m3)")
#fill.2dplot(flex3, para="conc", xlim=xlim, ylim=ylim, zlim=zlim, lev=level, log.scale=TRUE,
#        borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE,
#        nlevels=30, pnts=release, key.title="C (ng/m3)")

#	load icon output variable 'icon.var' as a 'reg.grid' object. 'reg.grid' is a class defined in 
#	'myRtools'. The 3D field is stored as icon$zz.  

#	single plot for both models with one color scale
png(paste0("/project/ivme/MCH-1/icon-art-BRM-CDO/3interpolations_comparison_25058_13h_level_",level,".png"), height=12, width=24, units="cm", res=800)
layout(matrix(1:3, ncol=3), widths=c(1,1,lcm(2)))
par(mar=c(4,4,2,1)+.1)
# my_rmse = rmse(icon$zz, flex$conc, na.rm = FALSE)
col.scale= fill.2dplot(flex, para="conc", xlim=xlim, ylim=ylim, zlim=zlim,lev= level, log.scale=TRUE,
    borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE, 
    nlevels=30, pnts=release, plot.key=FALSE, main="NNB")
col.scale= fill.2dplot(flex2, xlim=xlim, ylim=ylim, zlim=zlim, lev=level, log.scale=TRUE,
    borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE, 
    nlevels=30, pnts=release, plot.key=FALSE, main="RBF")
#col.scale= fill.2dplot(flex3, xlim=xlim, ylim=ylim, zlim=zlim, lev=level, log.scale=TRUE,
#    borders=TRUE, type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE,
#    nlevels=30, pnts=release, plot.key=FALSE, main="BYC")
plot.colorpalette(col.scale)
dev.off()
 }

