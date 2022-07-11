require(myRtools)
require(myRplots)
require(Rflexpart)
require(ncdf4)


#####################################################################
#	SETTINGS
#####################################################################
case = "20181221" # "150418"
interpolation = "nnb" #"rbf" #"byc" #"rbf" # "nnb"

icon.dir = file.path("/project/ivme/MCH-1/icon-art-BRM/icon_scripts_output_code/icon_art_output", 
	case,"structured") # "structured")

icon.pattern = "ICON-ART-OEM_211218CDO_DOM01_*"

icon.var = "testtr212"

flex.dir = file.path("/project/ivme/MCH-1/icon-art-BRM-CDO", 
	paste0("fp_output_", interpolation, "_interpolation"), case)

flex.pattern = "grid_conc_20181221110000.nc"

tidx = 3	#	use tidx to get output for specific time 1: 1 hour after start, etc
rel.com = "BEROMUENSTER 212"

out.dir ="/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_nnb_interpolation/20181221/images/plume_compare"

zlim = c(1, 1E5)
mu.tracer = 28.97	#	in the xml file for the icon run it says 6.4E-2. I assume this means kg/mol. g/mol would not make sense. No such molecule exists. 

levels = c(5,15)		#	level to be plotted, counted from t:he bottom

#	constants
mu.air = 28.97
R.air = 287

#####################################################################
#	END OF SETTINGS
#####################################################################


options(error=recover)

if (!dir.exists(out.dir)) dir.create(out.dir)
if (!dir.exists(file.path(out.dir, interpolation))) dir.create(file.path(out.dir, interpolation))

flex.hdr = read.header(file.path(flex.dir, "header"))
outlev.mids = flex.hdr$outheight - diff(c(0, flex.hdr$outheight))/2 
xlim = c(-1, 5) *   1 + flex.hdr$rel.lng1[1]
ylim = c(-1, 5) * 0.5 + flex.hdr$rel.lat1[1]
flex.out.dtm = flex.hdr$dtm + 
	seq(flex.hdr$loutstep, 23400, by=flex.hdr$loutstep)/86400
flex.out.dtm = round(flex.out.dtm, "mins")
nn.t = length(flex.out.dtm)

flex.fn = dir(flex.dir, flex.pattern, full.names=TRUE)

icon.fls = dir(icon.dir, icon.pattern)
case.dtm = chron(case, "00:00:00", format=c("ymd", "h:m:s"), out.format=c("y-m-d", "h:m:s"))

tmp = sapply(strsplit(icon.fls, "[_.]"), FUN=function(x) x[length(x)-1])
ss = as.numeric(substr(tmp, 7,8))
mm =  as.numeric(substr(tmp, 5,6))
hh =  as.numeric(substr(tmp, 3,4)) 
dd = as.numeric(substr(tmp, 1,2))
icon.dtm = round(case.dtm + dd + hh/24 + mm/24/60 + ss/86400, "min")

#	release location 
release = list(x=flex.hdr$rel.lng1, y=flex.hdr$rel.lat1, pch=4, cex=2, col="black")

plume.stats = data.frame(
	dtm = flex.out.dtm,
	rmse = rep(NA, nn.t),
	r = rep(NA, nn.t),
	dist.hor =  rep(NA, nn.t),
	dist.vert =  rep(NA, nn.t),
	flex.max.conc = rep(NA, nn.t),
	icon.max.conc = rep(NA, nn.t),
	flex.cm.lon = rep(NA, nn.t),
	flex.cm.lat = rep(NA, nn.t),
	flex.cm.height = rep(NA, nn.t),
	icon.cm.lon = rep(NA, nn.t),
	icon.cm.lat = rep(NA, nn.t),
	icon.cm.height = rep(NA, nn.t),
	flex.sd.lon = rep(NA, nn.t),
	flex.sd.lat = rep(NA, nn.t),
	flex.sd.height = rep(NA, nn.t),
	icon.sd.lon = rep(NA, nn.t),
	icon.sd.lat = rep(NA, nn.t),
	icon.sd.height = rep(NA, nn.t)
	)


for (ii in 1:nn.t){

	#	data is in 'conc' as a 5D array. Dimensions 4 and 5 only contain one entry. 
	flex = read.grid.output(flex.fn, release=rel.com, tidx=ii)
	#	drop not used dimension of conc
	flex$zz = array(flex$conc, dim=dim(flex$conc)[1:3])
	flex$conc = NULL
	

	#	load icon output variable 'icon.var' as a 'reg.grid' object. 
	#		'reg.grid' is a class defined in 'myRtools'. The 3D field is stored as icon$zz.  
	idx = which(icon.dtm == flex.out.dtm[ii])
	icon = read.reg.grid.from.ncdf(file.path(icon.dir, icon.fls[idx]), sn=icon.var)
	#	flip vertical order
	icon$zz = icon$zz[,,dim(icon$zz)[3]:1]
	
	#	convert units. ICON units are actually mole fractions (volume mixing ratios). That's not the same as in FLEXPART (mass mixing ratio). Conversion requires molar weight of the tracer in ICON and that of air. 
	icon$zz = icon$zz * mu.tracer/mu.air * 1E12 # mol/mol -> g/g -> ng/kg

	#	compare the total mass in both simulations. 
	###################################################################
	if (ii==1){
		#	1) Get air density from ICON output
		icon.rho = read.reg.grid.from.ncdf(file.path(icon.dir, icon.fls[idx]), sn="rho")
		icon.rho$zz = icon.rho$zz[,,dim(icon.rho$zz)[3]:1]
		#	2) get half-level heights 
		icon.hhl = read.reg.grid.from.ncdf(file.path(icon.dir, icon.fls[idx]), sn="z_ifc")
		icon.hhl$zz = icon.hhl$zz[,,dim(icon.hhl$zz)[3]:1]

		#	calculate full level heights above ground
		icon.hfl = icon.hhl
		icon.hfl$zz = (icon.hhl$zz[,,-1] + icon.hhl$zz[,, -dim(icon.hhl$zz)[3]])/2
		icon.hfl$zz = icon.hfl$zz - rep(icon.hhl$zz[,,1], dim(icon.hfl$zz)[3])

		#	level thickness
		dz = aperm(apply(icon.hhl$zz, 1:2, diff), c(2,3,1))
	
		#	surface are of grid
		A = get.surface.area(grid=icon, units="m^2")
		#	volume of grid cells
		V = rep(A$zz, dim(dz)[3]) * dz
		#	total air mass of grid cells
		mass.air = icon.rho
		mass.air$zz = V * icon.rho$zz

		icon.rho$zz = icon.rho$zz[,,1:33]
		icon.hfl$zz = icon.hfl$zz[,,1:33]
		mass.air$zz = mass.air$zz[,,1:33]

		icon.rho = crop.grid(icon.rho, xrng=xlim, yrng=ylim)
		icon.hfl = crop.grid(icon.hfl, xrng=xlim, yrng=ylim)
		mass.air = crop.grid(mass.air, xrng=xlim, yrng=ylim)
	
		#	for calculation of center of mass 
		lon = array(rep(seq(icon.rho$xrng[1], icon.rho$xrng[2], icon.rho$dx), dim(icon.rho$zz)[2]), 
			dim=dim(icon.rho$zz)[1:3])
		lat = array(rep(seq(icon.rho$yrng[1], icon.rho$yrng[2], icon.rho$dy), 
			each=dim(icon.rho$zz)[1]), dim=dim(icon.rho$zz)[1:3])
		flex.hfl = array(rep(outlev.mids, each=prod(dim(icon.rho$zz)[1:2])), 
			dim=dim(icon.rho$zz)[1:3])
	}


	#	crop everything to common area 
	#########################################################
	icon$zz = icon$zz[,,1:33]

	icon = crop.grid(icon, xrng=xlim, yrng=ylim)
	class(icon) = "reg.grid"
	flex = crop.grid(flex, xrng=xlim, yrng=ylim)
	class(flex) = "reg.grid"

	#	total mass 
	mass.icon = sum(mass.air$zz * icon$zz) * 1E-12	# ng -> kg
	mass.flex = sum(mass.air$zz * flex$zz) * 1E-12
	
	cat("Total mass ICON:", mass.icon, "(kg)\n")
	cat("Total mass FLEXPART:", mass.flex, "(kg)\n")

	#	calculate center of mass
	###################################################################################

	icon$cntr.lon = weighted.mean(lon, mass.air$zz*icon$zz)
	icon$cntr.lat = weighted.mean(lat, mass.air$zz*icon$zz)
	icon$cntr.height = weighted.mean(icon.hfl$zz, mass.air$zz*icon$zz)

	flex$cntr.lon = weighted.mean(lon, mass.air$zz*flex$zz)
	flex$cntr.lat = weighted.mean(lat, mass.air$zz*flex$zz)
	flex$cntr.height = weighted.mean(flex.hfl, mass.air$zz*flex$zz)

	#	fill plume statistics	
	msk = which(flex$zz>0 & icon$zz>0)
	plume.stats$rmse[ii] = rmse(flex$zz[msk], icon$zz[msk])
	plume.stats$r[ii] = cor(flex$zz[msk], icon$zz[msk])
	plume.stats$flex.max.conc[ii] = max(flex$zz)
	plume.stats$icon.max.conc[ii] = max(icon$zz)
	plume.stats$dist.hor[ii] = sphere.distance(icon$cntr.lat, icon$cntr.lon, flex$cntr.lat, 
		flex$cntr.lon)
	plume.stats$dist.vert[ii] = flex$cntr.height - icon$cntr.height 
	plume.stats$flex.cm.lon[ii] = flex$cntr.lon
	plume.stats$flex.cm.lat[ii] = flex$cntr.lat
	plume.stats$flex.cm.height[ii] = flex$cntr.height
	plume.stats$icon.cm.lon[ii] = icon$cntr.lon
	plume.stats$icon.cm.lat[ii] = icon$cntr.lat
	plume.stats$icon.cm.height[ii] = icon$cntr.height

	plume.stats$flex.sd.lon[ii] = weighted.sd(lon, mass.air$zz*flex$zz)
	plume.stats$flex.sd.lat[ii] = weighted.sd(lat, mass.air$zz*flex$zz)
	plume.stats$flex.sd.height[ii] = weighted.sd(flex.hfl, mass.air$zz*flex$zz)
	plume.stats$icon.sd.lon[ii] = weighted.sd(lon, mass.air$zz*icon$zz)
	plume.stats$icon.sd.lat[ii] = weighted.sd(lat, mass.air$zz*icon$zz)
	plume.stats$icon.sd.height[ii] = weighted.sd(icon.hfl$zz, mass.air$zz*icon$zz)


	#	concentration plot for both models with same color scale and different levels
	###################################################################################
	fn.png = file.path(out.dir, interpolation, 
		paste0("multiLevel_concentration_", interpolation, "_", 
			chron.2.string(flex.out.dtm[ii], "%Y%m%d_%H")) )
	open.graphics.device(fn.png, width=15, height=5*length(levels), pointsize=10)
	par(lwd=0.5)
	lay = cbind(matrix(1:(length(levels)*2), ncol=2, byrow=TRUE), 
			rep(length(levels)*2+1, length(levels)))
	layout(lay, widths=c(1,1,lcm(2)))
	par(mar=c(3,3,1,1)+.1)
	for (lev in levels){
		col.scale= fill.2dplot(flex, xlim=xlim, ylim=ylim, zlim=zlim, lev=lev, 
			log.scale=TRUE, borders=TRUE, type="grid", 
			color.palette=brewer.ramp, palette.par=list(name="YlOrRd"), 
			map.db="swiss.cantons", map.source="myRplots", borders.col="gray", plot.grid=TRUE, 
		    nlevels=30, pnts=release, plot.key=FALSE)
		#	center of mass
		points(flex$cntr.lon, flex$cntr.lat, pch=3, col="gray80", cex=1.5, lwd=2*par("lwd"))
		if (lev==levels[1]) textbox(par("usr")[1], par("usr")[4], "FLEXPART", adj=c(0,1), 
			fill="white")
		if (lev==levels[1]) {
			textbox(par("usr")[2], par("usr")[4], chron.2.string(flex.out.dtm[ii]), adj=c(1,1), 
				fill="white")
		}
		textbox(par("usr")[1], par("usr")[3], paste0(outlev.mids[lev], " (m a.g.l.)"), adj=c(0,0))
		col.scale= fill.2dplot(icon, xlim=xlim, ylim=ylim, zlim=zlim, lev=lev, log.scale=TRUE,
	    	borders=TRUE, type="grid", 
			color.palette=brewer.ramp, palette.par=list(name="YlOrRd"), 
			map.db="swiss.cantons", map.source="myRplots", borders.col="gray", plot.grid=TRUE, 
			nlevels=30, pnts=release, plot.key=FALSE, key.title = "C (ng/kg)")
		#	center of mass
		points(icon$cntr.lon, icon$cntr.lat, pch=4, col="gray80", cex=1.5, lwd=2*par("lwd"))
		if (lev==levels[1]) textbox(par("usr")[1], par("usr")[4], "ICON", adj=c(0,1), 
			fill="white")
		textbox(par("usr")[1], par("usr")[3], paste0(outlev.mids[lev], " (m a.g.l.)"), adj=c(0,0),
			fill="white")

	}
	plot.colorpalette(col.scale)

	dev.off()


	#	concentration difference plot 
	###################################################################################
	fn.png = file.path(out.dir, interpolation, 
		paste0("multiLevel_concentrationDifferences_", interpolation, "_", 
			chron.2.string(flex.out.dtm[ii], "%Y%m%d_%H")) )
	open.graphics.device(fn.png, width=8.5, height=5*length(levels), pointsize=7.5)
	par(lwd=0.5)
	lay = cbind(matrix(1:(length(levels)), ncol=1, byrow=TRUE), 
			rep(length(levels)+1, length(levels)))
	layout(lay, widths=c(1,lcm(2)))
	par(mar=c(3,3,1,1)+.1)
	for (lev in levels){
		col.scale= fill.2dplot(flex-icon, xlim=xlim, ylim=ylim, zlim=c(-1,1)*zlim[2], lev=lev, 
			log.scale=TRUE, borders=TRUE, type="grid", 
			color.palette=brewer.ramp, palette.par=list(name="RdBu", col0=NULL), 
			map.db="swiss.cantons", map.source="myRplots", borders.col="gray", plot.grid=TRUE, 
		    nlevels=30, pnts=release, plot.key=FALSE, key.title = "C (ng/kg)")
		points(flex$cntr.lon, flex$cntr.lat, pch=3, col="gray80", cex=1.5, lwd=2*par("lwd"))
		points(icon$cntr.lon, icon$cntr.lat, pch=4, col="gray80", cex=1.5, lwd=2*par("lwd"))
		if (lev==levels[1]) {
			textbox(par("usr")[2], par("usr")[4], chron.2.string(flex.out.dtm[ii]), adj=c(1,1),
				fill="white")
		}
		textbox(par("usr")[1], par("usr")[3], paste0(outlev.mids[lev], " (m a.g.l.)"), adj=c(0,0),
			fill="white")
	}
	plot.colorpalette(col.scale)
	dev.off()
	

}


#	save stats
fn = file.path(out.dir, interpolation, paste0("plume_stats_", interpolation, ".csv"))
save.file(fn, plume.stats)

stop()




#	default plot function for class 'reg.grid' translates to fill.2dplot. By default field 'zz' is plotted. If 3D array then level 1 is displayed. Since there is not flip of levels in 'read.reg.grid.from.ncdf', we need to specify level=80 to plot surface level.
fill.2dplot(icon, xlim=xlim, ylim=ylim, zlim=zlim, lev=lev, log.scale=TRUE, borders=TRUE, 
	type="grid", color.palette=special.rainbow, map.db="worldHires", plot.grid=TRUE, nlevels=30)




#	compare the total mass in both simulations. 
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


