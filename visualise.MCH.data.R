#	PURPOSE
#		visualise wind field taken from COSMO/ICON simulations
#		Currently ignores staggered character of grid
#		Vertical wind is in absolute coordintes.

#	PRELIMINARIES

require(myRplots)
require(Reccodes)		#	make sure to load grib_api module before starting R
require(sp)


#	SETTINGS
#######################################################################################

grib.path = "/input/ICON/MCH-1/REGGRID_STAG/regrid_icon_rbf/"
fn = "disp00230000"
dx = 0.55
dy = 0.55
plt.path = "/project/ivme/MCH-1/REGGRID_STAG/images" #/nnb_interpolation/images"
zlim=c(-1., 1.)
levels = c(1, 5, 10, 20, 25)	#	levels to visualise (flipped vertical order: 0 at surface)

#	location of Beromünster in geographic coordinates
pts = SpatialPoints(matrix(c(8.175482, 47.189569), byrow=TRUE, ncol=2), 
	proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

#	MAIN
#######################################################################################
#	get listing of all grib fields in file
grib.info = get.grib.file.info(file.path(grib.path, fn))

#	get vertical wind speed; vertical levels get flipped when reading -> level 0 at surface
ww = grib.get.field(file.path(grib.path, fn), sn="W", typeOfLevel="generalVertical")
uu = grib.get.field(file.path(grib.path, fn), sn="U", typeOfLevel="generalVerticalLayer")
vv = grib.get.field(file.path(grib.path, fn), sn="V", typeOfLevel="generalVerticalLayer")
tt = grib.get.field(file.path(grib.path, fn), sn="T", typeOfLevel="generalVerticalLayer")
hhl = grib.get.field(file.path(grib.path, fn), sn="HHL", typeOfLevel="generalVertical")

#	convert points to rotated coordinates
pts.rot = spTransform(pts, CRS(ww$proj4string), use_ob_tran=TRUE)

xlim = c(-1,1)*dx+pts.rot@coords[1,1]
ylim = c(-1,1)*dy+pts.rot@coords[1,2]
ww=crop.grid(ww,xrng=xlim,yrng=ylim)
uu = crop.grid(uu, xrng=xlim, yrng=ylim)
vv = crop.grid(vv, xrng=xlim, yrng=ylim)
hhl = crop.grid(hhl, xrng=xlim, yrng=ylim)
#       SETTINGS
#######################################################################################

grib.path1 = "/input/ICON/MCH-1/REGGRID_STAG/regrid_icon_nnb/21101500/"

#       MAIN
#######################################################################################


#       get listing of all grib fields in file
grib.info1 = get.grib.file.info(file.path(grib.path1, fn))

#       get vertical wind speed; vertical levels get flipped when reading -> level 0 at surface
ww1 = grib.get.field(file.path(grib.path1, fn), sn="W", typeOfLevel="generalVertical")
hhl1 = grib.get.field(file.path(grib.path1, fn), sn="HHL", typeOfLevel="generalVertical")
tt1 = grib.get.field(file.path(grib.path1, fn), sn="T", typeOfLevel="generalVertical")
ww1=crop.grid(ww1,xrng=xlim,yrng=ylim)
#	construct fields containing uu, vv as individual components. Used to plot vector components
names(uu)[names(uu)=="zz"] = "uu"
uu$vv = vv$zz
uu$vv[] = 0.

names(vv)[names(vv)=="zz"] = "vv"
vv$uu = uu$uu
vv$uu[] = 0.
class(ww) = 'reg.grid'
class(ww1) = 'reg.grid'

for (lev in levels){

	open.graphics.device(file.path(plt.path, paste0(fn, "_wind_rbf_", sprintf("%02d", lev))), 
		width=12, height=10)
	par(mar=c(4,4,1,1)+.1, lwd=0.5)
	
	#	plot level 10 of vertical wind speed
	fill.grid(ww, borders=TRUE, lev=lev, 
		key.title="w (m/s)", xlab="Rot. longitude (°E)", ylab="Rot. latitude (°N)",
		map.db=c("swiss.cantons"), borders.col=c(1), map.fill=c(FALSE), 
		map.source="myRplots", map.lwd=2*par("lwd"), grid.borders="gray30")
	points(pts.rot, pch=4, lwd=2, cex=1.5, col="red")
	
	fill.contour(hhl, add=TRUE, contour.lines=TRUE, nlevels=10, col="transparent",
		cntr.drawlabels=FALSE, cntr.col="gray30")
	
	#	add vectors of horizontal wind. Considering staggered grid!
	#fill.2dplot(uu, type="vector", lev=lev, add=TRUE, vec.x.skip=0, vec.y.skip=0,
	#	vec.col="darkorange4", vec.ref=2, vec.len=0.01)
#	fill.2dplot(vv, type="vector", lev=lev, add=TRUE, vec.x.skip=0, vec.y.skip=0,
      	 #   vec.col="darkorange4", vec.ref=2, vec.len=0.01)
	
	dev.off()

}



