require(myRtools)
require(myRplots)

plot.plume.stats = function(plume.stats, para, ylab, incl.icon=FALSE, log=""){

	if (incl.icon) {
		para.icon = paste0("icon.", para)
		para = paste0("flex.", para)
	}
	xlim = range(plume.stats[[1]]$dtm)
	ylim = range(sapply(plume.stats, function(x, para) range(x[[para]]), para=para))
	if (incl.icon) {
		ylim = range(c(ylim, plume.stats[[1]][[para.icon]]))
	}
	plot.new()
	plot.window(xlim, ylim, log=log)
	axis(2, tck=1, col="gray", lty=2)
	axis.chron(1, tck=1, col="gray", lty=2)
	axis(2)
	axis.chron(1)
	box()
	title(ylab=ylab)
	for (inter in interpolations){
		lines(plume.stats[[inter]]$dtm, plume.stats[[inter]][[para]], col=cols[inter], 
			lwd=2*par("lwd"))
	}
	if (incl.icon){
		lines(plume.stats[[inter]]$dtm, plume.stats[[inter]][[para.icon]], col="gray20", 
			lwd=2*par("lwd"))
		legend("topright", legend=c("ICON", interpolations), col=c("gray20", cols), lty=1)
	} else {
		legend("topright", legend=interpolations, col=cols, lty=1)
	}

}

out.dir = "/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_nnb_interpolation/20181221/images/plume_compare"

interpolations = c("nnb", "rbf", "byc")
cols = brewer.pal(length(interpolations), "Dark2")
names(cols) = interpolations

plume.stats = vector("list", length(interpolations))
names(plume.stats) = interpolations

for (inter in interpolations){
	fn = file.path(out.dir, inter, paste0("plume_stats_", inter, ".csv"))

	plume.stats[[inter]] = readfile.csv(fn)
}

para = "rmse"
ylab = "RMSE (ng/kg)"
incl.icon = FALSE
log=""
fn.png = file.path(out.dir, paste0("plumeCompare_", para))
open.graphics.device(fn.png, height=4, width=8, pointsize=6)
par(mar=c(2,4,1,1)+.1, lwd=0.5)
plot.plume.stats(plume.stats, para=para, ylab=ylab, incl.icon=incl.icon, log=log)
dev.off()

para = "r"
ylab = "Correlation coefficient"
incl.icon = FALSE
log=""
fn.png = file.path(out.dir, paste0("plumeCompare_", para))
open.graphics.device(fn.png, height=4, width=8, pointsize=6)
par(mar=c(2,4,1,1)+.1, lwd=0.5)
plot.plume.stats(plume.stats, para=para, ylab=ylab, incl.icon=incl.icon, log=log)
dev.off()

para = "dist.hor"
ylab = "Distance horizontal (km)"
incl.icon = FALSE
log=""
fn.png = file.path(out.dir, paste0("plumeCompare_", para))
open.graphics.device(fn.png, height=4, width=8, pointsize=6)
par(mar=c(2,4,1,1)+.1, lwd=0.5)
plot.plume.stats(plume.stats, para=para, ylab=ylab, incl.icon=incl.icon, log=log)
dev.off()

para = "dist.vert"
ylab = "Distance vertical (m)"
incl.icon = FALSE
log=""
fn.png = file.path(out.dir, paste0("plumeCompare_", para))
open.graphics.device(fn.png, height=4, width=8, pointsize=6)
par(mar=c(2,4,1,1)+.1, lwd=0.5)
plot.plume.stats(plume.stats, para=para, ylab=ylab, incl.icon=incl.icon, log=log)
dev.off()

para = "sd.lon"
ylab = "Longitude SD (°)"
incl.icon = TRUE
log=""
fn.png = file.path(out.dir, paste0("plumeCompare_", para))
open.graphics.device(fn.png, height=4, width=8, pointsize=6)
par(mar=c(2,4,1,1)+.1, lwd=0.5)
plot.plume.stats(plume.stats, para=para, ylab=ylab, incl.icon=incl.icon, log=log)
dev.off()


para = "sd.lat"
ylab = "Latitude SD (°)"
incl.icon = TRUE
log=""
fn.png = file.path(out.dir, paste0("plumeCompare_", para))
open.graphics.device(fn.png, height=4, width=8, pointsize=6)
par(mar=c(2,4,1,1)+.1, lwd=0.5)
plot.plume.stats(plume.stats, para=para, ylab=ylab, incl.icon=incl.icon, log=log)
dev.off()


para = "sd.height"
ylab = "Height SD (°)"
incl.icon = TRUE
log=""
fn.png = file.path(out.dir, paste0("plumeCompare_", para))
open.graphics.device(fn.png, height=4, width=8, pointsize=6)
par(mar=c(2,4,1,1)+.1, lwd=0.5)
plot.plume.stats(plume.stats, para=para, ylab=ylab, incl.icon=incl.icon, log=log)
dev.off()


para = "max.conc"
ylab = "Max. concentration (ng/kg)"
incl.icon = TRUE
log="y"
fn.png = file.path(out.dir, paste0("plumeCompare_", para))
open.graphics.device(fn.png, height=4, width=8, pointsize=6)
par(mar=c(2,4,1,1)+.1, lwd=0.5)
plot.plume.stats(plume.stats, para=para, ylab=ylab, incl.icon=incl.icon, log=log)
dev.off()

