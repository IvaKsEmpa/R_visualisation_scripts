require(myRplots)
require(Reccodes)


icon.dir = "/input/ICON/ivme/rbf/250518"
icon.fn.pattern = "IconPntSr00%H0000"
icon.dir = "/input/ICON/MCH-1/REGGRID_STAG"
icon.fn.pattern = "disp00%H0000"
icon.dt = 60

cosmo.dir = "/input/COSMO/kenda1"
cosmo.fn.pattern = "laf%Y%m%d%H"
cosmo.dt = 1

res.dir = "/project/hes134/projects/IG3IS/icon_vs_cosmo/mch_nnb"


if (!dir.exists(res.dir)) dir.create(res.dir, recursive=TRUE)

#		 sensible   momentum
paras = c("ASHFL_S", "AUMFL_S", "AVMFL_S", "T", "U")
cosmo.typeOfLevel = list("surface", "surface", "surface", "generalVerticalLayer",
	"generalVerticalLayer")
icon.typeOfLevel = list("surface", "surface", "surface", "generalVerticalLayer",
	"generalVerticalLayer")
zlim = list(c(-200, 200), c(-0.5,0.5), c(-0.5,0.5), c(270,300), c(-10,10))

dtm.start = string.2.chron("2018-05-25 02:00:00")
dtm.end = string.2.chron("2018-05-25 23:00:00")
dtm.start = string.2.chron("2021-10-15 02:00:00")
dtm.end = string.2.chron("2021-10-15 23:00:00")
dtm = seq.dates(dtm.start, dtm.end, by=1/24)
for (jj in 1:length(dtm)){

	icon.fn = chron.2.string(dtm[jj], icon.fn.pattern)
	icon.prev = chron.2.string(dtm[jj]-1/24, icon.fn.pattern)
	icon.info = get.grib.file.info(file.path(icon.dir, icon.fn))

	cosmo.fn = chron.2.string(dtm[jj], cosmo.fn.pattern)
	cosmo.prev = chron.2.string(dtm[jj]-1/24, cosmo.fn.pattern)
	cosmo.info = get.grib.file.info(file.path(cosmo.dir, cosmo.fn))
	
	icon.flds = vector("list", length(paras))
	cosmo.flds = vector("list", length(paras))
	
	for (ii in 1:length(paras)){
	
		#	ICON
		icon.flds[[ii]] = grib.get.field(file.path(icon.dir, icon.fn), sn=paras[ii], 
			grib.info=icon.info, typeOfLevel=icon.typeOfLevel[[ii]])
		icon.flds[[ii]]$zz[icon.flds[[ii]]$zz==9999] = NA
		
		if (icon.flds[[ii]]$end.step - icon.flds[[ii]]$start.step > icon.dt){
			cat(icon.flds[[ii]]$end.step, icon.flds[[ii]]$start.step, "\n")
			#	read previous field for deaccumulation
			tmp = grib.get.field(file.path(icon.dir, icon.prev), sn=paras[ii], 
				typeOfLevel=icon.typeOfLevel[[ii]])
	
			#	deaccumulate
			dt = icon.flds[[ii]]$end.step - tmp$end.step
			dt1 = tmp$end.step - tmp$start.step
			dt2 = icon.flds[[ii]]$end.step - icon.flds[[ii]]$start.step
	
			icon.flds[[ii]]$zz = (icon.flds[[ii]]$zz*dt2 - tmp$zz*dt1)/dt
		}
	
		#	COSMO
		cosmo.flds[[ii]] = grib.get.field(file.path(cosmo.dir, cosmo.fn), sn=paras[ii], 
			grib.info=cosmo.info, typeOfLevel=cosmo.typeOfLevel[[ii]])
		cosmo.flds[[ii]]$zz[cosmo.flds[[ii]]$zz==9999] = NA
		
		if (cosmo.flds[[ii]]$end.step - cosmo.flds[[ii]]$start.step > cosmo.dt){
			cat(cosmo.flds[[ii]]$end.step, cosmo.flds[[ii]]$start.step, "\n")
			#	read previous field for deaccumulation
			tmp = grib.get.field(file.path(cosmo.dir, cosmo.prev), sn=paras[ii], 
				typeOfLevel=cosmo.typeOfLevel[[ii]])
		
			#	deaccumulate
			dt = cosmo.flds[[ii]]$end.step - tmp$end.step
			dt1 = tmp$end.step - tmp$start.step
			dt2 = cosmo.flds[[ii]]$end.step - cosmo.flds[[ii]]$start.step
		
			cosmo.flds[[ii]]$zz = (cosmo.flds[[ii]]$zz*dt2 - tmp$zz*dt1)/dt
		}

	}
	
	fn = file.path(res.dir, chron.2.string(dtm[jj], "icon_vs_cosmo_%Y%m%d_%H"))
	open.graphics.device(fn, width=17, height=5*length(paras), pointsize=10)
	par(mar=c(3,3,1,1), lwd=0.5)
	layout(matrix(1:(length(paras)*3), ncol=3, byrow=TRUE), widths=c(1,1,lcm(2)))
	for (ii in 1:length(paras)){
		fill.grid(icon.flds[[ii]], borders=TRUE, zlim=zlim[[ii]], xlim=c(-3, 1), ylim=c(-2,1), 
			map.db="worldHires", plot.key=FALSE)
		ckey = fill.grid(cosmo.flds[[ii]], borders=TRUE, zlim=zlim[[ii]], xlim=c(-3, 1), 
			ylim=c(-2,1), map.db="worldHires", plot.key=FALSE, key.title=paras[ii])
		plot.colorpalette(ckey)
	}
	dev.off()


}

