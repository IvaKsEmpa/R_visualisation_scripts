require(myRplots)
require(Reccodes)

icon.dir = "/input/ICON/ivme/icon_art_pntSrc/Shweta_input/rbf/"
icon.fn.pattern = "IconCdoMergeTest%H"
icon.dt = 60

res.dir = "/project/ivme/pictures_plots_images/"


if (!dir.exists(res.dir)) dir.create(res.dir, recursive=TRUE)

#		 sensible   momentum
paras = c("ASHFL_S", "AUMFL_S", "AVMFL_S", "ASOB_S", "RAIN_GSP")
sprintf("%s i am here in paras", paras)
icon.typeOfLevel = list("surface", "surface", "surface", "surface",
        "surface")
zlim = list(c(-200, 200), c(-200, 200), c(-200, 200), c(-200, 200), c(-200, 200))
#zlim = list(c(-200, 200), c(-0.5,0.5), c(-0.5,0.5), c(270,300), c(-10,10))
dtm.start = string.2.chron("2021-08-07 01:00:00")
dtm.end = string.2.chron("2021-08-07 04:00:00")
dtm = seq.dates(dtm.start, dtm.end)
sprintf("%f dtm", dtm)
sprintf("%d length(dtm)", length(dtm))
for (jj in 1:length(dtm)){
        cat("jj loop, jj=", jj)
        print(jj) 
	icon.fn = chron.2.string(dtm[jj], icon.fn.pattern)
	cat("icon.fn", icon.fn)
	icon.prev = chron.2.string(dtm[jj]-1/24, icon.fn.pattern)
	icon.info = get.grib.file.info(file.path(icon.dir, icon.fn))

        	
	icon.flds = vector("list", length(paras))
        	
	for (ii in 1:length(paras)){
                cat("ii loop", ii)	
		#	ICON
		icon.flds[[ii]] = grib.get.field(file.path(icon.dir, icon.fn), sn=paras[ii], 
			grib.info=icon.info, typeOfLevel=icon.typeOfLevel[[ii]])
#		icon.flds[[ii]]$zz[icon.flds[[ii]]$zz==9999] = NA
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
	

	}
	
	fn = file.path(res.dir, chron.2.string(dtm[jj], "remappedWithFieldextraIconNetcdfOutput_%Y%m%d_%H"))
	open.graphics.device(fn, width=17, height=5*length(paras), pointsize=10)
	par(mar=c(3,3,1,1), lwd=0.5)
	layout(matrix(1:(length(paras)*2), ncol=2, byrow=TRUE), widths=c(1,lcm(2)))
	for (ii in 1:length(paras)){
	
		ckey=fill.grid(icon.flds[[ii]], borders=TRUE, zlim=zlim[[ii]], xlim=c(-3, 1), ylim=c(-2,1), 
			map.db="worldHires", plot.key=FALSE, key.title=paras[ii])
		plot.colorpalette(ckey)
	}
	dev.off()


}

