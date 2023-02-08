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

icon.pattern = "ICON-ART-OEM*"

icon.var = "testtr212"

flex.dir = file.path("/project/ivme/MCH-1/icon-art-BRM-CDO", 
	paste0("fp_output_", interpolation, "_interpolation"), case)

flex.pattern = "grid_conc_20181221110000.nc"

tidx = 3	#	use tidx to get output for specific time 1: 1 hour after start, etc
rel.com = "BEROMUENSTER 212"

out.dir ="/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_nnb_interpolation/20181221/images/plume_compare"

zlim = c(1, 1E5)
zlim.vert = c(1E1,1E4)

mu.tracer = 28.97	#	in the xml file for the icon run it says 6.4E-2. I assume this means kg/mol. g/mol would not make sense. No such molecule exists. 

levels = c(1:33)		#	level to be plotted, counted from t:he bottom

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

#plume.stats = data.frame(
#	dtm=c(1:33),
#	flex.vertical.mean=rep(NA, 33),
#	icon.vertical.mean=rep(NA, 33)
#	)
flex.vertical.mean = array(NA, dim=c(33, nn.t))
icon.vertical.mean = array(NA, dim=c(33, nn.t))
#ii =3
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

	#########################################################
	icon$zz = icon$zz[,,1:33]

	icon = crop.grid(icon, xrng=xlim, yrng=ylim)
	class(icon) = "reg.grid"
	flex = crop.grid(flex, xrng=xlim, yrng=ylim)
	class(flex) = "reg.grid"

	#	fill plume statistics	
	msk = which(flex$zz>0 & icon$zz>0)
	
# 	plume.stats$flex.vertical.mean = apply(flex$zz,c(1,2),mean)
# 	plume.stats$icon.vertical.mean = apply(icon$zz,c(1,2),mean)


	flex.vertical.mean[,ii] <-apply(flex$zz,3,mean)
	icon.vertical.mean[,ii] <-apply(icon$zz,3,mean)
 
 #       mean vertical cross section with longitude and lattitude 
        #       plot them
 #	xlim=range(c(icon.vertical.mean,flex.vertical.mean), na.rm=TRUE)
        xlim=max(c(icon.vertical.mean,flex.vertical.mean), na.rm=TRUE)*c(1E-2,1)
	ylim=range(outlev.mids)
plot(flex.vertical.mean[,ii], outlev.mids, xlim=xlim, ylim=ylim, log="x")
lines(icon.vertical.mean[,ii],outlev.mids, col="red")      


}

#	save stats
# fn = file.path(out.dir, interpolation, paste0("plume_mean_stats_", interpolation, ".csv"))
# save.file(fn, plume.stats)

fn = file.path(out.dir, interpolation, paste0("flex_vertical_mean_", interpolation, ".rda"))
save(file=fn, flex.vertical.mean)

fn = file.path(out.dir, interpolation, paste0("icon_vertical_mean_", interpolation, ".rda"))
save(file=fn, icon.vertical.mean)
#	load(fn) ->  flex.vertical.mean

fn = file.path(out.dir, interpolation, paste0("flex_vertical_mean_", interpolation, ".csv"))
write.csv(flex.vertical.mean, fn)

fn = file.path(out.dir, interpolation, paste0("icon_vertical_mean_", interpolation, ".csv"))
write.csv(icon.vertical.mean, fn)

# 	flex.vertical.mean = read.csv() or read.table 
stop()

# i have to convert data frame to list of arrays to take in account time also. 
