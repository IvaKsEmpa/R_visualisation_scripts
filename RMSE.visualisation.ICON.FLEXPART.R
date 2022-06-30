require(myRtools)
require(myRplots)
require(Rflexpart)
require(ncdf4)
require(plyr)
require(ggplot2)
#	SETTINGS

icon.dir = "/project/ivme/MCH-1/icon-art-BRM/icon_scripts_output_code/icon_art_output/211218/structured/"
icon.flist <- list.files(icon.dir,pattern='*.nc',full.names=FALSE)
cat("seq_along(icon.flist)", seq_along(icon.flist))
my_rmse=0
my_rmse2=0
my_rmse3=0
flex.dir = "/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_rbf_interpolation/211218"
flex.fn = "grid_conc_20181221110000.nc"

flex2.dir = "/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_rbf_interpolation/211218"

flex3.dir = "/project/ivme/MCH-1/icon-art-BRM-CDO/fp_output_byc_interpolation/211218"

xlim = c(7.,9.)
ylim = c(46.75, 47.5)
zlim = c(1E1, 1E4)
mu.tracer = 28.97       #       in the xml file for the icon run it says 6.4E-2. I assume this means kg/mol. g/mol would not make sense. No such molecule exists. 


levels= c(5,10,15,20)           #       level to be plotted, counted from the bottom


#       constants
mu.air = 28.97
R.air = 287

for(k in 13:length(icon.flist)) {
	cat("i am in the nc loop")
cat("length(icon.flist)",length(icon.flist))
icon.var = "testtr12"

#	END OF SETTINGS
tidx =k-12

cat("k",k)
#	data is in 'conc' as a 5D array. Dimensions 4 and 5 only contain one entry. 
flex = read.grid.output(file.path(flex.dir, flex.fn), release="BEROMUENSTER 12", tidx=tidx)
flex2 = read.grid.output(file.path(flex.dir, flex.fn), release="BEROMUENSTER 12", tidx=tidx)
flex3 = read.grid.output(file.path(flex.dir, flex.fn), release="BEROMUENSTER 12", tidx=tidx)
#	drop not used dimension of conc
flex$conc = array(flex$conc, dim=dim(flex$conc)[1:3])
flex2$conc = array(flex2$conc, dim=dim(flex2$conc)[1:3])
flex3$conc = array(flex3$conc, dim=dim(flex3$conc)[1:3])
str(crop.grid(flex, xrng=c(7,10), yrng=c(46, 48), para="conc"))
str(crop.grid(flex2, xrng=c(7,10), yrng=c(46, 48), para="conc"))
str(crop.grid(flex3, xrng=c(7,10), yrng=c(46, 48), para="conc"))
#	release location 
release = list(x=flex$hdr$rel.lng1, y=flex$hdr$rel.lat1, pch=4)

cat("icon.flist[k]",icon.flist[k])
icon = read.reg.grid.from.ncdf(file.path(icon.dir, icon.flist[k]), sn=icon.var)
#	flip vertical order
icon$zz = icon$zz[,,dim(icon$zz)[3]:1]
icon$zz = icon$zz[,,1:33]
str(crop.grid(icon, xrng=c(7,10), yrng=c(46, 48)))

#	convert units. ICON units are actually mole fractions (volume mixing ratios). That's not the same as in FLEXPART (mass mixing ratio). Conversion requires molar weight of the tracer in ICON and that of air. 
icon$zz = icon$zz * mu.tracer/mu.air * 1E12 # mol/mol -> g/g -> ng/kg

icon = crop.grid(icon, xrng=xlim, yrng=ylim)
flex = crop.grid(flex, para="conc", xrng=xlim, yrng=ylim)
flex2 = crop.grid(flex2, para="conc", xrng=xlim, yrng=ylim)
flex3 = crop.grid(flex3, para="conc", xrng=xlim, yrng=ylim)
layout(matrix(1:3, ncol=3), widths=c(1,1,lcm(2)))
par(mar=c(4,4,2,1)+.1)

my_rmse[k-12] <- rmse(icon$zz, flex$conc, na.rm = FALSE)
my_rmse2[k-12] <- rmse(icon$zz, flex2$conc, na.rm = FALSE)
my_rmse3[k-12] <- rmse(icon$zz, flex3$conc, na.rm = FALSE)
}

j<-13:length(icon.flist)
df <- data.frame(j,my_rmse, my_rmse2, my_rmse3)
png(paste0("/project/ivme/MCH-1/icon-art-BRM-CDO/scripts/rmse_3plots211218rbf.png"), height=12, width=24, units="cm", res=800)
plot(j, my_rmse,xlab="Time", ylab="RMSE", col = "red")
#g <- ggplot(df, aes(j))                    # basic graphical object
#g <- g + geom_line(aes(y=my_rmse), colour="red") # first layer
#g <- g + geom_line(aes(y=my_rmse2), colour="green")  # second layer
#g <- g + geom_line(aes(y=my_rmse3), colour="blue")  # third layer
#g <- g + ylab("RMSE") + xlab("Time")
#print(g)
dev.off()

