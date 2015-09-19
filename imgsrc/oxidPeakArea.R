
pdf('../imgbin/oxidPeakArea.pdf', width=7, height=3.3)
par(mar = c(1/32, 1/32, 1.8+1/8, 1/32), oma = c(0, 0, 0, 0), mgp=c(2.2,1,0))

plot(1, xlim = c(-0.25, 1), ylim = c(-0.05, 1.0), type='n', axes=F, xlab="", ylab="", xaxt='n', yaxt='n', bty="n")
title("3D chromatogram constructed from a run of LC-MS 
      that analyzes the mono-oxidized forms and unoxidized form of a peptide", cex.main=1)

polygon(c(0, -0.25, 0.7, 0.95), c(0, 0.5, 0.5, 0), border=NA, col = "#DDDDDD")

RT = seq(0,1,0.001) 

RT1  = RT
TIC1 = 0.90*exp(-15^2*(RT1-0.75)^2*1.5)
RT2  = RT
TIC2 = (
  0.05*exp(-16^2*(RT2-0.60)^2) 
+ 0.07*exp(-25^2*(RT2-0.44)^2)
+ 0.12*exp(-15^2*(RT2-0.40)^2)
+ 0.06*exp(-10^2*(RT2-0.32)^2)
+ 0.09*exp(-12^2*(RT2-0.20)^2)
)

for(i in 2:length(RT2)) {rect(RT2[i-1]-0.2, 0.4, RT2[i]-0.2, 0.4+TIC2[i], col="#FF8888", border=NA)}
text(0.40-0.2-0.01/3, 0.40+0.08,  " P' \n peak-area of \n mono-oxidized forms ", cex=0.75)
for(i in 2:length(RT1)) {rect(RT1[i-1]-0.05, 0.1, RT1[i]-0.05, 0.1+TIC1[i], col="#8888FF", border=NA)}
text(0.75-0.05-0.01/3, 0.10+0.08, " P  \n peak-area of \n unoxidized form ", cex=0.75)

Y = sum(TIC1-0.2); X = sum(TIC2-0.6)

arrows(0, 0,  1, 0, length=0.1)
arrows(0, 0,  0, 0.8, length=0.1)
arrows(0, 0, -0.28, 0.56, length=0.1)
text( 0.5,      -0.05,    "Retention time (RT)",    srt=0,         cex=0.99)
text(-0.02,     0.4,      "Absolute intensity",     srt=90,        cex=0.99)
text(-0.14-0.02, 0.28-0.02, "m/z", srt=104+asin(1/2)*(360/2/pi), cex=0.99)

text(0.35, 0.75, expression(paste("Relative frequency of oxidation is estimated to be ", "P'" %/% "(P'+P)")), cex=0.75)

dev.off()