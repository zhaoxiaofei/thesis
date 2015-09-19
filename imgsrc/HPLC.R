#library("Hmisc")

N1 = 6
N2 = 4
N3 = 8
N4 = 6

rt.end <- 9
rt <- 0:rt.end

x1 = rt 
y1 = pmax(0, rt)^1.5  /5.5
y1[y1>1] = -10
for (i in 2:(length(y1))){
    if(y1[i]>1 && y1[i-1] > 1) {
      y1[i] = 9
    }
}
x1 = rep(x1, N1) + (runif(length(x1)*N1)-1/2)/2  
y1 = rep(y1, N1) + rnorm(length(x1), sd=sqrt(pmax(0,x1-0))/2)/30 

x2 = rt 
y2 = pmax(0, rt-1)^1.5 /8.5
y2[y2>1] = -10
for (i in 2:(length(y2))){
  if(y2[i]>1 && y2[i-1] > 1) {
    y2[i] = 9
  }
}
x2 = rep(x2, N2) + (runif(length(x2)*N2)-1/2)/2  
y2 = rep(y2, N2) + rnorm(length(x2), sd=sqrt(pmax(0,x2-1))/2)/30 

x3 = rt 
y3 = pmax(0, rt-2)^1.6 /14
y3[y3>1] = -10
for (i in 2:(length(y3))){
  if(y3[i]>1 && y3[i-1] > 1) {
    y3[i] = 9
  }
}
x3 = rep(x3, N3) + (runif(length(x3)*N3)-1/2)/2  
y3 = rep(y3, N3) + rnorm(length(x3), sd=sqrt(pmax(0,x3-3))/2)/30 

x4 = rt 
y4 = pmax(0, rt-3)^1.6 /14
y4[y4>1] = -10
for (i in 2:(length(y4))){
  if(y4[i]>1 && y4[i-1] > 1) {
    y4[i] = 9
  }
}
x4 = rep(x4, N4) + (runif(length(x4)*N4)-1/2)/2  
y4 = rep(y4, N4) + rnorm(length(x4), sd=sqrt(pmax(0,x4-3))/2)/30 

pdf('../imgbin/RPLC.pdf', width=8.5, height=4)
par(mar = c(3.5, 1/2, 1/2, 1/2), oma = c(0, 0, 0, 0), mgp=c(2.2,1,0))
plot(1, 1, xlim = c(0-0.3,rt.end+0.3), ylim = c(1.4, -0.17), type='n', 
      xaxt='n', 
     yaxt='n', xlab='Retention time (RT) in minute' #, ylab='Position in Column'
     )
axis(side=1, at=rt, labels=rt*10)
cols <- rainbow(rt.end+4, alpha=0.4)
cols <- gray.colors(rt.end+1, start = 0.9, end = 0.5, gamma = 2.2, alpha = 1)
for (i in rt) {
  rect(i-0.33,0.9,i+0.33,-0.1*0,col=cols[i+1])
  arrows(i,0.5-(i+2)*0.013,i,0.5+(i+2)*0.013, length=0.1)
  #arrows(i,0.7-(i+2)*0.013,i,0.7+(i+2)*0.013, length=0.1)
}
points(x4, y4, col=adjustcolor('blue',  alpha.f = 0.9), pch=4)
points(x3, y3, col=adjustcolor('cyan',  alpha.f = 0.8), pch=3)
points(x2, y2, col=adjustcolor('green', alpha.f = 0.7), pch=2)
points(x1, y1, col=adjustcolor('red',   alpha.f = 0.6), pch=1)
legend('topright',legend=c('Very hydrophilic molecules', 
                                'Hydrophilic molecules', 
                                'Hydrophobic molecules', 
                           'Very hydrophobic molecules'),
       pch=c(1,2,3,4),col=c('red', 'green', 'cyan', 'blue'), pt.bg = cols[5], #cex=c(1,1,1,1,2),
       bg='#FFFFFF'
       )
text(2.7, -0.13, cex=0.66, labels=paste0(
    'Each rectangle is the column containing both the stationary phase and the mobile phase\n', 
    'Darkness in each rectangle is the proportion of high-elution-strength components in the mobile phase\n', 
    'Each downward arrow is the direction of and strength of elution in mobile phase'))
abline(1.05,0)
RT = seq(0,9,0.01) 
TIC1 = 0.3*exp(-5^2*(RT-3)^2/2)*5/5
TIC2 = 0.2*exp(-4^2*(RT-5)^2/2)*4/5
TIC3 = 0.4*exp(-3^2*(RT-7)^2/2)*3/5
TIC4 = 0.3*exp(-3^2*(RT-8)^2/2)*3/5
lines(RT,          1.45-(TIC1+TIC2+TIC3+TIC4), lwd=1)
lines(RT[200:400], 1.45-TIC1[200:400],col='red',   lty=2, lwd=2)
lines(RT[400:600], 1.45-TIC2[400:600],col='green', lty=2, lwd=2)
lines(RT[600:800], 1.45-TIC3[600:800],col='cyan',  lty=2, lwd=2)
lines(RT[700:900], 1.45-TIC4[700:900],col='blue',  lty=2, lwd=2)
text(4.5, 1.1, cex=0.88, labels=
       'Chromatogram: detected quantity of molecules that exit the column as a function of RT')
dev.off()