bty="n"

get.mass <- function (R) {
  switch(R,
         A = 71.037114  ,
         R = 156.101111 ,
         N = 114.042927 ,
         D = 115.026943 ,
         C = 103.009185 ,
         E = 129.042593 ,
         Q = 128.058578 ,
         G = 57.021464  ,
         H = 137.058912 ,
         I = 113.084064 ,
         L = 113.084064 ,
         K = 128.094963 ,
         M = 131.040485 ,
         F = 147.068414 ,
         P = 97.052764  ,
         S = 87.032028  ,
         T = 101.047679 ,
         W = 186.079313 ,
         Y = 163.06332  , 
         V = 99.068414
         ) 
}

get.y.height <- function (R) {
  switch(R,
         A = 0.2,
         D = 0.2,
         E = 0.2,
         F = 0.267,
         G = 0.333,
         H = 0.333,
         I = 0.2,
         K = 0.133,
         L = 0.2,
         M = 0.2,
         N = 0.267,
         P = 0.53,
         Q = 0.2,
         R = 0.133,
         S = 0.267,
         T = 0.267,
         V = 0.2,
         Y = 0.267
         ) 
}

get.b.height <- function (R) {
  switch(R,
         A = 0.2,
         D = 0.2,
         E = 0.2,
         F = 0.267,
         G = 0.267,
         H = 0.2,
         I = 0.267,
         K = 0.067,
         L = 0.267,
         M = 0.267,
         N = 0.267,
         P = 0.53,
         Q = 0.2,
         R = 0.067,
         S = 0.267,
         T = 0.267,
         V = 0.267,
         Y = 0.267
  ) 
}

discretize <- function(num) {
  rpois(1,lambda=num)
}
discretize.many <- function(vec) {
  rpois(length(vec), lambda=vec)
}

get.theoretical <- function(peptide) {
  chars = strsplit(substring(peptide, 1,nchar(peptide)),split='')[[1]] #19.0184
  prec.mass = Reduce(function(a,b){a+get.mass(b)}, chars, 18.01057+1.00782, accumulate=F)
  prec.peak = rpois(1,150)
  prec.mass = c(prec.mass, (prec.mass + 1.00782)/2)
  prec.peak = c(prec.peak, rpois(1,50))
  mz.b = NULL
  mz.bt = 1.0078
  pk.b = NULL
  rs.b = strsplit(substring(peptide, 1,nchar(peptide)-1),split='')[[1]]
  NH3.b.idx = NULL
  NH3.b.pk = NULL
  H2O.b.idx = NULL
  H2O.b.pk = NULL
  i = 0
  for (R in rs.b) {
    i = i + 1
    mz.bt = mz.bt + get.mass(R)
    mz.b = c(mz.b, mz.bt)
    pk.b = c(pk.b, discretize(get.b.height(R)*1000))
    if (R %in% c('R', 'N', 'Q')) {
      NH3.b.idx = c(NH3.b.idx, i); NH3.b.pk = c(NH3.b.pk, rbeta(1,1,1)*pk.b[length(pk.b)]);pk.b[length(pk.b)]=rbeta(1,1,1)*pk.b[length(pk.b)]
    }
    if (R %in% c('S', 'T')) {
      H2O.b.idx = c(H2O.b.idx, i); H2O.b.pk = c(H2O.b.pk, rbeta(1,1,1)*pk.b[length(pk.b)]);pk.b[length(pk.b)]=rbeta(1,1,1)*pk.b[length(pk.b)]
    }
  }
  mz.y = NULL
  mz.yt = 19.0184
  pk.y = NULL
  rs.y = rev(strsplit(substring(peptide, 2,nchar(peptide)),split='')[[1]])
  NH3.y.idx = NULL
  NH3.y.pk = NULL
  H2O.y.idx = NULL
  H2O.y.pk = NULL
  i = 0
  print(NH3.y.idx)
  for (R in rs.y) {
    i = i + 1
    mz.yt = mz.yt + get.mass(R)
    mz.y = c(mz.y, mz.yt)
    pk.y = c(pk.y, discretize(get.y.height(R)*1000))
    if (R %in% c('R', 'N', 'Q')) {
      NH3.y.idx = c(NH3.y.idx, i);NH3.y.pk = c(NH3.y.pk, rbeta(1,2,2)*pk.y[length(pk.y)]);pk.y[length(pk.y)]=rbeta(1,2,2)*pk.y[length(pk.y)]
    }
    if (R %in% c('S', 'T', 'E')) {
      H2O.y.idx = c(H2O.y.idx, i);H2O.y.pk = c(H2O.y.pk, rbeta(1,2,2)*pk.y[length(pk.y)]);pk.y[length(pk.y)]=rbeta(1,2,2)*pk.y[length(pk.y)]
    }
  }
  print(NH3.y.idx)
  max.height = max(pk.b,pk.y,H2O.b.pk,H2O.y.pk,NH3.b.pk,NH3.y.pk)
  #max.height = 1
  pk.b = pk.b / max.height
  pk.y = pk.y / max.height
  H2O.b.pk = H2O.b.pk / max.height
  H2O.y.pk = H2O.y.pk / max.height
  NH3.b.pk = NH3.b.pk / max.height
  NH3.y.pk = NH3.y.pk / max.height
  #pk.y = pk.y/max.height
  #pk.b = pk.b/max.height
  pdf(paste0('../imgbin/',peptide,'.pdf'), width=10, height=5)
  par(mar = c(3.5, 3.5, 1/2, 1/2), oma = c(0, 0, 0, 0), mgp=c(2.2,1,0))
  plot(0,0, xlim=c(0,700),ylim=c(-0.00, 1.15), xlab='m/z', xaxt='n', yaxt='n', ylab='Relative intensity', type='n')
  axis(side=1, at=50*(0:20), labels=50*(0:20))
  axis(side=2, at=c(0, 0.5, 1), labels=c('0%', '50%', '100%'))
  rect(-999, -999, 9999, 0, col='grey')
  lines(50+rchisq(1000,2)^2*1000, rpois(1000,10)/max.height, type='h')
  lines(50+rchisq(100,2)^2*1000,  discretize.many(20*rchisq(100,1))/max.height, type='h')
  
  if (length(NH3.b.idx)>0)
  for (i in 1:length(NH3.b.idx)) {
    lines(mz.b[NH3.b.idx[i]]-17.031,   NH3.b.pk[i], type='h', col='blue')
    text(mz.b[NH3.b.idx[i]]-17.031,    NH3.b.pk[i]+0.08, labels=substitute(atop(b[i]-N*H[3], mz), list(i = NH3.b.idx[i], mz = round(mz.b[NH3.b.idx[i]]-17.031,   2))), col='blue')
  }
  if (length(H2O.b.idx)>0)
  for (i in 1:length(H2O.b.idx)) {
    lines(mz.b[H2O.b.idx[i]]-18.01057, H2O.b.pk[i], type='h', col='blue')
    text(mz.b[H2O.b.idx[i]]-18.01057,  H2O.b.pk[i]+0.08, labels=substitute(atop(b[i]-H[2]*O, mz), list(i = H2O.b.idx[i], mz = round(mz.b[H2O.b.idx[i]]-18.01057, 2))), col='blue')
  }
  if (length(NH3.y.idx)>0)
  for (i in 1:length(NH3.y.idx)) {
    lines(mz.y[NH3.y.idx[i]]-17.031,   NH3.y.pk[i], type='h', col='red')
    text(mz.y[NH3.y.idx[i]]-17.031,    NH3.y.pk[i]+0.08, labels=substitute(atop(y[i]-N*H[3], mz), list(i = NH3.y.idx[i], mz = round(mz.y[NH3.y.idx[i]]-17.031,   2))), col='red')
  }
  if (length(H2O.y.idx)>0)
  for (i in 1:length(H2O.y.idx)) {
    lines(mz.y[H2O.y.idx[i]]-18.01057, H2O.y.pk[i], type='h', col='red')
    text(mz.y[H2O.y.idx[i]]-18.01057,  H2O.y.pk[i]+0.08, labels=substitute(atop(y[i]-H[2]*O, mz), list(i = H2O.y.idx[i], mz = round(mz.y[H2O.y.idx[i]]-18.01057, 2))), col='red')
  }
  
  lines(mz.b, pk.b, type='h', col='blue')
  lines(mz.y, pk.y, type='h', col='red')
  for (i in 1:length(rs.b)) {
    text(mz.b[i], pk.b[i]+0.08, labels=substitute(atop(b[i],mz), list(i = i, mz = round(mz.b[i],2))), col='blue')
  }
  for (i in 1:length(rs.y)) {
    text(mz.y[i], pk.y[i]+0.08, labels=substitute(atop(y[i],mz), list(i = i, mz = round(mz.y[i],2))), col='red')
  }
  points(prec.mass, prec.peak/max.height, type='h', col='green')
  for (prec.idx in 1:length(prec.mass)){
    print('AAA')
    #z = paste(rep('+',prec.idx), sep='', collapse='')
    z = if(prec.idx==1){''}else{prec.idx}
    text(prec.mass[prec.idx], prec.peak[prec.idx]/max.height+0.08, col='green',
         labels=substitute(atop((M+z*H)^{z*'+'}, m), list(z=z, m=round(prec.mass[prec.idx], 2)))
         )
  }
  dev.off()
}
get.theoretical('LTAKGR')
#get.theoretical('ANELLLNVK')