plot.map <- function(theme, poly, color = NULL, ncl = 9, main = NULL) {
  require(RColorBrewer); 
  require(maptools);
  require(classInt)
  int<- classIntervals(theme)$brks
  pal<-brewer.pal(ncl, color)
  cols <- pal[findInterval(theme, int,rightmost.closed = T)]
  plot(plot,col = cols)
  title(main = main)
}

library(spdep)
library(raster)
col.poly <- shapefile("SIMDgreen.shp")
col.data <- slot(col.poly, "data")
attach(col.data)
colnames(col.data)

plot.map(as.numeric(HlthRank), col.poly, main = "HlthRank")


col.nb <- poly2nb(col.poly, queen = F)
col.listw <- nb2listw(col.nb, style = "C")
plot(col.poly, col = "gray", border = "white")
coords<- coordinates(col.poly)
plot(col.nb, coords, add = T)

lm.ols<-lm(as.numeric(HlthDprs) ~ as.numeric(green300pe))
summary(lm.ols)
coef(lm.ols)
lm.ols.res <- residuals(lm.ols)
par(mfrow = c(2,2))
plot.map(Rankv2, col.poly, main = "HlthRank")
plot.map(lm.ols.res, col.poly,main = "Residuals")
par(mfrow = c(1,1))


moran.test(as.numeric(HlthRank), col.listw)
moran.test(lm.ols.res, col.listw)
shapiro.test(lm.ols.res)
library(lmtest)
bptest(lm.ols)


n<-length(col.nb)
C<-listw2mat(col.listw)
M<-diag(1,n)-1/n
MCM <-M%*%C%*%M
E<-eige(MCM)$vectors

X<-cbind(1,INC,HOVAL)
M<-diag(1,n)-tcrossprod(X%*%qr.solve(crossprod(X)),X)
cbind(M%*%HlthRank,lm.ols.res)
MCM<-M%*%C%*%M
eig<-eigen(MCM)
E<-eig$vector

ones<-rep(1,n)
mi<-eig$values*n/crossprod(ones,C%*%ones)
plot(mi,ylim = c(-1,1), pch = 20, xlab = "Eigenvector Spatial Pattern", ylab = "Moran's I")
abline(0,0,lty=3)


library(sm)
for(i in 1:n){
  plot.map(E[,i],col.poly,main = paste("EV",i,":Moran's I = ",round(mi[i],3)));
  pausee()
}

sf.err<- SpatialFiltering(as.numeric(HlthSMR)~as.numeric(green300pe),
                          data = col.poly, nb = col.nb)
sf.err<- SpatialFiltering(lm.ols, nb = col.nb, style = "C", alpha = 0.25, ExactEV = T)
E.sel <-fitted(sf.err)
dim(E.sel)
par(mfrow = c(2,2))
for(i in 1:4)
  plot.map(E.sel[,i],col.poly,main = colnames(E.sel)[i])
par(mfrow = c(1,1))


lm.sf<-lm(HlthSMR~green300pe+E.sel)
summary(lm.sf)
lm.sf.res <- residuals(lm.sf)
par(mfrow = c(2,2))
plot.map(HlthRank, col.poly, main="HlthRank")
plot.map(lm.sf.res, col.poly, main="SF Residuals")
par(mfrow = c(1,1))

moran.test(as.numeric(HlthDprs), col.listw)
moran.test(lm.sf.res, col.listw)
shapiro.test(lm.sf.res)
bptest(lm.sf)
################

library(spdep)
library(raster)
col.poly <- shapefile("SIMDgreen.shp")
col.data <- slot(col.poly, "data")
attach(col.data)
colnames(col.data)

col.nb <- poly2nb(col.poly, queen = F)
col.listw <- nb2listw(col.nb, style = "C")


sf.err<- SpatialFiltering(as.numeric(HlthSMR)~IncomePer + as.numeric(CrimeRate) , data = col.poly, nb = col.nb)
E.sel <-fitted(sf.err)
data <- data.frame(OBJECTID,E.sel)
newdata <- merge(col.poly, data)
raster::shapefile(newdata, "aaaaa.shp")

lm.sf<-lm(as.numeric(HlthSMR)~IncomePer + as.numeric(CrimeRate)  + E.sel)
summary(lm.sf)
lm.sf.res <- residuals(lm.sf)
moran.test(lm.sf.res, col.listw)

sf.err<- SpatialFiltering(as.numeric(HlthSMR)~IncomePer + as.numeric(CrimeRate) + green300pe , data = col.poly, nb = col.nb)
E.sel1 <-fitted(sf.err)
data1 <- data.frame(OBJECTID,E.sel1)
newdata1 <- merge(col.poly, data1)
raster::shapefile(newdata1, "bbbbb.shp")
lm.sf<-lm(as.numeric(HlthSMR)~IncomePer + as.numeric(CrimeRate) + green300pe  + E.sel)
summary(lm.sf)

sf.err<- SpatialFiltering(as.numeric(HlthSMR)~ green300pe , data = col.poly, nb = col.nb)
E.sel2 <-fitted(sf.err)
lm.sf<-lm(as.numeric(HlthSMR)~ green300pe   + E.sel2)
summary(lm.sf)
lm.sf.res <- residuals(lm.sf)
moran.test(lm.sf.res, col.listw)
