########################################################################
### Code for the analysis in the JSS paper for the R package spTest  ###
########################################################################

################################
### install/load packages    ###
################################

# install.packages("geoR")
# install.packages("spTest")
# install.packages("splines")

library("spTest")
library("fields")
library("geoR")
library("splines")
library("MASS")
library("rgdal")

#####################################
###  Code for plots in Figure 1   ###
#####################################
library("mvtnorm")
matern.2d <- function(h, phi, alpha, nu) {
    a <- (pi * phi)/(2^(nu - 1) * gamma(nu + 1) * alpha^(2 * nu))
    d <- h * alpha
    d[d == 0] <- 1e-10
    b <- (d)^nu
    c <- besselK(d, nu)
    
    cov.val <- a * b * c
    return(cov.val)
}

# Function to transform coordinates to create anisotropy
coords.aniso <- function (coords, aniso.pars, reverse = FALSE) {
    coords <- as.matrix(coords)
    n <- nrow(coords)
    if (length(aniso.pars) != 2) 
        stop("argument aniso.pars must be a vector with 2 elements: the anisotropy angle and anisotropy ratio, respectively")
    psiA <- aniso.pars[1]
    psiR <- aniso.pars[2]
    if (psiR < 1) {
        psiR <- round(psiR, digits = 8)
        if (psiR < 1) 
            stop("anisotropy ratio must be greater than 1")
    }
    rm <- matrix(c(cos(psiA), -sin(psiA), sin(psiA), cos(psiA)), 
                 ncol = 2)
    tm <- diag(c(1, 1/psiR))
    if (reverse) 
        coords.mod <- coords %*% solve(rm %*% tm)
    else coords.mod <- coords %*% rm %*% tm
    return(coords.mod)
}

set.seed(2017)

x <- seq(0, 12, by = 0.25)
y <- seq(0, 12, by = 0.25)

coords <- expand.grid(x, y)
n <- dim(coords)[1]

aniso.angle <- 3 * pi/4
aniso.ratio <- 2
coordsA <- coords.aniso(coords, c(aniso.angle, aniso.ratio))
Da <- as.matrix(dist(coordsA))

sigma.sq <- 1.0
alpha <- 1/2
sm <- 1
tau.sq <- 0

d <- seq(0, 20, by = 0.1)
phi <- sigma.sq * alpha^(2 * sm) * gamma(sm + 1)/ (pi * gamma(sm))
mcov2 <- matern.2d(d, phi, alpha, sm)
er.index <- which(round(mcov2, 2) == 0.05)
er <- d[er.index]
er
er <- er[1]

D <- as.matrix(dist(coords))
R <- matern.2d(D, phi, alpha, sm)
R <- R + diag(tau.sq, nrow = n, ncol = n)

Ra <- matern.2d(Da, phi, alpha, sm)
Ra <- Ra + diag(tau.sq, nrow = n, ncol = n)

z <- rmvnorm(1, rep(0, n), R, method = "chol")
z <- z - mean(z)
zmat <- matrix(z, nrow = length(y), ncol = length(x), byrow = TRUE)
image.plot(x, y, zmat)

x.shift <- -102
y.shift <- 37
image.plot(x + x.shift, y + y.shift, zmat,
           ylab = "", xlab = "",
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"),
           cex.main = 1.80, pty = "s", legend.cex = 0)
map("state", add = TRUE, lwd = 2)
mtext("Stationary and Isotropic Spatial Process", side = 3, line = 0.5, cex = 1.25, font = 2)
mtext("Longitude", side = 1, line = 2)
mtext("Latitude", side = 2, line = 2)

zmat.iso <- zmat

# pdf(file = "isotropic.pdf")
x.shift <- -102
y.shift <- 37
par(mai = c(0.8, 0.4, 0.8, 0.2), pty = "s")
image.plot(x + x.shift, y + y.shift, zmat,
           ylab = "", xlab = "",
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"),
           cex.main = 1.80, pty = "s", legend.cex = 0)
map("state", add = TRUE, lwd = 2)
mtext("Stationary and Isotropic Spatial Process", side = 3, line = 0.5, cex = 1.25, font = 2)
mtext("Longitude", side = 1, line = 2)
mtext("Latitude", side = 2, line = 2)
# dev.off()

########################
# Anisotropic Process
########################

set.seed(2)
z <- rmvnorm(1, rep(0, n), Ra, method = "chol")
z <- z - mean(z)
zmat <- matrix(z, nrow = length(y), ncol = length(x), byrow = TRUE)
image.plot(x, y, zmat)

zmat.anis <- zmat

x.shift <- -102
y.shift <- 37
par(mai = c(0.8, 0.4, 0.8, 0.2), pty = "s")
image.plot(x + x.shift, y + y.shift, zmat,
           ylab = "", xlab = "",
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"),
           cex.main = 1.80, pty = "s", legend.cex = 0)
map("state", add = TRUE, lwd = 2)
mtext("Stationary and Anisotropic Spatial Process", side = 3, line = 0.5, cex = 1.25, font = 2)
mtext("Longitude", side = 1, line = 2)
mtext("Latitude", side = 2, line = 2)

zmat.iso <- zmat

# pdf(file = "anisotropic.pdf")
x.shift <- -102
y.shift <- 37
par(mai = c(0.8, 0.4, 0.8, 0.2), pty = "s")
image.plot(x + x.shift, y + y.shift, zmat,
           ylab = "", xlab = "",
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"),
           cex.main = 1.80, pty = "s", legend.cex = 0)
map("state", add = TRUE, lwd = 2)
mtext("Stationary and Anisotropic Spatial Process", side = 3, line = 0.5, cex = 1.25, font = 2)
mtext("Longitude", side = 1, line = 2)
mtext("Latitude", side = 2, line = 2)
# dev.off()

###############################
###  gridded data example   ###
###############################

data("WRFG", package = "spTest")

### not included in JSS paper ###
image.plot(WRFG$lon-360, WRFG$lat, WRFG$WRFG.NCEP.tas, 
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"),
           legend.lab = "Temp (K)",legend.cex = 0.8,legend.line = 2.2,
           xlab = "Longitude", ylab = "Latitude", 
           main = "Mean WRFG-NCEP Temperatures")
world(add = TRUE)
#################################

coords <- expand.grid(WRFG$xc, WRFG$yc)
sub <- which(coords[, 1] > 2.90e6 & coords[, 1] < 3.95e6 & 
             coords[, 2] > 1.2e6 & coords[, 2] < 2.25e6)
coords.ll <- cbind((WRFG$lon-360)[sub], WRFG$lat[sub])
image.plot(WRFG$lon-360, WRFG$lat, WRFG$WRFG.NCEP.tas, 
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"), 
           legend.lab = "Temp (K)", legend.cex = 0.8, legend.line = 2.2, 
           xlab = "Longitude", ylab = "Latitude", 
           main = "Mean WRFG-NCEP Temperatures")
world(add = TRUE)
left <- seq(1, 400, by = 20)
right <- seq(20, 400, by = 20)
for (i in 2:20) {
    segments(coords.ll[i-1, 1], coords.ll[i-1, 2], coords.ll[i, 1], coords.ll[i, 2], lwd = 2)
    segments(coords.ll[left[i-1], 1], coords.ll[left[i-1], 2],
             coords.ll[left[i], 1], coords.ll[left[i], 2], lwd = 2)
    segments(coords.ll[right[i-1], 1], coords.ll[right[i-1], 2],
             coords.ll[right[i], 1], coords.ll[right[i], 2], lwd = 2)
    j <- i + 380
    segments(coords.ll[j-1, 1], coords.ll[j-1, 2], coords.ll[j, 1], coords.ll[j, 2], lwd = 2)
}

### not included in JSS paper ###
image.plot(WRFG$xc, WRFG$yc, WRFG$WRFG.NCEP.tas, 
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"), 
           legend.lab = "Temp (K)", legend.cex = 0.8, legend.line = 2.2, 
           xlab = "Easting", ylab = "Northing", 
           main = "Map Projection Coordinates")
my.coords <- coords[sub, ]
rect(min(my.coords[, 1]), min(my.coords[, 2]), 
     max(my.coords[, 1]), max(my.coords[, 2]), lwd = 1.5)
#################################

tas <- c(WRFG$WRFG.NCEP.tas)[sub]
x.coord <- unique(my.coords[, 1])
y.coord <- unique(my.coords[, 2])
nx <- length(x.coord)
ny <- length(y.coord)
tas.mat <- matrix(tas, nrow = nx, ncol = ny, byrow = FALSE)
image.plot(x.coord, y.coord, tas.mat, 
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"), 
           legend.lab = "Temp (K)", legend.cex = 0.8, legend.line = 2.2, 
           ylab = "Northing", xlab = "Easting", 
           main = "Subset of Temperatures")

### not included in JSS paper ###
hist(tas, freq = FALSE, xlab = "Temperature, K", main = "Histogram of Temperature")
z.tas <- (tas - mean(tas))/sd(tas)
qqnorm(z.tas)
abline(0, 1)
#################################

tas.geodat <- as.geodata(cbind(my.coords[, 1], my.coords[, 2], tas))
plot( variog4(tas.geodat), ylab = "estimated semivariogram")
title("Directional Sample Semivariograms")

# temp.data <- as.matrix(cbind(my.coords, tas))
my.delta <- 50000
mylags <-  rbind(c(1, 0), c(0, 1), c(1, 1), c(-1, 1))
myA <-  rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))
tr <- GuanTestGrid(spdata = tas.geodat, delta = my.delta, lagmat = mylags, A = myA, 
                   window.dims = c(4, 4), pt.est.edge = TRUE, sig.est.edge = TRUE, 
                   sig.est.finite = TRUE, df = 2)
tr$alternative <- NULL
tr$sigma.hat <- NULL
tr

### not included in JSS paper ###
par(mfrow = c(1, 2), oma = c(0, 0, 1, 0))
plot(my.coords[, 1], tas, xlab = "Easting", ylab = "Temperature")
plot(tas, my.coords[, 2], xlab = "Temperature", ylab = "Northing")
title(main = "Scatters of Temperature vs. Coordinates" , line = -1, outer = TRUE)
#################################

m1 <- lm(tas ~ ns(my.coords[, 1], df = 3) + ns(my.coords[, 2], df = 3))
summary(m1)

### not included in JSS paper ###
stres <- studres(m1)
hist(stres)
qqnorm(stres)
abline(0, 1)
#################################

resid.mat <- matrix(studres(m1), nrow = nx, ncol = ny, byrow = FALSE)
image.plot(x.coord, y.coord, resid.mat, 
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"), 
           xlab = "Easting", ylab = "Northing")
title("Heat Map of Residuals")

resid.geo <- as.geodata(as.matrix(cbind(my.coords[, 1:2], studres(m1))))
plot(variog4(resid.geo), ylab = "estimated semivariogram")
title("Directional Sample Semivariograms")

tr <- GuanTestGrid(spdata = resid.geo, delta = my.delta, lagmat = mylags, A = myA, 
                   window.dims = c(4, 4), pt.est.edge = TRUE, sig.est.edge = TRUE, 
                   sig.est.finite = TRUE, df = 2)
tr$p.value.finite

################################
### non-gridded data example ###
################################
data("COmonthlyMet", package = "fields")
sub30 <- CO.ppt[74:103, , ]
nstations <- 376
years <- 1968:1997
nyears <- length(years)
yr.avg <- matrix(data = NA_real_, nrow = nstations, ncol = nyears)
for (i in 1:nyears) {
    yr.dat <- sub30[i, , ]
    yr.avg[, i] <- apply(yr.dat, 2 , mean, na.rm = TRUE)
}
avg30 <- apply(yr.avg, 1, mean, na.rm = TRUE)
CO.loc.utm <- project(as.matrix(CO.loc), "+proj=utm +zone=13 ellps=WGS84") / 1e+05
quilt.plot(CO.loc.utm, log(avg30), 
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"), 
           legend.lab = "Precip (log mm)", legend.cex = 0.8, legend.line = 2.2, 
           xlab = "Easting", ylab = "Northing", 
           main = "Quilt Plot of log(precip)")
mp <- map("state", region = c("colorado", "wyoming", "nebraska", "utah", "new mexico", "oklahoma"),
          plot = FALSE)
states <- project(cbind(mp$x, mp$y), "+proj=utm +zone=13 ellps=WGS84") / 1e+05
points(states[, 1], states[, 2], type = "l", lwd = 1.5)

quilt.plot(CO.loc.utm, CO.elev, 
           col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"), 
           legend.lab = "Elevation (meters)", legend.cex = 0.8, legend.line = 2.7, 
           xlab = "Easting", ylab = "Northing", 
           main = "Quilt Plot of Elevation")
points(states[, 1], states[, 2], type = "l", lwd = 1.5)

### not included in JSS paper ###
par(mfrow = c(1, 2))
plot(CO.loc.utm[, 1], log(avg30), xlab = "Easting", ylab = "log(precip)", 
     main = "Scatter of log(precip) vs. Easting")
plot( log(avg30), CO.loc.utm[, 2], xlab = "log(precip)", ylab = "Northing", 
     main = "Scatter of Northing vs log(precip)")
#################################

plot(CO.elev, log(avg30) , xlab = "Elevation", ylab = "log(precip)", 
     main = "Scatter of log(precip) vs. Elevation")
m1 <- lm(log(avg30) ~ ns(CO.elev, df = 3))
summary(m1)
fits <- m1$fitted.values
bad <- is.na(avg30)
x <- CO.elev[which(!bad)]
lines(sort(x), fits[order(x)], lwd = 3, col = "red")

qqnorm(studres(m1))
abline(0, 1)

precip.resid <- cbind(CO.loc.utm[which(!bad), ][, 1], CO.loc.utm[which(!bad), ][, 2] , studres(m1))
precip.geo <- as.geodata(precip.resid) 
plot(variog4(precip.geo), xlab = "distance (100s of kilometers)",
     ylab = "estimated semivariogram", legend = FALSE)
legend("bottomright",
       legend = c(expression(0 * degree), expression(45 * degree),
                  expression(90 * degree), expression(135 * degree)), col = 1:4, lty = 1:4)
title("Directional Sample Semivariograms")

### not included in JSS paper ###
###################
# Test CSR
###################
plot(precip.resid[, 1:2])
library("spatstat")
min(precip.resid[, 1])
max(precip.resid[, 1])
min(precip.resid[, 2])
max(precip.resid[, 2])
my.xlims <-  c(min(precip.resid[, 1]), max(precip.resid[, 1]))
my.ylims <-  c(min(precip.resid[, 2]), max(precip.resid[, 2]))
mywin <- owin(xrange = my.xlims, yrange = my.ylims)
pp <- ppp(precip.resid[, 1], precip.resid[, 2], window = mywin)

# A test for complete spatial randomness
hopskel.test(pp)
###############################

mylags <- rbind(c(0.60, 0), c(0, 0.60), c(0.45, 0.45), c(-0.45, 0.45))
myA <- rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))

quilt.plot(precip.resid[, 1:2], precip.resid[, 3], 
col = two.colors(n = 256, start = "blue3", end = "red3", middle = "gray60"), 
xlab = "Longitude", ylab = "Latitude", xlim = c(0.82, 8.75), ylim = c(40.0, 46.3))
title("Quilt Plot of Residuals and Grid Used for Subsampling")
c(min(precip.resid[, 1]), max(precip.resid[, 1]))
c( min(precip.resid[, 2]), max(precip.resid[, 2]))

### not included in JSS paper ###
xdist <- max(precip.resid[, 1]) - min(precip.resid[, 1])
ydist <- max(precip.resid[, 2]) - min(precip.resid[, 2])
xdist
ydist
area <- xdist * ydist
area
n <- dim(precip.resid)[1]
pts.area <- n / area
pts.area
target.points.per.window <- sqrt(n)
target.points.per.window
target.window.area <- target.points.per.window/pts.area
# change these parameters to change the proposed window/block size
x.div <- 16
y.div <- 12
window.dim <- c(4, 3)
proposed.window.area <- (window.dim[1] * xdist / x.div) * (window.dim[2] * ydist / y.div) 
proposed.window.area
target.window.area

proposed.points.per.window <- proposed.window.area * pts.area
proposed.points.per.window
target.points.per.window
#################################

# The limits of the sampling region should be 
# slightly larger than the min/max coordinates.
tol <- 0.02
my.xlims <- c(min(precip.resid[, 1]) - tol, max(precip.resid[, 1]) + tol)
my.ylims <- c( min(precip.resid[, 2]) - tol, max(precip.resid[, 2]) + tol)
xlen <-  my.xlims[2] - my.xlims[1]
ylen <-  my.ylims[2] - my.ylims[1]
my.grid.spacing <-  c(xlen / 16, ylen / 12)
xgrid <- seq(my.xlims[1], my.xlims[2], by = my.grid.spacing[1])
ygrid <- seq(my.ylims[1], my.ylims[2], by = my.grid.spacing[2])
segments(x0 = xgrid, y0 = min(my.ylims), y1 = max(my.ylims), lty = 2)
segments(x0 = min(my.xlims), y0 = ygrid, x1 = max(my.xlims), lty = 2)

myh <-  0.70
myh.sb <-  0.85

tr.guan <- GuanTestUnif(spdata = precip.geo, lagmat = mylags, A = myA, df = 2, h = myh, 
                        kernel = "norm", truncation = 1.5, xlims = my.xlims, ylims = my.ylims, 
                        grid.spacing = my.grid.spacing, window.dims = c(4, 3), subblock.h = myh.sb)
tr.guan$estimate
tr.guan$p.value.finite

### not included in JSS paper ###
# compare estimates of the semivariance to those
# given by geoR
tr.guan$estimate
directions <-  c(0, 90, 45, 135)
dv <- variog4(precip.geo)
names(dv)
# variog gives angles in azimuth, so 0 degrees in variogg4 is north-south, 90 degress is east-west
dists <- c(dv$"90"$u[2], dv$"0"$u[2], dv$"45"$u[2], dv$"135"$u[2])
distance <- round(dists, 3)
semivario <- c(dv$"90"$v[2], dv$"0"$v[2], dv$"45"$v[2], dv$"135"$v[2])
sv <- round(semivario, 3)
rbind(directions, distance, sv)
tr.guan$estimate

plot(dv, legend = FALSE)
abline(v = 0.60)

# Check sensitivity to bandwidth parameters
myh <- 0.6
myh.sb <- 0.7
tr.guan <- GuanTestUnif(spdata = precip.geo, lagmat = mylags, A = myA, df = 2, h = myh, 
                        kernel = "norm", truncation = 1.5, xlims = my.xlims, ylims = my.ylims, 
                        grid.spacing = my.grid.spacing, window.dims = c(4, 3), subblock.h = myh.sb)
tr.guan$p.value.finite

myh <-  0.8
myh.sb <-  0.9
tr.guan <- GuanTestUnif(spdata = precip.geo, lagmat = mylags, A = myA, df = 2, h = myh, 
                        kernel = "norm", truncation = 1.5, xlims = my.xlims, ylims = my.ylims, 
                        grid.spacing = my.grid.spacing, window.dims = c(4, 3), subblock.h = myh.sb)
tr.guan$p.value.finite
###############################

set.seed(2016)
tr.maity <- MaityTest(spdata = precip.geo, lagmat = mylags, A = myA, df = 2, 
                      xlims = my.xlims, ylims = my.ylims, grid.spacing = my.grid.spacing, 
                      block.dims = c(4, 3), nBoot = 100, kappa = 1)
tr.maity$p.value

### not included in JSS paper ###
# check sensitivity to the number of bootstrap samples
tr.maity <- MaityTest(spdata = precip.geo, lagmat = mylags, A = myA, df = 2, 
                      xlims = my.xlims, ylims = my.ylims, grid.spacing = my.grid.spacing, 
                      block.dims = c(4, 3), nBoot = 200, kappa = 1)
tr.maity$p.value
#################################

