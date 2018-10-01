
## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(sf)
sf_extSoftVersion()

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(lwgeom)
lwgeom_extSoftVersion()

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
(pt = st_point(c(2,4)))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
(pt_bin = st_as_binary(pt))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
st_as_sfc(list(pt_bin))[[1]]

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
str(pt)
str(st_linestring(rbind(c(0,0), c(0,1), c(1,1))))
str(st_polygon(list(rbind(c(0,0), c(0,1), c(1,1), 
  c(0,0)))))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
st_crs("+proj=longlat +datum=WGS84")  # "Proj.4 string"
st_crs(3857)                          # EPSG code
st_crs(3857)$units                    # reveal units
st_crs(NA)                            # unknown

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library("tmap") # required version 1.11-1 or later
library("tmaptools") # required version 1.2-3 or later
library(sf)
data("World", "metro", package = "tmap")
metro$growth <- (metro$pop2020 - metro$pop2010) / (metro$pop2010 * 10) * 100

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
m1 <- tm_shape(World) +
    tm_polygons("income_grp", palette = "-Blues", 
      title = "Income class", contrast = 0.7, border.col = "grey30", id = "name") +
    tm_text("iso_a3", size = "AREA", col = "grey30", root = 3) +
  tm_shape(metro) +
    tm_bubbles("pop2010", col = "growth", border.col = "black",
      border.alpha = 0.5,
      breaks = c(-Inf, 0, 2, 4, 6, Inf) ,
      palette = "-RdYlGn",
      title.size = "Metro population (2010)", 
      title.col = "Annual growth rate (%)",
      id = "name",
      popup.vars = c("pop2010", "pop2020", "growth")) + 
  tm_format_World() + tm_style_gray(frame.lwd = 2)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
m1 <- tm_shape(World) +
    tm_polygons("income_grp", palette = "-Blues", 
      title = "Income class", contrast = 0.7, border.col = "grey30", id = "name") +
    tm_text("iso_a3", size = "AREA", col = "grey30", root = 3) +
  tm_shape(metro) +
    tm_bubbles("pop2010", col = "growth", border.col = "black",
      border.alpha = 0.5,
      breaks = c(-Inf, 0, 2, 4, 6, Inf) ,
      palette = "-RdYlGn",
      title.size = "Metro population (2010)", 
      title.col = "Annual growth rate (%)",
      id = "name",
      popup.vars = c("pop2010", "pop2020", "growth")) + 
  tm_format("World") + tm_style("gray", frame.lwd = 2)

## ----plot1, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
print(m1)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
m1 <- tm_shape(World) + tm_style("gray", frame.lwd = 2) +
    tm_polygons("income_grp", palette = "-Blues", 
      title = "Income class", contrast = 0.7, border.col = "grey30", id = "name") +
    tm_text("iso_a3", size = "AREA", col = "grey30", root = 3) +
  tm_shape(metro) +
    tm_bubbles("pop2010", col = "growth", border.col = "black",
      border.alpha = 0.5,
      breaks = c(-Inf, 0, 2, 4, 6, Inf) ,
      palette = "-RdYlGn",
      title.size = "Metro population (2010)", 
      title.col = "Annual growth rate (%)",
      id = "name",
      popup.vars = c("pop2010", "pop2020", "growth")) + 
  tm_format("World")

## ----plot2, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
print(m1)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
m0 <- tm_shape(metro) + 
  tm_bubbles(size = "pop2030") +
  tm_format_World() +
  tm_style_cobalt()

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
m0 <- tm_shape(metro) + 
  tm_bubbles(size = "pop2030") +
  tm_format("World") +
  tm_style("cobalt")

## ----plot3, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
print(m0)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
m0 <- tm_shape(metro) +
  tm_style("cobalt") + 
  tm_bubbles(size = "pop2030") +
  tm_format("World")

## ----plot4, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
print(m0)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
m21 <- tm_shape(World) + tm_polygons(c("blue", "red")) + tm_layout(frame.lwd = 1.5)

## ----plot5, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
print(m21)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
m21 <- m21 + tm_layout(bg.color="transparent")

## ----plot6, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
print(m21)

## ---- echo = TRUE, eval = FALSE, mysize=TRUE, size='\\tiny'--------------
## tmap_mode("view")
## m1

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(tmap)
library(spData)
library(spDataLarge)
library(magrittr)
# Add fill layer to nz shape
m1 = tm_shape(nz) + tm_fill() +
  tm_layout(title = "tm_shape(nz) +\n  tm_fill()", title.size = 0.7) + tm_layout(bg.color="transparent")
# Add border layer to nz shape
m2 = tm_shape(nz) + tm_borders() +
  tm_layout("tm_shape(nz) +\n  tm_borders()", title.size = 0.7) + tm_layout(bg.color="transparent")
# Add fill and border layers to nz shape
m3 = tm_shape(nz) + tm_fill() + tm_borders() +
  tm_layout("tm_shape(nz) +\n  tm_fill() +\n  tm_borders()", title.size = 0.7) + tm_layout(bg.color="transparent")
m123 <- tmap_arrange(m1, m2, m3, nrow = 1)

## ----plot7, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
print(m123)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
map_nz = tm_shape(nz) + tm_polygons() + tm_layout(bg.color="transparent")
map_nz1 = map_nz +
  tm_shape(nz_elev) + tm_raster(alpha = 0.7)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny', cache = TRUE--------------
library(elevatr)
nz_elev_new <- get_elev_raster(as(nz, "Spatial"), z=5, src="aws")
nz_water = st_union(nz) %>% st_buffer(22200) %>% 
  st_cast(to = "LINESTRING")
is.na(nz_elev_new) <- nz_elev_new < 0
nz_elev_new <- crop_shape(nz_elev_new, nz_water)
map_nz1 = map_nz +
  tm_shape(nz_elev_new) + tm_raster(alpha = 0.7)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
map_nz2 = map_nz1 +
  tm_shape(nz_water) + tm_lines()
map_nz3 = map_nz2 +
  tm_shape(nz_height) + tm_dots()
map_nz123 <- tmap_arrange(map_nz1, map_nz2, map_nz3)

## ----plot8, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
print(map_nz123)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
nz_region = st_bbox(c(xmin = 1340000, xmax = 1450000,
                      ymin = 5130000, ymax = 5210000)) %>% 
  st_as_sfc() %>% 
  st_set_crs(st_crs(nz_height))
nz_region

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
nz_height_map = tm_shape(nz_elev, bbox = tmaptools::bb(nz_region)) +
  tm_raster(style = "cont", palette = "YlGn", legend.show = TRUE) +
  tm_shape(nz_height) + tm_symbols(shape = 2, col = "red", size = 1) +
  tm_scale_bar(position = c("left", "bottom"))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
nz_map = tm_shape(nz) + tm_polygons() +
  tm_shape(nz_height) + tm_symbols(shape = 2, col = "red", size = 0.1) + 
  tm_shape(nz_region) + tm_borders(lwd = 3) 

## ----plot9, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
library(grid)
nz_height_map
print(nz_map, vp = viewport(0.8, 0.27, width = 0.5, height = 0.5))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
nz_map = nz_map + tm_layout(bg.color="transparent") 

## ----plot10, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
nz_height_map
print(nz_map, vp = viewport(0.8, 0.27, width = 0.5, height = 0.5))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny', message=FALSE, warning=FALSE----
library(spatialprobit)
data(Katrina)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(sf)
sf_katrina <- st_as_sf(Katrina.raw, coords=c("long", "lat"))
st_crs(sf_katrina) <- "+proj=longlat +datum=WGS84"

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
is.na(sf_katrina$days) <- sf_katrina$days == 99999
sf_katrina$Days <- as.Date("2005-08-29") + sf_katrina$days
sf_katrina$y0_90 <- as.numeric(sf_katrina$days < 90)
diff(as.Date(names(table(sf_katrina$Days))))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(car)
sf_katrina$f_owntype <- recode(sf_katrina$owntype, "1='sole'; 2='local'; 3='national'", as.factor=TRUE)
sf_katrina$f_sesstatus <- recode(sf_katrina$sesstatus, "c(1,2)='lower'; c(4,5)='upper'; else='av'", as.factor=TRUE)
sf_katrina$f_sizeemp <- recode(sf_katrina$sizeemp, "1='small'; 2='av'; 3='large'", as.factor=TRUE)
sf_katrina$f_street1 <- recode(sf_katrina$street1, "1='Magazine Street'; 2='Carrollton Avenue'; 3='St. Claude Avenue'", as.factor=TRUE)

## ---- echo = TRUE, eval=FALSE, mysize=TRUE, size='\\tiny'----------------
## library(mapview)
## mapview(sf_katrina)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny', results='hide'------------
library(sf)
SPol2_RT38 <- st_read("SPol2_RT38.gpkg")
#library(RColorBrewer)
rds <- colorRampPalette(c("grey95", "#EB811B"))
dive <- colorRampPalette(c("#0E302F", "grey95", "#EB811B"))

## ----plot27a, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=6, fig.width=6, dev.args=list(family="Fira Sans")----
library(tmap)
tm_shape(SPol2_RT38) + tm_polygons("farms", pal=rds(20), breaks=seq(0, 36.001, length.out=21)) + tm_layout("Farms eligible", legend.outside = TRUE) + tm_layout(bg.color="transparent")

## ----plot27b, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=6, fig.width=6, dev.args=list(family="Fira Sans")----
tm_shape(SPol2_RT38) + tm_polygons(names(SPol2_RT38)[2:5], pal=rds(20), breaks=seq(0, 36.001, length.out=21)) + tm_layout(c("1928-1929", "1930", "1931", "1932"), legend.show = FALSE) + tm_layout(bg.color="transparent")

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
x <- y <- seq(-12.5, 12.5, 1)
x5 <- y5 <- seq(-12.5, 12.5, 5)
xy <- expand.grid(x=x, y=y)
xy$z <- sqrt(xy$x^2 + xy$y^2)
xy$f <- exp(0.7966 + -1.585*log(xy$z))
m <- matrix(xy$f, length(x), length(y))
m1 <- matrix(NA, 5, 5)
s <- c(1,6,11,16,21,26)
for (i in 1:5)
  for (j in 1:5)
    m1[i,j] <- sum(m[s[i]:s[i+1], s[j]:s[j+1]])
m1[3,3] <- 57.4
MIF <- as.matrix(read.csv("asbyMIF.csv", header=FALSE))
MIFgrd <- expand.grid(x=1:5, y=1:5)

## ----plot28, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
layout(matrix(c(1,1,2,1,1,3), 2, 3, byrow = TRUE))
curve(exp(0.7966 + -1.585*log(x)), sqrt(0.5^2 + 0.5^2), 12.5, xlim=c(0, 12.5), xlab="distance km", ylab="function of distance", main="1953 fitted model", col="#EB811B")
plot(c(0.5,5.5), c(0.5,5.5), type="n", axes=FALSE, xaxs="i",
  yaxs="i", main="Original MIF", xlab="", ylab="")
box(col="#EB811B")
text(MIFgrd, label=format(c(MIF), digits=2), cex=0.9)
abline(v=seq(1.5,4.5,1), col="#EB811B")
abline(h=seq(1.5,4.5,1), col="#EB811B")
plot(c(0.5,5.5), c(0.5,5.5), type="n", axes=FALSE, xaxs="i",
  yaxs="i", main="Generated MIF", xlab="", ylab="")
box(col="#EB811B")
text(MIFgrd, label=format(c(m1/sum(c(m1))), digits=2), cex=0.9)
abline(v=seq(1.5,4.5,1), col="#EB811B")
abline(h=seq(1.5,4.5,1), col="#EB811B")
layout(1)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny', results='hide'------------
library(sp)
MIF <- as.matrix(read.csv("asbyMIF.csv", header=FALSE))
GT <- GridTopology(c(-2*4972.927, -2*4972.927), c(4972.927, 4972.927), c(5, 5))
adopts <- names(SPol2_RT38)[grep("adopt", names(SPol2_RT38))]

pre_sim_lut <- function(GT, SPols) {
  require(maptools)
  SPix_MIF <- as(SpatialGrid(GT), "SpatialPixels")
  crds <- coordinates(SPols)
  lut <- matrix(as.integer(NA), ncol=length(SPix_MIF), nrow=length(SPols))
  for (i in 1:length(SPols)) {
    m0 <- elide(SPix_MIF, shift=crds[i,])
    proj4string(m0) <- CRS(proj4string(SPols))
    lut[i,] <- over(m0, as(SPols, "SpatialPolygons"))
  }
  lut
}

## ---- echo = TRUE, mysize=TRUE, size='\\tiny', results='hide'------------
do_sim <- function(SPols, MIF, lut, adopts, farms, t=1) {
  cMIF <- cumsum(c(MIF))
  if (cMIF[length(cMIF)] <  1) cMIF[length(cMIF)] <- 1
  T <- length(adopts)
  xx <- sapply(adopts, function(x) sum(SPols[[x]]))
  y <- SPols[[adopts[t]]]
  yy <- xx[t]
  while (t < T) {
    t_cands <- rep(x=1:length(SPols), times=SPols[[adopts[t]]])
    while (yy < xx[(t+1)]) {
      sender <- sample(t_cands, 1)
      hits <- which(runif(1) > cMIF)
      if (length(hits) == 0) hits <- 1
      hit <- hits[length(hits)]
      ihit <- lut[sender, hit]
      if (!is.na(ihit)) {
        if (y[ihit] < farms[ihit]) {
          y[ihit] <- y[ihit] + 1
          yy <- yy + 1
        }
      }
    }
#    cat(t, yy, "\n")
    t <- t + 1
  }
  y
}

## ---- echo = TRUE, mysize=TRUE, size='\\tiny', results='hide', cache=TRUE----
lut <- pre_sim_lut(GT, as(SPol2_RT38, "Spatial"))
nsim <- 500
res <- matrix(0, nrow=length(as(SPol2_RT38, "Spatial")), ncol=nsim)
set.seed(1)
for (i in 1:nsim) res[, i] <- do_sim(SPols=as(SPol2_RT38, "Spatial"), MIF=MIF, lut=lut, adopts=adopts, farms=SPol2_RT38$farms, t=1)
res1 <- cbind(res, SPol2_RT38$adopt32)
mm1 <- apply(res1, 1, mean)
hope1 <- apply(res1, 2, function(x) cor(x, mm1))
for (i in 1:nsim) res[, i] <- do_sim(SPols=as(SPol2_RT38, "Spatial"), MIF=MIF, lut=lut, adopts=adopts, farms=SPol2_RT38$farms, t=2)
res2 <- cbind(res, SPol2_RT38$adopt32)
mm2 <- apply(res2, 1, mean)
hope2 <- apply(res2, 2, function(x) cor(x, mm2))
for (i in 1:nsim) res[, i] <- do_sim(SPols=as(SPol2_RT38, "Spatial"), MIF=MIF, lut=lut, adopts=adopts, farms=SPol2_RT38$farms, t=3)
res3 <- cbind(res, SPol2_RT38$adopt32)
mm3 <- apply(res3, 1, mean)
hope3 <- apply(res3, 2, function(x) cor(x, mm3))

## ----plot29, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=6, fig.width=6, dev.args=list(family="Fira Sans")----
plot(density(hope1), lty=1, main="Hope-type test (correlation)", xlim=c(0.6,1), ylim=c(0, 47), xlab="")
abline(v=hope1[501], lwd=2, lty=1, col="#EB811B")
lines(density(hope2), lty=2)
abline(v=hope2[501], lwd=2, lty=2, col="#EB811B")
lines(density(hope3), lty=3)
abline(v=hope3[501], lwd=2, lty=3, col="#EB811B")
legend("top", title="start year", legend=c("1928-29", "1930", "1931"), lty=1:3, bty="n", title.col="#EB811B")

## ----plot30, mysize=TRUE, size='\\tiny', fig.show='hide', fig.height=4, fig.width=8, dev.args=list(family="Fira Sans")----
SPol2_RT38$mm1 <- mm1
SPol2_RT38$mm2 <- mm2
SPol2_RT38$mm3 <- mm3
tm_shape(SPol2_RT38) + tm_polygons(c("adopt32", "mm1", "mm2"), pal=rds(20), breaks=seq(0, 25.001, length.out=21)) + tm_layout(c("observed", "mean (start 1928-29)", "mean (start 1930)"), legend.show = FALSE) + tm_layout(bg.color="transparent")

## ---- echo = TRUE, eval=FALSE, mysize=TRUE, size='\\tiny'----------------
## TASSIM <- st_read("TASSIM.gpkg")
## mapview(TASSIM)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(sf)
library(spData)
b506 <- st_read(system.file("shapes/boston_tracts.shp", package="spData")[1])

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
b489 <- b506[b506$censored == "no",]
t0 <- aggregate(b489, list(ids = b489$NOX_ID), head, n = 1)
b94 <- t0[, c("ids", attr(t0, "sf_column"))]

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
st_queen <- function(a, b = a) st_relate(a, b, pattern = "F***T****")
qm1 <- st_queen(b94)
any(sapply(qm1, length) == 0)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
NOX_ID_no_neighs <- b94$ids[which(sapply(qm1, length) == 0)]
b487 <- b489[is.na(match(b489$NOX_ID, NOX_ID_no_neighs)),]
t0 <- aggregate(b487, list(ids = b487$NOX_ID), head, n = 1)
b93 <- t0[, c("ids", attr(t0, "sf_column"))]

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
qm_93 <- st_queen(b93)
class(qm_93) <- "nb"
attr(qm_93, "region.id") <- as.character(b93$ids)

## ---- echo = TRUE, eval = FALSE, mysize=TRUE, size='\\tiny'--------------
## library(ggplot2)
## ggplot(b487) + geom_sf(aes(fill=NOX))

## ----fig26, fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
library(ggplot2)
ggplot(b487) + geom_sf(aes(fill=NOX)) + theme(plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(colour = NA, fill = "transparent"))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(classInt)
set.seed(1)
cI <- classIntervals(b487$NOX, n=7L, style="bclust", verbose=FALSE)
cI

## ----fig25, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
library(RColorBrewer)
ybrpal <- brewer.pal(6, "YlOrBr")
fC <- findColours(cI, ybrpal)
pal <- attr(fC, "palette")
p <- ggplot(b487, aes(NOX)) + stat_ecdf() + ylab("") +
  geom_vline(xintercept=cI$brks, colour="darkgrey", linetype=2) +
  geom_hline(yintercept=c(0, 1), colour="darkgrey", linetype=2) + ylim(c(0, 1))
for (i in seq(along=pal)) {
  p <- p + geom_rect(ymin=-0.05, ymax=0, xmin=cI$brks[i],
    xmax=cI$brks[i+1], fill=pal[i])
}
p + ggtitle("NOX class intervals") + theme(plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(colour = NA, fill = "transparent"))

## ----fig24, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
b487$NOX_ <- factor(findCols(cI), levels=length(pal):1, labels=rev(names(attr(fC, "table"))))
ggplot(b487) + geom_sf(aes(fill=NOX_)) + scale_fill_manual(values=rev(pal)) + theme(plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(colour = NA, fill = "transparent"))

## ---- echo = TRUE, eval = FALSE, mysize=TRUE, size='\\tiny'--------------
## library(tmap)
## qtm(b487, fill="NOX", fill.n=7L, fill.style="bclust")

## ----fig23, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
library(tmap)
null <- capture.output(qtm(b487, fill="NOX", fill.n=7L, fill.style="bclust", layout.bg.color="transparent"))

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
form <- formula(log(median) ~ CRIM + ZN + INDUS + CHAS + I((NOX*10)^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + I(BB/100) + log(I(LSTAT/100)))

## ----sI, echo = TRUE, mysize=TRUE, size='\\tiny'-------------------------
sessionInfo()

