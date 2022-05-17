#

# evolutionary allometry analysis:
##############
R.size.pgls <- procD.pgls(f1 = coords~csize, phy = tree, data = R.df)
summary(R.size.pgls)

R.size.pgls.plot <- plot(R.size.pgls, type = "regression", predictor = R.size,reg.type ="RegScore",
                         bg = pts.col, pch = pts.pch, cex = rescale.numeric(R.size, to = c(0.5,2.25)),lwd = 3,
                         main = "Evolutionary Allometry (PGLS)", xlab = "log-transformed centroid size", ylab = "cranial shape")

R.size.pgls.preds <- shape.predictor(R.size.pgls$GM$pgls.fitted, x= R.size.pgls.plot$RegScore, Intercept = FALSE,
                        predmin = min(R.size.pgls.plot$RegScore),
                        predmax = max(R.size.pgls.plot$RegScore))


R.size.subfam <- procD.pgls(f1 = coords ~ csize*subfamily, phy = tree, data=R.df)
anova(R.size.subfam)
# different evolutionary allometric slopes per tribe

# does the interaction affect the fit of the full model?
R.size.subfam2 <- procD.pgls(f1 = coords ~ csize+subfamily, phy = tree, data=R.df)
anova(R.size.subfam2)

anova(R.size.subfam2,R.size.subfam)
# ya
##############

# interspecific morphospace
# figure 1
###############
as.proc <- gpagen(ruminants)
allspecs <- as.proc$coords
allspecs.size <- as.proc$Csize

allspecs.scores <- matrix(nrow=dim(allspecs)[3],ncol=dim(R.pca$rotation)[2])
for(i in 1:dim(allspecs)[3]){
  x <- rotonto(mshape(R),allspecs[,,i],scale=F)$yrot - mshape(R)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  allspecs.scores[i,] <- a.preds[1,] %*% R.pca$rotation
}

family.cols <- c()
family.pch <- c()
 for(i in 1:NROW(allspecs.scores)){
   name <- name.id[which(name.id[,1]==si.FV[i]),3]#specimenInfo$species.FV[i]),3]
   nf <- full.inventory1[which(rownames(full.inventory1)==name),2]
   if(nf=="Cervidae"){family.cols[[i]] <- addTrans(family.colors[1],255*0.15) ; family.pch[[i]] <- 21}else{
   if(nf=="Bovidae"){family.cols[[i]] <- addTrans(family.colors[2],255*0.15); family.pch[[i]] <- 22}else{
   if(nf=="Tragulidae"){family.cols[[i]] <- addTrans(family.colors[3],255*0.15); family.pch[[i]] <- 24}else{
   if(nf=="Moschidae"){family.cols[[i]] <- addTrans(family.colors[4],255*0.15); family.pch[[i]] <- 25}}}}}
 family.cols <- unlist(family.cols)
  family.pch <- unlist(family.pch)

 cervid.scores <- allspecs.scores[which(family.cols == addTrans(family.colors[1],255*0.15)),]
 cervid.chull <- chull(cervid.scores[,1:2])
 cervid.chull <- c(cervid.chull,cervid.chull[1])
 bovid.scores <- allspecs.scores[which(family.cols == addTrans(family.colors[2],255*0.15)),]
 bovid.chull <- chull(bovid.scores[,1:2])
 bovid.chull <- c(bovid.chull,bovid.chull[1])
 moschid.scores <- allspecs.scores[which(family.cols == addTrans(family.colors[4],255*0.15)),]
 moschid.chull <- chull(moschid.scores)
 moschid.chull <- c(moschid.chull,moschid.chull[1])
 tragulid.scores <- allspecs.scores[which(family.cols == addTrans(family.colors[3],255*0.15)),]
 tragulid.chull <- chull(tragulid.scores[,1:2])
 tragulid.chull <- c(tragulid.chull,tragulid.chull[1])


# R.bm <- mvBM(R.df$tree, data= R.pca$x[,1:2])

 R.d <- distRoot(R.df$tree)[1]



# Interspecific morphospace - PC1 & 2

# with densities
par(fig=c(0,.9,0,.9),mar=c(4,4,1,1))
# ADD ALLOMETRY VECTOR
   plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=1,cex.lab=1.2,
       xlab = paste("PC1: ",toString(round(R.pca$d[1]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("PC2: ",toString(round(R.pca$d[2]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

points(allspecs.scores[,1],allspecs.scores[,2],cex = rescale(log(allspecs.size),c(0.5,2.5)),col=family.cols,bg=family.cols,pch=family.pch)
  lines(cervid.scores[cervid.chull,1:2],col=family.colors[1],lwd=3)
  lines(bovid.scores[bovid.chull,1:2],col=family.colors[2],lwd=3)
  lines(moschid.scores[moschid.chull,1:2],col=family.colors[4],lwd=3)
  lines(tragulid.scores[tragulid.chull,1:2],col=family.colors[3],lwd=3)

  mixtools::ellipse(mu = as.numeric(R.bm$theta)[1:2], sigma = R.bm$sigma[1:2,1:2]*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
  mixtools::ellipse(mu = as.numeric(R.bm$theta)[1:2], sigma = R.bm$sigma[1:2,1:2]*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)

  addtree2(R.pca$x[,1:2],tree = R.df$tree)
  points(R.pca$x[,1],R.pca$x[,2],pch=pts.pch,bg = pts.col,lwd = 3, cex = rescale(log(R.size),c(0.5,2.25)))
  points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 2, lwd = 2)
  legend("topleft", c("Tragulidae","Cervidae","Moschidae","Bovidae","Common Ancestor"), col = 'black',pt.lwd = 2, pt.cex = 2,pt.bg = c(family.colors[c(3,1,4,2)],"goldenrod1"), pch = c(24,21,25,22,23),title.adj = -1, cex = 1.25)


# adding the density plots:
w <- table(si.FV)

# make vector of the relative weights of each specimen, based on the amount of specimens representing the species in this dataset
wts <- 1:length(si.FV)
names(wts) <- si.FV
for(i in 1:length(si.FV)){
  wts[i] <- 1/w[which(names(w)==si.FV[i])]
}

all.fams <- c()
for(i in 1:length(family.cols)){
  a <- name.id[which(name.id[,1] == si.FV[i]),3]
  all.fams[[i]] <- R.df$family[which(R.df.n == a)]
}
all.fams <- unlist(all.fams)


library(sm)

uaf <- unique(all.fams) ; uaf <- uaf[c(2,3,4,1)]

wts2 <- wts
wts2 <- (wts2/sum(wts2))


# if you're in Rstudio, you'll have to fidget around with the console size to get the density curves to match
par(fig=c(0,0.9,0.78,1),new=T)
plot(NA,xlim=c(-0.2,0.2),ylim=c(0,6),xlab=NA,ylab=NA,
     frame = F, axes = F)
#axis(4,labels=FALSE)
for(i in 1:4){
famid <- which(all.fams == uaf[i])
wtsf <- wts2[famid]
dens <- density(allspecs.scores[famid,1],weights=wtsf/sum(wtsf))
dens$y <- dens$y * sum(wtsf)
polygon(dens,col = addTrans(family.colors[i],255*0.1), border = family.colors[i], lwd =2)
}


par(fig=c(0.785,1,0,0.9),new=T)
plot(NA,ylim=c(-0.2,0.2),xlim=c(0,6),xlab=NA,ylab=NA,
     frame = F, axes = F)
#axis(3, labels = FALSE)
for(i in 1:4){
famid <- which(all.fams == uaf[i])
wtsf <- wts2[famid]
dens <- density(allspecs.scores[famid,2],weights=wtsf/sum(wtsf))
dens$y <- dens$y * sum(wtsf)
dens2 <- dens
dens$x <- dens2$y
dens$y <- dens2$x
polygon(dens,col = addTrans(family.colors[i],255*0.1), border = family.colors[i], lwd =2)
}


par(mfrow = c(1,1), mar = c(5,5,2,2))
# THE PHYLOGENETIC DISTRIBUTION SCHEMATIC
###
library(mvtnorm)
x.points <- seq(-0.5,0.5,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
mu <- R.bm$theta
sigma <- R.bm$sigma*R.d
for (i in 1:100) {
for (j in 1:100) {
z[i,j] <- mvtnorm::dmvnorm(c(x.points[i],y.points[j]),
mean=mu,sigma=sigma)
}
}
sigma.null <- sigma ; sigma.null[1,2] <- 0 ; sigma.null[2,1] <- 0


par(mar=c(5,5,2,2))
  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=1,cex.lab=1.2,
       xlab = "PC1",ylab="PC2")
.filled.contour(x.points,y.points,z, levels = seq(0,max(z),length.out = 50), col = colorRampPalette(c("gray10","gray85"))(50))
    mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
   mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)
   addtree2(R.pca$x[,1:2], tree = R.df$tree, col = 'white')
  points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 2, lwd = 2)
shape::colorlegend(col= colorRampPalette(c("gray10","gray85"))(50), zlim=c(0,max(z)), zlevels = NULL,
            posx = c(0.81,0.91),posy=c(0.06,0.26), main = "expected density\nof species", main.col = "white", main.cex = 1.2)

# 95 & 99% labels put in in illustrator because I'm lazy


# SHAPE VARIATION ALONG THE MORPHOSPACE
####
density.2d <- kde2d.weighted(allspecs.scores[,1], allspecs.scores[,2], w = wts,
                             lims = c(range(allspecs.scores[,1]),range(allspecs.scores[,2])))

par(mar=c(5,5,2,2))
  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1.2,
       xlab = "PC1",ylab="PC2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 2, lwd = 2)
    mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
   mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)


xseg <- c(-0.15,0,0.15)   ; yseg <- c(-0.15,0,0.15)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 2, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.pca$rotation[,1]*scores[i,1]+R.pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  if(i==5){SHADE<-F}else{SHADE<-T}
  scrs <- as.matrix(scores[i,])
  up <- c(scrs[1],scrs[2]+0.0375)
  down <- c(scrs[1]-0.0071,scrs[2]-0.025)
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(R)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.165, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = up,
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)
  tlm <- shapes1[dor.drop,-2,i] ; tlm[,2] <- -tlm[,2]
  rlm <- mshape(R)[dor.drop,-2] ; rlm[,2] <- -rlm[,2]
  tps(target.lm = tlm,
      reference.lm = rlm,
      n.grid.col = 24,
      scale.lm = 0.16, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = down,
      links = links.dor,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)
}

text(c(-0.025,-0.025), y=c(0.069,-0.004),labels=c("lateral","dorsal"))


# DIET SURFACE
# BROWSING-GRAZING LANDSCAPE
####
plot(NA,xlab = "PC1", ylab = "PC2", xlim = range(R.pca$x[,1]), ylim = range(R.pca$x[,2]))
points(R.pca$x[,1:2],cex = rescale.numeric(R.size,c(0.75,2.25)))
points(R.pca$x[R2.df.n,1:2],cex = rescale.numeric(R.size,c(0.75,2.25))[R2.df.n],
       col = make.color.vec(100-R2.df$percentgrass,11,brewer_pal(palette = "BrBG")), pch = 19)

pgrass.surf <- surf.ls(2,R.pca$x[R2.df.n,1],R.pca$x[R2.df.n,2],100-R2.df$percentgrass)
pgrass.trmat <- trmat(pgrass.surf,-.3,.3,-.3,.3,n = 24)
anova(pgrass.surf)

#plot(NA,xlab = "PC1", ylab = "PC2", xlim = range(R.pca$x[,1]), ylim = range(R.pca$x[,2]), asp = T, cex.lab = 1.5)
plot(NA,xlab = "", ylab = "", xlim = range(R.pca$x[,1]), ylim = range(R.pca$x[,2]), asp = 1, cex.lab = 1.5, axes = F)
grid(lty=1) ;
box()
.filled.contour(x=pgrass.trmat$x,y=pgrass.trmat$y,z = pgrass.trmat$z,col = brewer_pal(palette = "BrBG")(11),levels = seq(0,100,length.out = 11))
#contour(x=pgrass.trmat$x,y=pgrass.trmat$y,z = pgrass.trmat$z,col = brewer_pal(palette = "BrBG")(11),levels = seq(100,0,length.out = 11),lwd = 5, labels = NULL, add = T)
#mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
#mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)
 addtree2(R.pca$x[,1:2],tree = R.df$tree, col = "black", lwd = 1)
points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 3, lwd = 2)
#points(R.pca$x[R2.df.n,1:2],cex = rescale.numeric(R.size,c(0.75,2.25))[R2.df.n],
 #      col = make.color.vec(100-R2.df$percentgrass,11,brewer_pal(palette = "BrBG")), pch = 19)
shape::colorlegend(col = brewer_pal(palette = "BrBG",direction=-1)(11), zlim = c(0,100), zlevels = seq(0,100,length.out = 11),
            posx = c(0.1,0.3), posy = c(0.7,0.9), main = "", lab.col = "transparent")

#########################

# FAMILY MORPHOSPACES
# figure 2
#########################


plot.contour <- FALSE
f2.CL <- 1.25

# cervid phylomorphospace

# predicting PC scores of minimum and maximum shapes in cervidae
allo.Cervidae <- procD.pgls(f1 = cervid.lm~cervid.size, phy = cervid.tree, data = NULL); anova(allo.Cervidae)
#             Df       SS        MS     Rsq      F      Z Pr(>F)
# cervid.size  1 0.016565 0.0165655 0.35361 17.506 4.2427  0.001 **
# Residuals   32 0.030281 0.0009463 0.64639
# Total       33 0.046847
#write.csv(summary(allo.Cervidae)$table,file="results/CervidAlloPGLS.csv")
allo.Cervidae1 <- procD.lm(f1 = cervid.lm~cervid.size, data = NULL); anova(allo.Cervidae1)


# without Alces, just to see
allo.Cervidae2 <- procD.pgls(f1 = cervid.lm[,,-c(31,32)]~cervid.size[-c(31,32)], phy = keep.tip(cervid.tree,cervid.tree$tip.label[-c(31,32)]), data = NULL); anova(allo.Cervidae2)


cervid.fit <- mvBM(tree = cervid.tree, data = cervid.pca$x[,1:2], model = "BM1")
cervid.dist <- as.numeric(distRoot(cervid.tree)[1])
buf <- 1.33 ; buf2 <- buf*1.5

cervid.x.points <- seq((range(cervid.pca$x[,1])*buf2)[1],(range(cervid.pca$x[,1])*buf2)[2],length.out=100)
cervid.y.points <- seq((range(cervid.pca$x[,2])*buf2)[1],(range(cervid.pca$x[,2])*buf2)[2],length.out=100)
cervid.z <- matrix(0,nrow=100,ncol=100)
cervid.mu <-  as.numeric(cervid.fit$theta)

for (i in 1:100) {
for (j in 1:100) {
cervid.z[i,j] <- dmvnorm(c(cervid.x.points[i],cervid.y.points[j]), mu=cervid.mu,sigma=cervid.fit$sigma * cervid.dist)
}
}
#contour(cervid.x.points,cervid.y.points,cervid.z )
#plot(expand.grid(cervid.x.points,cervid.y.points),cex = rescale.numeric(cervid.z ), pch=19)
NLEVEL <- 24
cervid.levels <- seq(0,max(range(cervid.z)), length.out = NLEVEL)


bovid.mean.scores <- matrix(nrow=dim(bovid.lm)[3],ncol=4)
for(i in 1:dim(bovid.lm)[3]){
  x <- rotonto(mshape(cervid.lm),bovid.lm[,,i],scale=FALSE)$yrot - mshape(cervid.lm)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  bovid.mean.scores[i,] <- a.preds[1,] %*% (cervid.pca$rotation)[,1:4]
}
rownames(bovid.mean.scores) <- dimnames(bovid.lm)[[3]]


cervidPC1.d <- round(cervid.pca$d[1]/sum(cervid.pca$d),3)*100
cervidPC2.d <- round(cervid.pca$d[2]/sum(cervid.pca$d),3)*100

cervid.allo.plot <- plot(allo.Cervidae, type = "regression", predictor = cervid.size,reg.type ="RegScore", pch = 21, lwd = 2,
                         bg = family.colors[1], cex = rescale.numeric(cervid.size,to = c(0.5,2)), xlab = "log-transformed centroid size")
cervid.preds <- shape.predictor(allo.Cervidae$GM$pgls.fitted, x= cervid.allo.plot$RegScore, Intercept = FALSE,
                        predmin = max(cervid.allo.plot$RegScore),
                        predmax = min(cervid.allo.plot$RegScore))
cervid.allo.scores2 <- matrix(data=NA,nrow=2,ncol=27)
x <- bindArr(rotonto(mshape(cervid.lm),cervid.preds$predmin,scale=TRUE)$yrot,rotonto(mshape(cervid.lm),cervid.preds$predmax,scale=TRUE)$yrot,along=3)
x <- sweep(x, 1:2, mshape(cervid.lm))
c.preds <- two.d.array(x)
cervid.allo.scores2 <- c.preds %*% (cervid.pca$rotation)[,1:27]


# plotting the phylomorphospace
plot(NA,xlim = range(cervid.pca$x[,1])*buf, ylim =range(cervid.pca$x[,2])*buf, asp = 1,
     xlab = paste("PC1: ",cervidPC1.d,"% of variation",sep=""),
     ylab = paste("PC2: ",cervidPC2.d,"% of variation",sep=""), cex.lab = 1.25)
title("Cervidae",adj=0,line=0.1,cex.main=4,col.main=family.colors[5],font=2)
grid(lty=1,col=addTrans('lightgray',255*0.5))
plot.contour = FALSE
if(plot.contour){.filled.contour(cervid.x.points,cervid.y.points,cervid.z,
          levels = cervid.levels, col = addTrans(colorRampPalette(c('white',family.colors[3],family.colors[1],family.colors[5]))(NLEVEL), 255*0.33))}
cervid.ellipse95 <- mixtools::ellipse(mu = as.numeric(cervid.fit$theta), sigma = cervid.fit$sigma*cervid.dist,
        alpha = 0.05, col = family.colors[1], lwd = 3)
cervid.ellipse99 <- mixtools::ellipse(mu = as.numeric(cervid.fit$theta), sigma = cervid.fit$sigma*cervid.dist,
        alpha = 0.01, col = family.colors[1], lwd = 3, lty = 2)
addtree2(bovid.mean.scores[,1:2],bovid.tree,"darkgray",1, lwd = 1)
addtree2(cervid.pca$x[,1:2],cervid.tree,"black",1, lwd = 1)
points(x=bovid.mean.scores[,1:2],pch= 15, col = addTrans(family.colors[2],255*0.5), cex = rescale.numeric(bovid.size,c(1,2.5)))
points(x = cervid.pca$x[,1:2], pch=21,bg=family.colors[1],cex = rescale.numeric(cervid.size,c(1,2.5)))
arrows(cervid.allo.scores2[1,1],cervid.allo.scores2[1,2],cervid.allo.scores2[2,1],cervid.allo.scores2[2,2],
        lwd = 5)
arrows(cervid.allo.scores2[1,1],cervid.allo.scores2[1,2],cervid.allo.scores2[2,1],cervid.allo.scores2[2,2],
       col=family.colors[1], lwd = 3)
points(cervid.fit$theta, bg = "goldenrod1", pch = 23, cex =2.5, lwd = 3, col = family.colors[1])
text(x = c(0.16,0.2025), y = c(0.05,0.05), labels = c('95%','99%'), cex = , col = family.colors[1])

cervid.pts95.1 <- matrix( c(0,.14,.165,-.105,0.0965,.09,0,-0.111), nrow = 4, ncol=2)
    points(x = cervid.pts95.1,lwd=3,pch=21,bg='white',cex=2)
cervid.pts95 <- cervid.pts95.1
cervid.pts95[1,1] <-  -0.05; cervid.pts95[2,1] <- 0.21
cervid.pts95[3,] <- c(.18,-.025) ; cervid.pts95[4,] <- c(-0.05,-0.0875)
scores <- cervid.pts95.1
shapes <- matrix(NA, ncol = dim(cervid.lm)[1]*dim(cervid.lm)[2], nrow = nrow(scores))
#PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(cervid.lm)))+cervid.pca$rotation[,1]*scores[i,1]+cervid.pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)

# to add the TPS figs
for(i in 1:NROW(scores)){
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(cervid.lm)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = .11, mag = 1,
      shade = T, shade.trans = 0.7,
      add = T, at = as.matrix(cervid.pts95)[i,],
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)
}
tps(target.lm = R[lat.drop,-3,"Alces_alces"],
      reference.lm = mshape(cervid.lm)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = .11, mag = 1,
      shade = T, shade.trans = 0.7,
      add = T, at = c(.2,-.12),
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)




###




# Bovidae phylomorphospace
allo.Bovidae <- procD.pgls(f1 = bovid.lm~bovid.size, phy = bovid.tree, data = NULL); anova(allo.Bovidae)$table
#write.csv(anova(allo.Bovidae)$table,file="results/BovidAlloPGLS.csv")
allo.Bovidae1 <- procD.lm(f1 = bovid.lm~bovid.size, data = NULL); anova(allo.Bovidae1)$table

bovid.fit <- mvBM(tree = bovid.tree, data = bovid.pca$x[,1:2], model = "BM1")
bovid.dist <- as.numeric(distRoot(bovid.tree)[1])
buf <- 1.2 ; buf2 <- buf*1.5

bovid.x.points <- seq((range(bovid.pca$x[,1])*buf2)[1],(range(bovid.pca$x[,1])*buf2)[2],length.out=100)
bovid.y.points <- seq((range(bovid.pca$x[,2])*buf2)[1],(range(bovid.pca$x[,2])*buf2)[2],length.out=100)
bovid.z <- matrix(0,nrow=100,ncol=100)
bovid.mu <-  as.numeric(bovid.fit$theta)

for (i in 1:100) {
for (j in 1:100) {
bovid.z[i,j] <- dmvnorm(c(bovid.x.points[i],bovid.y.points[j]), mu=bovid.mu,sigma=bovid.fit$sigma * bovid.dist)
}
}
#contour(bovid.x.points,bovid.y.points,bovid.z )
#plot(expand.grid(bovid.x.points,bovid.y.points),cex = rescale.numeric(bovid.z ), pch=19)

bovid.levels <- seq(0,max(range(bovid.z)), length.out = NLEVEL)










cervid.mean.scores <- matrix(nrow=dim(cervid.lm)[3],ncol=4)
for(i in 1:dim(cervid.lm)[3]){
  x <- rotonto(mshape(bovid.lm),cervid.lm[,,i],scale=FALSE)$yrot - mshape(bovid.lm)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  cervid.mean.scores[i,] <- a.preds[1,] %*% (bovid.pca$rotation)[,1:4]
}
rownames(cervid.mean.scores) <- dimnames(cervid.lm)[[3]]

bovidPC1.d <- round(bovid.pca$d[1]/sum(bovid.pca$d),3)*100
bovidPC2.d <- round(bovid.pca$d[2]/sum(bovid.pca$d),3)*100

bovid.allo.plot <- plot(allo.Bovidae, type = "regression", predictor = bovid.size ,reg.type ="RegScore", xlab = "log-transformed centroid size", pch =21, lwd = 2, bg = family.colors[2], cex = rescale.numeric(bovid.size,to=c(0.5,2)))
bovid.preds <- shape.predictor(allo.Bovidae$GM$pgls.fitted, x= bovid.allo.plot$RegScore, Intercept = FALSE,
                        predmin = min(bovid.allo.plot$RegScore),
                        predmax = max(bovid.allo.plot$RegScore))

bovid.allo.scores2 <- matrix(data=NA,nrow=2,ncol=27)
x <- bindArr(rotonto(mshape(bovid.lm),bovid.preds$predmin,scale=TRUE)$yrot,rotonto(mshape(bovid.lm),bovid.preds$predmax,scale=TRUE)$yrot,along=3)
x <- sweep(x, 1:2, mshape(bovid.lm))
c.preds <- two.d.array(x)
bovid.allo.scores2 <- c.preds %*% (bovid.pca$rotation)[,1:27]




plot(NA,xlim = range(bovid.pca$x[,1])*buf, ylim =range(bovid.pca$x[,2])*buf, asp = 1,
     xlab = paste("PC1: ",bovidPC1.d,"% of variation",sep=""),
     ylab = paste("PC2: ",bovidPC2.d,"% of variation",sep=""), cex.lab = f2.CL)
grid(lty=1,col=addTrans('lightgray',255*0.5))
title("Bovidae",adj=0,line=0.1,cex.main=4,col.main=family.colors[6],font=2)
if(plot.contour){.filled.contour(bovid.x.points,bovid.y.points,bovid.z,
          levels = bovid.levels, col = addTrans(colorRampPalette(c('white',family.colors[4],family.colors[2],family.colors[6]))(NLEVEL), 255*0.33))}
bovid.ellipse95 <- mixtools::ellipse(mu = as.numeric(bovid.fit$theta), sigma = bovid.fit$sigma*bovid.dist,
        alpha = 0.05, col = family.colors[2], lwd = 3)
bovid.ellipse99 <- mixtools::ellipse(mu = as.numeric(bovid.fit$theta), sigma = bovid.fit$sigma*bovid.dist,
        alpha = 0.01, col = family.colors[2], lwd = 3, lty = 2)
addtree2(cervid.mean.scores[,1:2],cervid.tree,"darkgray",1, lwd = 1)
addtree2(bovid.pca$x[,1:2],bovid.tree,"black",1, lwd = 1)

points(x=cervid.mean.scores[,1:2],pch= 16, col = addTrans(family.colors[1],255*0.5), cex = rescale.numeric(cervid.size,c(1,2.5)))
points(x = bovid.pca$x[,1:2], pch=22,bg=family.colors[2],cex = rescale.numeric(bovid.size,c(1,2.5)))
arrows(bovid.allo.scores2[1,1],bovid.allo.scores2[1,2],bovid.allo.scores2[2,1],bovid.allo.scores2[2,2],
        lwd = 5)
arrows(bovid.allo.scores2[1,1],bovid.allo.scores2[1,2],bovid.allo.scores2[2,1],bovid.allo.scores2[2,2],
       col=family.colors[2], lwd = 3)
points(bovid.fit$theta, bg = "goldenrod1", pch = 23, cex =2.5, lwd = 3, col = family.colors[2])

text(x = c(0.11,0.15), y = c(0.05,0.05), labels = c('95%','99%'), cex = f2.CL, col = family.colors[2])

bovid.pts95.1 <- rbind(matrix( c(0,0,0.1185,-.0975), nrow = 2, ncol = 2), matrix( c(0.137,-0.148,0,0), nrow = 2, ncol = 2))
  points(x = bovid.pts95.1,lwd=3,pch=21,bg='white',cex=2)


bovid.pts95 <- bovid.pts95.1
bovid.pts95[1,2] <- 0.14 ; bovid.pts95[2,] <- c(0.03,-.13)
bovid.pts95[3,1] <- .2 ; bovid.pts95[4,] <- c(-0.155,0.025)
scores <- bovid.pts95.1
shapes <- matrix(NA, ncol = dim(bovid.lm)[1]*dim(bovid.lm)[2], nrow = nrow(scores))
PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(bovid.lm)))+bovid.pca$rotation[,1]*scores[i,1]+bovid.pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(bovid.lm)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.135, mag = 1,
      shade = T, shade.trans = 0.7,
      add = T, at = as.matrix(bovid.pts95)[i,],
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)
}





# ruminant populations
####

# specimen scores

# all cervid specimens
cervid.spec.n <- name.id[name.id[,3] %in% names(f.fam)[f.fam=="Cervidae"],1] ; names(cervid.spec.n) <- name.id[name.id[,3] %in% names(f.fam)[f.fam=="Cervidae"],3]

cervid.specs <- specimenInfo$species.FV %in% name.id[name.id[,3] %in% names(f.fam)[f.fam=="Cervidae"],1]
cervid.specs2 <- specimenInfo$species.FV %in% c("Odo_he_he1","Odo_he_ca2","Odo_he_co4","Odo_vi_bo6","Odo_vi_le3","Odo_vi_co7")
cervid.spec <- ruminant.lm[1:25,,cervid.specs|cervid.specs2]
cervid.size <- log(apply(cervid.spec,3,cSize))

bovid.spec.n <- name.id[name.id[,3] %in% names(f.fam)[f.fam=="Bovidae"],1] ; names(bovid.spec.n) <- name.id[name.id[,3] %in% names(f.fam)[f.fam=="Bovidae"],3]

bovid.specs <- specimenInfo$species.FV %in% name.id[name.id[,3] %in% names(f.fam)[f.fam=="Bovidae"],1]
bovid.spec <- ruminant.lm[1:25,,bovid.specs]
bovid.size <- log(apply(bovid.spec,3,cSize))


cervid.spec.n2 <- name.id[name.id[,3] %in% names(fam)[fam=="Cervidae"],1] ; names(cervid.spec.n2) <- name.id[name.id[,3] %in% names(fam)[fam=="Cervidae"],3]
bovid.spec.n2 <- name.id[name.id[,3] %in% names(fam)[fam=="Bovidae"],1] ; names(bovid.spec.n2) <- name.id[name.id[,3] %in% names(fam)[fam=="Bovidae"],3]
# loop that calculates specimen PC scores for each specimen (i) for each species (j)

ALSPSIZE <- rescale.numeric(log(allspecs.size),c(0.5,2.5))

cervidspecs.scores <- c()
cervid.SIZES <- c()
for(j in 1:length(cervid.spec.n)){
  SPECS <- allspecs[,,which(dimnames(allspecs)[[3]] == cervid.spec.n[[j]])]
  cervid.SIZES[[j]] <- ALSPSIZE[which(dimnames(allspecs)[[3]] == cervid.spec.n[[j]])]
  SPECSscores <- matrix(nrow=dim(SPECS)[3],ncol=dim(cervid.pca$rotation)[2])

for(i in 1:dim(SPECS)[3]){
  x <- rotonto(mshape(cervid.lm),SPECS[,,i],scale=FALSE, reflection = T)$yrot - mshape(cervid.lm)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  SPECSscores[i,] <- a.preds[1,] %*% cervid.pca$rotation
}
  cervidspecs.scores[[j]] <- SPECSscores
  names(cervidspecs.scores)[[j]] <- names(cervid.spec.n)[[j]]
}

bovidspecs.scores <- c()
bovid.SIZES <- c()
for(j in 1:length(bovid.spec.n)){
  SPECS <- allspecs[,,which(dimnames(allspecs)[[3]] == bovid.spec.n[[j]])]
  SPECSscores <- matrix(nrow=dim(SPECS)[3],ncol=dim(bovid.pca$rotation)[2])
  bovid.SIZES[[j]] <- ALSPSIZE[which(dimnames(allspecs)[[3]] == bovid.spec.n[[j]])]

for(i in 1:dim(SPECS)[3]){
  x <- rotonto(mshape(bovid.lm),SPECS[,,i],scale=F)$yrot - mshape(bovid.lm)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  SPECSscores[i,] <- a.preds[1,] %*% bovid.pca$rotation
}
  bovidspecs.scores[[j]] <- SPECSscores
  names(bovidspecs.scores)[[j]] <- names(bovid.spec.n)[[j]]
}


Vproj.pop.col <- make.color.vec(c(0,Vproj),n=100, col.pal.sunset)
Vproj.pop.col <- Vproj.pop.col[-1] # scaled so that Vproj=0 is dark blue

#cervid.pop.col <- na.omit(Vproj.pop.col[cervid.tree$tip.label])
#bovid.pop.col <- na.omit(Vproj.pop.col[bovid.tree$tip.label])

cervid.pop.col <- make.color.vec(c(0,na.omit(Vproj[cervid.tree$tip.label])), 100, colorRampPalette(c('white',family.colors[1])))
cervid.pop.col <- cervid.pop.col[-1]
bovid.pop.col <- make.color.vec(c(0,na.omit(Vproj[bovid.tree$tip.label])), 100, colorRampPalette(c('white',family.colors[2])))
bovid.pop.col <- bovid.pop.col[-1]

cervid.specs.size <- log(allspecs.size[which(names(allspecs.size) %in% cervid.spec.n)])
bovid.specs.size <- log(allspecs.size[which(names(allspecs.size) %in% bovid.spec.n)])


# the plots:
plot(NA,xlim = range(cervid.pca$x[,1])*buf, ylim =range(cervid.pca$x[,2])*buf, asp = T,
     xlab = paste("PC1: ",cervidPC1.d,"% of variation",sep=""),
     ylab = paste("PC2: ",cervidPC2.d,"% of variation",sep=""), cex.lab = f2.CL)
grid(lty=1,col=addTrans('lightgray',255*0.5))
addtree2(x = cervid.pca$x[,1:2],tree = cervid.tree, col = "lightgray")

for(i in 1:length(cervidspecs.scores)){
  xy <- cervidspecs.scores[[i]][,1:2]
  points(xy, pch = 19, col = addTrans(family.colors[1],255*0.1), cex=cervid.SIZES[[i]])
  mixtools::ellipse(mu = colMeans(xy), sigma = cov(xy), alpha = 0.2, col = 'black', lwd = 4, lty = 1)
  mixtools::ellipse(mu = colMeans(xy), sigma = cov(xy), alpha = 0.2, col = cervid.pop.col[i], lwd = 2, lty = 1)
}
arrows(cervid.allo.scores2[1,1],cervid.allo.scores2[1,2],cervid.allo.scores2[2,1],cervid.allo.scores2[2,2],
        lwd = 5)
arrows(cervid.allo.scores2[1,1],cervid.allo.scores2[1,2],cervid.allo.scores2[2,1],cervid.allo.scores2[2,2],
       col=family.colors[1], lwd = 3)
points(cervid.fit$theta, bg = "goldenrod1", pch = 23, cex =2.5, lwd = 3, col = family.colors[1])


plot(NA,xlim = range(bovid.pca$x[,1])*buf, ylim =range(bovid.pca$x[,2])*buf, asp = T,
     xlab = paste("PC1: ",bovidPC1.d,"% of variation",sep=""),
     ylab = paste("PC2: ",bovidPC2.d,"% of variation",sep=""), cex.lab = f2.CL)
grid(lty=1,col=addTrans('lightgray',255*0.5))
addtree2(x = bovid.pca$x[,1:2],tree = bovid.tree, col = "lightgray")

for(i in 1:length(bovidspecs.scores)){
  xy <- bovidspecs.scores[[i]][,1:2]
  points(xy, pch = 15, col = addTrans(family.colors[2],255*0.1), cex = bovid.SIZES[[i]])
  mixtools::ellipse(mu = colMeans(xy), sigma = cov(xy), alpha = 0.2, col = 'black', lwd = 4, lty = 1)
  mixtools::ellipse(mu = colMeans(xy), sigma = cov(xy), alpha = 0.2, col = bovid.pop.col[i], lwd = 2, lty = 1)
}
arrows(bovid.allo.scores2[1,1],bovid.allo.scores2[1,2],bovid.allo.scores2[2,1],bovid.allo.scores2[2,2],
        lwd = 5)
arrows(bovid.allo.scores2[1,1],bovid.allo.scores2[1,2],bovid.allo.scores2[2,1],bovid.allo.scores2[2,2],
       col=family.colors[2], lwd = 3)
points(bovid.fit$theta, bg = "goldenrod1", pch = 23, cex =2.5, lwd = 3, col = family.colors[2])



# contours of face length and nasal retraction
cervid.grid <- pca.grid(cervid.pca, PC1.segments = 10, PC2.segments = 10)

#cervid.grid.x.points <- unlist(lapply(cervid.grid$score.grid,FUN = function(x)x[1]))
#cervid.grid.y.points <- unlist(lapply(cervid.grid$score.grid,FUN = function(x)x[2]))
cervid.z.nasalretraction <- apply(cervid.grid$shape.grid,3,nasal.retraction)
cervid.z.facelength <- apply(cervid.grid$shape.grid,3,face.length)

CNRsurf <- surf.ls(np = 3,cervid.grid$score.grid[,1], cervid.grid$score.grid[,2],cervid.z.nasalretraction)
CNRmat <- trmat(CNRsurf,range(cervid.grid$score.grid[,1])[1],range(cervid.grid$score.grid[,1])[2],
                range(cervid.grid$score.grid[,2])[1],range(cervid.grid$score.grid[,2])[2], n = 100)

# plot(NA, xlim = range(cervid.grid$score.grid[,1])*0.875,
#         ylim = range(cervid.grid$score.grid[,2])*0.875, xlab = "", ylab = "", axes = F, frame = T,
#      main = "nasal retraction", cex.main = 2.5)
#
# .filled.contour(x=CNRmat$x,y=CNRmat$y,z=CNRmat$z, levels = seq(min(CNRmat$z), max(CNRmat$z), length.out = 24),
#                 col = colorRampPalette(c('white', family.colors[1]))(24)) ; box()
# contour(CNRmat, add=T, nlevels = 23, drawlabels =  F)
contour(CNRmat, xlim = range(cervid.grid$score.grid[,1]),
        ylim = range(cervid.grid$score.grid[,2]),
        zlim = range(cervid.z.nasalretraction), add = F, nlevels = 40, drawlabels = F,
        col = colorRampPalette(c('white',family.colors[1],family.colors[5]))(40), lwd = 4, axes = F, frame = T,
        main= "nasal retraction", font.main = 1, cex.main = 2)


CFLsurf <- surf.ls(np = 3,cervid.grid$score.grid[,1], cervid.grid$score.grid[,2],cervid.z.facelength)
CFLmat <- trmat(CFLsurf,range(cervid.grid$score.grid[,1])[1],range(cervid.grid$score.grid[,1])[2],
                range(cervid.grid$score.grid[,2])[1],range(cervid.grid$score.grid[,2])[2], n = 100)


contour(CFLmat, xlim = range(cervid.grid$score.grid[,1]),
        ylim = range(cervid.grid$score.grid[,2]),
        zlim = range(cervid.z.facelength), add = F, nlevels = 40, drawlabels = F,
        col = colorRampPalette(c('white',family.colors[1],family.colors[5]))(40), lwd = 4, axes = F, frame = T,
        main= "face length", font.main = 1, cex.main = 2)


bovid.grid <- pca.grid(bovid.pca, PC1.segments = 10, PC2.segments = 10)

#cervid.grid.x.points <- unlist(lapply(cervid.grid$score.grid,FUN = function(x)x[1]))
#cervid.grid.y.points <- unlist(lapply(cervid.grid$score.grid,FUN = function(x)x[2]))
bovid.z.nasalretraction <- apply(bovid.grid$shape.grid,3,nasal.retraction)
bovid.z.facelength <- apply(bovid.grid$shape.grid,3,face.length)

BNRsurf <- surf.ls(np = 3,bovid.grid$score.grid[,1], bovid.grid$score.grid[,2],bovid.z.nasalretraction)
BNRmat <- trmat(BNRsurf,range(bovid.grid$score.grid[,1])[1],range(bovid.grid$score.grid[,1])[2],
                range(bovid.grid$score.grid[,2])[1],range(bovid.grid$score.grid[,2])[2], n = 100)

contour(BNRmat, xlim = range(bovid.grid$score.grid[,1]),
        ylim = range(bovid.grid$score.grid[,2]),
        zlim = range(bovid.z.nasalretraction), add = F, nlevels = 48, drawlabels = F,
        col = colorRampPalette(c('white',family.colors[2], family.colors[6]))(64), lwd = 4, axes = F, frame = T,
        main= "nasal retraction", font.main = 1, cex.main = 2)

BFLsurf <- surf.ls(np = 3,bovid.grid$score.grid[,1], bovid.grid$score.grid[,2],bovid.z.facelength)
BFLmat <- trmat(BFLsurf,range(bovid.grid$score.grid[,1])[1],range(bovid.grid$score.grid[,1])[2],
                range(bovid.grid$score.grid[,2])[1],range(bovid.grid$score.grid[,2])[2], n = 100)

contour(BFLmat, xlim = range(bovid.grid$score.grid[,1]),
        ylim = range(bovid.grid$score.grid[,2]),
        zlim = range(bovid.z.facelength), add = F, nlevels = 48, drawlabels = F,
        col = colorRampPalette(c('white',family.colors[2],family.colors[6]))(64), lwd = 4, axes = F, frame = T,
        main= "face length", font.main = 1, cex.main = 2)



#########################

# projected variances
# figure 3
############


# B - PROJECTED VARIANCE DISTRIBUTIONS

# testing against a null distribution of evolutionary divergences
nsim <- 999
Nmin <- 27 # minimum of 27 specimens
Vproj.boot <- c()
for(i in 1:nsim){
  x <- numeric(length(proc.list.20))
    for(j in 1:length(x)){
      l <- dim(proc.list.20[[j]]$coords)[3]
        vcv <- var(two.d.array(A=proc.list.20[[j]]$coords[,,sample(1:l, Nmin)]))
        random.Z <- scalar1(sample(seq(-1,1,length.out=200),size=75))
        x[j] <- projected.variance(vcv,random.Z)
    }
  Vproj.boot[[i]] <- x
}

# RANDOM DIRECTIONS DRAWN FROM THE
random.phylo.Z <- mvrnorm(Sigma = R.bm$sigma, n = nsim, mu = rep(0,2))


Vproj.boot.phylo <- c()
for(i in 1:nsim){
  x <- numeric(length(proc.list.20))
    for(j in 1:length(x)){
            l <- dim(proc.list.20[[j]]$coords)[3]
        vcv <- var(two.d.array(A=proc.list.20[[j]]$coords[,,sample(1:l, Nmin)]))
        random.Z <- flatten(M-restoreShapes(scores = random.phylo.Z[i,], PC = R.pca$rotation[,1:2], mshape = M))
        x[j] <- projected.variance(vcv,random.Z)
    }
Vproj.boot.phylo[[i]] <- x
}

# interspecific PCs
PC1 <- R.pca$rotation[,1] ; PC2 <- R.pca$rotation[,2]
pPC1 <- R.Ppca$rotation[,1] ; pPC2 <- R.Ppca$rotation[,2]
PaC1 <- R.paca$rotation[,1] ; PaC2 <- R.paca$rotation[,2]

Vproj.allo.distribution <- c()
var.dist <- c()
var.dist2 <- c()
pc1.d <- c()
pc2.d <- c()
PPC1.d <- c()
PPC2.d <- c()
paca1.d <- c()
paca2.d <- c()
#var.dist3 <- c()
#intra.allo.dist <- c()
for(i in 1:nsim){
  ad <- numeric(length(proc.list.20))
  PC1D <- numeric(length(ad))
  PC2D <- numeric(length(ad))
  pPC1D <- numeric(length(ad))
  pPC2D <- numeric(length(ad))
  PAC1D <- numeric(length(ad))
  PAC2D <- numeric(length(ad))
  v <- numeric(length(ad))
  v2 <- numeric(length(ad))
  v3 <- numeric(length(ad))

  for(j in 1:length(proc.list.20)){
      l <- dim(proc.list.20[[j]]$coords)[3]
      sam <- sample(1:l, Nmin, replace = T)
      A2d <- two.d.array(A=proc.list.20[[j]]$coords[,,sam])
      vcv <- var(A2d)
      svd1 <- svd(vcv)
      ##v[j] <- svd1$d[1]/sum(svd1$d)
      ##v2[j] <- svd1$d[2]/sum(svd1$d)
      ##v3[j] <- svd1$d[3]/sum(svd1$d)
      #v[j] <- projected.variance(vcv,svd1$u[,1], scale=T)
      #v2[j] <- projected.variance(vcv,svd1$u[,2], scale=T)
      ad[j] <- projected.variance(vcv,allo.vector, scale=T)
      pPC1D[j] <- projected.variance(vcv,pPC1, scale=T)
      pPC2D[j] <- projected.variance(vcv,pPC2, scale=T)
      PAC1D[j] <- projected.variance(vcv,PaC1, scale=T)
      PAC2D[j] <- projected.variance(vcv,PaC2, scale=T)
      PC1D[j] <- projected.variance(vcv,PC1, scale=T)
      PC2D[j] <- projected.variance(vcv,PC2, scale=T)
  }

  Vproj.allo.distribution[[i]] <- ad
  var.dist[[i]] <- v
  var.dist2[[i]] <- v2
  pc1.d[[i]] <- PC1D
  pc2.d[[i]] <- PC2D
  PPC1.d[[i]] <- pPC1D
  PPC2.d[[i]] <- pPC2D
  paca1.d[[i]] <- PAC1D
  paca2.d[[i]] <- PAC2D
}


PC.cols <- c("dodgerblue4", "dodgerblue",
             "dodgerblue4", "dodgerblue",
             "dodgerblue4", "dodgerblue",
             "black", "firebrick2")

lwd.ins <- 4; lwd.outs <- 6
adjust <- 1
par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,0.2),ylim=c(0,100), xlab = "bootstrapped projected variances", ylab="density")
#for(i in 1:nsim){lines(density(Vproj.boot[[i]], adjust = adjust), lwd = 1, col = addTrans(PC.cols[7], trans=255*0.05))}
#for(i in 1:nsim){lines(density(Vproj.allo.distribution[[i]], adjust = adjust), lwd = 1, col = addTrans(PC.cols[8], trans=255*0.05))}
#for(i in 1:nsim){lines(density(pc1.d[[i]], adjust = adjust), lwd = 1, col = addTrans(PC.cols[1], trans=255*0.05))}
lines(density(unlist(pc1.d), adjust = adjust), lwd = lwd.ins, col = addTrans(PC.cols[1], trans=255*1), lty=1)
#for(i in 1:nsim){lines(density(pc2.d[[i]], adjust = adjust), lwd = 1, col = addTrans(PC.cols[2], trans=255*0.05))}
lines(density(unlist(pc2.d), adjust = adjust), lwd = lwd.ins, col = addTrans(PC.cols[2], trans=255*1), lty=1)
#for(i in 1:nsim){lines(density(PPC1.d[[i]], adjust = adjust), lwd = 1, col = addTrans(PC.cols[3], trans=255*0.05))}
lines(density(unlist(PPC1.d), adjust = adjust), lwd = lwd.ins, col = addTrans(PC.cols[3], trans=255*1), lty=5)
#for(i in 1:nsim){lines(density(PPC2.d[[i]], adjust = adjust), lwd = 1, col = addTrans(PC.cols[4], trans=255*0.05))}
lines(density(unlist(PPC2.d), adjust = adjust), lwd = lwd.ins, col = addTrans(PC.cols[4], trans=255*1), lty=5)
#for(i in 1:nsim){lines(density(paca1.d[[i]], adjust = adjust), lwd = 1, col = addTrans(PC.cols[5], trans=255*0.05))}
lines(density(unlist(paca1.d), adjust = adjust), lwd = lwd.ins, col = addTrans(PC.cols[5], trans=255*1), lty=3)
#for(i in 1:nsim){lines(density(paca2.d[[i]], adjust = adjust), lwd = 1, col = addTrans(PC.cols[6], trans=255*0.05))}
lines(density(unlist(paca2.d), adjust = adjust), lwd = lwd.ins, col = addTrans(PC.cols[6], trans=255*1), lty=3)

lines(density(unlist(Vproj.boot), adjust = adjust), lwd = lwd.outs, col = addTrans(PC.cols[7], trans=255*1))
#lines(density(unlist(Vproj.boot.phylo), adjust = adjust), lwd = lwd.outs, col = "gray")
lines(density(unlist(Vproj.allo.distribution), adjust = adjust), lwd = lwd.outs, col = addTrans(PC.cols[8], trans=255*1))

#legend("topright", legend = c("random directions", "directions pulled from W", "evolutionary allometry (CREA)",
  #                            "PC1", "PC2", "pPC1", "pPC2", "PaC1", "PaC2"),
#              col = c(PC.cols[7],'gray',PC.cols[8], PC.cols[1:6]), lty = c(1,1,5,5,3,3,1,1), lwd = c(3,3,2,2,3,3,3,3))
legend("topright", legend = c("random directions",  "CREA",
                              "PC1", "PC2", "pPC1", "pPC2", "PaC1", "PaC2"),
              col = c(PC.cols[7],PC.cols[8], PC.cols[1:6]), lty = c(1,1,5,5,3,3,1,1), lwd = c(3,3,2,2,3,3,3,3))


#############


# integration ~ disparity regressions
# figure 4, supplemental figure 3
#############
intra.N <- unlist(lapply(proc.list.20,function(x)dim(x$coords)[3]))
names(intra.N) <- tree.20$tip.label

intra.df <- data.frame(disparity1 = dis, disparity2 = maxD, integration1 = evd, integration2 = rdi.R, Vproj = Vproj, row.names = tree.20$tip.label)
intra.df2 <- intra.df ; colnames(intra.df2) <- c("A","B","C","D","E")

intra.tree <- ladderize(tree.20)

intra.tree <- paintSubTree(intra.tree, node = 50, state = "Ruminantia")
intra.tree <- paintSubTree(intra.tree, node = 51, state = "Tragulidae")
intra.tree <- paintSubTree(intra.tree, node = 53, state = "Cervidae")
intra.tree <- paintSubTree(intra.tree, node = 65, state = "Bovidae")
intra.tree.cols <- setNames(c("black",family.colors[c(3,1,2)]),c("Ruminantia","Tragulidae","Cervidae","Bovidae"))

plotSimmap(intra.tree,colors=intra.tree.cols,lwd = 4)

phylo.heatmap(intra.tree,intra.df2,colors = col.pal.sunset(100), standardize = T)

panel.cor <- function(x, y, ncomp = 10){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    P <- anova(gls(x~y, correlation = corPagel(0,phy=tree.20)))$`p-value`[2]
    #P <- anova(lm(x~y))$`Pr(>F)`[1]
    P <- round(P, digits=3)
    corr <- cor(x,y)
    txt <- paste0("p = ", P,"\nr = ",round(corr, digits=3))
    if(P<(0.05/ncomp)){ft<-2 ; cex <- 2}else{if(P<0.05){ft<-3}else{ft<-1} ; cex <- 1.5}
    text(0.5, 0.5, txt, cex = cex, font = ft)
}
TREEE <- keep.tip(phy = tree.20, tip = rownames(intra.df)[which(intra.N > 39)])
panel.cor2 <- function(x, y, ncomp = 10){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    P <- anova(gls(x~y, correlation = corPagel(0,phy=TREEE)))$`p-value`[2]
    #P <- anova(lm(x~y))$`Pr(>F)`[1]
    P <- round(P, digits=3)
    corr <- cor(x,y)
    txt <- paste0("p = ", P,"\nr = ",round(corr, digits=3))
    if(P<(0.05/ncomp)){ft<-2 ; cex <- 2}else{if(P<0.05){ft<-3}else{ft<-1} ; cex <- 1.5}
    text(0.5, 0.5, txt, cex = cex, font = ft)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch=pts.pch[tree.20$tip.label],
     bg=pts.col[tree.20$tip.label],)
  abline(lm(y~x),col="lightgray")
}
upper.panel2<-function(x, y){
  points(x,y, pch=pts.pch[TREEE$tip.label],
     bg=pts.col[TREEE$tip.label],)
  abline(lm(y~x),col="lightgray")
}

panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    D1 <- density(x) ; D1$y <- D1$y/max(D1$y)
    lines(D1,lwd=3,col="lightgray")
    fam.porp <- rev(table(f.fam))
    uf <- unique(f.fam)
    for(i in 2:3){
    D <- density(x[which(f.fam == uf[i])])
    D$y <- D$y / max(D$y)
    D$y <- D$y * (fam.porp[i] / sum(fam.porp))
    lines(D, lwd = 3, col = family.colors[c(3,1,2)][i])
    }
}
panel.hist2 <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    D1 <- density(x) ; D1$y <- D1$y/max(D1$y)
    lines(D1,lwd=3,col="lightgray")
    FFAM <-  f.fam[rownames(intra.df)[which(intra.N > 39)]]
    fam.porp <- rev(table(FFAM))
    uf <- unique(FFAM)
    for(i in 2:3){
    D <- density(x[which(FFAM == uf[i])])
    D$y <- D$y / max(D$y)
    D$y <- D$y * (fam.porp[i] / sum(fam.porp))
    lines(D, lwd = 3, col = family.colors[c(3,1,2)][i])
    }
}

# FIGURE 4:
pairs(intra.df,
      lower.panel = panel.cor,
      upper.panel = upper.panel,diag.panel = panel.hist,
      labels = c("A - disparity 1", "B - disparity 2","C - integration 1", "D - integration 2", "E - Vproj"))


# only including the specimens with at least 40 individuals
pairs(intra.df[which(intra.N > 39),],
      lower.panel = panel.cor2,
      upper.panel = upper.panel2,
      diag.panel = panel.hist2,
      labels = c("A - disparity 1", "B - disparity 2","C - integration 1", "D - integration 2", "E - Vproj"))



# with allometry removed:
par(mar=c(5,5,2,2),mfrow=c(1,1))

proc.list.MF <- c()
dimorphism.procD <- c()
dimorphism.allo.procD <- c()
for(i in 1:length(proc.list.20)){
proc.list.MF[[i]] <- specimenInfo$mf[proc.list.20.ID[[i]]] ; MF <- proc.list.MF[[i]]

procfoo <- proc.list.20[[i]]
A <- procfoo$coords ; dimnames(A)[[3]] <- proc.list.20.ID[[i]]
dimorphism.allo.procD[[i]] <- procD.lm(A~log(procfoo$Csize)*MF)
dimorphism.procD[[i]] <- procD.lm(A~MF)
}


dimorphism.Z <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$Z[2])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.Z, outline = FALSE, lwd = 6)

dimorphism.frat <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$F[2])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.frat, outline = FALSE, lwd = 6)

dimorphism.pval <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$P[2])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.pval, outline = FALSE, lwd = 6)


allo.Z <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$Z[1])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.Z, outline = FALSE, lwd = 6)

allo.frat <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$F[1])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.frat, outline = FALSE, lwd = 6)

allo.pval <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$P[1])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.pval, outline = FALSE, lwd = 6)


dimorphism.allo.Z <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$Z[3])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.Z, outline = FALSE, lwd = 6)

dimorphism.allo.frat <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$F[3])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.frat, outline = FALSE, lwd = 6)

dimorphism.allo.pval <- setNames(unlist(lapply(dimorphism.allo.procD, function(x)x$aov.table$P[3])), tree.20$tip.label)
#contMap(tree = ladderize(tree.20),  x = dimorphism.pval, outline = FALSE, lwd = 6)

# the effect sizes
zmat <- matrix(nrow = length(allo.Z), ncol = 3)
zmat[,1] <- allo.Z ; zmat[,2] <- dimorphism.Z ; zmat[,3] <- dimorphism.allo.Z
rownames(zmat) <- names(allo.Z)
colnames(zmat) <- c("static allometry", "sexual dimorphism", 'allometry/dimorphism interaction')
phylo.heatmap(tree = ladderize(tree.20), X = zmat, colors = col.pal.sunset(100), standardize = F)

# figure of the P-values
pmatstat <- matrix(nrow=length(allo.Z),ncol=3)
pmatstat[,1] <- allo.pval ; pmatstat[,2] <- dimorphism.pval ; pmatstat[,3] <- dimorphism.allo.pval
pmatstat2 <- pmatstat
for(i in 1:dim(pmatstat)[1]){
  for(j in 1:3){
    x <- pmatstat2[i,j]
    if(x > 0.05){pmatstat[i,j] <- 1}else{
    if(0.05 > x && x > 0.0011){pmatstat[i,j] <- 0.5}else{
    if(x < 0.0011){pmatstat[i,j] <- 0}}}
  }
}
rownames(pmatstat) <- names(allo.Z)
colnames(pmatstat) <- c("static allometry", "sexual dimorphism", 'allometry/dimorphism interaction')

# P-values lower than 0.05 are gray, p-values below 0.001 are black
phylo.heatmap(tree = ladderize(tree.20), X = pmatstat, colors = colorRampPalette(c('gray1','gray99'))(25), legend=FALSE)


# the f-ratios
fmat <- matrix(nrow = length(allo.Z), ncol = 3)
fmat[,1] <- allo.frat ; fmat[,2] <- dimorphism.frat ; fmat[,3] <- dimorphism.allo.frat
rownames(fmat) <- names(allo.Z)
colnames(fmat) <- c("static allometry", "sexual dimorphism", 'allometry/dimorphism interaction')
#phylo.heatmap(tree = ladderize(tree.20), X = fmat, colors = col.pal.sunset(100), standardize = F)
# log-transformed \/ \/ \/ \/
phylo.heatmap(tree = ladderize(tree.20), X = log(fmat), colors = col.pal.sunset(100), standardize = F)


par(mar=c(5,5,2,2),mfrow=c(1,1))
# ENCOURAGINGLY, AND EXPECTEDLY, PROJECTED VARIANCE POSITIVELY CORRELATES WITH STATIC ALLOMETRY EFFECT SIZE
# this means that when there is strong static allometry, its in the direction of CREA
plot(Vproj, allo.Z)
addtree2(cbind(Vproj, allo.Z), tree = tree.20, col = "gray")
points(Vproj, allo.Z, pch = pts.pch[tree.20$tip.label], bg = pts.col[tree.20$tip.label],
     cex = rescale.numeric(R.size[tree.20$tip.label], c(0.5,2)))

par(mar=c(5,5,2,2),mfrow=c(1,1))




allo.free <- c()
for(i in 1:length(proc.list.20)){
  M <- proc.list.20[[i]]$consensus
  N <- dim(proc.list.20[[i]]$coords)[3]
  X <- matrix(NA, nrow = N, ncol = ncol(dimorphism.allo.procD[[i]]$residuals))

  for(j in 1:N){
  X[j,] <-  flatten(M) + dimorphism.allo.procD[[i]]$residuals[j,]}

  allo.free[[i]] <- arrayspecs(X,25,3)
}

#for(i in 1:dim(allo.free[[1]])[3]){spheres3d(allo.free[[1]][,,i], radius = 0.001)}
# okay coolio, just making sure

# B - PROJECTED VARIANCE DISTRIBUTIONS

# testing against a null distribution of evolutionary divergences
nsim <- 999
Nmin <- 27 # minimum of 27 specimens
Vproj.boot <- c()
for(i in 1:nsim){
  x <- numeric(length(allo.free))
    for(j in 1:length(x)){
      l <- dim(allo.free[[j]])[3]
        vcv <- var(two.d.array(A=allo.free[[j]][,,sample(1:l, Nmin)]))
        random.Z <- scalar1(sample(seq(-1,1,length.out=200),size=75))
        x[j] <- projected.variance(vcv,random.Z)
    }
  Vproj.boot[[i]] <- x
}

# interspecific PCs
PC1 <- R.pca$rotation[,1] ; PC2 <- R.pca$rotation[,2]
pPC1 <- R.Ppca$rotation[,1] ; pPC2 <- R.Ppca$rotation[,2]
PaC1 <- R.paca$rotation[,1] ; PaC2 <- R.paca$rotation[,2]

Vproj.allo.distribution2 <- c()
pc1.d2 <- c()
pc2.d2 <- c()
PPC1.d2 <- c()
PPC2.d2 <- c()
paca1.d2 <- c()
paca2.d2 <- c()
#var.dist3 <- c()
#intra.allo.dist <- c()
for(i in 1:nsim){
  ad <- numeric(length(proc.list.20))
  PC1D <- numeric(length(ad))
  PC2D <- numeric(length(ad))
  pPC1D <- numeric(length(ad))
  pPC2D <- numeric(length(ad))
  PAC1D <- numeric(length(ad))
  PAC2D <- numeric(length(ad))
  v <- numeric(length(ad))
  v2 <- numeric(length(ad))
  v3 <- numeric(length(ad))

  for(j in 1:length(proc.list.20)){
      l <- dim(allo.free[[j]])[3]
      sam <- sample(1:l, Nmin, replace = T)
      A2d <- two.d.array(A=allo.free[[j]][,,sam])
      vcv <- var(A2d)
      svd1 <- svd(vcv)
      ##v[j] <- svd1$d[1]/sum(svd1$d)
      ##v2[j] <- svd1$d[2]/sum(svd1$d)
      ##v3[j] <- svd1$d[3]/sum(svd1$d)
      #v[j] <- projected.variance(vcv,svd1$u[,1], scale=T)
      #v2[j] <- projected.variance(vcv,svd1$u[,2], scale=T)
       ad[j] <- projected.variance(vcv,allo.vector, scale=T)
      # pPC1D[j] <- projected.variance(vcv,pPC1, scale=T)
      # pPC2D[j] <- projected.variance(vcv,pPC2, scale=T)
      # PAC1D[j] <- projected.variance(vcv,PaC1, scale=T)
      # PAC2D[j] <- projected.variance(vcv,PaC2, scale=T)
      # PC1D[j] <- projected.variance(vcv,PC1, scale=T)
      # PC2D[j] <- projected.variance(vcv,PC2, scale=T)
  }

  Vproj.allo.distribution2[[i]] <- ad
  # pc1.d2[[i]] <- PC1D
  # pc2.d2[[i]] <- PC2D
  # PPC1.d2[[i]] <- pPC1D
  # PPC2.d2[[i]] <- pPC2D
  # paca1.d2[[i]] <- PAC1D
  # paca2.d2[[i]] <- PAC2D
}



# intraspecific metrics
evd2 <- c() # integration 1, eigenvalue dispersion (deviation of the expected relative SD of eigenvalues)
dis2 <- c() # disparity 1, Procrustes variance
dmats2 <- c()
for(i in 1:length(proc.list.20)){
  A <- allo.free[[i]]
#  dimnames(A)[[3]] <- 1:length(dimnames(A)[[3]])
  M <- mshape(A)
  dis2[[i]]<-morphol.disparity(A~1,print.progress = FALSE)
  evd2[[i]] <- eigen.var(A, sd = TRUE, rel = TRUE, sample = dim(A)[3])
  dmats2[[i]] <- gpagen(A,verbose = TRUE, print.progress = FALSE)$procD
  }
dis2 <- unlist(dis2)
names(dis2) <- names(specs20)
evd2 <- unlist(evd2)
names(evd2) <- names(specs20)
maxD2 <- unlist(lapply(dmats2, FUN = function(X)max(X)))
names(maxD2) <- names(specs20)

# integration 2, effective eigenvalue dispersion
rdi.R2 <- numeric(length(allo.free))
for(i in 1:length(proc.list.20)){
    A <- allo.free[[i]]
    dimnames(A)[[3]] <- 1:dim(A)[3]
    rdi.R2[i] <- rdi(A, replace = TRUE, nsim = 99, nreps = Nmin)$EffectiveDispersion
    }
names(rdi.R2) <- tree.20$tip.label

allo.vector <- R.size.pgls$pgls.coefficients[2,]

Vproj2 <- numeric(length(proc.list.20))
Mr <- mshape(R.df$coords)
for(i in 1:length(proc.list.20)){
  A1 <- allo.free[[i]]
  vcv <- var(two.d.array(A=A1))
  Vproj2[i] <- projected.variance(vcv,allo.vector)
}
names(Vproj2)<-tree.20$tip.label


par(mar=c(5,5,2,2),mfrow=c(1,1))
# C - PHYLOHEATMAP
intra.df.allo.free <- data.frame(disparity1 = dis2, disparity2 = maxD2, integration1 = evd2, integration2 = rdi.R2, Vproj = Vproj2, row.names = tree.20$tip.label)
intra.df2.allo.free <- intra.df.allo.free ; colnames(intra.df2.allo.free) <- c("A","B","C","D","E")

phylo.heatmap(intra.tree,intra.df2.allo.free,colors = col.pal.sunset(100), standardize = T)



# D - PAIRWISE LINEAR REGRESSIONS
panel.cor <- function(x, y, ncomp = 10){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    P <- anova(gls(x~y, correlation = corPagel(0,phy=tree.20)))$`p-value`[2]
    #P <- anova(lm(x~y))$`Pr(>F)`[1]
    P <- round(P, digits=3)
    corr <- cor(x,y)
    txt <- paste0("p = ", P,"\nr = ",round(corr, digits=3))
    if(P<(0.05/ncomp)){ft<-2 ; cex <- 2}else{if(P<0.05){ft<-3}else{ft<-1} ; cex <- 1.5}
    text(0.5, 0.5, txt, cex = cex, font = ft)
}
TREEE <- keep.tip(phy = tree.20, tip = rownames(intra.df)[which(intra.N > 39)])
panel.cor2 <- function(x, y, ncomp = 10){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    P <- anova(gls(x~y, correlation = corPagel(0,phy=TREEE)))$`p-value`[2]
    #P <- anova(lm(x~y))$`Pr(>F)`[1]
    P <- round(P, digits=3)
    corr <- cor(x,y)
    txt <- paste0("p = ", P,"\nr = ",round(corr, digits=3))
    if(P<(0.05/ncomp)){ft<-2 ; cex <- 2}else{if(P<0.05){ft<-3}else{ft<-1} ; cex <- 1.5}
    text(0.5, 0.5, txt, cex = cex, font = ft)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch=pts.pch[tree.20$tip.label],
     bg=pts.col[tree.20$tip.label],)
  abline(lm(y~x),col="lightgray")
}
upper.panel2<-function(x, y){
  points(x,y, pch=pts.pch[TREEE$tip.label],
     bg=pts.col[TREEE$tip.label],)
  abline(lm(y~x),col="lightgray")
}

panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    D1 <- density(x) ; D1$y <- D1$y/max(D1$y)
    lines(D1,lwd=3,col="lightgray")
    fam.porp <- rev(table(f.fam))
    uf <- unique(f.fam)
    for(i in 2:3){
    D <- density(x[which(f.fam == uf[i])])
    D$y <- D$y / max(D$y)
    D$y <- D$y * (fam.porp[i] / sum(fam.porp))
    lines(D, lwd = 3, col = family.colors[c(3,1,2)][i])
    }
}
panel.hist2 <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    D1 <- density(x) ; D1$y <- D1$y/max(D1$y)
    lines(D1,lwd=3,col="lightgray")
    FFAM <-  f.fam[rownames(intra.df)[which(intra.N > 39)]]
    fam.porp <- rev(table(FFAM))
    uf <- unique(FFAM)
    for(i in 2:3){
    D <- density(x[which(FFAM == uf[i])])
    D$y <- D$y / max(D$y)
    D$y <- D$y * (fam.porp[i] / sum(fam.porp))
    lines(D, lwd = 3, col = family.colors[c(3,1,2)][i])
    }
}

pairs(intra.df.allo.free,
      lower.panel = panel.cor,
      upper.panel = upper.panel,diag.panel = panel.hist,
      labels = c("A - disparity 1", "B - disparity 2","C - integration 1", "D - integration 2", "E - Vproj"))


# only including the specimens with at least 40 individuals
pairs(intra.df.allo.free[which(intra.N > 39),],
      lower.panel = panel.cor2,
      upper.panel = upper.panel2,
      diag.panel = panel.hist2,
      labels = c("A - disparity 1", "B - disparity 2","C - integration 1", "D - integration 2", "E - Vproj"))

#############

# Allometry by subfamily
# supplemental figure 1
############
pcafit <- gm.prcomp(R.size.subfam$GM$pgls.fitted)
preds.fit <- pcafit$x[,1]
  preds.resids <- gm.prcomp(R.size.tribe$GM$pgls.residuals)$x[,1]


par(mfrow=c(1,2))
plot(NA,xlim=c(5,6.67),ylim=c(-.15,.155),xlab = "log(size)", ylab = "regression score")
grid(lty=1)
  ut <- unique(R.df$subfamily) ; utcols <- colourWheel(length(ut))
for(i in 1:length(ut)){
  index <- which(R.df$subfamily == ut[i])
    points(R.df$csize[index],preds.fit[index],col=utcols[i],type='l',lwd=2)
  }
legend('topleft',fill=utcols, legend = ut, cex = 0.7)
plot(NA,xlim=c(-0.5,0.5),ylim=c(-1,1),axes=F,frame=F,xlab=NA,ylab=NA,asp=1)
tps(pcafit$shapes$shapes.comp1$min[lat.drop,-3],mshape(R.df$coords)[lat.drop,-3],
    scale.lm=1.3,add=T,at=c(0,-0.83),shade = T, links = links.lat, plot.lm = F, plot.ref.links = T)
tps(mshape(R.df$coords)[lat.drop,-3],mshape(R.df$coords)[lat.drop,-3],
    scale.lm=1.3,add=T,at=c(0,0),shade = F, links = links.lat, plot.lm = F, plot.ref.links = F)
tps(pcafit$shapes$shapes.comp1$max[lat.drop,-3],mshape(R.df$coords)[lat.drop,-3],
    scale.lm=1.3,add=T,at=c(0,0.85),shade = T, links = links.lat, plot.lm = F, plot.ref.links = T)
#############

# all ordination plots
# supplemental figure 2
##############



R.Ppca <- gm.prcomp(A = R, phy = R.df$tree, GLS=T,transform = T)
R.paca <- gm.prcomp(A = R, phy = R.df$tree, align.to.phy = T)

cex <- 1

par(mfrow=c(3,2),mar=c(4,4,1,1))


# ADD ALLOMETRY VECTOR
   plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=1,cex.lab=1.2,
       xlab = paste("PC1: ",toString(round(R.pca$d[1]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("PC2: ",toString(round(R.pca$d[2]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.pca$x[,1:2],tree = R.df$tree)
  points(R.pca$x[,1],R.pca$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)
  legend("topleft", c("Tragulidae","Cervidae","Moschidae","Bovidae","Common Ancestor"), col = 'black',pt.lwd = 1, pt.cex = 1,pt.bg = c(family.colors[c(3,1,4,2)],"goldenrod1"), pch = c(24,21,25,22,23),title.adj = -1, cex = 0.5)



  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=1,cex.lab=1,
       xlab = "PC1",ylab="PC2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.pca$x[,1:2],tree = R.df$tree)
points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 2)
   #  mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
   # mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)


xseg <- c(-0.15,0,0.15)   ; yseg <- c(-0.15,0,0.15)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 2, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.pca$rotation[,1]*scores[i,1]+R.pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  if(i==5){SHADE<-F}else{SHADE<-T}
  scrs <- as.matrix(scores[i,])
  up <- c(scrs[1],scrs[2]+0.0375)
  down <- c(scrs[1]-0.0071,scrs[2]-0.025)
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(R)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.165, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = up,
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
  tlm <- shapes1[dor.drop,-2,i] ; tlm[,2] <- -tlm[,2]
  rlm <- mshape(R)[dor.drop,-2] ; rlm[,2] <- -rlm[,2]
  tps(target.lm = tlm,
      reference.lm = rlm,
      n.grid.col = 24,
      scale.lm = 0.16, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = down,
      links = links.dor,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
}









   plot(NA,xlim=c(-0.05,0.05),ylim=c(-.05,.05),asp=T,cex.lab=1,
       xlab = paste("pPC1: ",toString(round(R.Ppca$d[1]/sum(R.Ppca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("pPC2: ",toString(round(R.Ppca$d[2]/sum(R.Ppca$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.Ppca$x[,1:2],tree = R.df$tree)
  points(R.Ppca$x[,1],R.Ppca$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.Ppca$anc.x[1,1:2],R.Ppca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)







  plot(NA,xlim=c(-0.05,0.05),ylim=c(-.05,.05),asp=T,cex.lab=1,
       xlab = "pPC1",ylab="pPC2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.Ppca$x[,1:2],tree = R.df$tree)
points(rbind(R.Ppca$anc.x[1,1:2],R.Ppca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)

xseg <- c(-0.0375,0,0.0375)   ; yseg <- c(-0.0375,0,0.0375)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 2, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.Ppca$rotation[,1]*scores[i,1]+R.Ppca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  if(i==5){SHADE<-F}else{SHADE<-T}
  scrs <- as.matrix(scores[i,])
  up <- c(scrs[1],scrs[2]+0.0075)
  down <- c(scrs[1]-0.00142,scrs[2]-0.005)
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(R)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.0375, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = up,
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
  tlm <- shapes1[dor.drop,-2,i] ; tlm[,2] <- -tlm[,2]
  rlm <- mshape(R)[dor.drop,-2] ; rlm[,2] <- -rlm[,2]
  tps(target.lm = tlm,
      reference.lm = rlm,
      n.grid.col = 24,
      scale.lm = 0.0375, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = down,
      links = links.dor,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
}












   plot(NA,xlim=c(-0.15,0.15),ylim=c(-.15,.15),asp=T,cex.lab=1,
       xlab = paste("C1: ",toString(round(R.paca$d[1]/sum(R.paca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("C2: ",toString(round(R.paca$d[2]/sum(R.paca$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.paca$x[,1:2],tree = R.df$tree)
  points(R.paca$x[,1],R.paca$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.paca$anc.x[1,1:2],R.paca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)



  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1,
       xlab = "C1",ylab="C2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.paca$x[,1:2],tree = R.df$tree)
points(rbind(R.paca$anc.x[1,1:2],R.paca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)


xseg <- c(-0.15,0,0.15)   ; yseg <- c(-0.15,0,0.15)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 1, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.paca$rotation[,1]*scores[i,1]+R.paca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  if(i==5){SHADE<-F}else{SHADE<-T}
  scrs <- as.matrix(scores[i,])
  up <- c(scrs[1],scrs[2]+0.0375)
  down <- c(scrs[1]-0.0071,scrs[2]-0.025)
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(R)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.165, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = up,
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
  tlm <- shapes1[dor.drop,-2,i] ; tlm[,2] <- -tlm[,2]
  rlm <- mshape(R)[dor.drop,-2] ; rlm[,2] <- -rlm[,2]
  tps(target.lm = tlm,
      reference.lm = rlm,
      n.grid.col = 24,
      scale.lm = 0.16, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = down,
      links = links.dor,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
}



# with interspecific allometry removed

R.mean <- mshape(R.df$coords)

R.noallo <- array(NA,dim=dim(R.df$coords))
for(i in 1:dim(R.noallo)[3]){R.noallo[,,i] <- unflatten(flatten(R.mean) + R.size.pgls$pgls.residuals[i,])}
dimnames(R.noallo)[[3]] <- R.df$tree$tip.label

R.pca.noallo <- gm.prcomp(R.noallo, phy = R.df$tree)



R.Ppca <- gm.prcomp(A = R.noallo, phy = R.df$tree, GLS=T,transform = T)
R.paca <- gm.prcomp(A = R.noallo, phy = R.df$tree, align.to.phy = T)
R.pca <- gm.prcomp(R.noallo,R.df$tree)

R.bm12.noallo <- mvBM(tree=R.df$tree, data = R.pca$x[,1:2])


cex <- 1

par(mfrow=c(3,2),mar=c(4,4,1,1))


# ADD ALLOMETRY VECTOR
   plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=1,cex.lab=1.2,
       xlab = paste("PC1: ",toString(round(R.pca$d[1]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("PC2: ",toString(round(R.pca$d[2]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.pca$x[,1:2],tree = R.df$tree)
  points(R.pca$x[,1],R.pca$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)
  legend("topleft", c("Tragulidae","Cervidae","Moschidae","Bovidae","Common Ancestor"), col = 'black',pt.lwd = 1, pt.cex = 1,pt.bg = c(family.colors[c(3,1,4,2)],"goldenrod1"), pch = c(24,21,25,22,23),title.adj = -1, cex = 0.5)



  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=1,cex.lab=1,
       xlab = "PC1",ylab="PC2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.pca$x[,1:2],tree = R.df$tree)
points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 2)
#     mixtools::ellipse(mu = as.numeric(R.bm12.noallo$theta), sigma = R.bm12.noallo$sigma*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
#    mixtools::ellipse(mu = as.numeric(R.bm12.noallo$theta), sigma = R.bm12.noallo$sigma*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)


xseg <- c(-0.15,0,0.15)   ; yseg <- c(-0.15,0,0.15)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 2, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.pca$rotation[,1]*scores[i,1]+R.pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  if(i==5){SHADE<-F}else{SHADE<-T}
  scrs <- as.matrix(scores[i,])
  up <- c(scrs[1],scrs[2]+0.0375)
  down <- c(scrs[1]-0.0071,scrs[2]-0.025)
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(R)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.165, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = up,
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
  tlm <- shapes1[dor.drop,-2,i] ; tlm[,2] <- -tlm[,2]
  rlm <- mshape(R)[dor.drop,-2] ; rlm[,2] <- -rlm[,2]
  tps(target.lm = tlm,
      reference.lm = rlm,
      n.grid.col = 24,
      scale.lm = 0.16, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = down,
      links = links.dor,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
}









   plot(NA,xlim=c(-0.05,0.05),ylim=c(-.05,.05),asp=1,cex.lab=1,
       xlab = paste("pPC1: ",toString(round(R.Ppca$d[1]/sum(R.Ppca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("pPC2: ",toString(round(R.Ppca$d[2]/sum(R.Ppca$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.Ppca$x[,1:2],tree = R.df$tree)
  points(R.Ppca$x[,1],R.Ppca$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.Ppca$anc.x[1,1:2],R.Ppca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)







  plot(NA,xlim=c(-0.05,0.05),ylim=c(-.05,.05),asp=1,cex.lab=1,
       xlab = "pPC1",ylab="pPC2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.Ppca$x[,1:2],tree = R.df$tree)
points(rbind(R.Ppca$anc.x[1,1:2],R.Ppca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)

xseg <- c(-0.0375,0,0.0375)   ; yseg <- c(-0.0375,0,0.0375)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 2, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.Ppca$rotation[,1]*scores[i,1]+R.Ppca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  if(i==5){SHADE<-F}else{SHADE<-T}
  scrs <- as.matrix(scores[i,])
  up <- c(scrs[1],scrs[2]+0.0075)
  down <- c(scrs[1]-0.00142,scrs[2]-0.005)
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(R)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.0375, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = up,
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
  tlm <- shapes1[dor.drop,-2,i] ; tlm[,2] <- -tlm[,2]
  rlm <- mshape(R)[dor.drop,-2] ; rlm[,2] <- -rlm[,2]
  tps(target.lm = tlm,
      reference.lm = rlm,
      n.grid.col = 24,
      scale.lm = 0.0375, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = down,
      links = links.dor,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
}












   plot(NA,xlim=c(-0.15,0.15),ylim=c(-.15,.15),asp=1,cex.lab=1,
       xlab = paste("C1: ",toString(round(R.paca$d[1]/sum(R.paca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("C2: ",toString(round(R.paca$d[2]/sum(R.paca$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.paca$x[,1:2],tree = R.df$tree)
  points(R.paca$x[,1],R.paca$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.paca$anc.x[1,1:2],R.paca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)



  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1,
       xlab = "C1",ylab="C2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.paca$x[,1:2],tree = R.df$tree)
points(rbind(R.paca$anc.x[1,1:2],R.paca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)


xseg <- c(-0.15,0,0.15)   ; yseg <- c(-0.15,0,0.15)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 1, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.paca$rotation[,1]*scores[i,1]+R.paca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  if(i==5){SHADE<-F}else{SHADE<-T}
  scrs <- as.matrix(scores[i,])
  up <- c(scrs[1],scrs[2]+0.0375)
  down <- c(scrs[1]-0.0071,scrs[2]-0.025)
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(R)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.165, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = up,
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
  tlm <- shapes1[dor.drop,-2,i] ; tlm[,2] <- -tlm[,2]
  rlm <- mshape(R)[dor.drop,-2] ; rlm[,2] <- -rlm[,2]
  tps(target.lm = tlm,
      reference.lm = rlm,
      n.grid.col = 24,
      scale.lm = 0.16, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = down,
      links = links.dor,plot.ref.links = T,
      ref.link.aes = list(lwd = 1, col = "darkgray"),
      link.aes = list(lwd = 1, col = "black"),
      grid.aes = list(col = 'gray', trans = 1, lwd = 0.25),
      plot.lm = FALSE)
}




###############

