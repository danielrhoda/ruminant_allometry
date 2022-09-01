# Supplemental Figure 1
# Allometric slopes by subfamily
########
R.size.subfam <- procD.pgls(f1 = coords ~ csize*subfamily, phy = tree, data=R.df)
anova(R.size.subfam)
# DIFFERENT EVOLUTIONARY ALLOMETRIC TRAJECTORIES BY TRIBE

# DOES THE INTERACTION AFFECT THE FIT OF THE MODEL?
R.size.subfam2 <- procD.pgls(f1 = coords ~ csize, phy = tree, data=R.df)
anova(R.size.subfam2)

anova(R.size.subfam2,R.size.subfam)
# it does

#write.csv(x=R.size.tribe$aov.table,file = "results/R.size.tribe.PGLS.csv")


pcafit <- gm.prcomp(R.size.subfam$GM$pgls.fitted)
preds.fit <- pcafit$x[,1]
  preds.resids <- gm.prcomp(R.size.subfam$GM$pgls.residuals)$x[,1]


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
##########




# Supplemental Figure 2 - ALL ORDINATIONS (PCA, PPCA, PACA)
#######

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



# the same ordinations but with evolutinoary allometry removed
R.mean <- mshape(R.df$coords)

R.noallo <- array(NA,dim=dim(R.df$coords))
for(i in 1:dim(R.noallo)[3]){R.noallo[,,i] <- unflatten(flatten(R.mean) + R.size.pgls$pgls.residuals[i,])}
dimnames(R.noallo)[[3]] <- R.df$tree$tip.label

R.pca.noallo <- gm.prcomp(R.noallo, phy = R.df$tree)


R.pca2 <- gm.prcomp(R.noallo, phy = R.df$tree)
R.Ppca2 <- gm.prcomp(A = R.noallo, phy = R.df$tree, GLS=T,transform = T)
R.paca2 <- gm.prcomp(A = R.noallo, phy = R.df$tree, align.to.phy = T)

R.bm12.noallo <- mvBM(tree=R.df$tree, data = R.pca.noallo$x[,1:5])


cex <- 1

par(mfrow=c(3,2),mar=c(4,4,1,1))


# ADD ALLOMETRY VECTOR
   plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1.2,
       xlab = paste("PC1: ",toString(round(R.pca2$d[1]/sum(R.pca2$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("PC2: ",toString(round(R.pca2$d[2]/sum(R.pca2$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.pca2$x[,1:2],tree = R.df$tree)
  points(R.pca2$x[,1],R.pca2$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.pca2$anc.x[1,1:2],R.pca2$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)
  legend("topleft", c("Tragulidae","Cervidae","Moschidae","Bovidae","Common Ancestor"), col = 'black',pt.lwd = 1, pt.cex = 1,pt.bg = c(family.colors[c(3,1,4,2)],"goldenrod1"), pch = c(24,21,25,22,23),title.adj = -1, cex = 0.5)



  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1,
       xlab = "PC1",ylab="PC2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.pca2$x[,1:2],tree = R.df$tree)
points(rbind(R.pca2$anc.x[1,1:2],R.pca2$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 2)
#     mixtools::ellipse(mu = as.numeric(R.bm12.noallo$theta), sigma = R.bm12.noallo$sigma*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
#    mixtools::ellipse(mu = as.numeric(R.bm12.noallo$theta), sigma = R.bm12.noallo$sigma*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)


xseg <- c(-0.15,0,0.15)   ; yseg <- c(-0.15,0,0.15)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 2, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.pca2$rotation[,1]*scores[i,1]+R.pca2$rotation[,2]*scores[i,2]
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
       xlab = paste("pPC1: ",toString(round(R.Ppca2$d[1]/sum(R.Ppca2$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("pPC2: ",toString(round(R.Ppca2$d[2]/sum(R.Ppca2$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.Ppca2$x[,1:2],tree = R.df$tree)
  points(R.Ppca2$x[,1],R.Ppca2$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.Ppca2$anc.x[1,1:2],R.Ppca2$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)







  plot(NA,xlim=c(-0.05,0.05),ylim=c(-.05,.05),asp=T,cex.lab=1,
       xlab = "pPC1",ylab="pPC2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.Ppca2$x[,1:2],tree = R.df$tree)
points(rbind(R.Ppca2$anc.x[1,1:2],R.Ppca2$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)

xseg <- c(-0.0375,0,0.0375)   ; yseg <- c(-0.0375,0,0.0375)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 2, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.Ppca2$rotation[,1]*scores[i,1]+R.Ppca2$rotation[,2]*scores[i,2]
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
       xlab = paste("C1: ",toString(round(R.paca2$d[1]/sum(R.paca2$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("C2: ",toString(round(R.paca2$d[2]/sum(R.paca2$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

  addtree2(R.paca2$x[,1:2],tree = R.df$tree)
  points(R.paca2$x[,1],R.paca2$x[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale(log(R.size),c(0.5,2.25)))

  points(rbind(R.paca2$anc.x[1,1:2],R.paca2$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)



  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1,
       xlab = "C1",ylab="C2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
  addtree2(R.paca2$x[,1:2],tree = R.df$tree)
points(rbind(R.paca2$anc.x[1,1:2],R.paca2$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 1, lwd = 1)


xseg <- c(-0.15,0,0.15)   ; yseg <- c(-0.15,0,0.15)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 1, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.paca2$rotation[,1]*scores[i,1]+R.paca2$rotation[,2]*scores[i,2]
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



#######

# Supplemental Figure 3 - Raw %Grass in diet data
#########

plot(NA, xlim = range(R.pca$x[,1]), ylim = range(R.pca$x[,2]), asp = 1,
     xlab = paste("PC1: ",toString(round(R.pca$d[1]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("PC2: ",toString(round(R.pca$d[2]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""),
     main = "Percent Grass in Diet Data")
grid(lty=1,col="lightgray")
points(R.pca$x[,1:2],cex = rescale.numeric(R.size,c(1.5,3)))
points(R.pca$x[R2.df.n,1:2],cex = rescale.numeric(R.size,c(1.5,3))[R2.df.n],
       col = make.color.vec(100-R2.df$percentgrass,11,brewer_pal(palette = "BrBG")), pch = 19)
text(x = c(-0.2,-0.2),
     y = c(0.08, 0.18),
     labels = c("0%", "100%"))
shape::colorlegend(col = brewer_pal(palette = "BrBG",direction=-1)(11), zlim = c(0,100), zlevels = seq(0,100,length.out = 11),
            posx = c(0.1,0.2), posy = c(0.7,0.9), main = "", lab.col = "transparent")


#########

