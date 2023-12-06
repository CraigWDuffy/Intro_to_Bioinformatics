
tajimad1 <- function(haps,lege,Gene,msnp,minn)
{
keep1 <- lege[,"gene"] == Gene
haps1 <- haps[keep1,]
lege1 <- lege[keep1,]
results <- matrix(NA,nrow=1,ncol=11)
colnames(results) <- c("Gene","H","S","theta","khat","tajimasd","vars","Uniq","subtelo","newpos","colors")
results[1,1] <- Gene
ss <- sum(1*( keep1));
results[1,3] <- ss
if(ss > 1)
{
vars <-  mean(lege1[,"vars"],na.rm=T)
subtelo <- mean(lege1[,"subtelo"],na.rm=T)
Uniq <- mean(lege1[,"uniq"],na.rm=T)
newpos <- mean(lege1[,"newpos"],na.rm=T)
colors <- mean(lege1[,"colors"],na.rm=T)
temp <- haps1

temps <- apply(temp,2,collap); count.temp <- table(temps)

n.ind <- sum(count.temp); freqs.temp <- count.temp / n.ind
results[1,2] <- length(count.temp);
results[1,7] <- vars; results[1,8] <- Uniq; results[1,9] <- subtelo; results[1,10] <- newpos
results[1,11] <- colors
if(as.numeric(results[1,3]) >= msnp)
{
all  <- as.matrix( dist(t(temp), "manhattan",diag=T,upper=T))
khat <- mean(c((all)[lower.tri(all)]),na.rm=T)
a1 <- sum( 1/(1:(n.ind-1)))
a2 <- sum( 1/(1:(n.ind-1))^2);
b1 <- (n.ind + 1)/3/(n.ind-1)
b2 <- 2*(n.ind^2 + n.ind +3)/9/n.ind/(n.ind-1); c1 <- b1 - 1/a1
c2 <- b2 - (n.ind + 2)/ a1 / n.ind + a2/a1/a1; e1 <- c1/a1
e2 <- c2 / (a1^2 + a2)
results[1,4] <- ss / a1
results[1,5] <- khat
results[1,6] <- (khat - ss/a1) / sqrt ( e1 * ss + e2 * ss * (ss-1))
}
}
results[!is.na(results[1,"tajimasd"]),]
}

collap <- function(x)  {  paste(x,collapse="") }

plottd <- function(resgene,thres)
{
plot(resgene[,"newpos"],resgene[,"tajimasd"],pch=16,cex=1, xlab="Genes sorted by genomic position", ylab="Per gene Tajima's D score",ylim=c(-3,4))
points(resgene[,"newpos"],resgene[,"tajimasd"],pch=16,cex=1, col=c("black","red")[as.numeric(resgene[,"colors"])])
abline(h=0,lty=2)
hits <- which(as.numeric(resgene[,"tajimasd"]) > thres)
colblue <- resgene[hits,]
points(colblue[,"newpos"],colblue[,"tajimasd"],pch=16,cex=1,col="blue")
}



converts <- function(x,ref)
{
y <- matrix(0.5,nrow=nrow(x),ncol=ncol(x))
colnames(y) <- colnames(x)

refs  <- which(ref=="A");
yy <- as.matrix(x[refs,]); yy[yy == "A"] <- 0; yy[yy == "C"] <- 1;
yy[yy == "G"] <- 1; yy[yy == "T"] <- 1; yy[yy == "N"] <- NA; yy[yy == "-"] <- NA;
yy[yy == "AG"] <- 0.5; yy[yy == "AC"] <- 0.5; yy[yy == "AT"] <- 0.5;

y[refs,] <- yy

refs  <- which(ref=="C")
yy <- as.matrix(x[refs,]); yy[yy == "C"] <- 0; yy[yy == "A"] <- 1;
yy[yy == "G"] <- 1; yy[yy == "T"] <- 1; yy[yy == "N"] <- NA; yy[yy == "-"] <- NA;
yy[yy == "AC"] <- 0.5; yy[yy == "CG"] <- 0.5; yy[yy == "CT"] <- 0.5;
y[refs,] <- yy

refs  <- which(ref=="G")
yy <- as.matrix(x[refs,]); yy[yy == "G"] <- 0; yy[yy == "C"] <- 1;
yy[yy == "A"] <- 1; yy[yy == "T"] <- 1; yy[yy == "N"] <- NA; yy[yy == "-"] <- NA;
yy[yy == "AG"] <- 0.5; yy[yy == "CG"] <- 0.5; yy[yy == "GT"] <- 0.5;
y[refs,] <- yy

refs  <- which(ref=="T")
yy <- as.matrix(x[refs,]); yy[yy == "T"] <- 0; yy[yy == "C"] <- 1;
yy[yy == "G"] <- 1; yy[yy == "A"] <- 1; yy[yy == "N"] <- NA; yy[yy == "-"] <- NA;
yy[yy == "AT"] <- 0.5; yy[yy == "CT"] <- 0.5; yy[yy == "GT"] <- 0.5;
y[refs,] <- yy

y
}

