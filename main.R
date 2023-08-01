library(MFPCA);library(tidyverse);library(fda);library(survAUC);library(Triangulation);library(BPST); library(readr); options(warn=-1)

source("ancillary.R")

train_dat.id = tdat.id = read_rds("training")
test_dat.id = read_rds("testing")

ncr = n1= n2 =40;
npix=n1*n2
u1=seq(0,1,length.out=n1)
v1=seq(0,1,length.out=n2)
uu=rep(u1,each=n2)
vv=rep(v1,times=n1)
Z=as.matrix(cbind(uu,vv))

load("V2.rda");load("Tr2.rda")
ind.inside=inVT(Brain.V2, Brain.Tr2, Z[,1], Z[,2])$ind.inside
V.est<-Brain.V2; Tr.est<-Brain.Tr2


idx = which(df.left$id %in% train_dat.id$id)

est1 = bernstein(Y = train_dat.id$Y1, V.est = V.est, Tr.est = Tr.est, d.est = 2, 
                 r = 1, Z = Z, lambda = c(0, 1, 10, 10^2, 10^3, 10^6))
B.est = est1[[1]]; S1.est = t(est1[[2]]) 

est2 = bernstein(Y = train_dat.id$Y2, V.est = V.est, Tr.est = Tr.est, d.est = 2, 
                 r = 1, Z = Z, lambda = c(0, 1, 10, 10^2, 10^3, 10^6))
S2.est = t(est2[[2]]) 

Emat = cbind(S1.est, S2.est)

npc = c(ncol(S1.est), ncol(S2.est)); M = 3
Z.hat = Emat / (nrow(S1.est)-1)
tmpSVD = irlba::irlba(as.matrix(Z.hat), nv = M)
vectors <- tmpSVD$v
values <- tmpSVD$d^2

normFactors <- 1/sqrt(diag(as.matrix(Matrix::crossprod(Z.hat %*% vectors))))
Mscore <- crossprod(t(Z.hat), vectors) * sqrt(nrow(S1.est)-1)
Mscore <- as.matrix( crossprod(t(Mscore), diag(sqrt(values)) * normFactors, nrow = M, ncol = M)) # normalization

tmpWeights <- as.matrix(Matrix::crossprod(Z.hat, Z.hat %*%vectors))
npcCum <- cumsum(c(0, npc)) # indices for blocks (-1)

score1_in = 1/sqrt(values) * normFactors * t(tmpWeights[npcCum[1]+seq_len(npc[1]), , drop = FALSE])
score2_in = 1/sqrt(values) * normFactors * t(tmpWeights[npcCum[2]+seq_len(npc[2]), , drop = FALSE])

org.idx = c(1:ncr*ncr)
desMat = matrix(0, nrow = ncr*ncr, ncol = ncol(score1_in))

for(zz in 1:ncol(desMat)){
  desMat[ind.inside,zz] = B.est[,zz]
}

psi1 = t(tcrossprod(score1_in, desMat))
psi2 = t(tcrossprod(score2_in, desMat))

psi1 = psi1[ind.inside,]
psi2 = psi2[ind.inside,]


#------- score training ---------#

score_fpca = as.matrix(Mscore[,c(1:M)])
score_names = c()
for(q in 1:(ncol(score_fpca))){
  tname = paste("score", as.character(q), sep = "")
  score_names = c(score_names, tname)
}
tdat.id = train_dat.id[,c(1:6)]
tdat.id = cbind(tdat.id, score_fpca)

colnames(tdat.id)[(ncol(tdat.id) - ncol(score_fpca) + 1) : ncol(tdat.id)] = score_names

fmla = as.formula(paste("Surv(time,event) ~ ", paste(score_names, collapse= "+")))

fitted_obj = coxph(fmla, data = tdat.id, x = T, y = T) # this is our survival training obj with saved parameters

#----------- testing -------------#

est1 = bernstein(Y = test_dat.id$Y1, V.est = V.est, Tr.est = Tr.est, d.est = 2,
                 r = 1, Z = Z, lambda = c(0, 1, 10, 10^2, 10^3, 10^6))
est2 = bernstein(Y = test_dat.id$Y2, V.est = V.est, Tr.est = Tr.est, d.est = 2,
                 r = 1, Z = Z, lambda = c(0, 1, 10, 10^2, 10^3, 10^6))
b.basis.test1 = est1[[1]]; b.scores.test1 = t(est1[[2]]) 
b.basis.test2 = est2[[1]]; b.scores.test2 = t(est2[[2]]) 

Z.hat.test = cbind(b.scores.test1, b.scores.test2)/(nrow(b.scores.test1)-1)

normFactors.test <- 1/sqrt(diag(as.matrix(Matrix::crossprod(Z.hat.test %*% vectors))))
Mscore.test <- crossprod(t(Z.hat.test), vectors) * sqrt(nrow(b.scores.test1)-1) 
Mscore.test <- as.matrix( crossprod(t(Mscore.test), diag(sqrt(values)) * normFactors.test, nrow = M, ncol = M)) # normalization

score_names = c()
for(q in 1:(ncol( Mscore.test))){
  tname = paste("score", as.character(q), sep = "")
  score_names = c(score_names, tname)
}

testd.id = test_dat.id[,c(1:4)]

testd.id = cbind(testd.id, Mscore.test)

colnames(testd.id)[(ncol(testd.id) - ncol(Mscore.test) + 1) : ncol(testd.id)] = score_names

train.prob = predict(fitted_obj)
surv.prob = predict(fitted_obj, newdata = testd.id)

Surv.train = Surv(tdat.id$time, tdat.id$event)
Surv.test = Surv(testd.id$time, testd.id$event)

times = seq(0, 15, by = 0.1)

# iAUC 
AUC.uno(Surv.train, Surv.test, surv.prob, times)$iauc

