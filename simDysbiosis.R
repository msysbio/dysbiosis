# Simulate the effect of bloomers that grow better in perturbed but not in healthy samples.
# Bloomers do not change their interactions in the dysbiotic case, only their growth rate,
# simulating an environment that is more favorable to them (e.g. with higher levels of oxygen).
# Note that a global interaction matrix and growth rate vector are sub-sampled to create local communities. 
# Local communities are created separately for the dysbiotic and healthy case.
# Different rounds simulate data sets with different global interaction matrices and growth rates.
# With the current parameters, this simulation results in higher evenness and beta diversity and 
# tends to result in lower positive edge percentage (PEP) for healthy than dysbiotic communities. 
# Note that PEP in inferred networks is much higher than defined PEP in the known interaction matrix.

# for gLV with stochasticity consider function simulateGLV in miaSim
library(seqtime)
library(seqgroup)
library(beemStatic)
library(ggplot2)
library(dplyr)
sparseiCov <- SpiecEasi:::sparseiCov # needs to be exported because no longer present in huge library

# Parameters
numBloomers=6 # the set of global taxa that bloom in IBD (note that there can be IBD samples without bloomers if none makes it into the local species pool)
randBloom=TRUE # a third of the bloomers are randomly selected in each IBD sample to bloom instead of all of them (increases beta diversity in IBD samples)
bloomerFactor=4 # bloomers get boosted by multiplying their growth rates with this factor (decreases evenness in IBD samples)
numOtherTaxa=24 # non-blooming taxa that grow worse in IBD
nonBloomerPenalty=0.1 # non-bloomers get penalized by multiplying their growth rates with this factor
numSubtaxa=20 # the number of taxa ending up in sub-communities; should be <= numBloomers+numOtherTaxa 
samples=50 # number of samples in healthy and IBD matrix
pep=20 # pre-defined PEP in global interaction matrix; do not choose too high to avoid explosions in gLV simulations (20 or below)
minocc=10 # only keep rows with at least the given number of non-zero values for network inference; applied after rarefaction and after evenness and beta diversity are computed; CoNet and BEEM require a prevalence filter: 5 for CoNet is enough, BEEM needs a much more stringent one, e.g. 20
method="mb" # network inference method (beem, mb, glasso, pearson/spearman/bray/kld and others, see details in buildNetwork in seqgroup), mb and glasso are provided by SpiecEasi
rounds=10 # data generation rounds with different global interaction matrices
testmode=FALSE # do not re-generate A 

# Functions

# modify A to avoid explosion
stabiliseA<-function(A,b,Atweak="tweak",c=0.2,pep="20",maxIter=20){
  iter=1
  divisor = 3 # num positive edges/divisor are converted into negative edges during tweaking with tinker method tweak
  stability.method="ricker"
  K=b
  y=b
  Atype="klemm"
  N=nrow(A)
  explosion.bound = 10^8
  stable = testStability(A,method=stability.method, K=K, explosion.bound=explosion.bound)
  # try to generate a stable A
  while(stable == FALSE && iter < maxIter){
    print(paste("Generating A, iter",iter))
    iter = iter + 1
    A=generateA(N=N,type=Atype,c=c,pep=pep)
    stable=testStability(A,method=stability.method, K=K, y=y, explosion.bound=explosion.bound)
    if(stable == FALSE){
      if(Atweak == "tweak"){
        print("Tweaking A")
        # tweak matrix until it is stable or too many arcs have become positive
        pos.num=length(A[A>0])
        # number of tweaking trials
        trials = round(pos.num/divisor)
        for(arc.index in 1:trials){
          A=modifyA(A=A,mode=Atweak)
          stable = testStability(A,method=stability.method, K=K, y=y, explosion.bound=explosion.bound)
          if(stable == TRUE){
            print(paste("It took",arc.index,"positive arcs converted into negative ones to generate a stable A."))
            break
          }
        } # end loop over positive arcs
      }else if(Atweak == "schur"){
        # tweak matrix using Schur decomposition
        print("Applying Schur decomposition...")
        A=modifyA(A=A,mode=Atweak)
        stable = testStability(A,method=stability.method, K=K, y=y, explosion.bound=explosion.bound)
        if(stable == TRUE){
          print(paste("Schur decomposition made interaction matrix stable."))
          break
        }
      }else{
        stop(paste("Interaction matrix modification method",Atweak,"is not supported"))
      }
    } # not stable
  } # iterations
  return(A)
}

# prevalence filter
prevalenceFilter<-function(x, minocc=0, keepSum=TRUE){
  toFilter=c()
  xcopy=x
  # convert into presence/absence matrix
  xcopy[xcopy>0]=1
  # sum for each taxon = number of occurrences across samples
  rowsums=apply(xcopy,1,sum)
  toFilter=which(rowsums<minocc)
  if(keepSum) print(paste("Summing",length(toFilter),"taxa that do not pass the prevalence filter."))
  else print(paste("Removing",length(toFilter),"taxa that do not pass the prevalence filter."))
  indices.tokeep=setdiff(c(1:nrow(x)),toFilter)
  # if any taxa were filtered: sum filtered taxa and add sum
  if(keepSum==TRUE && length(toFilter)>0){
    filtered=x[toFilter,]
    x=x[indices.tokeep,]
    # single taxon is being removed
    if(is.null(nrow(filtered))){
      sums.filtered=filtered
    }else{
      sums.filtered=apply(filtered,2,sum)
    }
    x=rbind(x,sums.filtered)
  }else{
    x=x[indices.tokeep,]
  }
  return(x)
}


# check parameters
if(numSubtaxa > (numBloomers+numOtherTaxa)){
  stop("numSubTaxa needs to be smaller or equal to the sum of numBloomers and numOtherTaxa!")
}

if(rounds>1 && testmode){
  print("In test mode, only one round is carried out. Disable test mode for more rounds.")
  rounds=1
}

# prep result vectors
healthy.mediantotalcount=c()
ibd.mediantotalcount=c()
healthy.sheldon=c()
ibd.sheldon=c()
healthy.betadiv=c()
ibd.betadiv=c()
healthy.edgenum=c()
ibd.edgenum=c()
healthy.pep=c()
ibd.pep=c()
simulated.pep=c()

# pick bloomers - these are the first numBloomers taxa, because we later assign growth rates and initial abundances
# with a geometric series and the first ones grow the fastest
bloomer.indices=1:numBloomers

# generate rounds data sets
for(round in 1:rounds){
  print(paste("Round",round))
  
  # initialize matrices with 0 abundances
  healthyMatrix=matrix(nrow=(numBloomers+numOtherTaxa),ncol=samples,0)
  ibdMatrix=matrix(nrow=(numBloomers+numOtherTaxa),ncol=samples,0)
  
  # generate global growth rates with a geometric series that is tuned to be rather even (unevenness comes from bloomers)
  b=generateAbundances(count=(numBloomers+numOtherTaxa)/2,mode=6, k=0.1,N=numBloomers+numOtherTaxa)
  
  # generate growth rates for the IBD case: bloomers get boosted, everyone else gets suppressed
  b.ibd=b
  if(!randBloom){
    b.ibd[bloomer.indices]=b[bloomer.indices]*bloomerFactor
  }
  b.ibd[(numBloomers+1):length(b.ibd)]=b[(numBloomers+1):length(b.ibd)]*nonBloomerPenalty
  
  if(!testmode){
    # generate global interaction matrix (can take some time for larger taxon numbers)
    A=generateA(N=(numBloomers+numOtherTaxa),type="klemm", c=0.2, pep=pep) 
    #plotA(A, method="ggplot",removeOrphans = FALSE)
    # tweak A if necessary to avoid explosions (by replacing positive by negative arcs one by one until stable)
    A=stabiliseA(A, b=b, c=0.2, pep=pep) 
    a.pep=getPep(A)
    print(paste("Final positive edge percentage of A:",a.pep))
    simulated.pep=c(simulated.pep, a.pep)
  }
  
  # generate the requested number of samples using the global interaction matrix
  for(iter in 1:samples){
    b.ibd.local=b.ibd
    
    # randomly sub-sample A twice to create 2 local communities, one for the healthy and the other for the IBD case
    local.taxa.indices=sample(1:nrow(A))[1:numSubtaxa]
    ibd.local.taxa.indices=sample(1:nrow(A))[1:numSubtaxa]
    
    # local realization of global interaction matrix
    Asub=A[local.taxa.indices,local.taxa.indices]
    Asub.ibd=A[ibd.local.taxa.indices,ibd.local.taxa.indices]
    
    # run simulation for healthy case
    healthyCounts =  glv(N=length(local.taxa.indices), b=b[local.taxa.indices], y=b[local.taxa.indices],A=Asub)
    
    # sample last time point as steady-state abundance
    healthyMatrix[local.taxa.indices,iter]=healthyCounts[,ncol(healthyCounts)]
    
    if(randBloom){
      # select a sub-set of the bloomers randomly and raise their growth rates with the bloomer factor
      randBloomerNum=round(numBloomers/3)
      randBloomers=sample(1:numBloomers)[1:randBloomerNum]
      b.ibd.local[randBloomers]=b.ibd.local[randBloomers]*bloomerFactor
    }
    
    # run simulation for IBD case
    ibdCounts = glv(N=length(ibd.local.taxa.indices), b=b.ibd.local[ibd.local.taxa.indices], y=b.ibd.local[ibd.local.taxa.indices],A=Asub.ibd)
    #tsplot(ibdCounts, legend=TRUE)
    # sample last time point as steady-state abundance
    ibdMatrix[ibd.local.taxa.indices,iter]=ibdCounts[,ncol(ibdCounts)]
  }
  
  print("Simulations done")
  
  # compute data properties
  
  # deal with gLV simulation failures: remove columns with missing values
  naColumns.healthy=which(colSums(is.na(healthyMatrix)) > 0)
  if(length(naColumns.healthy)>0){
    print(paste("Discarding ",length(naColumns.healthy),"with missing values."))
    healthyMatrix=healthyMatrix[,setdiff(1:ncol(healthyMatrix),naColumns.healthy)]
  }
  naColumns.ibd=which(colSums(is.na(ibdMatrix)) > 0)
  if(length(naColumns.ibd)>0){
    print(paste("Discarding ",length(naColumns.ibd),"with missing values."))
    ibdMatrix=ibdMatrix[,setdiff(1:ncol(ibdMatrix),naColumns.ibd)]
  }
  
  # set small negative values introduced by simulation to zero
  healthyMatrix[healthyMatrix<0]=0
  ibdMatrix[ibdMatrix<0]=0
  
  # multiply by a factor and round to simulate read counts
  healthyMatrix=round(healthyMatrix*100)
  ibdMatrix=round(ibdMatrix*100)
  
  # take out differences in "sequencing depth" (here these are absolute abundances, but for sequencing data, we have to rarefy or normalize)
  print(paste("Median total abundance in healthy data: ",median(colSums(healthyMatrix))))
  print(paste("Median total abundance in IBD data: ",median(colSums(ibdMatrix))))
  print(paste("Rarefying healthy data to: ",median(colSums(healthyMatrix))/2))
  print(paste("Rarefying IBD data to: ",median(colSums(ibdMatrix))/2))
  healthyMatrix=rarefyFilter(healthyMatrix,min = median(colSums(healthyMatrix))/2)$rar
  ibdMatrix=rarefyFilter(ibdMatrix,min = median(colSums(ibdMatrix))/2)$rar
  print(paste("Zero-sum taxon number after rarefaction in healthy: ",length(which(rowSums(healthyMatrix)==0))))
  print(paste("Zero-sum taxon number after rarefaction in IBD: ",length(which(rowSums(ibdMatrix)==0))))
  healthy.mediantotalcount=c(healthy.mediantotalcount,median(colSums(healthyMatrix)))
  ibd.mediantotalcount=c(ibd.mediantotalcount,median(colSums(ibdMatrix)))
  
  # compute evenness
  healthy.values=apply(healthyMatrix,2,sheldon)
  ibd.values=apply(ibdMatrix,2,sheldon)
  print(paste("mean Sheldon healthy",mean(healthy.values)))
  print(paste("mean Sheldon IBD",mean(ibd.values)))
  p <- wilcox.test(healthy.values,ibd.values)$p.value
  print(paste("p-value Wilcoxon",p))
  
  # compute beta diversity
  healthy.bc = as.matrix(vegdist(t(healthyMatrix), method = "bray"))
  healthy.bc.values=healthy.bc[upper.tri(healthy.bc)]
  print(paste("Sample-wise BC distribution healthy mean:",mean(healthy.bc.values),"sd:",sd(healthy.bc.values)))
  
  ibd.bc = as.matrix(vegdist(t(ibdMatrix), method = "bray"))
  ibd.bc.values=ibd.bc[upper.tri(ibd.bc)]
  print(paste("Sample-wise BC distribution IBD mean:",mean(ibd.bc.values),"sd:",sd(ibd.bc.values)))
  
  healthy.sheldon=c(healthy.sheldon,mean(healthy.values))
  ibd.sheldon=c(ibd.sheldon,mean(ibd.values))
  healthy.betadiv=c(healthy.betadiv,mean(healthy.bc.values))
  ibd.betadiv=c(ibd.betadiv,mean(ibd.bc.values))
  
  # apply a prevalence filter
  if(minocc>0){
    # sum of filtered taxa is kept, so as not to change total counts
    healthyMatrix=prevalenceFilter(healthyMatrix,minocc=minocc,keepSum = TRUE)
    ibdMatrix=prevalenceFilter(ibdMatrix,minocc=minocc,keepSum = TRUE)
  }
  
  # compute networks
  # SPIECEASI needs lineages, assign fake ones (also to garbage taxon if sum of filtered taxa is kept)
  data(ibd_lineages)
  rownames(healthyMatrix)=rownames(ibd_lineages)[1:nrow(healthyMatrix)]
  rownames(ibdMatrix)=rownames(ibd_lineages)[1:nrow(ibdMatrix)]
  
  N=nrow(healthyMatrix)
  if(method=="beem"){
    beemOK=TRUE
    res <- tryCatch(
      func.EM(healthyMatrix),
      error = function(e) {
        beemOK<<-FALSE
        message("BEEM failed: ", conditionMessage(e))
        return(NULL)
      }
    )
    if(beemOK){
      est <- beem2param(res)
      pep.healthy=getPep(est$b.est)
    }else{
      pep.healthy=NA
    }
    edgeNum=NA
    numNeg=NA
  }else{
    healthyNetwork=buildNetwork(healthyMatrix,method=method,lineages=ibd_lineages[1:nrow(healthyMatrix),],initEdgeNum = N*(N-1)/4, legend=FALSE)
    edgeNum=length(E(healthyNetwork))
    numNeg=length(which(E(healthyNetwork)$color=="red"))
    numPos=edgeNum-numNeg
    pep.healthy=numPos/(edgeNum/100)
  }
  print(paste("Healthy network has PEP",pep.healthy,",",numNeg," negative edges"," and ",edgeNum,"total edges"))
  healthy.edgenum=c(healthy.edgenum,edgeNum)
  healthy.pep=c(healthy.pep,pep.healthy)
  
  N=nrow(ibdMatrix)
  if(method=="beem"){
    beemOK=TRUE
    res <- tryCatch(
      func.EM(ibdMatrix),
      error = function(e) {
        beemOK<<-FALSE
        message("BEEM failed: ", conditionMessage(e))
        return(NULL)
      }
    )
    if(beemOK){
      est <- beem2param(res)
      pep.ibd=getPep(est$b.est)
    }else{
      pep.ibd=NA
    }
    edgeNum=NA
    numNeg=NA
  }else{
    ibdNetwork=buildNetwork(ibdMatrix,method=method,lineages=ibd_lineages[1:nrow(ibdMatrix),],initEdgeNum = N*(N-1)/4, legend=FALSE)
    edgeNum=length(E(ibdNetwork))
    numNeg=length(which(E(ibdNetwork)$color=="red"))
    numPos=edgeNum-numNeg
    pep.ibd=numPos/(edgeNum/100)
  }
  print(paste("IBD network has PEP",pep.ibd,",",numNeg," negative edges"," and ",edgeNum,"total edges"))
  ibd.edgenum=c(ibd.edgenum,edgeNum)
  ibd.pep=c(ibd.pep,pep.ibd)  
  
} # end rounds 

print("Results")
print("PEP in simulations:")
print(simulated.pep)
print("PEP in healthy networks:")
print(healthy.pep)
print("PEP in IBD networks:")
print(ibd.pep)

if(rounds > 1){
  # PEP difference
  box.df <- data.frame(
    PEP = c(healthy.pep, ibd.pep),
    Group = factor(c(rep("Healthy", length(healthy.pep)),
                     rep("IBD", length(ibd.pep))),
                   levels = c("Healthy", "IBD"))
  )
  g.pep<-ggplot(box.df, aes(x = Group, y = PEP, fill = Group)) +
    geom_violin(alpha = 0.4, trim = FALSE) +
    geom_jitter(aes(color = Group), width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = c("Healthy" = "#0078B9", "IBD" = "#EA0017")) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "IBD" = "#EA0017")) +
    annotate("text",
             x = 1.5, y = max(box.df$PEP, na.rm = TRUE),
             label = paste0("p = ", formatC(p.val, format = "e", digits = 1)),
             size = 4, vjust = -0.5) +
    labs(title = "Positive edge percentage in simulated data",
         y = "Positive edge percentage", x = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
  plot(g.pep)
  
  # Correlation of PEP to evenness
  healthy.res=cbind(healthy.sheldon,healthy.pep,rep("healthy",length(healthy.pep)))
  ibd.res=cbind(ibd.sheldon,ibd.pep,rep("IBD",length(ibd.pep)))
  res.df=as.data.frame(rbind(healthy.res,ibd.res))
  colnames(res.df)=c("Sheldon","PEP","Group")
  res.df$Sheldon <- as.numeric(res.df$Sheldon)
  res.df$PEP <- as.numeric(res.df$PEP)
  
  # compute per-group stats
  stats.df <- res.df %>%
    group_by(Group) %>%
    summarise(
      r2 = summary(lm(PEP ~ Sheldon))$r.squared,
      pval = summary(lm(PEP ~ Sheldon))$coefficients[2, 4],
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("R² = ", round(r2, 3), ", p = ", signif(pval, 3))
    )
  
  overall.fit <- lm(PEP ~ Sheldon, data = res.df)
  overall.label <- paste0(
    "Global: R² = ", round(summary(overall.fit)$r.squared, 3),
    ", p = ", signif(summary(overall.fit)$coefficients[2, 4], 3)
  )
  
  gout<-ggplot(res.df, aes(x = Sheldon, y = PEP, color = Group, fill = Group)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_smooth(aes(x = Sheldon, y = PEP), method = "lm",
                se = TRUE, inherit.aes = FALSE,
                color = "black", linetype = "dashed") +
    scale_color_manual(values = c("healthy" = "#0078B9", "IBD" = "#EA0017")) +
    scale_fill_manual(values = c("healthy" = "#0078B9", "IBD" = "#EA0017")) +
    ggtitle("Evenness versus positive edge percentage")+
    geom_text(
      data = stats.df,
      aes(label = label, color = Group),
      x = Inf, y = Inf,
      hjust = 1.1,
      vjust = c(1.5, 3.0),
      size = 4,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.1, vjust = 4.5,
             label = overall.label,
             size = 4, color = "black") +
    labs(x = "Sheldon evenness", y = "PEP") +
    theme_minimal()
  plot(gout)
  
  # diagnostic plots for data properties
  box.df <- data.frame(
    BC = c(healthy.betadiv, ibd.betadiv),
    Group = factor(c(rep("Healthy", length(healthy.betadiv)),
                     rep("IBD", length(ibd.betadiv))),
                   levels = c("Healthy", "IBD"))
  )
  g.bc<-ggplot(box.df, aes(x = Group, y = BC, fill = Group)) +
    geom_violin(alpha = 0.4, trim = FALSE) +
    geom_jitter(aes(color = Group), width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = c("Healthy" = "#0078B9", "IBD" = "#EA0017")) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "IBD" = "#EA0017")) +
    annotate("text",
             x = 1.5, y = max(box.df$BC, na.rm = TRUE),
             label = paste0("p = ", formatC(p.val, format = "e", digits = 1)),
             size = 4, vjust = -0.5) +
    labs(title = "Beta diversity in simulated data",
         y = "Mean Bray Curtis", x = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
  plot(g.bc)
  
  box.df <- data.frame(
    Sheldon = c(healthy.sheldon, ibd.sheldon),
    Group = factor(c(rep("Healthy", length(healthy.sheldon)),
                     rep("IBD", length(ibd.sheldon))),
                   levels = c("Healthy", "IBD"))
  )
  g.sheldon<-ggplot(box.df, aes(x = Group, y = Sheldon, fill = Group)) +
    geom_violin(alpha = 0.4, trim = FALSE) +
    geom_jitter(aes(color = Group), width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = c("Healthy" = "#0078B9", "IBD" = "#EA0017")) +
    scale_color_manual(values = c("Healthy" = "#0078B9", "IBD" = "#EA0017")) +
    annotate("text",
             x = 1.5, y = max(box.df$Sheldon, na.rm = TRUE),
             label = paste0("p = ", formatC(p.val, format = "e", digits = 1)),
             size = 4, vjust = -0.5) +
    labs(title = "Evenness in simulated data",
         y = "Sheldon evenness", x = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
  plot(g.sheldon)
  
  # library(patchwork)
  # combined <- (g.pep + gout) / (g.sheldon + g.bc) 
  # ggsave("combined.pdf",combined,wdith=12,height=10)
  
  boxplot(list(Healthy = healthy.mediantotalcount, IBD = ibd.mediantotalcount), col = c("#0078B9", "#EA0017"),main="Total abundance in simulated data after rarefaction",ylab="Abundance")
}
