# Qst - Fst calculations for unbalanced breeding designs 
#	by Kimberly J. Gilbert (kgilbert@zoology.ubc.ca) and Michael C. Whitlock


#  Based on Whitlock and Guillaume (2009)
#  Code edited and expanded to function with AFLP data and with additional breeding designs




#'
#' Compares the QST of a single phenotypic trait to the mean FST of series of marker loci. 
#' It calculates the distribution of QST - FST under a model assuming neutrality of both 
#' the phenotypic trait and the genetic markers from which FST is estimated.
#' Returns the simulated estimates of Qst - Fst under neutrality following the procedure
#' described in Gilbert and Whitlock (2014) and Whitlock & Guillaume (2009). Also returns 
#' the simulated estimates of Fst and Qst used to compute the null distribution.
#'
#' @title Compare Qst to Fst
#'
#' @param fst.dat  A data frame with the first column indicating population of origin and
#'          the following columns representing genotypes at loci; see the 
#'          README \url{https://github.com/kjgilbert/QstFstComp/blob/master/README.md} for further description.
#'          If using AFLPs, this is a data frame of q-hat values, with pops in columns, loci in rows
#'          and the corresponding q-hat variances in the following columns, and \code{AFLP=TRUE} must be designated.
#'
#' @param qst.dat  the input table of the breeding data
#'  \itemize{
#'      \item    For \code{breeding.design = "half.sib.sire"}:  qst.dat should have four
#'           columns in this order: population, sire, dam, and the trait value of
#'           the individual. Each population, sire, and dam should have unique
#'           names or numbers.
#'      \item    For \code{breeding.design = "half.sib.dam"}:  qst.dat should have three
#'           columns in this order: population, dam, and the trait value of
#'           the individual. Each population and dam should have unique
#'           names or numbers.
#'  }
#'
#' @param numpops  number of populations in the sample
#'
#' @param nsim     number of simulation replicates to perform to create the null
#'          distributions and bootstraps
#' 
#'  @param AFLP whether or not to use AFLP data
#'  @param breeding.design	the breeding design used when collecting the trait data
#'  OPTIONS for breeding design:
#'  \enumerate{
#'  \item 	"half.sib.sire" is a half sib design with dam nested within sire nested within
#'  		population which works for either balanced or unbalanced sampling designs
#'  \item 	"half.sib.dam" is a half sib design with dam nested within population
#'  		which works for either balanced or unbalanced sampling designs
#'  }
#'
#'  @param dam.offspring.relatedness relatedness between offspring in the dam model, default is 1/4, i.e. half-sib
#'
#'  @param output whether to output full or concise results, see details below
#'
#'  @return
#'
#' Returns either a concise list of a subset of results or a full list with all possible results:
#'
#'  Concise list returns (default)
#'  \itemize{
#'  \item 	the calculated difference between Qst and Fst with 95\% critical values,
#'  \item 	one- and two- tailed p-values for this difference, and
#'  \item 	the additive genetic variance for the trait with 95\% confidence intervals
#'  }
#'  Full list returns 
#'  \itemize{
#'  		\item  the calculated difference between Qst and Fst with 95\% critical values,
#'  		\item  one- and two- tailed p-values for this difference,
#'  		\item  the Fst as calculated from the genetic data provided, with 95\% confidence intervals,
#'  		\item  the resampled Fst as calculated from bootstrapping across simulations, with standard deviation and 95\% confidence intervals,
#'  		\item  the resampled neutral Qst as calculated from bootstrapping across simulations, with standard deviation and 95\% critical values,
#'  		\item  the ANOVA table for the calculated means squared, n coefficients, and degrees of freedom,
#'  		\item  the additive genetic variance for the trait with 95\% confidence intervals, and
#'  		\item  the coefficient of additive genetic variance for the trait with 95\% confidence intervals
#'  }
#'
#' @author Kimberly J Gilbert & Michael C Whitlock
#'
#' @import hierfstat
#'
#' @references Gilbert KJ and MC Whitlock (\emph{In prep.}) \emph{Qst} \emph{Fst} comparisons with half-sib designs.
#'
#' @references Whitlock MC and F Guillaume (2009) \href{http://www.genetics.org/content/183/3/1055}{Testing for spatially divergent selection: Comparing \emph{Qst} to \emph{Fst}.} \emph{Genetics}, 183:1055-1063.
#'
#'
#'
#' @examples
#'
#' ## using balanced half-sib sire trait data and biallelic marker data 
#' data(hssire) # trait data
#' data(biallelic) # marker data
#' QstFstComp(biallelic, hssire, numpops=15, nsim=100, breeding.design="half.sib.sire", output="full")
#' 
#' data(hsdam)
#' data(aflp)
#' QstFstComp(aflp, hsdam, numpops=15, nsim=100, AFLP=TRUE, breeding.design="half.sib.dam", output="concise")
#' 
#' 
#' @export




QstFstComp <- function(fst.dat, qst.dat, numpops, nsim=1000, AFLP=FALSE, breeding.design, dam.offspring.relatedness=0.25, output="concise")
{	
  if(missing(fst.dat)) stop("Genotypic data must be provided.")
  if(missing(qst.dat)) stop("Phenotypic data must be provided.")
  if(missing(numpops)) stop("Number of populations must be defined.")
  if(missing(breeding.design)) stop("Breeding design must be defined.")

  	if(breeding.design=="half.sib.sire"){
  		qst.MS <- MeanSq.unbalanced.sire(qst.dat)
  		qst.temp <- QSTfromSireModel(qst.MS$MSpops, qst.MS$MSsires, qst.MS$MSdams, qst.MS$MSwithin, 
  				qst.MS$n0primeprime, qst.MS$n0prime, qst.MS$n0, qst.MS$nc0prime, qst.MS$nc0, qst.MS$ncb0)
   		qst.obs <- qst.temp$Qst
   		mean.trait.value <- mean(qst.dat[,4], na.rm=TRUE) # this takes the mean of all trait values across all ind.s, no bootstrapping
		Va <- qst.temp$Va
		if(Va < 0){Va <- 0}	# if Va includes negative values, the distribution is truncated to zero
      	CVa <- sqrt(Va)/mean.trait.value * 100
  	}
  	if(breeding.design=="half.sib.dam"){
  		qst.MS <- MeanSq.unbalanced.dam(qst.dat)
  		qst.obs <- QSTfromDamModel(qst.MS$MSpops,qst.MS$MSdams,qst.MS$MSwithin,qst.MS$n0prime,qst.MS$n0,qst.MS$nb0,
  			dam.offspring.relatedness)
   		mean.trait.value <- mean(qst.dat[,3], na.rm=TRUE) # this takes the mean of all trait values across all ind.s, no bootstrapping  
		Va <- 1/dam.offspring.relatedness*(qst.MS$MSdams-qst.MS$MSwithin)/qst.MS$n0
		if(Va < 0){Va <- 0}	# if Va includes negative values, the distribution is truncated to zero
      	CVa <- sqrt(Va)/mean.trait.value * 100
  	}

 # for non-AFLP markers, i.e. non-dominant markers  
  if(AFLP==FALSE){
  	nloci <- ncol(fst.dat)-1
  	abc.mat <- wc.calc(fst.dat, diploid=TRUE)
  	nalleles <- nrow(abc.mat)
	#the observed Fst:
	fst.obs <- sum(abc.mat[,1])/sum(abc.mat[,1] + abc.mat[,2] + abc.mat[,3])
  }

 # for AFLP markers, i.e. dominant markers  
  if(AFLP==TRUE){
  	# get the data for the Fst calculation and resampling
 	full.fst.dat <- read.fst.input(fst.dat, num.pops=numpops)
	fst.data <- as.matrix(full.fst.dat [[1]])	# this is a matrix of the q_hat values
 	var.dat <- as.matrix(full.fst.dat [[2]])	# this is a matrix of the variances
	nloci <- nrow(fst.data)
	#the observed Fst:
	fst.obs <- mean.fst(dat=fst.data, vardat=var.dat, num.loci=nloci, num.pops=numpops)
  }

  # values to be filled in across replicate simulations
  sim.est <- vector(length = nsim)
  qst.neut  <- vector(length = nsim)
  fst.est  <- vector(length = nsim)
  qstForCI.est  <- vector(length = nsim)
  Va.est  <- vector(length = nsim)
  CVa.est <- vector(length = nsim)

  #simulating
  for(i in 1:nsim) {
   #1. Fst simulation replicate; sample nloci from the neutral markers, with replacement
    if(AFLP==FALSE){ fst.repl <- fst.sample(fst.dat, nalleles) }
	if(AFLP==TRUE){	fst.repl <- fst.sample.aflp(fst.data, var.dat, nloci) }
   #2. get a simulated replicate of Qst by sampling the null distribution 					 
  	if(breeding.design=="half.sib.sire") qst.repl <- qst.parboot.siremodel(qst.MS, fst.obs)[[1]] # this function now spits out Qst and Va, just need the Qst here
  	if(breeding.design=="half.sib.dam") qst.repl <- qst.parboot.dammodel(qst.MS, fst.obs)
   #3. get the simulated replicate of Qst - Fst, and store it
    sim.est[i] <- qst.repl - fst.repl  
    #store values
    fst.est[i] <- fst.repl
    qst.neut[i] <- qst.repl

    #bootstraps for Qst and Va
 	 if(breeding.design=="half.sib.sire") { temp <- qstVa.parbootForCI.siremodel(qst.MS) }
	 if(breeding.design=="half.sib.dam") { temp <- qstVa.parbootForCI.dammodel(qst.MS) }
     qstForCI.est[i] <- temp[1]
     Va.est[i] <- temp[2]
     if(Va.est[i] < 0){Va.est[i] <- 0}	# if Va includes negative values, the distribution is truncated to zero
     CVa.est[i] <- sqrt(Va.est[i])/mean.trait.value * 100
  }
  

  # calculate the p-value for the difference between the neutral qst and fst
  diff.repl <- qst.neut - fst.est
  Q.obsMinusF.obs <- qst.obs-fst.obs

	right.one.tailed.p <- sum(Q.obsMinusF.obs < diff.repl)/length(diff.repl)  # right hand 1-tailed p value
	left.one.tailed.p <- sum(Q.obsMinusF.obs > diff.repl)/length(diff.repl)	# left hand 1-tailed p value
	two.tailed.p <- 2*min(right.one.tailed.p, left.one.tailed.p)

## Create return items for either the concise or the full list of outputs, as specified from input parameters
  if(output=="concise"){return(list(
		QminusF	<- c(
				"Calculated Qst-Fst" = Q.obsMinusF.obs, 
				"Lower Bound crit. value" = quantile(sim.est,0.025,na.rm=TRUE), 
				"Upper bound crit. value" = quantile(sim.est,0.975,na.rm=TRUE)),
	    QminusF.p.values <- c( 
    			"Lower one-tailed p value" = left.one.tailed.p, 
       			"Upper one-tailed p value" = right.one.tailed.p, 
    			"Two-tailed p value" = two.tailed.p),
    	Va  <- c("Va"=Va, 
    			"Lower bound CI" = quantile(Va.est,0.025,na.rm=TRUE), 
    			"Upper bound CI" = quantile(Va.est,0.975,na.rm=TRUE))
	))
  } # end return of concise list
  if(output=="full"){
  	if(breeding.design=="half.sib.sire"){ # this full output has the ANOVA table for the half-sib sire breeding design
  		return(list(
		QminusF	<- c(
				"Calculated Qst-Fst" = Q.obsMinusF.obs, 
				"Lower Bound crit. value" = quantile(sim.est,0.025,na.rm=TRUE), 
				"Upper bound crit. value" = quantile(sim.est,0.975,na.rm=TRUE)),
	    QminusF.p.values <- c(  
    			"Lower one-tailed p value" = left.one.tailed.p, 
       			"Upper one-tailed p value" = right.one.tailed.p, 
    			"Two-tailed p value" = two.tailed.p),
		Fst	<- c("Estimated Fst" = fst.obs, 
				"Lower Bound CI" = quantile(fst.est,0.025,na.rm=TRUE), 
				"Upper bound CI" = quantile(fst.est,0.975,na.rm=TRUE)),
    	Fst.resampled <- c(
    			"Fst Resampled" = mean(fst.est, na.rm=TRUE), 
    			"Fst std. dev." = sd(fst.est, na.rm=TRUE), 
    			"Fst 95% CI" = quantile(fst.est,c(0.025,0.975), na.rm=TRUE), 
    			"Fst 99% CI" = quantile(fst.est,c(0.005,0.995), na.rm=TRUE)),
		Qst.neutral	<- c(
				"Qst Resampled" = mean(qst.neut, na.rm=TRUE), 
				"Qst std. dev." = sd(qst.neut, na.rm=TRUE), 
				"Qst 95% crit. values" = quantile(qst.neut,c(0.025,0.975), na.rm=TRUE), 
				"Qst 99% crit. values" = quantile(qst.neut,c(0.005,0.995), na.rm=TRUE)),
    	ANOVA.table	<- c(
					"MS Pop" = qst.MS[[1]],
					"MS Sire" = qst.MS[[2]],
					"MS Dam" = qst.MS[[3]],
					"MS Within" = qst.MS[[4]],
					"n0primeprime" = qst.MS[[5]],
					"n0prime" = qst.MS[[6]],
					"n0" = qst.MS[[7]],
					"nc0prime" = qst.MS[[8]],
					"nc0" = qst.MS[[9]],
					"ncb0" = qst.MS[[10]],
					"df Pop" = qst.MS[[11]],
					"df Sire" = qst.MS[[12]],
					"df Dam" = qst.MS[[13]],
					"df Within" = qst.MS[[14]],
					"Sigma^2 sires" = qst.MS[[15]]),				
		Va <- c("Va" = Va,
				"Lower bound CI" = quantile(Va.est,0.025,na.rm=TRUE),
				"Upper bound CI" = quantile(Va.est,0.975,na.rm=TRUE)),
    	CVa <- c("CVa" <- CVa,
    			"Lower bound CI" = quantile(CVa.est,0.025,na.rm=TRUE),
    			"Upper bound CI" = quantile(CVa.est,0.975,na.rm=TRUE))
	))} # end return of half.sib.sire full output
	  	if(breeding.design=="half.sib.dam"){ # this full output has the ANOVA table for the half-sib dam breeding design
  		return(list(
		QminusF	<- c(
				"Calculated Qst-Fst" = Q.obsMinusF.obs, 
				"Lower Bound crit. value" = quantile(sim.est,0.025,na.rm=TRUE), 
				"Upper bound crit. value" = quantile(sim.est,0.975,na.rm=TRUE)),
	    QminusF.p.values <- c(  
    			"Lower one-tailed p value" = left.one.tailed.p, 
       			"Upper one-tailed p value" = right.one.tailed.p, 
    			"Two-tailed p value" = two.tailed.p),
		Fst	<- c("Estimated Fst" = fst.obs, 
				"Lower Bound CI" = quantile(fst.est,0.025,na.rm=TRUE), 
				"Upper bound CI" = quantile(fst.est,0.975,na.rm=TRUE)),
    	Fst.resampled <- c(
    			"Fst Resampled" = mean(fst.est, na.rm=TRUE), 
    			"Fst std. dev." = sd(fst.est, na.rm=TRUE), 
    			"Fst 95% CI" = quantile(fst.est,c(0.025,0.975), na.rm=TRUE), 
    			"Fst 99% CI" = quantile(fst.est,c(0.005,0.995), na.rm=TRUE)),
		Qst.neutral	<- c(
				"Qst Resampled" = mean(qst.neut, na.rm=TRUE), 
				"Qst std. dev." = sd(qst.neut, na.rm=TRUE), 
				"Qst 95% crit. values" = quantile(qst.neut,c(0.025,0.975), na.rm=TRUE), 
				"Qst 99% crit. values" = quantile(qst.neut,c(0.005,0.995), na.rm=TRUE)),
    	ANOVA.table	<- c(
					"MS Pop" = qst.MS[[1]],
					"MS Dam" = qst.MS[[2]],
					"MS Within" = qst.MS[[3]],
					"n0prime" = qst.MS[[4]],
					"n0" = qst.MS[[5]],
					"nb0" = qst.MS[[6]],
					"df Pop" = qst.MS[[7]],
					"df Dam" = qst.MS[[8]],
					"df Within" = qst.MS[[9]],
					"Sigma^2 dams" = qst.MS[[10]]),				
		Va <- c("Va" = Va,
				"Lower bound CI" = quantile(Va.est,0.025,na.rm=TRUE),
				"Upper bound CI" = quantile(Va.est,0.975,na.rm=TRUE)),
    	CVa <- c("CVa" = CVa,
    			"Lower bound CI" = quantile(CVa.est,0.025,na.rm=TRUE),
    			"Upper bound CI" = quantile(CVa.est,0.975,na.rm=TRUE))
	))} # end return of half.sib.dam full output
  } # end return of full list
} # end main function
###########################################################################





###########################################################################
#######                                                             #######
#######    additional functions required for main function above    #######
#######                                                             #######
###########################################################################






################       MeanSq unbalanced sire model      #####################
#
# Finding no and noprime values, as per Sokal and Rohlf p297, and also returning
#	key elements of the anova, including MS. 
# This will be expecting to receive a data frame with three columns of numbers:
#	the first column labelled 'Population',
#	the second giving dam names,
#	the third giving sire names, and
#	the fourth column is the trait values, with NA for any missing values.
# ** Sires and dams should be labelled distinctly (i.e. the same name refers to the same
#	individual dam, without forcing reference to the population name)

MeanSq.unbalanced.sire <- function(dat){
	poplist <- dat[,1]
	sirelist <- dat[,2]
	damlist <- dat[,3]
	traitlist <- dat[,4]
	
  #Remove all lines for which there is no data for the trait, to simplify calculations later.
  pops <- poplist[!is.na(traitlist)]
  sires <- sirelist[!is.na(traitlist)]
  dams <- damlist[!is.na(traitlist)]
  sires <- sirelist[!is.na(traitlist)]
  nTable <- table(as.data.frame(cbind(pops,sires,dams)))
  
  traits <- as.numeric(traitlist[!is.na(traitlist)])
  dataTable <- cbind(as.data.frame(cbind(pops,sires,dams)),traits=I(traits))
  
  nPops <- length(table(poplist))
  nSiresTot <-length(table(sirelist))
  nDamsTot <- length(table(damlist))
  dfsires <- (nSiresTot-nPops)
  dfpops <- nPops-1
  dfdams <- nDamsTot-nSiresTot
  dfwithin <- length(traits)-dfsires-dfdams-dfpops-1  
  
  ##Calculating the n values (based on Sokal and Rohlf pp302ff)
  
  #Quantity 1
  q1 <- sum(nTable^2)
  
  #Quantity 2
  nTablePopSire <- table(as.data.frame(cbind(pops,sires)))
  q2 <- sum(nTablePopSire^2)
  
  #Quantity 3
  nTablePop <- table(as.data.frame(cbind(pops)))
  q3 <- sum(nTablePop^2)
  
  #Quantity 4 
  nTot <- length(traits)
  q4 <- nTot
  
  #Quantity 5
  nTableSireDams <- table(as.data.frame(cbind(sires,dams)))
  q5 <-sum(rowSums(nTableSireDams^2)/rowSums(nTableSireDams))

  #Quantity 6
  nTablePopsDams <- table(as.data.frame(cbind(pops,dams)))
  q6 <-sum(rowSums(nTablePopsDams^2)/rowSums(nTablePopsDams))
  
  #Quantity 7
  q7 <- sum(rowSums(nTablePopSire^2)/rowSums(nTablePopSire))
  
  
  ###Calculating n's
  
  n0primeprime <- (q6-q1/q4)/dfpops
  
  n0prime <- (q5-q6)/dfsires
  
  n0 <- (q4-q5)/dfdams
  
  nc0prime <- (q7 - (q2/q4))/dfpops
  
  nc0 <- (q4-q7)/dfsires
  
  ncb0 <- (q4-(q3/q4))/dfpops
  
  
  #####Calculating SS
  
  grandTotal <- sum(traits)
  sumSqObs <- sum(traits^2)
  
  #Sires
	siresq <- tapply(dataTable$traits,dataTable$sires,FUN=sum)^2
	siren <- table(dataTable$sires)
	qsire <- sum(siresq/siren)

  #Dams
	damsq <- tapply(dataTable$traits,dataTable$dams,FUN=sum)^2
	damn <- table(dataTable$dams)
	qdam <- sum(damsq/damn)

  #Pops
	popsq <- tapply(dataTable$traits,dataTable$pops,sum)^2
	popn <- table(dataTable$pops)
	qpop <- sum(popsq/popn)

  #CT
	CT <- grandTotal^2/length(dataTable[,1])

  #SS
	SStotal <- sumSqObs-CT
	SSpops <- qpop-CT
	SSsires <- qsire-qpop
	SSdams <- qdam-qsire
	SSwithin <- sumSqObs-qdam

  #MS
	MSpops <- SSpops/dfpops
	MSsires <- SSsires/dfsires
	MSdams <- SSdams/dfdams
	MSwithin <- SSwithin/dfwithin

	sigma2dams <- (MSdams-MSwithin)/n0
	sigma2sires <- (MSsires - MSwithin - n0prime*sigma2dams)/nc0


return(list(
  		MSpops = MSpops,
  		MSsires = MSsires,
  		MSdams = MSdams,
  		MSwithin = MSwithin,
  		n0primeprime = n0primeprime,
  		n0prime = n0prime,
  		n0 = n0,
  		nc0prime = nc0prime,
  		nc0 = nc0,
  		ncb0 = ncb0,
  		dfpops = dfpops,
  		dfsires = dfsires,
  		dfdams = dfdams,
  		dfwithin = dfwithin,
  		sigma2sires = sigma2sires,
  		sigma2dams = sigma2dams
  	))
}
###########################################################################




################       MeanSq unbalanced dam model      ##################### 
#    
# Finding n0 and n0prime values, as per Sokal and Rohlf p297, and also returning
#	key elements of the anova, including MS. 
# This will be expecting to receive a data frame with three columns of numbers:
#	the first column labelled 'Population',
#	the second giving dam names, and 
#	the third column is the trait values, with NA for any missing values.
# ** Dams should be labelled distinctly (i.e. the same name refers to the same
#	individual dam, without forcing reference to the population name)
MeanSq.unbalanced.dam <- function(dat){
	poplist <- dat[,1]
	damlist <- dat[,2]
	traitlist <- dat[,3]
	  
  #Remove all lines for which there is no data for the trait, to simplify calculations later.
  pops <- poplist[!is.na(traitlist)]
  dams <- damlist[!is.na(traitlist)]
  nTable <- table(as.data.frame(cbind(pops,dams)))
  
  traits <- as.numeric(traitlist[!is.na(traitlist)])
  dataTable <- cbind(as.data.frame(cbind(pops,dams)),traits=I(traits))
  
  nPops <- length(table(poplist))
  nDamsTot <- length(table(damlist))
  dfdam <- (nDamsTot-nPops)
  dfpops <- nPops-1
  dfwithin <- length(traits)-dfdam-dfpops-1      
 
  ##Calculating the n values
  #Quantity 1
  nTot <- length(traits)
  #Quantity 2
  nsqTot <- sum(nTable^2) 
  #Quantity 3
  sumnsq <- sum(rowSums(nTable)^2)
  #Quantity 4
  eq4 <- sum(rowSums(nTable^2)/rowSums(nTable))
 
  n0prime <- (eq4-nsqTot/nTot)/(nPops-1)
  n0 <- (nTot - eq4)/dfdam
  nb0 <- (nTot-sumnsq/nTot)/(nPops-1)
 
  #Calculating SS
  grandTotal <- sum(traits)
  sumSqObs <- sum(traits^2)
  
  #Quantity 3 for SS
  damsq <- tapply(dataTable$traits,dataTable$dams,FUN=sum)^2
  damn <- table(dataTable$dams)
  q3 <- sum(damsq/damn)
  
  #Quantity 4 for SS
  popsq <- tapply(dataTable$traits,dataTable$pops,sum)^2
  popn <- table(dataTable$pops)
  q4 <- sum(popsq/popn)
  
  #CT
  CT <- grandTotal^2/length(dataTable[,1])
  
  #SS
  SStotal <- sumSqObs-CT
  SSpops <- q4-CT
  SSdams <- q3-q4
  SSwithin <- sumSqObs-q3
  
  #MS
  MSpops <- SSpops/dfpops
  MSdams <- SSdams/dfdam
  MSwithin <- SSwithin/dfwithin
 
  return(list(
  		MSpops = MSpops,
  		MSdams = MSdams,
  		MSwithin = MSwithin,
  		n0prime = n0prime,
  		n0 = n0,
  		nb0 = nb0,
  		dfpops = dfpops,
  		dfdam = dfdam,
  		dfwithin = dfwithin,
  		sigma2dams = (MSdams-MSwithin)/n0
  	))
}
###########################################################################




########################     Qst.sireModel     #############################
#
# Calculating Qst from the sire halfsib model, with arbitrary numbers of dams per sire and offspring per dam

QSTfromSireModel <- function(MSpops, MSsires, MSdams, MSwithin, n0primeprime, n0prime, n0, nc0prime, nc0, ncb0){
  sigma2dams <- (MSdams-MSwithin)/n0
  sigma2sires <- (MSsires - MSwithin - n0prime*sigma2dams)/nc0
  sigma2pops <- (MSpops - MSwithin - n0primeprime*sigma2dams - nc0prime*sigma2sires)/ncb0
  
  Va <- 4*sigma2sires
  
  Qst <- sigma2pops/(sigma2pops+2*Va)
  
  return(list(Qst = Qst, Va = Va))
}	
###########################################################################




########################     Qst.damModel     #############################
#
# Calculating Qst from the dam halfsib model, assuming that the offpsring of dams are all by unique fathers

QSTfromDamModel <- function(MSpops,MSdams,MSwithin,n0prime,n0,nb0, dam.offspring.relatedness){
  sigma2dams <- (MSdams-MSwithin)/n0
  sigma2pops <- (MSpops-sigma2dams*n0prime-MSwithin)/nb0
  VA <- 1/dam.offspring.relatedness*sigma2dams
  
  Qst <- sigma2pops/(sigma2pops+2*VA)
  
  return(Qst)
}	
###########################################################################




##############################		Qst parboot sire model		######################################
#
# New qst.parboot for unbalanced model (collapses to balanced model as well) returns a Qst pseudovalue for that simulated null (neutral)
# distribution based on half-sib sire model. This requires MS for sires, dams and residuals,
# along with the various n's from the nested ANOVA and dfs for each level (These are all encapsulated in the results list from MeanSq.unbalanced.sire.)
	
qst.parboot.siremodel <- function(MSdflist, meanFst){
  MSsiresResample <- MSdflist$MSsires * rchisq(1,MSdflist$dfsires) / MSdflist$dfsires
  MSdamsResample <- MSdflist$MSdams * rchisq(1,MSdflist$dfdam) / MSdflist$dfdam
  MSwithinResample <- MSdflist$MSwithin * rchisq(1,MSdflist$dfwithin) / MSdflist$dfwithin
 
  VpopNeut <- 8*meanFst / (1-meanFst) * MSdflist$sigma2sires
  
  MSpopNeut <- MSdflist$MSwithin + MSdflist$n0primeprime * MSdflist$sigma2dams + MSdflist$nc0prime * MSdflist$sigma2sires + MSdflist$ncb0 * VpopNeut
 
  MSpopNeutResample <- MSpopNeut * rchisq(1,MSdflist$dfpops) / MSdflist$dfpops

  sim.qst <- QSTfromSireModel(MSpopNeutResample,MSsiresResample,MSdamsResample,MSwithinResample,MSdflist$n0primeprime,MSdflist$n0prime,MSdflist$n0,MSdflist$nc0prime,MSdflist$nc0,MSdflist$ncb0)  
  return(sim.qst)
}
###########################################################################




##############################		Qst parboot dam model		######################################
# New qst.parboot for dam model returns a Qst pseudovalue for that simulated null
# distribution based on half-sib dam model. This requires MS for dams and residuals,
# along with the n0 value from the nested ANOVA and dfs for each level	

qst.parboot.dammodel <- function(MSdflist, meanFst, dam.offspring.relatedness){
  MSdamsResample <- MSdflist$MSdams * rchisq(1,MSdflist$dfdam) / MSdflist$dfdam
  MSwithinResample <- MSdflist$MSwithin * rchisq(1,MSdflist$dfwithin) / MSdflist$dfwithin
 
  VpopNeut <- 2/dam.offspring.relatedness*meanFst / (1-meanFst) * MSdflist$sigma2dams
  MSpopNeut <- MSdflist$MSwithin + MSdflist$n0prime * MSdflist$sigma2dams + MSdflist$nb0 * VpopNeut
 
  MSpopNeutResample <- MSpopNeut * rchisq(1,MSdflist$dfpops) / MSdflist$dfpops

  sim.qst <- QSTfromDamModel(MSpopNeutResample,MSdamsResample,MSwithinResample,MSdflist$n0prime,MSdflist$n0,MSdflist$nb0)  
  return(sim.qst)
}
###########################################################################




#########################       QstVa Parboot for sire model      ###############################
#
#New qstVa.parbootForCI.siremodel for sire model returns a Qst pseudovalue
#resampled by a paremetric bootstrap from the sire model. This follows the logic
#of O'Hara and Merila. It also returns a parametric bootrap value for Va. This
#requires MS for sires, dams and residuals, along with the n's and dfs for each level

qstVa.parbootForCI.siremodel <- function(MSdflist){
  MSsiresResample <- MSdflist$MSsires * rchisq(1,MSdflist$dfsires) / MSdflist$dfsires
  MSdamsResample <- MSdflist$MSdams * rchisq(1,MSdflist$dfdam) / MSdflist$dfdam
  MSwithinResample <- MSdflist$MSwithin * rchisq(1,MSdflist$dfwithin) / MSdflist$dfwithin
  
  MSpopResample <- MSdflist$MSpops * rchisq(1,MSdflist$dfpops) / MSdflist$dfpops
  
  sim.qstForCI <- QSTfromSireModel(MSpopResample,MSsiresResample,MSdamsResample,MSwithinResample,MSdflist$n0primeprime,MSdflist$n0prime,MSdflist$n0,MSdflist$nc0prime,MSdflist$nc0,MSdflist$ncb0)
  sim.qstForCI <- sim.qstForCI[[1]] # only want to return the Qst here, not also the Va  
  
  sigma2damsResample <- (MSdamsResample-MSwithinResample)/MSdflist$n0

  sigma2siresResample <- (MSsiresResample - MSwithinResample - MSdflist$n0prime*sigma2damsResample)/MSdflist$nc0
  sim.Va <- 4*sigma2siresResample
  
  return(c(sim.qstForCI, sim.Va))
}
###########################################################################



  
#########################       QstVa Parboot for dam model      ###############################
#
# New qstVa.parbootForCI.dammodel for dam model returns a Qst pseudovalue
# resampled by a paremetric bootstrap from the dam model. This follows the logic
# of O'Hara and Merila. It also returns a parametric bootrap value for Va. This
# requires MS for dams and residuals, along with the n0 value from the nested
# ANOVA and dfs for each level

qstVa.parbootForCI.dammodel <- function(MSdflist, dam.offspring.relatedness){
  MSdamsResample <- MSdflist$MSdams * rchisq(1,MSdflist$dfdam) / MSdflist$dfdam
  MSwithinResample <- MSdflist$MSwithin * rchisq(1,MSdflist$dfwithin) / MSdflist$dfwithin
  MSpopResample <- MSdflist$MSpops * rchisq(1,MSdflist$dfpops) / MSdflist$dfpops
  
  sim.qstForCI <- QSTfromDamModel(MSpopResample,MSdamsResample,MSwithinResample,MSdflist$n0prime,MSdflist$n0,MSdflist$nb0)  
  sim.Va <- 1/dam.offspring.relatedness*(MSdamsResample-MSwithinResample)/MSdflist$n0
  
  return(c(sim.qstForCI, sim.Va))
}
###########################################################################



#########################         wc.calc        #########################       
#
# returns the a, b, and c values required to compute Fst according to Weir and Cockerham 1984
#	modified from code previously implemented in hierfstat package 
#	(Jerome Goudet, http://cran.r-project.org/web/packages/hierfstat/index.html)
#
#	ndat		data frame with first column indicating population of origin and following representing loci
#	diploid 	Whether data are diploid

wc.calc <- function(ndat, diploid = TRUE){
	    if (!diploid){
        dum <- ndat[, -1]
        nd <- max(dum, na.rm = TRUE)
        modu <- 1000
        if(nd < 10)		modu <- 10
        if(nd < 100)	modu <- 100
        dum <- dum * modu + dum
        ndat <- data.frame(ndat[, 1], dum)
    }
    pop <- ndat[, 1]
    ni <- length(pop)
    dat <- ndat
    loc.names <- names(dat)[-1]
    n <- t(ind.count(dat)) ## NEED in.count to give number of inds genotyped per locus and per pop
    nt <- apply(n, 1, sum, na.rm = TRUE)
    untyped.loc <- which(nt == 0)
    typed.loc <- which(nt != 0)
    if(length(untyped.loc) > 0){
        dat <- dat[, -(untyped.loc + 1)]
        n <- t(ind.count(dat))
        nt <- apply(n, 1, sum, na.rm = TRUE)
    }
    alploc <- nb.alleles(cbind(rep(1, ni), dat[, -1]))
    np <- dim(n)[2]
    npl <- apply(n, 1, tempfun <- function(x) sum(!is.na(x)))
    nl <- dim(n)[1]
    p <- pop.freq(dat, diploid)
    pb <- pop.freq(cbind(rep(1, length(pop)), dat[, -1]), diploid)
    n <- matrix(unlist(n), ncol = np)
    nal <- n[rep(1:nl, alploc), ]
    nc <- (nt - apply(n^2, 1, sum, na.rm = TRUE)/nt)/(npl - 1)
    ntal <- rep(nt, alploc)
    ncal <- rep(nc, alploc)
    p <- matrix(unlist(lapply(p, t)), ncol = np, byrow = TRUE)
    pb <- matrix(unlist(pb), ncol = 1)
    if(diploid){
        dum <- getal.b(dat[, -1])
        all.loc <- apply(dum, 2, tempfun1 <- function(y) as.numeric(dimnames(table(y))[[1]]))
        hetpl <- apply(dum, 2, fun <- function(z) {
            lapply(as.numeric(dimnames(table(z))[[1]]), who.is.het <- function(y) apply(z == 
                y, 1, ind.is.het <- function(x) xor(x[1], x[2])))
        })
        mho <- lapply(hetpl, tempfun2 <- function(x) matrix(unlist(lapply(x, 
            tempfun3 <- function(y) tapply(y, pop, sum, na.rm = TRUE))), 
            ncol = np))
        mho <- matrix(unlist(mho), ncol = np, byrow = TRUE)
        mhom <- (2 * nal * p - mho)/2
    }else{mhom <- nal * p}
    SSG <- apply(nal * p - mhom, 1, sum, na.rm = TRUE)
    dum <- nal * (p - 2 * p^2) + mhom
    SSi <- apply(dum, 1, sum, na.rm = TRUE)
    dum1 <- nal * (sweep(p, 1, pb))^2
    SSP <- 2 * apply(dum1, 1, sum, na.rm = TRUE)
    ntalb <- rep(npl, alploc)
    MSG <- SSG/ntal
    MSP <- SSP/(ntalb - 1)
    MSI <- SSi/(ntal - ntalb)
    sigw <- MSG
    sigb <- 0.5 * (MSI - MSG)
    siga <- 1/2/ncal * (MSP - MSI)

    abc.mat <- cbind(siga, sigb, sigw)
    # this returns a matrix of a, b, and c values in columns with one allele per row
    return(abc.mat) 
}
###########################################################################



#########################       fst.sample      #########################       
#
# returns the Fst value computed on a sample of the alleles in the data
# @param:
#  - obs:  a table containing the components of variance for each locus
#          the table must have one line per allele and at least 3 columns corresponding
#          to the three coefficient a, b, and c as defined in Weir&Cockerham 1984
#
#  - nalleles: the size of the sample (i.e. num of alleles to draw from the table)

fst.sample <- function(obs, nalleles) {
  # choose the alleles to randomly sample:
  allele.smpl <- sample(1: nalleles,size= nalleles,replace=TRUE)
  #select the sampled alleles from the input table:
  dat <- obs[allele.smpl,]
  # Fst = a/(a+b+c); from Weir & Cockerham 1984
  return( sum(dat[,1])/sum(dat[,1]+dat[,2]+dat[,3]) )  
}
###########################################################################



#########################       read.fst.input      #########################
#
#function takes from user: filename, number of pops, the number of any 
#	extra columns in front (default is 1), and if there is a header row (default is yes)
#data must be in the form: n columns, columns of all q_hat values by pop and in order, columns 
#	of all q_hat variances by pop in order, any additional columns all saved in a .csv

read.fst.input <- function(q_hat.dat, num.pops, num.extra.columns=1){
	num.pops <- num.pops
	# which columns are the q_hats per population
	q.columns <- num.extra.columns + (1:num.pops)
	# which columns are the corresponding variances in q_hat per population
	last.q.column <- tail(q.columns, n=1)
	first.var.column <- last.q.column+1
	var.columns <- first.var.column:(first.var.column+num.pops-1)	
	
	# extract the data from input file
	dat <- q_hat.dat
	q_hat.matrix <- dat[,q.columns[1]:tail(q.columns, n=1)]
	var.q_hat.matrix <- dat[,var.columns[1]:tail(var.columns, n=1)]
	
	return(list(
	q_hat.matrix, 
	var.q_hat.matrix))
	# this returns a list containing the 2 matrices
	# assign results of the function to an object, then query [1] or [2]
	# 		for the q_hat.matrix and variance.matrix respectively
}
###########################################################################



#########################       mean.fst       #########################   
#
# calculate the mean fst from q_hat values for one sample
# this function will then be used in the fst.sample function
# necessary calculations come from Lynch and Milligan 1994, as implemented in AFLP-SURV (Vekemans et al. 2002)
#
#  - dat: a MATRIX of q_hat values per population in columns, per locus in rows (can be produced with AFLP-SURV)
#
#  - vardat: a MATRIX of the variance in q_hat, corresponding to the values in dat (can be produced with AFLP-SURV)
#
#  - num.loci: number of loci individuals are genotyped at
#
#  - num.pops: number of populations sampled

mean.fst <- function(dat, vardat, num.loci, num.pops){
	
##### calculate H_j(i)
	H_j_i.matrix <- matrix(NA, nrow=num.loci, ncol=num.pops)
	#	H_j(i) = 2q_j(i)[1-q_j(i)] + 2Var[q_j(i)]
	# j is one population, i is one locus
	for(j in 1:num.pops){ # loop through population columns
		for(i in 1:num.loci){ # loop through loci rows
			q <- dat[i,j]	# q value of one locus at one population
			#x <- xdat[i,j]	# x value of one locus at one population
			#N <- samps[i,j]	# sample size for one locus at one population
			var.q <- vardat[i,j] # variance in q for that locus across all populations
			H_j_i <- 2*q*(1-q) + 2*var.q # measure of gene diversity at one locus
			H_j_i.matrix[i,j] <- H_j_i # put it into a matrix			
		}	
	}
	
##### calculate H_j
	# H_j = average of H_j(i) across all loci for that population
	H_j.matrix <- matrix(NA, nrow=1, ncol=num.pops)
	for(j in 1:num.pops){ # get the H_j for each pop
		H_j.temp <- mean(H_j_i.matrix[,j]) # sum over all loci for one pop and divide by number of loci
			# done by taking mean per pop
		H_j.matrix[1,j] <- H_j.temp
	}	
	H_j <- H_j.matrix # this is the mean observed gene diversity in each pop, j
		
##### calculate H_jk
	H_jk <- matrix(0, nrow=num.pops, ncol=num.pops) 
		# make a matrix of num.pops x num.pops (j by k)
	for(j in 1:num.pops){
		for(k in 1:num.pops){
			sum.H_jki <- 0
			for(i in 1:num.loci){
				# first need H'_jk(i) where j and k are a pair of populations,  eqn 9a
				Hprime_jki <- dat[i,j] + dat[i,k] - 2*dat[i,j]*dat[i,k]
				# then calculate H_jk(i)  (i.e. not 'prime'),  eqn 10a
				H_jki <- Hprime_jki - (H_j_i.matrix[i,j] + H_j_i.matrix[i,k])/2
				# sum this across all loci in this loop
				sum.H_jki <- sum.H_jki + H_jki
			}
			# then divide by number of loci to get the H_jk for that pair of pops (eqn 12), and store it in a matrix
			H_jk[j,k] <- sum.H_jki/num.loci
		}
	}
	
##### calculate H_B, mean between-pop gene diversity
	# H_B = 2/(n(n-1)) * sum(H_jk)  where n is the number of populations  
	first.half.distinct.list <- NULL
	second.half.distinct.list <- NULL
	H_jk.distinct <- NULL
	for(a in 1:(num.pops-1)){
		b <- a+1
		pair <- expand.grid(a, b:num.pops)
		list5 <- pair[,1]
		list6 <- pair[,2]
		first.half.distinct.list <- c(first.half.distinct.list, list5)
		second.half.distinct.list <- c(second.half.distinct.list, list6)
	} # this made 2 lists that when lined up are all possible distinct pairs of populations
	# extract all the H_jk values of distinct population pairs using these lists
	for(z in 1:length(first.half.distinct.list)){
		H_jk.distinct.temp <- H_jk[first.half.distinct.list[z], second.half.distinct.list[z]]
		# print(c(first.half.distinct.list[z], second.half.distinct.list[z])) # check only getting distinct pop pairs, yes.
		H_jk.distinct <- c(H_jk.distinct, H_jk.distinct.temp)
	} # now do eqn 13a
	H_B <- (2 / (num.pops*(num.pops - 1))) * sum(H_jk.distinct)
	
##### calculate H_W, mean within-pop diversity
	# H_W = 1/n * sum(H_j)  where n is the number of populations
	H_W <- (1/num.pops)*sum(H_j)
	
##### calculate H_T
	# H_T = H_B + H_W
	H_T <- H_B + H_W
	
##### calculate Var(H_W)
	var.H_W <- (1 / (num.pops*(num.pops - 1))) * sum( (H_j - H_W)^2 )
	
##### calculate V_B
	# modified from Xavier Vekemans C code in AFLP-SURV, found in the Isomain.c file
	# V_B is the variance among diversity measures in non-overlapping pairs of populations - the variance in H_jk
	count <- 0
	sum.H_jk <- 0
	for(k in 1:(num.pops/2)*2){ # a sequence to a decimal only goes to the rounded down decimal, so if odd num.pops still works
		j <- k-1
		temp.H_jk <- H_jk[j,k]
		# sum the non-overlapping pop pairs H_jk's
		sum.H_jk <- sum.H_jk + temp.H_jk
		count <- count+1
		# print(c(j,k))
	}
	mean.H_jk <- sum.H_jk / count # mean non-overlapping pop H_jk
	sum.diffs.squared <- 0
	for(k in 1:(num.pops/2)*2){ # a sequence to a decimal only goes to the rounded down decimal, so if odd num.pops still works
		j <- k-1
		temp.diff.squared <- (H_jk[j,k] - mean.H_jk)^2
		# sum the squared differences to get the variance
		sum.diffs.squared <- sum.diffs.squared + temp.diff.squared
	}
	V_B <- sum.diffs.squared / ((count) * (count-1)) 
	if(num.pops <=3){V_B <- 0} # for 2 and 3 deme cases, V_B and C_B are set to zero 
	
##### calculate C_B
	# modified from Xavier Vekemans C code in AFLP-SURV, found in the Isomain.c file
	# C_B is the covariance among diversity measures in overlapping population pairs - the covariance in H_jk
	# covariance is sum over all x minus x_bar times y minus y_bar all divided by n-1
	# overlapping pops are those that share either of the pops in the pair
	count <- 0
	sum.H_jk <- 0
	sum.H_j2k2 <- 0
	for(j in 1:(num.pops-1)){ # loop through pop 1 of X pair
		for(k in (j+1):num.pops){ # loop through pop 2 of X pair
			# get first set of Y pair of populations, those which share pop 1 of X
			j2 <- j # this is first pop of Y pair, shared with first pop of X pair
			k2 <- k+1 # this is second pair of Y pop
			while(k2 <= num.pops){ # loop through first set, second pair of Y pop
				X.H_jk <- H_jk[j, k]
				sum.H_jk <- sum.H_jk + X.H_jk # get the mean of X pairs
				Y.H_jk <- H_jk[j2, k2]
				sum.H_j2k2 <- sum.H_j2k2 + Y.H_jk # get the mean of Y pairs
				# print(c(j,k,"  ",j2,k2)) 
				k2 <- k2+1
				count <- count + 1
			}
			# get second set of Y pair of populations, those which share pop 2 of X
			j2 <- k # this is first pop of Y pair, shared with second pop of X pair
			k2 <- k+1 # this is second pair of Y pop
			while(k2 <= num.pops){ # loop through second set, second pair of Y pop
				X.H_jk <- H_jk[j, k]
				sum.H_jk <- sum.H_jk + X.H_jk
				Y.H_jk <- H_jk[j2, k2]
				sum.H_j2k2 <- sum.H_j2k2 + Y.H_jk
				# print(c(j,k,"  ",j2,k2)) 
				k2 <- k2+1
				count <- count + 1
			}
		}
	}
	mean.H_jk <- sum.H_jk/count
	mean.H_j2k2 <- sum.H_j2k2/count
	cov <- 0
	# do the same as above but now get the covariance, i.e subtracting the mean from each and multiplying together
	for(j in 1:(num.pops-1)){
		for(k in (j+1):num.pops){
			j2 <- j
			k2 <- k+1
			while(k2 <= num.pops){
				cov.temp <- (H_jk[j, k] - mean.H_jk) * (H_jk[j2, k2] - mean.H_j2k2)
				cov <- cov + cov.temp
				k2 <- k2+1
			}
			j2 <- k
			k2 <- k+1
			while(k2 <= num.pops){
				cov.temp <- (H_jk[j, k] - mean.H_jk) * (H_jk[j2, k2] - mean.H_j2k2)
				cov <- cov + cov.temp
				k2 <- k2+1
			}
		}
	}
	C_B <- cov / ((count) * (count-1))
	if(num.pops <=3){C_B <- 0} # for 2 and 3 deme cases, V_B and C_B are set to zero 
	
##### calculate Var(H_B)
	var.H_B <- (2 * (V_B + 2*(num.pops-2)*C_B)) / (num.pops * (num.pops - 1))
	
##### calculate Cov(H_B, H_W)
	# need the sum of the sums term in there - sum of all H_jk per population j (except where j=k) times H_j
	# multiply the sum of each row of the H_jk matrix by its corresponding H_j from the H_j.matrix excluding pops where j=k in that sum
	sum.H_j.H_jk <- 0
	for(j in 1:num.pops){
		# sum the row of the H_jk matrix, then subtract off the H_jk where j=k
		temp.sum.H_jk <- 0
		for(k in 1:num.pops){
			if(j==k) same.jk <- H_jk[j,k]  #this will be the value to subtract off the sum of H_jk for that row
			temp.sum.H_jk <- temp.sum.H_jk + H_jk[j,k]
			# print(c(j,k))
		}
		sum.H_jk <- temp.sum.H_jk - same.jk
		temp.H_j.H_jk <- H_j[j] * sum.H_jk
		sum.H_j.H_jk <- sum.H_j.H_jk + temp.H_j.H_jk
	}
	cov.H_B.H_W <- ( (sum.H_j.H_jk/(num.pops*(num.pops-1)))-(H_W*H_B) )/num.pops
		
##### Fst!
	mean.Fst <- (H_B / H_T) * (1 + (H_B*var.H_W - H_W*var.H_B + (H_B - H_W)*cov.H_B.H_W) / (H_B*H_T^2 ))^(-1)
	
	return(mean.Fst)
}
###########################################################################



#########################       fst.sample.aflp       ########################
#
# returns the Fst value computed on a sample of the loci given in input as q_hat values
# @param:
#  - obs:  a table containing the q_hat values needed to calculate Fst
#
#  - qvar: the corresponding variances for given q_hat values
#
#  - nloci: the size of the sample (i.e. num of loci to draw from the table)

fst.sample.aflp <- function(obs, qvar, nloci) {
  # take a random sample of loci, with replacement, of the specified size:
  loc.smpl <- sample(1:nloci,size=nloci,replace=TRUE)
  # select the sampled loci from the input table of data:
  q_hats <- obs[loc.smpl,]
  # using this number of loci,
  num.loci <- length(obs[,1])
  # select the corresponding variances
  vars <- qvar[loc.smpl,]
  # and lastly need number of pops
  num.pops <- length(obs[1,])

  # calculate Fst based on function built in mean.fst
    return( mean.fst(q_hats, vars, num.loci, num.pops) )  
}

###########################################################################


NULL
#' AFLP marker dataset
#'
#' An example dataset created from simulations of 15 neutral populations under 
#' an island model. Individuals are diploid with genotypes at 100 dominant, AFLP loci.
#'
#' \describe{
#'  \item{pop}{population of origin}
#'  \item{q_hat}{columns 2 through 16 are the q_hat values for each locus in each population}
#'	\item{q_var}{columns 17 through 31 are the q_hat variances for each locus in each population}
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 100 observations on the following 31 variables.
#' @name aflp

NULL
#' Biallelic marker dataset
#'
#' An example dataset created from simulations of 15 neutral populations under 
#' an island model. Individuals are diploid with genotypes at 25 biallelic loci.
#'
#' \describe{
#'  \item{pop}{population of origin}
#'  \item{loc}{loc1 through loc25 are the 25 biallelic loci for which each diploid individual is genotyped}
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 375 observations on the following 26 variables.
#' @name biallelic

NULL
#' Half-sib dam trait dataset
#'
#' An example dataset of a neutral trait created from simulations of 15 
#' populations under an island model. Offspring per dam are unbalanced across 
#' the breeding design with 6 dams per population, an average of 10 offspring per 
#' dam, and 900 individuals total.
#'
#' \describe{
#'  \item{pop}{a column of population identifiers}
#'  \item{dam}{a column of unique dam identifiers}
#'  \item{trait}{column of trait values per individual}
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 900 observations on the following 3 variables.
#' @name hsdam

NULL
#' Half-sib sire trait dataset
#'
#' An example dataset of a neutral trait created from simulations of 15 populations 
#' under an island model. 900 individuals total come from an average of 4 dams 
#' per sire and an average of 12 sires per population for an overall average of 60 
#' individuals per population.
#'
#' \describe{
#'  \item{pop}{a column of population identifiers}
#'  \item{sire}{a column of unique sire identifiers}
#'  \item{dam}{a column of unique dam identifiers}
#'  \item{trait}{column of trait values per individual}
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 900 observations on the following 4 variables.
#' @name hssire

NULL