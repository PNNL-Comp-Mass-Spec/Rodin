#' intact.binom()
#'
#' Performs binom.tests  on an <.intact> object column .
#'
#' @param X any column of a <Query.intact> object
#' @param Y any column of a <Universe.intact> object (has to be the same column as X)
#'
#' @details Peforms a binomial enrichment test on a column of a <Query.intact>  object compared to the same column of a <Universe.intact> object.
#' @details The universe file db.intact is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.binom(queryExample.intact$`Main class`, universeExample.intact$`Main class`)
#'
#' @author Geremy Clair
#' @export
#'
intact.binom<-function(X,Y){
  if(typeof(X)!="character"){stop("X(your query) is not a column from an intact object")}
  if(typeof(Y)!="character"){stop("Y(your universe) is not a column from an intact object")}

  query.factor <- factor(X)
  universe.factor <- factor(Y)
  query.counts<-as.data.frame(table(query.factor))
  universe.counts<-as.data.frame(table(universe.factor))
  classes<-as.character(query.counts[,1])

  param.vectors<- list()
  for (i in 1:length(classes))
  {
    param.vectors[[i]]<- vector(mode="numeric", length = 4)
    param.vectors[[i]][1]<- query.counts[classes[i]==query.counts[,1],2]
    param.vectors[[i]][2]<- length(X)
    param.vectors[[i]][3]<- universe.counts[classes[i]==universe.counts[,1],2]
    param.vectors[[i]][4]<- length(Y)
  }
  names(param.vectors) <- classes

  binom.results<- data.frame(matrix(NA, nrow=length(classes), ncol=8))
  colnames(binom.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")
  binom.results[,1] <- classes

  for (i in 1:length(classes))
  {
    test<-binom.test(query.counts[classes[1]==query.counts[,1],2],length(X),p = universe.counts[classes[i]==universe.counts[,1],2]/length(Y))
    binom.results$Count.query[i] <- paste(param.vectors[[i]][1],"/",param.vectors[[i]][2], sep="")
    binom.results$Count.universe[i] <- paste(param.vectors[[i]][3],"/",param.vectors[[i]][4], sep="")
    binom.results$`%.query`[i]<-param.vectors[[i]][1]/param.vectors[[i]][2]*100
    binom.results$`%.universe`[i]<-  param.vectors[[i]][3]/param.vectors[[i]][4]*100
    binom.results$Pvalue[i] <-test$p.value
    binom.results$BHadjustPvalue[i]<-p.adjust(binom.results$Pvalue[i], method = "BH", n=length(classes))
    binom.results$fold.change [i]<- binom.results$`%.query`[i]/binom.results$`%.universe`[i]
  }
  binom.results
}

#' allchains.binom()
#'
#' Performs binom.tests  on an <.allchains> object
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#'
#' @details Peforms a binomial enrichment test on a <Query.allchains> object compared to the same column of a <Universe.allchains> object.
#' @details The universe file db.intact is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.binom(queryExample.allchains, universeExample.allchains)
#'
#' @author Geremy Clair
#' @export
#'
allchains.binom<- function(X,Y){
  if(typeof(X)!="list"){stop("X(your query) is not a .allchains object")}
  if(typeof(Y)!="list"){stop("Y(your universe) is not a .allchains object")}

  allchains.query<-X[,2]
  allchains.universe<-Y[,2]
  unique.query<-unique(allchains.query)

  param.vectors<- list()
  for (i in 1:length(unique.query))
  {
    param.vectors[[i]]<- vector(mode="numeric", length = 4)
    param.vectors[[i]][1]<- sum(allchains.query==unique.query[i])
    param.vectors[[i]][2]<- length(allchains.query)
    param.vectors[[i]][3]<- sum(allchains.universe==unique.query[i])
    param.vectors[[i]][4]<- length(allchains.universe)
  }
  names(param.vectors) <- unique.query

  binom.results<-data.frame(matrix(NA, nrow=length(unique.query), ncol=8))
  colnames(binom.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")
  binom.results$Classifier <-unique.query

  for (i in 1:length(unique.query))
  {
    test<-binom.test(param.vectors[[i]][1],param.vectors[[i]][2],p=param.vectors[[i]][3]/param.vectors[[i]][4])
    binom.results$Count.query[i] <- paste(param.vectors[[i]][1],"/",param.vectors[[i]][2], sep="")
    binom.results$Count.universe[i] <- paste(param.vectors[[i]][3],"/",param.vectors[[i]][4], sep="")
    binom.results$`%.query`[i]<-param.vectors[[i]][1]/param.vectors[[i]][2]*100
    binom.results$`%.universe`[i]<-  param.vectors[[i]][3]/param.vectors[[i]][4]*100
    binom.results$Pvalue[i] <-test$p.value
    binom.results$BHadjustPvalue[i]<-p.adjust(binom.results$Pvalue[i], method = "BH", n=length(unique.query))
    binom.results$fold.change [i]<- binom.results$`%.query`[i]/binom.results$`%.universe`[i]
  }
  binom.results
}

#' chain.binom()
#'
#' Performs binom.tests  on an <.chain> object
#'
#' @param X any <Query.chain> object
#' @param Y any <Universe.chain> object
#'
#' @details Peforms a binomial enrichment test on a <Query.chain> object compared to the same column of a <Universe.chain> object.
#' @details The universe file db.intact is provided if you are lacking an universe file.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.binom(queryExample.chain, universeExample.chain)
#'
#' @author Geremy Clair
#' @export
#'
chain.binom<- function(X,Y){
  if(typeof(X)!="list"){stop("X(your query) is not a .chain object")}
  if(typeof(Y)!="list"){stop("Y(your universe) is not a .chain object")}

  for (i in 2:9){
    X[,i]<-as.numeric(X[,i])
    Y[,i]<-as.numeric(Y[,i])
  }

  param.vectors<- list()
  for (i in 2:9)
  {
    param.vectors[[i]]<- vector(mode="numeric", length = 4)
    param.vectors[[i]][1]<- sum(X[,i])
    param.vectors[[i]][2]<- sum(colSums(X[,2:5]))
    param.vectors[[i]][3]<- sum(Y[,i])
    param.vectors[[i]][4]<- sum(colSums(Y[,2:5]))
  }

  names(param.vectors)<- colnames(X[,1:9])
  binom.results<-data.frame(matrix(NA, nrow=9, ncol=8))
  colnames(binom.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")
  binom.results$Classifier<- colnames(X[,1:9])

  for (i in 2:9)
  {
    if(i!=1)
    {
      test<-binom.test(param.vectors[[i]][1],param.vectors[[i]][2],p=param.vectors[[i]][3]/param.vectors[[i]][4])
      binom.results$Count.query[i] <- paste(param.vectors[[i]][1],"/",param.vectors[[i]][2], sep="")
      binom.results$Count.universe[i] <- paste(param.vectors[[i]][3],"/",param.vectors[[i]][4], sep="")
      binom.results$`%.query`[i]<-param.vectors[[i]][1]/param.vectors[[i]][2]*100
      binom.results$`%.universe`[i]<-  param.vectors[[i]][3]/param.vectors[[i]][4]*100
      binom.results$Pvalue[i] <-test$p.value
      binom.results$BHadjustPvalue[i]<-p.adjust(binom.results$Pvalue[i], method = "BH", n=4)
      if(binom.results$`%.universe`[i]==0){binom.results$fold.change[i]<-0}
      else{binom.results$fold.change[i]<-  binom.results$`%.query`[i]/binom.results$`%.universe`[i]}

     }
  }
  binom.results[2:9,]
}
