#' intact.hyper()
#'
#' Performs hyper.tests  on an <.intact> object column .
#'
#' @param X any column of a <Query.intact> object
#' @param Y any column of a <Universe.intact> object (has to be the same column as X)
#'
#' @details Peforms a hypergeometric enrichment test on a column of a <Query.intact>  object compared to the same column of a <Universe.intact> object.
#' @details The universe file db.intact is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.hyper(queryExample.intact$`Main class`, universeExample.intact$`Main class`)
#'
#' @author Geremy Clair
#' @export
#'

intact.hyper<- function(X,Y){
  if(typeof(X)!="character"){stop("X(your query) is not a column from an intact object")}
  if(typeof(Y)!="character"){stop("Y(your universe) is not a column from an intact object")}

  query.factor <- factor(X)
  universe.factor <- factor(Y)
  query.counts<-as.data.frame(table(query.factor))
  universe.counts<-as.data.frame(table(universe.factor))
  classes<-as.character(query.counts[,1])

  #create contingency matrixes
  contingency.matrix<- list()
  for (i in 1:length(classes))
  {
    contingency.matrix[[i]]<- matrix(NA,ncol=2,nrow=2)
    contingency.matrix[[i]][1,1]<- query.counts[classes[i]==query.counts[,1],2]
    contingency.matrix[[i]][1,2]<- universe.counts[classes[i]==universe.counts[,1],2]
    contingency.matrix[[i]][2,1]<- sum(query.counts[,2])-contingency.matrix[[i]][1,1]
    contingency.matrix[[i]][2,2]<- sum(universe.counts[,2])-contingency.matrix[[i]][1,2]
  }

  names(contingency.matrix) <- classes

  hyper.results<-data.frame(matrix(NA, nrow=length(classes), ncol=8))
  colnames(hyper.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")
  hyper.results$Classifier <- classes

  for (i in 1:length(classes))
  {
    test<- min(1-cumsum(dhyper(0:(contingency.matrix[[i]][1,1]-1),contingency.matrix[[i]][1,2],contingency.matrix[[i]][2,2],contingency.matrix[[i]][2,1])))
    hyper.results$Pvalue[i]<-test
    hyper.results$BHadjustPvalue[i]<-p.adjust(hyper.results$Pvalue[i], method = "BH", n=length(classes))
    hyper.results$Count.query [i]<- paste(contingency.matrix[[i]][1,1],"/",sum(query.counts[,2]), sep="")
    hyper.results$Count.universe [i]<- paste(contingency.matrix[[i]][1,2],"/",sum(universe.counts[,2]), sep="")
    hyper.results$`%.query` [i]<- query.counts[classes[i]==query.counts[,1],2]/sum(query.counts[,2])*100
    hyper.results$`%.universe`[i]<- universe.counts[classes[i]==universe.counts[,1],2]/sum(universe.counts[,2])*100
    hyper.results$fold.change [i]<- hyper.results$`%.query`[i]/hyper.results$`%.universe`[i]
  }
  hyper.results
}

#' allchains.hyper()
#'
#' Performs hyper.tests  on an <.allchains> object
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#'
#' @details Peforms a hypergeometric enrichment test on a <Query.allchains> object compared to the same column of a <Universe.allchains> object.
#' @details The universe file db.intact is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.hyper(queryExample.allchains, universeExample.allchains)
#'
#' @author Geremy Clair
#' @export
#'

allchains.hyper<- function(X,Y){
  if(typeof(X)!="list"){stop("X(your query) is not a .allchains object")}
  if(typeof(Y)!="list"){stop("Y(your universe) is not a .allchains object")}

  allchains.query<-X[,2]
  allchains.universe<-Y[,2]
  query.totalchains<-length(allchains.query)
  universe.totalchains<-length(allchains.universe)

  unique.query<-unique(allchains.query)

  contingency.matrix<- list()
  for (i in 1:length(unique.query))
  {
    contingency.matrix[[i]]<- matrix(NA,ncol=2,nrow=2)
    contingency.matrix[[i]][1,1]<- sum(allchains.query==unique.query[i])
    contingency.matrix[[i]][1,2]<- sum(allchains.universe==unique.query[i])
    contingency.matrix[[i]][2,1]<- query.totalchains - contingency.matrix[[i]][1,1]
    contingency.matrix[[i]][2,2]<- universe.totalchains - contingency.matrix[[i]][1,2]
  }

  names(contingency.matrix)<- unique.query

  hyper.results<-data.frame(matrix(NA, nrow=length(unique.query), ncol=8))
  colnames(hyper.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")
  hyper.results$Classifier <- unique.query

  for (i in 1:length(unique.query))
  {
    hyper.results$Pvalue[i]<-min(1-cumsum(dhyper(0:(contingency.matrix[[i]][1,1]-1),contingency.matrix[[i]][1,2],contingency.matrix[[i]][2,2],contingency.matrix[[i]][2,1])))
    hyper.results$BHadjustPvalue[i]<-p.adjust(hyper.results$Pvalue[i], method = "BH", n=length(unique.query))
    hyper.results$Count.query[i]<- paste(contingency.matrix[[i]][1,1],"/",query.totalchains,sep="")
    hyper.results$Count.universe[i]<- paste(contingency.matrix[[i]][1,2],"/",universe.totalchains,sep="")
    hyper.results$`%.query`[i]<- (contingency.matrix[[i]][1,1])/query.totalchains*100
    hyper.results$`%.universe`[i]<- contingency.matrix[[i]][1,2]/universe.totalchains*100
    hyper.results$fold.change[i]<-  hyper.results$`%.query`[i]/hyper.results$`%.universe`[i]
  }

  hyper.results

}

#' chain.hyper()
#'
#' Performs hyper.tests  on an <.chain> object
#'
#' @param X any <Query.chain> object
#' @param Y any <Universe.chain> object
#'
#' @details Peforms a hypergeometric enrichment test on a <Query.chain> object compared to the same column of a <Universe.chain> object.
#' @details The universe file db.intact is provided if you are lacking an universe file.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.hyper(queryExample.chain, universeExample.chain)
#'
#' @author Geremy Clair
#' @export
#'

chain.hyper<- function(X,Y){
  if(typeof(X)!="list"){stop("X(your query) is not a .chain object")}
  if(typeof(Y)!="list"){stop("Y(your universe) is not a .chain object")}

  for (i in 2:9){
    X[,i]<-as.numeric(X[,i])
    Y[,i]<-as.numeric(Y[,i])
  }

  # Calculate number of chains in Univ and query
  query.totalchain<- sum(colSums(X[,2:5]))
  universe.totalchain<- sum(colSums(Y[,2:5]))

  contingency.matrix<- list()
  contingency.matrix[[1]]<-matrix(1,ncol=2,nrow=2)

  for (i in 2:9)
  {
    contingency.matrix[[i]]<- matrix(NA,ncol=2,nrow=2)
    contingency.matrix[[i]][1,1]<- sum(X[,i])
    contingency.matrix[[i]][1,2]<- sum(Y[,i])
    contingency.matrix[[i]][2,1]<- query.totalchain-sum(X[,i])
    contingency.matrix[[i]][2,2]<- universe.totalchain-sum(Y[,i])
  }

  names(contingency.matrix)<- colnames(X)
  hyper.results<-data.frame(matrix(NA, nrow=9, ncol=8))
  colnames(hyper.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")
  hyper.results$Classifier<- colnames(X[,1:9])

  for (i in 2:9)
  {
    hyper.results$Count.query[i]<- paste(sum(X[,i]),"/",query.totalchain,sep="")
    hyper.results$Count.universe[i]<- paste(sum(Y[,i]),"/",universe.totalchain,sep="")
    hyper.results$`%.query`[i]<- sum(X[,i])/query.totalchain*100
    hyper.results$`%.universe`[i]<- sum(Y[,i])/universe.totalchain*100
    #    hyper.results$fold.change[i]<-  hyper.results$`%.query`[i]/hyper.results$`%.universe`[i]
    if(hyper.results$`%.universe`[i]==0){hyper.results$fold.change[i]<-0}
    else{hyper.results$fold.change[i]<-  hyper.results$`%.query`[i]/hyper.results$`%.universe`[i]}

    hyper.results$Pvalue[i]<- min(1-cumsum(dhyper(0:(contingency.matrix[[i]][1,1]-1),contingency.matrix[[i]][1,2],contingency.matrix[[i]][2,2],contingency.matrix[[i]][2,1])))
    hyper.results$BHadjustPvalue[i]<-p.adjust(hyper.results$Pvalue[i], method = "BH", n=4)
  }
  hyper.results[2:9,]

}
