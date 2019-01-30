#' intact.EASE()
#'
#' Performs EASE.tests  on an <.intact> object column .
#'
#' @param X any column of a <Query.intact> object
#' @param Y any column of a <Universe.intact> object (has to be the same column as X)
#'
#' @details Peforms an EASE enrichment test (modified Fisher's Exact) on a column of a <Query.intact>  object compared to the same column of a <Universe.intact> object.
#' @details The universe file db.intact is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.EASE(queryExample.intact$`Main class`, universeExample.intact$`Main class`)
#'
#' @author Geremy Clair
#' @export
#'

intact.EASE<- function(X,Y){

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
    contingency.matrix[[i]][1,1]<- query.counts[classes[i]==query.counts[,1],2]-1
    contingency.matrix[[i]][1,2]<- universe.counts[classes[i]==universe.counts[,1],2]
    contingency.matrix[[i]][2,1]<- sum(query.counts[,2])-contingency.matrix[[i]][1,1]
    contingency.matrix[[i]][2,2]<- sum(universe.counts[,2])-contingency.matrix[[i]][1,2]
  }

  names(contingency.matrix) <- classes

  EASE.results<-data.frame(matrix(NA, nrow=length(classes), ncol=8))
  colnames(EASE.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")
  EASE.results$Classifier <- classes

  for (i in 1:length(classes))
  {
    test<-fisher.test(contingency.matrix[[i]])
    EASE.results$`p-value`[i]<-test$p.value
    EASE.results$Count.query [i]<- paste((contingency.matrix[[i]][1,1]+1),"/",sum(query.counts[,2]), sep="")
    EASE.results$Count.universe [i]<- paste(contingency.matrix[[i]][1,2],"/",sum(universe.counts[,2]), sep="")
    EASE.results$`%.query` [i]<- query.counts[classes[i]==query.counts[,1],2]/sum(query.counts[,2])*100
    EASE.results$`%.universe`[i]<- universe.counts[classes[i]==universe.counts[,1],2]/sum(universe.counts[,2])*100
    EASE.results$Fold.change [i]<- EASE.results$`%.query`[i]/EASE.results$`%.universe`[i]
  }
  EASE.results$`p-value`[is.na(EASE.results$`p-value`)|EASE.results$`p-value`==TRUE|EASE.results$`p-value`>1]<-1
  if(length(EASE.results$`p-value`)<=1){EASE.results$`FDR.q-value`<-1}else{EASE.results$`FDR.q-value`<- qvalue(EASE.results$`p-value`,lambda=0)$qvalues}

  EASE.results
}
#' allchains.EASE()
#'
#' Performs EASE.tests  on an <.allchains> object
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#'
#' @details Peforms an EASE enrichment test (modified Fisher's Exact) on a <Query.allchains> object compared to the same column of a <Universe.allchains> object.
#' @details The universe file db.intact is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.EASE(queryExample.allchains, universeExample.allchains)
#'
#' @author Geremy Clair
#' @export
#'
allchains.EASE<- function(X,Y){
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
    contingency.matrix[[i]][1,1]<- sum(allchains.query==unique.query[i])-1
    contingency.matrix[[i]][1,2]<- sum(allchains.universe==unique.query[i])
    contingency.matrix[[i]][2,1]<- query.totalchains - contingency.matrix[[i]][1,1]
    contingency.matrix[[i]][2,2]<- universe.totalchains - contingency.matrix[[i]][1,2]
  }

  names(contingency.matrix)<- unique.query

  EASE.results<-data.frame(matrix(NA, nrow=length(unique.query), ncol=8))
  colnames(EASE.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")
  EASE.results$Classifier <- unique.query

  for (i in 1:length(unique.query))
  {
    test<-fisher.test(contingency.matrix[[i]])
    EASE.results$`p-value`[i]<-test$p.value
    EASE.results$Count.query[i]<- paste(contingency.matrix[[i]][1,1]+1,"/",query.totalchains,sep="")
    EASE.results$Count.universe[i]<- paste(contingency.matrix[[i]][1,2],"/",universe.totalchains,sep="")
    EASE.results$`%.query`[i]<- (contingency.matrix[[i]][1,1]+1)/query.totalchains*100
    EASE.results$`%.universe`[i]<- contingency.matrix[[i]][1,2]/universe.totalchains*100
    EASE.results$Fold.change[i]<-  EASE.results$`%.query`[i]/EASE.results$`%.universe`[i]
  }
  EASE.results$`p-value`[is.na(EASE.results$`p-value`)|EASE.results$`p-value`==TRUE|EASE.results$`p-value`>1]<-1
  if(length(EASE.results$`p-value`)<=1){EASE.results$`FDR.q-value`<-1}else{EASE.results$`FDR.q-value`<- qvalue(EASE.results$`p-value`,lambda=0)$qvalues}
  EASE.results

}

#' chain.EASE()
#'
#' Performs EASE.tests  on an <.chain> object
#'
#' @param X any <Query.chain> object
#' @param Y any <Universe.chain> object
#'
#' @details Peforms an EASE enrichment test (modified Fisher's Exact) on a <Query.chain> object compared to the same column of a <Universe.chain> object.
#' @details The universe file db.intact is provided if you are lacking an universe file.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.EASE(queryExample.chain, universeExample.chain)
#'
#' @author Geremy Clair
#' @export
#'
chain.EASE<- function(X,Y){
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
    contingency.matrix[[i]][1,1]<- sum(X[,i])-1
    contingency.matrix[[i]][1,2]<- sum(Y[,i])
    contingency.matrix[[i]][2,1]<- query.totalchain-sum(X[,i])
    contingency.matrix[[i]][2,2]<- universe.totalchain-sum(Y[,i])
  }

  names(contingency.matrix)<- colnames(X)
  EASE.results<-data.frame(matrix(NA, nrow=9, ncol=8))
  colnames(EASE.results)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")
  EASE.results$Classifier<- colnames(X[,1:9])

  for (i in 2:9)
  {
    EASE.results$Count.query[i]<- paste(sum(X[,i]),"/",query.totalchain,sep="")
    EASE.results$Count.universe[i]<- paste(sum(Y[,i]),"/",universe.totalchain,sep="")
    EASE.results$`%.query`[i]<- sum(X[,i])/query.totalchain*100
    EASE.results$`%.universe`[i]<- sum(Y[,i])/universe.totalchain*100
    if(EASE.results$`%.universe`[i]==0){EASE.results$Fold.change[i]<-0}
    else{EASE.results$Fold.change[i]<-  EASE.results$`%.query`[i]/EASE.results$`%.universe`[i]}

    if(contingency.matrix[[i]][1,1]>-1){
      test<-fisher.test(contingency.matrix[[i]])
      EASE.results$`p-value`[i]<-test$p.value
    }else{
      EASE.results$`p-value`[i]<-1
    }
  }
  EASE.results<-EASE.results[2:9,]
  EASE.results$`p-value`[is.na(EASE.results$`p-value`)|EASE.results$`p-value`==TRUE|EASE.results$`p-value`>1]<-1
  if(length(EASE.results$`p-value`)<=1){EASE.results$`FDR.q-value`<-1}else{EASE.results$`FDR.q-value`<- qvalue(EASE.results$`p-value`,lambda=0)$qvalues}
  EASE.results
}


