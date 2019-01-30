#' plot_ranking()
#'
#' plots the ranking values of an object of type df with the lipid names as first column and the ranking values as second column
#'
#' @param df a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#'
#' @details plots the ranking values
#'
#' @return This function will return a plot generated through ggplot2
#'
#' @examples plot_ranking(cleaned.RTexample,order="ascending")
#'
#' @author Geremy Clair
#' @export
plot_ranking<-function(df,order="ascending"){
  df<-data.frame(df)
  if(ncol(df)!=2){stop("df should contain two columns")}
  df[df==""]<-NA
  if(sum(is.na(df))>0){stop("Your ranking table should not contain missing values")}

  if(missing(order)){order="ascending"}
  if(!order %in% c("ascending", "descending")){stop("order should either be 'ascending' or 'descending'")}

  #order based on the ranking values
  if(order=="ascending"){df<-df[order(df[,2]),]}else{df<-df[order(df[,2],decreasing = TRUE),]}
  df$Lipid_ranked<-1:nrow(df)

  #polot
  ggplot(df,aes(x=Lipid_ranked,y=df[,2]))+geom_point(color="darkslategray4")+theme_bw()+ylab(colnames(df)[2])

}

#' intact.KS()
#'
#' Performs weighted Kolmogorov-Smirnov test on a column of an <.intact> object based on Ranking values contained in a rankingTable
#'
#' @param X a <Rodin.intact> object
#' @param colname a colname from the <Rodin.intact> object
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#'
#' @details Peforms a ks test enrichment test on a column of a <Query.intact> based on ranking values contained in a ranking table
#' @details colname has to be in :"Lipid", "Category","Main class","Sub class","Total Number of Carbon","Double Bonds","Chain1","Chain2","Chain3","Chain4","NCarbonChain1","NCarbonChain2","NCarbonChain3","NCarbonChain4","DBChain1","DBChain2","DBChain3","DBChain4"
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.KS(RTexample.intact, colname="Main class", cleaned.RTexample)
#'
#' @author Geremy Clair
#' @export
intact.KS<- function(X,colname="Main class", rankingTable, order="ascending"){
  if(!ncol(X)==18 & !sum(colnames(X)[1:3]==c("Lipids","Category","Main class"))==3){stop("X should be an <.intact> object")}
  if(!is.data.frame(rankingTable)){stop("rankingTable should be a data frame of type <Ranking Table> containing the lipid list as first column and the ranking values as second column, we recommend that you use the function clean.RankingTable() to clean the table")}
  if(ncol(rankingTable)!=2){stop("rankingTable should contain two columns")}
  rankingTable[rankingTable==""]<-NA
  if(sum(is.na(rankingTable))>0){stop("rankingTable should not contain missing values")}

  if(!colname %in% colnames(X)){
    print("colname is not present in the object X, the columns available are the folling:")
    print(colnames(X))
    stop()
  }

  if(missing(order)){order=="ascending"}
  if(!order %in% c("ascending", "descending")){stop("order should either be 'ascending' or 'descending'")}


  #extract the column of interest
  X_vect<-X[,colnames(X)==colname]
  names(X_vect)<-X$Lipid

  #order the ranking table by the ranking values
  if(order=="ascending"){rankingTable<-rankingTable[order(rankingTable[,2]),]}else{rankingTable<-rankingTable[order(rankingTable[,2],decreasing=TRUE),]}
  rankingTable$KSRank<-1:nrow(rankingTable)

  #order X and rankingTable
  X<-X[order(X[,1]),]
  rankingTable<-rankingTable[order(rankingTable[,1]),]

  #verify that rankingTable contains the same lipids than the ones contained in X
  if(!sum(X[,1] %in% rankingTable[,1])>0){stop("The lipids contained in your <.intact> object were not identical as the ones in the rankingTable")}

  #extract the column of interest
  X_col<-X[,colname==colnames(X)]
  names(X_col)<-X$Lipid

  #identify the levels
  X_levels<-unique(X_col)

  #create a list of sub RankingTables
  RT_by_factor<-list()
  for (i in 1:length(X_levels)){
    Lipids<-X_col[X_col==X_levels[i]]
    Rank<-rankingTable[rankingTable[,1] %in% names(Lipids),]
    Rank<-Rank[order(Rank$KSRank),]
    RT_by_factor[[i]]<-Rank}

  names(RT_by_factor)<-X_levels

  #realize the KS.tests
  p<-numeric()
  for (i in 1:length(X_levels)){
    if(length(seq_len(nrow(X)))>length(RT_by_factor[[i]]$KSRank)){
      p[i]<-ks.test(RT_by_factor[[i]]$KSRank,seq_len(nrow(X))[-RT_by_factor[[i]]$KSRank],alternative="greater")$p.value}else{p[i]<-1}
  }
  #adjust pvalues using qvalue method
  if(length(p)<=1){q<-p}else{q <- qvalue(p,lambda=0)$qvalues}

  #create the final table
  final_table<-data.frame(matrix(ncol=4,nrow= length(X_levels)))
  colnames(final_table) <- c("Classifier","Count.in.list", "p-value", "FDR.q-value")
  final_table$Classifier<-X_levels
  for(i in 1:length(X_levels)){final_table$Count.in.list[i]<-nrow(RT_by_factor[[i]])}
  final_table$`p-value`<-p
  final_table$`FDR.q-value`<-q

  return(final_table)

}

#' chains.KS()
#'
#' Performs weighted Kolmogorov-Smirnov test on a <.chain> object based on Ranking values contained in a rankingTable
#'
#' @param X any <RT.chain> object
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#'
#' @details Peforms a ks test enrichment test on the chain parameters saved in a <.chain> object based on ranking values contained in a ranking table
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.KS(RTexample.chain, cleaned.RTexample)
#'
#' @author Geremy Clair
#' @export
#'
chain.KS<- function(X, rankingTable, order="ascending"){
  if(!ncol(X)==9 & !sum(colnames(X)[1:3]==c("Lipids","SCFA","MCFA"))==3){stop("X should be an <.chain> object")}

  if(!is.data.frame(rankingTable)){stop("rankingTable should be a data frame of type <Ranking Table> containing the lipid list as first column and the ranking values as second column, we recommend that you use the function clean.RankingTable() to clean the table")}
  if(ncol(rankingTable)!=2){stop("rankingTable should contain two columns")}
  rankingTable[rankingTable==""]<-NA
  if(sum(is.na(rankingTable))>0){stop("rankingTable should not contain missing values")}

  if(missing(order)){order=="ascending"}
  if(!order %in% c("ascending", "descending")){stop("order should either be 'ascending' or 'descending'")}

  #order X
  X<-X[order(X[,1]),]
  #keep only the values of the RankingTable that are in X (if TG rm is active this is necessary)
  rankingTable<-rankingTable[rankingTable[,1] %in% as.character(t(X[,1])),]
  #order ranking table
  rankingTable<-rankingTable[order(rankingTable[,1]),]

  #verify that rankingTable contains the same lipids than the ones contained in X
  if(sum(X[,1]!=rankingTable[,1])>0){stop("The lipids contained in your <.allchains> object were not identical as the ones in the rankingTable")}

  merged_table<-cbind(X,Ranking_column=rankingTable[,2])

  merged_table.list<-list()

  for (i in 1:8){
  merged_table.list[[i]] <- merged_table[rep(row.names(merged_table), merged_table[,i+1]), c(1,i+1,10)]
  merged_table.list[[i]]<-merged_table.list[[i]][,c(1,3)]
  merged_table.list[[i]]<-merged_table.list[[i]][order(merged_table.list[[i]]$Ranking_column),]
  merged_table.list[[i]]$Type<-rep(colnames(X)[1+i],nrow(merged_table.list[[i]]))
  }

  #calculate the ranking for the length
  all_list_length<-data.frame()
  all_list_length<-rbind(merged_table.list[[1]],merged_table.list[[2]],merged_table.list[[3]],merged_table.list[[4]])
  if(order=="ascending"){all_list_length<-all_list_length[order(all_list_length$Ranking_column),]}else{all_list_length<-all_list_length[order(all_list_length$Ranking_column,decreasing=TRUE),]}
  all_list_length$Rank<-1:nrow(all_list_length)

  #calculate the ranking for the unsat
  all_list_unsat<-data.frame()
  all_list_unsat<-rbind(merged_table.list[[5]],merged_table.list[[6]],merged_table.list[[7]],merged_table.list[[8]])
  if(order=="ascending"){all_list_unsat<-all_list_unsat[order(all_list_unsat$Ranking_column),]}else{all_list_unsat<-all_list_unsat[order(all_list_unsat$Ranking_column,decreasing=TRUE),]}
  all_list_unsat$Rank<-1:nrow(all_list_unsat)

  #Replace this in the merged_table.list (to add the rank there)
  merged_table.list[[1]]<- all_list_length[ all_list_length$Type=="SCFA",]
  merged_table.list[[2]]<- all_list_length[ all_list_length$Type=="MCFA",]
  merged_table.list[[3]]<- all_list_length[ all_list_length$Type=="LCFA",]
  merged_table.list[[4]]<- all_list_length[ all_list_length$Type=="VLCFA",]
  merged_table.list[[5]]<- all_list_unsat[ all_list_unsat$Type=="Saturated",]
  merged_table.list[[6]]<- all_list_unsat[ all_list_unsat$Type=="Monounsaturated",]
  merged_table.list[[7]]<- all_list_unsat[ all_list_unsat$Type=="Diunsaturated",]
  merged_table.list[[8]]<- all_list_unsat[ all_list_unsat$Type=="Polyunsaturated",]

  #realize the KS.tests
  p<-numeric()
  for (i in 1:4){
    if(nrow(merged_table.list[[i]])==0){p[i]<-1}else{
      if(length(seq_len(nrow(all_list_length)))>length(merged_table.list[[i]]$Rank)){
        p[i]<-ks.test(merged_table.list[[i]]$Rank,seq_len(nrow(all_list_length))[-merged_table.list[[i]]$Rank],alternative="greater")$p.value}else{p[i]<-1}
    }
  }
  for (i in 5:8){
    if(nrow(merged_table.list[[i]])==0){p[i]<-1}else{
      if(length(seq_len(nrow(all_list_length)))>length(merged_table.list[[i]]$Rank)){
      p[i]<-ks.test(merged_table.list[[i]]$Rank,seq_len(nrow(all_list_unsat))[-merged_table.list[[i]]$Rank],alternative="greater")$p.value}else{p[i]<-1}
    }
  }

  #adjust pvalues using qvalue method
  if(length(p)<=1){q<-p}else{q <- qvalue(p,lambda=0)$qvalues}

  #create the final table
  final_table<-data.frame(matrix(ncol=4,nrow=8))
  colnames(final_table) <- c("Classifier","Count.in.list", "p-value", "FDR.q-value")
  final_table$Classifier<-c("SCFA","MCFA","LCFA","VLCFA","Saturated","Monounsaturated","Diunsaturated","Polyunsaturated")
  for(i in 1:8){final_table$Count.in.list[i]<-nrow(merged_table.list[[i]])}
  final_table$`p-value`<-p
  final_table$`FDR.q-value`<-q

  return(final_table)


}

#' allchains.KS()
#'
#' Performs weighted Kolmogorov-Smirnov test on a <.allchains> object based on Ranking values contained in a rankingTable
#'
#' @param X any <RT.allchains> object
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#'
#' @details Peforms a ks test enrichment test on the chain parameters saved in a <.chain> object based on ranking values contained in a ranking table
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.KS(RTexample.allchains, cleaned.RTexample)
#'
#' @author Geremy Clair
#' @export
#'
allchains.KS<- function(X,rankingTable, order="ascending"){
  if(!ncol(X)==2 & !sum(colnames(X)[1:2]==c("Lipids","Chain"))==2){stop("X should be an <.allchains> object")}

  if(!is.data.frame(rankingTable)){stop("rankingTable should be a data frame of type <Ranking Table> containing the lipid list as first column and the ranking values as second column, we recommend that you use the function clean.RankingTable() to clean the table")}
  if(ncol(rankingTable)!=2){stop("rankingTable should contain two columns")}
  rankingTable[rankingTable==""]<-NA
  if(sum(is.na(rankingTable))>0){stop("rankingTable should not contain missing values")}

  if(missing(order)){order=="ascending"}
  if(!order %in% c("ascending", "descending")){stop("order should either be 'ascending' or 'descending'")}

  #remove the chain position from the LipidName
  X$Lipid<- gsub("_Chain.*","",X$Lipid)

  #Verify if all lipids are present in the rankingTable
  unique_lipid_X<-unique(X$Lipid)
  if(sum(unique_lipid_X %in% rankingTable[,1])!=length(unique_lipid_X)){stop("Some of the lipids contained in your <.allchains> object were not present in your ranking table")}

  #attach the rankingValues to X
  X$Ranking_values<- rankingTable[match(X$Lipid,rankingTable$Lipid),2]

  #order X by those values
  if(order=="ascending"){
  X<-X[order(X$Ranking_values),]
  }else{
  X<-X[order(X$Ranking_values,decreasing=TRUE),]
  }

  #add a rank
  X$KSRank<-1:nrow(X)

  #create an object containing the levels of the chains
  chain_levels<-unique(X$Chain)

  #create an object containing the
  list_by_chain<-list()
  for (i in 1:length(chain_levels)){
    list_by_chain[[i]]<-X[X$Chain==chain_levels[i],]
  }
  names(list_by_chain)<-chain_levels

  #realize the KS.tests
  p<-numeric()
  for (i in 1:length(chain_levels)){
    if(length(seq_len(nrow(X)))>length(list_by_chain[[i]]$KSRank)){
    p[i]<-ks.test(list_by_chain[[i]]$KSRank,seq_len(nrow(X))[-list_by_chain[[i]]$KSRank],alternative="greater")$p.value}else{p[i]<-1}
  }

  #adjust pvalues using qvalue method
  if(length(p)<=1){q<-p}else{q <- qvalue(p,lambda=0)$qvalues}

  #create the final table
  final_table<-data.frame(matrix(ncol=4,nrow= length(chain_levels)))
  colnames(final_table) <- c("Classifier","Count.in.list", "p-value", "FDR.q-value")
  final_table$Classifier<-chain_levels
  for(i in 1:length(chain_levels)){final_table$Count.in.list[i]<-nrow(list_by_chain[[i]])}
  final_table$`p-value`<-p
  final_table$`FDR.q-value`<-q

  return(final_table)
}


