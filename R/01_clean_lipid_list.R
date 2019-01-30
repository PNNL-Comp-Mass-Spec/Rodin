#' clean.lipid.list()
#'
#' Clean and prepare a character vector prior to input it in the lipid.miner function.
#'
#' @param X is a character vector containing the lipids names.
#'
#' @details We recommend to use the clean.lipid.list() prior to use the lipid.miner()
#' @details to maximize its compatibility.
#' @details This function removes missing values, usless space characters, and redundancy.
#' @details clean.lipid.list() informs you using warnings on the steps it has achieved.
#'
#' @return The R object outputed by this function is a character vector
#'
#' @examples clean.lipid.list(QueryExample)
#' @examples clean.lipid.list(UniverseExample)
#'
#' @author Geremy Clair
#' @export
#'

clean.lipid.list<-function(X){
  # if the data had more than one column -> remove the data in the other columns
  cleaned<- X
  if(ncol(cleaned)>0){
    cleaned<-as.character(unlist(cleaned[,1]))
    warning("Your lipid list format was inappropriate (it had more than one column)")
    warning("only the first column was kept as lipid identifier for the analysis")
  }

  #convert in character vector (if not already)
  if (typeof(cleaned)!="character"){
    cleaned<- as.character(unlist(cleaned[,1]))
  }

  #remove m
  if(sum(cleaned=="")>0){
    cleaned[cleaned==""]<-NA
    cleaned<- na.omit(cleaned)
    warning("Your lipid list was containing empty ids")
    warning("These ids were removed")
  }

  #Remove duplicates (original duplicates)
  if(sum(duplicated(cleaned))>0){
    n.duplicates<-sum(duplicated(cleaned))
    cleaned<-unique(cleaned)
    warning(paste("Your list was containing", n.duplicates, "duplicated lipid identifiers"))
    warning("the duplicated values were removed")
  }

  #Convert return delimited list into vector
  if(length(cleaned)==1&&grepl("\n",cleaned)){
    cleaned<- unlist(strsplit(cleaned,"\n"))
    warning("Your list was return delimited and was converted in an appropriate format")
    if(length(cleaned)>1&&grepl(";",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a semicolon or contains both semicolons and returns. Please ensure to use only one identifier for each lipid")
    }
    if(length(cleaned)>1&&grepl(",",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a comma or contains both comma and returns. Please ensure to use only one identifier for each lipid")
    }
  }

  #Convert return delimited list into vector
  if(length(cleaned)==1&&grepl("\t",cleaned)){
    cleaned<- unlist(strsplit(cleaned,"\t"))
    warning("Your list was tab delimited and was converted in an appropriate format")
    if(length(cleaned)>1&&grepl(";",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a semicolon or contains both semicolons and tab. Please ensure to use only one identifier for each lipid")
    }
    if(length(cleaned)>1&&grepl(",",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a comma or contains both comma and tab. Please ensure to use only one identifier for each lipid")
    }
  }

  #Convert comma delimited list into vector
  if(length(cleaned)==1&&grepl(",",cleaned)){
    cleaned<- unlist(strsplit(cleaned,","))
    warning("Your list was comma delimited and was converted in an appropriate format")
    if(length(cleaned)>1&&grepl(";",cleaned)){
      warning("Your list contain some elements with multiple identifiers separated by a semicolon or contains both semicolons and commas. Please ensure to use only one identifier for each lipid")
    }
  }

  #Convert semicolon delimited list into vector
  if(length(cleaned)==1&&grepl(";",cleaned)){
    cleaned<- unlist(strsplit(cleaned,";"))
    warning("Your list was semicolon delimited and was converted in an appropriate format")
  }

  #remove _CoE, _NEG, _POS, and spaces
  rm.list<- c("_CoE", "_Coe", "_coe", "_NEG", "_Neg", "_neg", "_POS", "_Pos", "_pos", " ")
  for (i in 1:length(rm.list)){
    cleaned<-gsub(rm.list[i],"",cleaned)
  }

  #Remove duplicates (final duplicates)
  if(sum(duplicated(cleaned))>0){
    n.duplicates<-sum(duplicated(cleaned))
    cleaned<-unique(cleaned)
    warning(paste("Your list was containing", n.duplicates, "duplicates"))
    warning("the duplicated values were removed")
  }

  #count identifiers
  warning(paste("After reformating, your", "lipid" , "list comprise", length(cleaned), "lipid(s)"))
  cleaned
}

#' clean.RankingTable()
#'
#' Clean and prepare a character vector prior to input it in the lipid.miner function.
#'
#' @param X is a ranking table
#'
#' @details We recommend to use the clean.RankingTable() prior to use the lipid.miner() on the first column of the rankingTable (containing the lipid list)
#' @details to maximize its compatibility.
#' @details This function removes missing values, usless space characters, and redundancy.
#' @details clean.lipid.list() informs you using warnings on the steps it has achieved.
#'
#' @return The R object outputed by this function is a data.frame with two columns the first one contains the lipids, the second one the ranked weights
#'
#' @examples cleaned.rankingTable<-clean.rankingTable(RankingTableExample)
#'
#' @author Geremy Clair
#' @export
#'

clean.rankingTable<-function(df){
  # if the data had more than one column -> remove the data in the other columns
    if(ncol(df)>2){
    df<-df[,1:2]
    warning("Your lipid list format was inappropriate (it had more than two columns)")
    warning("only the first column was kept as lipid identifier and the second columns as ranking values")
  }

  #ensure that the first column of df is of character type
  df[,1]<-as.character(t(df[,1]))

  #ensure that the secon column of df is of numeric type
  df[,2]<-as.numeric(t(df[,2]))

  #remove the rows containing NAs or missing values
  df[df==""]<-NA
  df[df=="NA"]<-NA
  if(sum(is.na(df))>0){
    warning("Your lipid list was containing empty ids or ranking values, these were removed from the list")
    df<-df[!rowSums(is.na(df))>0,]
  }

  #Remove ID duplicates
  if(sum(duplicated(df[,1]))>0){
    n.duplicates<-sum(duplicated(df))
    df<-unique(df)
    warning(paste("Your list was containing", n.duplicates, "only the first of the duplicated ID(s) were kept"))
  }

  #remove _CoE, _NEG, _POS, and spaces from the lipid names
  rm.list<- c("_CoE", "_Coe", "_coe", "_NEG", "_Neg", "_neg", "_POS", "_Pos", "_pos", " ")
  for (i in 1:length(rm.list)){
    cleaned<-gsub(rm.list[i],"",df)
  }

  #count identifiers
  warning(paste("After reformating, your", "lipid" , "list comprise", nrow(df), "lipid(s)"))

  #output the table
  return(df)

}

