#' lipid.miner()
#'
#' This function performs the text minning and extract lipid information from a <Query> or a <Universe> lipid list. It creates three new distinct R objects in the global environment: <name>.intact, <name>.chain and <name>.allchain
#'
#' @param X is character vector containing the lipids names
#' @param name is character string setting the name of the <name>.intact, <name>.chain and <name>.allchain object created
#' @param TGcollapse.rm is a boolean value. If TRUE, it removes the Triglycerides from the <name>.chain and <name>.allchains files when collapsed
#' @param output.list is a boolean value, which defaults to FALSE. If FALSE, the function returns objects as described below. If TRUE, the function instead returns a list of these same objects (minus the appended prefix '<name>.'), thus allowing the user to set the output of this function to a variable of their naming.
#'
#' @details We recommend to use the clean.lipid.list() function prior to use lipid.miner().
#' @details lipid.miner() mines the lipid names with the following format:
#' @details mainclass(chainLenght:doubleBond/chainLenght:doubleBond/etc.)
#' @details the mainclasses recognized by lipid.miner() are contained in the main.colors object
#' @details Any lipid that don't belong to one of recognized classes
#' @details will be considered as 'Uncategorized'.
#' @details Contact us if you need the inclusion of another lipid class
#'
#' @return This function return 3 objects named :
#' @return <name>.intact, <name>.chain and <name>.allchains
#' @return <name>.intact contains : category, main class, subclass, and various chain info
#' @return <name>.chain contains information regarding the length of the fatty acid chains
#' @return <name>.allchain contains the list of the individual fatty acid
#'
#' @examples lipid.miner(cleaned.queryExample, name= "QueryExample", TGcollapse.rm=TRUE)
#' @examples lipid.miner(cleaned.universeExample, name= "UniverseExample", TGcollapse.rm=TRUE)

#' @author Geremy Clair
#' @export

lipid.miner <- function(X,name="yourname", TGcollapse.rm=TRUE, output.list=FALSE){
  if(!is.character(X)){stop("X has to be a character vector")}
  if(!(is.character(name)&&length(name)==1)){"The output name is not properly defined"}
  if(!exists("TGcollapse.rm")){TGcollapse.rm<-TRUE}

  #remove the deuterated standards from the list
  X <- X[!grepl("\\(d5\\)", X)]

  #Create the output df
  df<-data.frame(matrix(NA,ncol = 27, nrow = length(X)))
  colnames(df)<- c("Lipid","Category", "Main class", "Sub class", "Total Number of Carbon", "Double Bonds", "Chain 1", "Chain 2", "Chain 3", "Chain 4", "N Carbon Chain1", "N Carbon Chain2", "N Carbon Chain3", "N Carbon Chain4", "DB Chain1", "DB Chain2","DB Chain3","DB Chain4","Lyso","SCFA","MCFA","LCFA","VLCFA", "Saturated", "Monounsaturated", "Diunsaturated", "Polyunsaturated" )
  df$Lipid <- X

  # Compute the main class
  X<- gsub("\\(3'-sulfo)Galbeta-Cer\\b","Sulfatide",X, ignore.case = TRUE)
  X<- gsub("M\\(IP\\)2C", "MIP2C", X) #rename M(IP)2C by M(IP)2C #PARENTHESIS PROBLEM FOR M(IP)2C
  X<- gsub(".*anandamide\\(", "EA(", X,ignore.case = TRUE)
  X<- gsub(".*CoenzymeQ.*", "CoQ(", X,ignore.case = TRUE)
  X<- gsub(".*CoQ.*", "CoQ(", X,ignore.case = TRUE)
  X<- gsub(".*Sulfatide", "Sulfatide(", X,ignore.case = TRUE)
  X<- gsub("Cholesterol", "Cholesterol(", X,ignore.case = TRUE)
  X<- gsub("Sulfatide\\(\\(","Sulfatide(",X,ignore.case = TRUE)
  X<- gsub("caritine\\(\\(","caritine(",X,ignore.case = TRUE)
  X<- gsub("PE-NMe2\\(","PE\\(NMe2\\-",X,ignore.case = TRUE)
  X<- gsub("PE-NMe\\(","PE\\(NMe\\-",X,ignore.case = TRUE)
  X<- gsub("\\,n\\-3","",X,ignore.case = TRUE)
  X<- gsub("\\,n\\-6","",X,ignore.case = TRUE)
  X<- gsub("\\,n\\-9","",X,ignore.case = TRUE)


  mainclass.list<- gsub("\\(.*", "",X)
  mainclass.list<-gsub("(.*)\\[.*", "\\1", mainclass.list)

  df$`Main class`<-mainclass.list

  # compute the lipid Category based on the main class
  #create the different lists
  FattyAcyls<-c("EA", "carnitine","FAHFA", "WE" )
  Prenol<- c("CoQ")
  Sterol<- c("Cholesterol","CE")
  Sphingolipids<- c("Cer","CerP","HexCer","SM","GM3","GD","PE-Cer","PI-Cer","Sulfatide","GM1","GM2","LacCer","GalCer","GD1", "GD3","GlcCer", "MIP2C","MIPC","GT3")
  Glycerophospholipids<-c("PE", "PC", "PS", "PI", "PG", "PA", "CL","PIP","PIP2", "PIP3")
  Glycerolipids<- c("TG","DG","MG","DGDG","MGDG","SQMG", "SQDG", "DGTS/A")


  #Fill the category section
  df$Category<- rep("Uncategorized", nrow(df))
  df$Category[df$`Main class` %in% Sphingolipids]<-"Sphingolipid"
  df$Category[df$`Main class` %in% Prenol]<-"Prenol"
  df$Category[df$`Main class` %in% Glycerolipids]<-"Glycerolipid"
  df$Category[df$`Main class` %in% Sterol]<-"Sterol"
  df$Category[df$`Main class` %in% FattyAcyls]<-"FattyAcyls"
  df$Category[df$`Main class` %in% Glycerophospholipids]<-"Glycerophospholipid"

  # Modify the input to remove unwanted strings and introduce the <oxidized> and <hydroxy>
  #remove text after ";"
  X<-gsub("\\;.*","",X)

  #Remove the text inside brackets from X
  X<-gsub("\\[[^\\]]*\\]", "", X)

  #transform the (CHO) (OH) and (COOH) into <oxidized>
  X<-gsub("\\(CHO\\)", "<oxidized>", X) #CHO
  X<-gsub("\\(COOH\\)", "<oxidized>", X) # (COOH)
  X<-gsub("\\(12OH\\)","<oxidized>",X) #(12OH[S])
  X<-gsub("\\(14OH\\)","<oxidized>",X) #(12OH[S])
  X<-gsub("\\(15OH\\)","<oxidized>",X) #(12OH[S])
  X<-gsub("\\(5OH\\)","<oxidized>",X) #(5OH[S])
  X<-gsub("\\(2OH\\)", "<hydroxy>", X) ##(2OH[S])

  #modify the CardioLipids for proper minning
  X<- gsub("1\\'\\-\\[","",X)#remove the 1'-[
  X<- gsub("\\,3\\'\\-\\[","\\/",X)#remove the 3'-[
  X<- gsub("\\]","",X)#remove the ]

  #compute the chains, the N of carbon, the N of DB, the total N of C and of Db by mining between the outermost parenthesis
  #extract the fatty acid chains
  chains.list<-X
  chains.list<- gsub("(.*)\\).*", "\\1", chains.list) #after last parenthesis
  chains.list<- sub(".*\\(.*?", "", chains.list) #before the first parenthesis
  chains.list<-gsub("NMe[2]\\-","",chains.list)
  chains.list<-gsub("NMe\\-","",chains.list)

  #split using '/' character
  chains.list<-as.list(strsplit(as.character(chains.list), '/'))#split into
  chains.list<-lapply(chains.list, `length<-`,4)
  chains.df <- data.frame(matrix(unlist(chains.list), nrow=length(X),byrow=T))


  df$`Chain 1`<- chains.df$X1
  df$`Chain 2`<- chains.df$X2
  df$`Chain 3`<- chains.df$X3
  df$`Chain 4`<- chains.df$X4

  #REMOVE 0:0 in DG and MG in df
  DGMG<-c("DG","MG")
  df$`Chain 1`[grepl("^0\\:0",df$`Chain 1`)&df$`Main class`%in% DGMG]<-NA
  df$`Chain 2`[grepl("^0\\:0",df$`Chain 2`)&df$`Main class`%in% DGMG]<-NA
  df$`Chain 3`[grepl("^0\\:0",df$`Chain 3`)&df$`Main class`%in% DGMG]<-NA
  df$`Chain 4`[grepl("^0\\:0",df$`Chain 3`)&df$`Main class`%in% DGMG]<-NA

  #Bring this back in the chain object
  chains<-paste(df$`Chain 1`,"/",df$`Chain 2`,"/",df$`Chain 3`,"/",df$`Chain 4`, sep="")

  #remove O- P- d t <ox> <hydroxyle>, etc
  chains<- gsub("<oxidized>", "",chains)
  chains<- gsub("<hydroxy>", "",chains)
  chains<- gsub("5-O-","O-",chains)
  chains<- gsub("6-O-","O-",chains)
  chains<- gsub("7-O-","O-",chains)
  chains<- gsub("8-O-","O-",chains)
  chains<- gsub("9-O-","O-",chains)
  chains<- gsub("10-O-","O-",chains)
  chains<- gsub("11-O-","O-",chains)
  chains<- gsub("12-O-","O-",chains)
  chains<- gsub("13-O-","O-",chains)
  chains<- gsub("O-", "",chains)
  chains<- gsub("P-", "",chains)
  chains<- gsub("d", "",chains)
  chains<- gsub("t", "",chains)
  chains<- gsub("h", "",chains)
  chains<- gsub("Me\\)", "",chains)
  chains<- gsub("m", "",chains)

  #chains<- gsub("\\)Me", "",chains)

  #remove internal parenthesis
  chains.list<- gsub("\\s*\\([^\\)]+\\)","",chains)
  chains.list<- gsub("\\(5E,9E","",chains.list)



  #split the string
  chains.list<-as.list(strsplit(as.character(chains.list), '/'))#split into list
  chains.list<-lapply(chains.list, `length<-`,4)
  chains.df <- data.frame(matrix(unlist(chains.list), nrow=length(X),byrow=T))

  #N carbon Chains 1->4
  df$`N Carbon Chain1` <- gsub(":.*","",chains.df$X1)
  df$`N Carbon Chain2` <- gsub(":.*","",chains.df$X2)
  df$`N Carbon Chain3` <- gsub(":.*","",chains.df$X3)
  df$`N Carbon Chain4` <- gsub(":.*","",chains.df$X4)

  #DB Chains 1->4
  df$`DB Chain1` <- gsub(".*:","",chains.df$X1)
  df$`DB Chain2` <- gsub(".*:","",chains.df$X2)
  df$`DB Chain3` <- gsub(".*:","",chains.df$X3)
  df$`DB Chain4` <- gsub(".*:","",chains.df$X4)

  #create a function to make real NAs
  make.true.NA <- function(x) if(is.character(x)||is.factor(x)){
    is.na(x) <- x=="NA"; x} else {
      x}
  df[] <- lapply(df, make.true.NA)

  #Total number of carbon
  suppressWarnings(chains.df$X1<-as.numeric(df$`N Carbon Chain1`))
  suppressWarnings(chains.df$X2<-as.numeric(df$`N Carbon Chain2`))
  suppressWarnings(chains.df$X3<-as.numeric(df$`N Carbon Chain3`))
  suppressWarnings(chains.df$X4<-as.numeric(df$`N Carbon Chain4`))

  df$`Total Number of Carbon`<- rowSums (chains.df, na.rm = T)

  #Total number of DB
  chains.df$X1<-as.numeric(df$`DB Chain1`)
  chains.df$X2<-as.numeric(df$`DB Chain2`)
  chains.df$X3<-as.numeric(df$`DB Chain3`)
  chains.df$X4<-as.numeric(df$`DB Chain4`)


  df$`Double Bonds`<- rowSums (chains.df, na.rm = T)

  remove(chains.df)

  #comupte the subclass by minning between the parenthesis and adding the subclass to the main class
  #O P d t DiEther ox and Lyso
  subclass.list<- gsub("<oxidized>", "<@>", X)
  subclass.list<- gsub("<hydroxy>", "<#>", subclass.list)
  subclass.list<- gsub("NMe[2]\\-", "\\|", subclass.list)
  subclass.list<- gsub("NMe\\-", "\\%", subclass.list)
  subclass.list<- gsub("(.*)\\).*", "\\1", subclass.list) #after last parenthesis
  subclass.list<- sub(".*\\(.*?", "", subclass.list) #before the first parenthesis
  subclass.list<- gsub("\\d", "",subclass.list) #remove numbers
  subclass.list<- gsub("\\:","",subclass.list)#remove the :
  subclass.list<- gsub("\\/","",subclass.list)#remove the /
  subclass.list<- gsub("\\(","",subclass.list)#remove the (
  subclass.list<- gsub("\\Z","",subclass.list)#remove the Z
  subclass.list<- gsub("\\E","",subclass.list)#remove the E
  subclass.list<- gsub("\\,","",subclass.list)#remove the ,
  subclass.list<- gsub("\\)","",subclass.list)#remove the )
  subclass.list<- gsub("\\A","",subclass.list)#remove the )
  subclass.list<- gsub("\\-O","O",subclass.list)#remove the -O replace by O
  subclass.list<- gsub("OO\\-","diester",subclass.list)#remove the -O replace by O
  subclass.list<-gsub( "Cholesterol","", subclass.list)
  subclass.list<-gsub( "CSulfatide","", subclass.list)
  subclass.list<-gsub( "CoQ","",subclass.list)

  #bring back oxidized and hydroxy
  subclass.list<-gsub( "<@>","<oxidized>", subclass.list)
  subclass.list<-gsub( "<#>","<hydroxy>", subclass.list)
  subclass.list<- gsub( "\\|","NMe2\\-", subclass.list)
  subclass.list<- gsub( "\\%","NMe\\-", subclass.list)


  #create the subclass list
  df$`Sub class`<- paste(df$`Main class`,"(",subclass.list,sep="")

  #Lyso (in the subclass and Lyso column)
  noNA<-data.frame(df[,11:18])
  noNA[is.na(noNA)]<- 123456789

  df[noNA$N.Carbon.Chain1==0 & noNA$`DB.Chain1`==0,4]<-paste("L",df[noNA$N.Carbon.Chain1==0 & noNA$`DB.Chain1`==0,4],sep="")
  df[noNA$N.Carbon.Chain2==0 & noNA$`DB.Chain2`==0,4]<-paste("L",df[noNA$N.Carbon.Chain2==0 & noNA$`DB.Chain2`==0,4],sep="")
  df[noNA$N.Carbon.Chain3==0 & noNA$`DB.Chain3`==0,4]<-paste("L",df[noNA$N.Carbon.Chain3==0 & noNA$`DB.Chain3`==0,4],sep="")
  df[noNA$N.Carbon.Chain4==0 & noNA$`DB.Chain4`==0,4]<-paste("L",df[noNA$N.Carbon.Chain4==0 & noNA$`DB.Chain4`==0,4],sep="")

  remove(noNA)

  # Chains, Calculate the chains properties
  #this can probably be simplified for faster calculation

  chains.df<-cbind(df$Lipid,df[,7:27])
  chains.df[,6:13]<-df[,11:18]
  chains.df[,6:13][is.na(chains.df[,6:13])]<-"-1"
  chains.df[,14:22]<-0
  colnames(chains.df)[1]<-"Lipid"

  for(i in 6:13){
    chains.df[,i]<-as.numeric(chains.df[,i])
  }

  for (i in 1:nrow(chains.df)){
    if(chains.df$`N Carbon Chain1`[i]==0){chains.df$Lyso[i]<-chains.df$Lyso[i]+1}
    if(chains.df$`N Carbon Chain2`[i]==0){chains.df$Lyso[i]<-chains.df$Lyso[i]+1}
    if(chains.df$`N Carbon Chain3`[i]==0){chains.df$Lyso[i]<-chains.df$Lyso[i]+1}
    if(chains.df$`N Carbon Chain4`[i]==0){chains.df$Lyso[i]<-chains.df$Lyso[i]+1}
    if(chains.df$`N Carbon Chain1`[i]>=1 && chains.df$`N Carbon Chain1`[i]<=5 ){chains.df$SCFA[i]<-chains.df$SCFA[i]+1}
    if(chains.df$`N Carbon Chain2`[i]>=1 && chains.df$`N Carbon Chain2`[i]<=5 ){chains.df$SCFA[i]<-chains.df$SCFA[i]+1}
    if(chains.df$`N Carbon Chain3`[i]>=1 && chains.df$`N Carbon Chain3`[i]<=5 ){chains.df$SCFA[i]<-chains.df$SCFA[i]+1}
    if(chains.df$`N Carbon Chain4`[i]>=1 && chains.df$`N Carbon Chain4`[i]<=5 ){chains.df$SCFA[i]<-chains.df$SCFA[i]+1}
    if(chains.df$`N Carbon Chain1`[i]>5 && chains.df$`N Carbon Chain1`[i]<=12 ){chains.df$MCFA[i]<-chains.df$MCFA[i]+1}
    if(chains.df$`N Carbon Chain2`[i]>5 && chains.df$`N Carbon Chain2`[i]<=12 ){chains.df$MCFA[i]<-chains.df$MCFA[i]+1}
    if(chains.df$`N Carbon Chain3`[i]>5 && chains.df$`N Carbon Chain3`[i]<=12 ){chains.df$MCFA[i]<-chains.df$MCFA[i]+1}
    if(chains.df$`N Carbon Chain4`[i]>5 && chains.df$`N Carbon Chain4`[i]<=12 ){chains.df$MCFA[i]<-chains.df$MCFA[i]+1}
    if(chains.df$`N Carbon Chain1`[i]>12 && chains.df$`N Carbon Chain1`[i]<=21 ){chains.df$LCFA[i]<-chains.df$LCFA[i]+1}
    if(chains.df$`N Carbon Chain2`[i]>12 && chains.df$`N Carbon Chain2`[i]<=21 ){chains.df$LCFA[i]<-chains.df$LCFA[i]+1}
    if(chains.df$`N Carbon Chain3`[i]>12 && chains.df$`N Carbon Chain3`[i]<=21 ){chains.df$LCFA[i]<-chains.df$LCFA[i]+1}
    if(chains.df$`N Carbon Chain4`[i]>12 && chains.df$`N Carbon Chain4`[i]<=21 ){chains.df$LCFA[i]<-chains.df$LCFA[i]+1}
    if(chains.df$`N Carbon Chain1`[i]>=22){chains.df$VLCFA[i]<-chains.df$VLCFA[i]+1}
    if(chains.df$`N Carbon Chain2`[i]>=22){chains.df$VLCFA[i]<-chains.df$VLCFA[i]+1}
    if(chains.df$`N Carbon Chain3`[i]>=22){chains.df$VLCFA[i]<-chains.df$VLCFA[i]+1}
    if(chains.df$`N Carbon Chain4`[i]>=22){chains.df$VLCFA[i]<-chains.df$VLCFA[i]+1}
    if(chains.df$`DB Chain1`[i]==0){chains.df$Saturated[i]<-chains.df$Saturated[i]+1}
    if(chains.df$`DB Chain2`[i]==0){chains.df$Saturated[i]<-chains.df$Saturated[i]+1}
    if(chains.df$`DB Chain3`[i]==0){chains.df$Saturated[i]<-chains.df$Saturated[i]+1}
    if(chains.df$`DB Chain4`[i]==0){chains.df$Saturated[i]<-chains.df$Saturated[i]+1}
    if(chains.df$`DB Chain1`[i]==1){chains.df$Monounsaturated[i]<-chains.df$Monounsaturated[i]+1}
    if(chains.df$`DB Chain2`[i]==1){chains.df$Monounsaturated[i]<-chains.df$Monounsaturated[i]+1}
    if(chains.df$`DB Chain3`[i]==1){chains.df$Monounsaturated[i]<-chains.df$Monounsaturated[i]+1}
    if(chains.df$`DB Chain4`[i]==1){chains.df$Monounsaturated[i]<-chains.df$Monounsaturated[i]+1}
    if(chains.df$`DB Chain1`[i]==2){chains.df$Diunsaturated[i]<-chains.df$Diunsaturated[i]+1}
    if(chains.df$`DB Chain2`[i]==2){chains.df$Diunsaturated[i]<-chains.df$Diunsaturated[i]+1}
    if(chains.df$`DB Chain3`[i]==2){chains.df$Diunsaturated[i]<-chains.df$Diunsaturated[i]+1}
    if(chains.df$`DB Chain4`[i]==2){chains.df$Diunsaturated[i]<-chains.df$Diunsaturated[i]+1}
    if(chains.df$`DB Chain1`[i]>=3){chains.df$Polyunsaturated[i]<-chains.df$Polyunsaturated[i]+1}
    if(chains.df$`DB Chain2`[i]>=3){chains.df$Polyunsaturated[i]<-chains.df$Polyunsaturated[i]+1}
    if(chains.df$`DB Chain3`[i]>=3){chains.df$Polyunsaturated[i]<-chains.df$Polyunsaturated[i]+1}
    if(chains.df$`DB Chain4`[i]>=3){chains.df$Polyunsaturated[i]<-chains.df$Polyunsaturated[i]+1}
  }
  df[,19:27]<-chains.df[,14:22]



  if(TGcollapse.rm==TRUE){
    #remove the chains that are from collapsed TGs
    ################################################################
    ###ADD AN OPTION BUTTON HERE
    ################################################################
    chains.df<- chains.df[!(chains.df$`N Carbon Chain2`==-1&df$`Main class`=="TG"),]
  }

  chains.df[chains.df==-1]<- NA

  ################################################################

  #create the last object that combine the existing chains in a vector
  allchains<- data.frame(matrix(nrow = 4*nrow(chains.df), ncol=2))
  allchains[,1]<- c(paste(chains.df$Lipid,"_Chain1",sep=""),paste(chains.df$Lipid,"_Chain2",sep=""),paste(chains.df$Lipid,"_Chain3",sep=""),paste(chains.df$Lipid,"_Chain4",sep=""))
  allchains[,2]<-c(as.character(chains.df$`Chain 1`),as.character(chains.df$`Chain 2`),as.character(chains.df$`Chain 3`),as.character(chains.df$`Chain 4`))
  colnames(allchains)<- c("Lipid", "Chain")
  allchains<- allchains[!is.na(allchains$Chain),]

  #finalize the object chains
  chains.df<- cbind(chains.df[,1],chains.df[,15:22])
  colnames(chains.df)[1]<- "Lipid"

  #if unknown main.class/subclass replace with Uncategorized
  df$`Main class`[!df$`Main class` %in% main.colors$MainClass]<-"Uncategorized"
  df$`Sub class`[!df$`Sub class` %in% sub.colors$SubClass]<-"Uncategorized"

  # Final  output
  if(output.list == FALSE){
    assign(paste(name,".intact",sep=""),df[,1:18], envir = globalenv())
    assign(paste(name,".chain",sep=""),data.frame(chains.df), envir = globalenv())
    assign(paste(name,".allchains",sep=""),data.frame(allchains),envir = globalenv())
  }else{
    # output.list == TRUE #
    return(list(intact = df[,1:18], chain = data.frame(chains.df), allchains = data.frame(allchains)))
  }

}
