#' chain.length.pie
#'
#' This function generate a pie chart indicating the distribution of the chain lengths within 4 chain lengths categories.
#' @param X is a ~.chain object
#' @details This function will place the chains in four categories based on their length: 1-SCFA (Short-Chain Fatty Acid, 5 or fewer carbons), 2-MCFA (Medium-Chain Fatty Acid, 6-12 carbons), 3-LCFA(Long Chain Fatty Acid, 13-21 carbons), and 4-VLCFA(Very Long Chain Fatty Acid, >22 carbons). This function require ggplot2 to generate the graphics.
#' @examples chain.length.pie(queryExample.chain)
#' @author Geremy Clair
#' @export

chain.length.pie<- function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  chain.length.percentage <- colSums(X[2:5])
  chain.length.percentage <- data.frame(chain.length.percentage/sum(chain.length.percentage)*100)
  colnames(chain.length.percentage)<-c("percentage")
  chain.length.percentage$len<- as.character(row.names(chain.length.percentage))
  chain.length.percentage$percentage<-round(chain.length.percentage$percentage, digits = 1)
  chain.length.percentage$tag<-""
  for (i in 1:4)
  {
    if(chain.length.percentage$percentage[i]!=0)
    {
      chain.length.percentage$tag[i]<- paste(chain.length.percentage$len[i],"(",chain.length.percentage$percentage[i],"%)",sep="")
    }
  }

  ggplot(chain.length.percentage, aes(x="", y=percentage,label = tag, fill=factor(len, levels = c("VLCFA","LCFA","MCFA","SCFA")))) +
    geom_bar(width = 1, stat = "identity") +
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    scale_fill_manual( values=c("#17aeae","#6bbfb9", "#a7dbd9","#d7f4f0" )) +
    theme(legend.title=element_blank()) +
    coord_polar("y", start=0)+
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(),legend.position="bottom")

  }

#' chain.unsat.pie
#'
#' This function generate a pie chart indicating the distribution of the chain based on their unsaturation count within 4 chain categories.
#' @param X is a ~.chain object
#' @details This function will place the chains in four categories based on their number of unsaturation: 1-Saturated, 2-Monounsaturated, 3-Diunsaturated, and 4-Polyunsaturated(i.e. >2 unsaturations). This function require ggplot2 to generate the graphics. This function require ggplot2 to generate the graphics.
#' @examples chain.unsat.pie(queryExample.chain)
#' @author Geremy Clair
#' @export

chain.unsat.pie<-function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  unsat.percentage <- colSums(X[6:9])
  unsat.percentage <- data.frame(unsat.percentage/sum(unsat.percentage)*100)
  colnames(unsat.percentage)<-c("percentage")
  unsat.percentage$len<- as.character(row.names(unsat.percentage))
  unsat.percentage$percentage<-round(unsat.percentage$percentage, digits = 1)
  unsat.percentage$tag<-""
  for (i in 1:4)
  {
    if(unsat.percentage$percentage[i]!=0)
    {
      unsat.percentage$tag[i]<- paste(unsat.percentage$len[i],"(",unsat.percentage$percentage[i],"%)",sep="")
    }
  }

  ggplot(unsat.percentage, aes(x="", y=percentage,label = tag, fill=factor(len, levels = c("Polyunsaturated","Diunsaturated","Monounsaturated","Saturated"))))+
    geom_bar(width = 1, stat = "identity") +
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    scale_fill_manual(values=c("#ff281d", "#ff6840", "#ff9750", "#ffc950")) +
    coord_polar("y", start=0)+
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(), legend.position= "bottom")

}

#' chain.length.stack
#'
#' This function generate a stacked barplot indicating the distribution of the chain lengths within 4 chain lengths categories.
#' @param X is a ~.chain object
#' @details This function will place the chains in four categories based on their length: 1-SCFA (Short-Chain Fatty Acid, 5 or fewer carbons), 2-MCFA (Medium-Chain Fatty Acid, 6-12 carbons), 3-LCFA(Long Chain Fatty Acid, 13-21 carbons), and 4-VLCFA(Very Long Chain Fatty Acid, >22 carbons). This function require ggplot2 to generate the graphics.
#' @examples chain.length.stack(queryExample.chain)
#' @author Geremy Clair
#' @export

chain.length.stack<-function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  chain.length.percentage <- colSums(X[2:5])
  chain.length.percentage <- data.frame(chain.length.percentage/sum(chain.length.percentage)*100)
  colnames(chain.length.percentage)<-c("percentage")
  chain.length.percentage$len<- as.character(row.names(chain.length.percentage))
  chain.length.percentage$percentage<-round(chain.length.percentage$percentage, digits = 1)
  chain.length.percentage$tag<-""
  for (i in 1:4)
  {
    if(chain.length.percentage$percentage[i]!=0)
    {
      chain.length.percentage$tag[i]<- paste(chain.length.percentage$len[i],"(",chain.length.percentage$percentage[i],"%)",sep="")
    }
  }
  ggplot(chain.length.percentage, aes(x="", y=percentage,label = tag, fill=factor(len, levels = c("VLCFA","LCFA","MCFA","SCFA"))))+
    geom_bar(width = 1, stat = "identity") +
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    scale_fill_manual( values=c("#17aeae","#6bbfb9", "#a7dbd9","#d7f4f0" )) +
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.position="bottom",legend.title=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank())

}

#' chain.unsat.stack
#'
#' This function generate a stacked bar plot indicating the distribution of the chain based on their unsaturation count within 4 chain categories.
#' @param X is a ~.chain object
#' @details This function will place the chains in four categories based on their number of unsaturation: 1-Saturated, 2-Monounsaturated, 3-Diunsaturated, and 4-Polyunsaturated(i.e. >2 unsaturations).This function require ggplot2 to generate the graphics.
#' @examples chain.unsat.stack(queryExample.chain)
#' @author Geremy Clair
#' @export

chain.unsat.stack<-function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  unsat.percentage <- colSums(X[6:9])
  unsat.percentage <- data.frame(unsat.percentage/sum(unsat.percentage)*100)
  colnames(unsat.percentage)<-c("percentage")
  unsat.percentage$len<- as.character(row.names(unsat.percentage))
  unsat.percentage$percentage<-round(unsat.percentage$percentage, digits = 1)
  unsat.percentage$tag<-""
  for (i in 1:4)
  {
    if(unsat.percentage$percentage[i]!=0)
    {
      unsat.percentage$tag[i]<- paste(unsat.percentage$len[i],"(",unsat.percentage$percentage[i],"%)",sep="")
    }
  }

  ggplot(unsat.percentage, aes(x="", y=percentage,label = tag, fill=factor(len, levels = c("Polyunsaturated","Diunsaturated","Monounsaturated","Saturated"))))+
  geom_bar(width = 1, stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values=c("#ff281d", "#ff6840", "#ff9750", "#ffc950")) +
  theme_minimal(base_size = 11, base_family = "")+
  theme(legend.position="bottom",legend.title=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank())

}

#' intact.cat.pie
#'
#' This script is meant to produce a pie chart depicting the distribution of the lipid of a given ~.intact object within the different lipid categories.
#' @param X is a ~.intact object
#' @details This function will use the Categories included in the ~.intact object to generate a pie chart of the lipid distribution in a given list. This function require ggplot2 to generate the graphics.
#' @examples intact.cat.pie(queryExample.intact)
#' @author Geremy Clair
#' @export

intact.cat.pie<- function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  cat<-data.frame(matrix(ncol=6,nrow=length(unique(X$Category))))
  colnames(cat)<- c("Index","Category","Count","Percentage","Color","TxtColor")
  cat$Category<- unique(X$Category)
  for (i in 1:nrow(cat)){
    cat$Count[i]<-sum( X$Category==cat$Category[i])
    cat$Percentage[i]<-round(cat$Count[i]/nrow(X)*100,2)
  }
  cat$Index<- cat.colors[match(cat$Category, cat.colors$Category), 2]
  cat$Color<- cat.colors[match(cat$Category, cat.colors$Category), 3]
  cat$TxtColor<- cat.colors[match(cat$Category, cat.colors$Category), 4]
  cat<-cat[order(cat$Index),]
  cat$tag<- paste(cat$Category,"(",cat$Percentage,"%)",sep="")

  cat$Category<- factor(cat$Category, levels = cat$Category[order(cat$Index)])
  ord<- rev(cat.colors$Category)

  ggplot(cat,aes("", Percentage, label= tag,fill=factor(Category, levels = ord)))+
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values=rev(as.character(cat$Color)))+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    coord_polar("y", start=0)+
    theme(legend.title=element_blank()) +
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(),legend.position="none")

}

#' intact.cat.stack
#'
#' This script is meant to produce a stacked bar plot depicting the distribution of the lipid of a given ~.intact object within the different lipid categories.
#' @param X is a ~.intact object
#' @details This function will use the Categories included in the ~.intact object to generate a stacked bar plot of the lipid distribution in a given list. This function require ggplot2 to generate the graphics.
#' @examples intact.cat.pie(queryExample.intact)
#' @author Geremy Clair
#' @export

intact.cat.stack<- function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  cat<-data.frame(matrix(ncol=6,nrow=length(unique(X$Category))))
  colnames(cat)<- c("Index","Category","Count","Percentage","Color","TxtColor")
  cat$Category<- unique(X$Category)
  for (i in 1:nrow(cat)){
    cat$Count[i]<-sum(cat$Category[i]==X$Category)
    cat$Percentage[i]<-round(cat$Count[i]/nrow(X)*100,2)
  }
  cat$Index<- cat.colors[match(cat$Category, cat.colors$Category), 2]
  cat$Color<- cat.colors[match(cat$Category, cat.colors$Category), 3]
  cat$TxtColor<- cat.colors[match(cat$Category, cat.colors$Category), 4]
  cat<-cat[order(cat$Index),]
  cat$tag<- paste(cat$Category,"(",cat$Percentage,"%)",sep="")

  cat$Category<- factor(cat$Category, levels = cat$Category[order(cat$Index)])
  ord<- cat.colors$Category

  ggplot(cat,aes("", Percentage, label= tag,fill=factor(Category, levels = ord)))+
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values=as.character(cat$Color))+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    theme(legend.title=element_blank()) +
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(),legend.position="none")

}

#' intact.main.pie
#'
#' This script is meant to produce a pie chart depicting the distribution of the lipid of a given ~.intact object within the different lipid main classes.
#' @param X is a ~.intact object
#' @details This function will use the Main classes included in the ~.intact object to generate a pie chart of the lipid distribution in a given list. This function require ggplot2 to generate the graphics.
#' @examples intact.main.pie(queryExample.intact)
#' @author Geremy Clair
#' @export

intact.main.pie<- function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  main<-data.frame(matrix(ncol=6,nrow=length(unique(X$`Main class`))))
  colnames(main)<- c("Index","MainClass","Count","Percentage","Color","TxtColor")
  main$MainClass<- unique(X$`Main class`)
  for (i in 1:nrow(main)){
    main$Count[i]<-sum(main$MainClass[i]==X$`Main class`)
    main$Percentage[i]<-round(main$Count[i]/nrow(X)*100,2)
  }
  main$Index<- main.colors[match(main$MainClass, main.colors$MainClass), 2]
  main$Color<- main.colors[match(main$MainClass, main.colors$MainClass), 3]
  main$TxtColor<- main.colors[match(main$MainClass, main.colors$MainClass), 4]
  main<-main[order(main$Index),]
  main$tag<- paste(main$MainClass,"(",main$Percentage,"%)",sep="")

  main$MainClass<- factor(main$MainClass, levels = main$MainClass[order(main$Index)])
  ord<- rev(main.colors$MainClass)

  ggplot(main,aes("", Percentage, label= tag,fill=factor(MainClass, levels = ord)))+
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values=rev(as.character(main$Color)))+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    coord_polar("y", start=0)+
    theme(legend.title=element_blank()) +
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(),legend.position="none")

}

#' intact.main.stack
#'
#' This script is meant to produce a stacked barplot depicting the distribution of the lipid of a given ~.intact object within the different lipid main classes.
#' @param X is a ~.intact object
#' @details This function will use the Main classes included in the ~.intact object to generate a stacked barplot  of the lipid distribution in a given list. This function require ggplot2 to generate the graphics.
#' @examples intact.main.stack(queryExample.intact)
#' @author Geremy Clair
#' @export

intact.main.stack<- function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  main<-data.frame(matrix(ncol=6,nrow=length(unique(X$`Main class`))))
  colnames(main)<- c("Index","MainClass","Count","Percentage","Color","TxtColor")
  main$MainClass<- unique(X$`Main class`)
  for (i in 1:nrow(main)){
    main$Count[i]<-sum(main$MainClass[i]==X$`Main class`)
    main$Percentage[i]<-round(main$Count[i]/nrow(X)*100,2)
  }
  main$Index<- main.colors[match(main$MainClass, main.colors$MainClass), 2]
  main$Color<- main.colors[match(main$MainClass, main.colors$MainClass), 3]
  main$TxtColor<- main.colors[match(main$MainClass, main.colors$MainClass), 4]
  main<-main[order(main$Index),]
  main$tag<- paste(main$MainClass,"(",main$Percentage,"%)",sep="")

  main$MainClass<- factor(main$MainClass, levels = main$MainClass[order(main$Index)])
  ord<- main.colors$MainClass

  ggplot(main,aes("", Percentage, label= tag,fill=factor(MainClass, levels = ord)))+
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values=as.character(main$Color))+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    theme(legend.title=element_blank()) +
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(),legend.position="none")

}

#' intact.sub.pie
#'
#' This script is meant to produce a pie chart depicting the distribution of the lipid of a given ~.intact object within the different lipid sub classes.
#' @param X is a ~.intact object
#' @details This function will use the Sub classes included in the ~.intact object to generate a pie chart of the lipid distribution in a given list. This function require ggplot2 to generate the graphics.
#' @examples intact.sub.pie(queryExample.intact)
#' @author Geremy Clair
#' @export

intact.sub.pie<- function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  sub<-data.frame(matrix(ncol=6,nrow=length(unique(X$`Sub class`))))
  colnames(sub)<- c("Index","SubClass","Count","Percentage","Color","TxtColor")
  sub$SubClass<- unique(X$`Sub class`)
  subspecial<- gsub("\\(","\\\\(",sub$SubClass)
  for (i in 1:nrow(sub)){
    sub$Count[i]<-length(grep(paste("^",subspecial[i],"$",sep=""), X$`Sub class`))
    sub$Percentage[i]<-round(sub$Count[i]/nrow(X)*100,2)
  }
  sub$Index<- sub.colors[match(sub$SubClass, sub.colors$SubClass), 2]
  sub$Color<- sub.colors[match(sub$SubClass, sub.colors$SubClass), 3]
  sub$TxtColor<- sub.colors[match(sub$SubClass, sub.colors$SubClass), 4]
  sub<-sub[order(sub$Index),]
  sub$tag<- paste(sub$SubClass,"_(",sub$Percentage,"%)",sep="")

  sub$SubClass<- factor(sub$SubClass, levels = sub$SubClass[order(sub$Index)])
  ord<- rev(sub.colors$SubClass)

  ggplot(sub,aes("", Percentage, label= tag,fill=factor(SubClass, levels = ord)))+
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values=rev(as.character(sub$Color)))+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    coord_polar("y", start=0)+
    theme(legend.title=element_blank()) +
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(),legend.position="none")

}

#' intact.sub.stack
#'
#' This script is meant to produce a stacked barplot depicting the distribution of the lipid of a given ~.intact object within the different lipid sub classes.
#' @param X is a ~.intact object
#' @details This function will use the Sub classes included in the ~.intact object to generate a stacked barplot of the lipid distribution in a given list. This function require ggplot2 to generate the graphics.
#' @examples intact.sub.stack(queryExample.intact)
#' @author Geremy Clair
#' @export

intact.sub.stack<- function(X){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }
  sub<-data.frame(matrix(ncol=6,nrow=length(unique(X$`Sub class`))))
  colnames(sub)<- c("Index","SubClass","Count","Percentage","Color","TxtColor")
  sub$SubClass<- unique(X$`Sub class`)
  subspecial<- gsub("\\(","\\\\(",sub$SubClass)
  for (i in 1:nrow(sub)){
    sub$Count[i]<-length(grep(paste("^",subspecial[i],"$",sep=""), X$`Sub class`))
    sub$Percentage[i]<-round(sub$Count[i]/nrow(X)*100,2)
  }
  sub$Index<- sub.colors[match(sub$SubClass, sub.colors$SubClass), 2]
  sub$Color<- sub.colors[match(sub$SubClass, sub.colors$SubClass), 3]
  sub$TxtColor<- sub.colors[match(sub$SubClass, sub.colors$SubClass), 4]
  sub<-sub[order(sub$Index),]
  sub$tag<- paste(sub$SubClass,"_(",sub$Percentage,"%)",sep="")

  sub$SubClass<- factor(sub$SubClass, levels = sub$SubClass[order(sub$Index)])
  ord<- sub.colors$SubClass

  ggplot(sub,aes("", Percentage, label= tag,fill=factor(SubClass, levels = ord)))+
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values=as.character(sub$Color))+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    theme(legend.title=element_blank()) +
    theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank(),legend.position="none")

}

#' allchains.barplot
#'
#' this script will generate a barplot to compare the proportion of the specific chains contained in a Query ~.allchains object and in a ~.allchains Universe object.
#' @param X is a <.allchains> object
#' @param Y is a <.allchains> object,
#' @details This barplot will show the percentage of the different unique chains contained in one or two <.allchains> objects. This function require ggplot2 to generate the graphics.
#' @examples allchains.barplot(queryExample.allchains,universeExample.allchains)
#' @author Geremy Clair
#' @export

allchains.barplot<- function(X,Y){
  if (!require('ggplot2')) {
    stop('The package ggplot2 is required to run this function')
  }

  if(!missing(Y)){
  unique.chain<- unique(Y$Chain)
  query.table<-data.frame(matrix(ncol=3,nrow = length(unique.chain)))
  colnames(query.table)<-c("Chain","Count","Percentage")
  query.table$Chain<-unique.chain
  univ.table<-query.table
  for (i in 1:length(unique.chain)){
  query.table$Count[i]<- sum(X$Chain==query.table$Chain[i])
  univ.table$Count[i]<- sum(Y$Chain==univ.table$Chain[i])
  query.table$Percentage[i]<- round(query.table$Count[i]/nrow(X)*100,2)
  univ.table$Percentage[i]<- round(univ.table$Count[i]/nrow(Y)*100,2)
  }
  query.table$list<-"Query"
  univ.table$list<-"Universe"
  final.table<-rbind(univ.table,query.table)
  ggplot(data=final.table, aes(x=Chain, y=Percentage, fill=factor(list, levels = c("Universe","Query")))) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("#999999","#E69F00"))+
    #geom_text(aes(label=Percentage), vjust=1.6, color="black",position = position_dodge(0.9), size=3.5)+
    #theme_minimal(base_size = 11, base_family = "")+
    theme(legend.title=element_blank(),axis.text.x = element_text(angle = 90),axis.ticks = element_blank())

  }else{

    unique.chain<- unique(X$Chain)
    query.table<-data.frame(matrix(ncol=3,nrow = length(unique.chain)))
    colnames(query.table)<-c("Chain","Count","Percentage")
    query.table$Chain<-unique.chain
    univ.table<-query.table
    for (i in 1:length(unique.chain)){
      query.table$Count[i]<- sum(X$Chain==query.table$Chain[i])
      query.table$Percentage[i]<- round(query.table$Count[i]/nrow(X)*100,2)
    }
    query.table$list<-"Lipids"
    ggplot(data=query.table, aes(x=Chain, y=Percentage)) +
      geom_bar(stat="identity", position=position_dodge()) +
      scale_fill_manual(values="#999999")+
      #geom_text(aes(label=Percentage), vjust=1.6, color="black",position = position_dodge(0.9), size=3.5)+
        #theme_minimal(base_size = 11, base_family = "")+
        theme(legend.title=element_blank(),axis.text.x = element_text(angle = 90),axis.ticks = element_blank())}

    }
