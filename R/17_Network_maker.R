#' lipid_edge_maker()
#'
#' Create a list of lipid classification edges from a cleaned lipid list
#'
#' @param X Any lipid character object, preferably cleaned using the clean.lipid.list() function
#'
#' @details Create a list of lipid classification edges from a cleaned lipid list
#'
#' @return This function will return a data frame with a list of edges
#'
#' @examples Query.edge<-lipid_edge_maker(cleaned.queryExample)
#'
#' @author Geremy Clair
#' @export
#'

lipid_edge_maker<-function(X){
  #create the required objects using lipid miner
  X.miner<- lipid.miner(X,name= "X.miner",output.list = T)
  X.intact<-data.frame(X.miner$intact)
  X.chain<- data.frame(X.miner$chain)


  X.results<- data.frame(matrix(ncol=3, nrow=0))
  colnames(X.results) <- c("Lipid name","Classifier","Class")
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Category","Class"=X.intact$Category))
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Main class","Class"=X.intact$Main.class))
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Sub class","Class"=X.intact$Sub.class))
X.allchains.results<- data.frame(matrix(ncol=6, nrow=0))
  colnames(X.allchains.results) <- c("Lipid name","Classifier","Class","Category","Main class","Sub class")
  if (sum(!is.na(X.intact$Chain.1))!=0){
  X.allchains.results<- rbind(X.allchains.results, data.frame("Lipid name"=X.intact$Lipid[!is.na(X.intact$Chain.1)],"Classifier"="Specific chain","Class"=X.intact$Chain.1[!is.na(X.intact$Chain.1)],"Category"=X.intact$Category[!is.na(X.intact$Chain.1)],"Main class"=X.intact$Main.class[!is.na(X.intact$Chain.1)],"Sub class"=X.intact$Sub.class[!is.na(X.intact$Chain.1)]))
  }
  if (sum(!is.na(X.intact$Chain.2))!=0){
    X.allchains.results<- rbind(X.allchains.results, data.frame("Lipid name"=X.intact$Lipid[!is.na(X.intact$Chain.2)],"Classifier"="Specific chain","Class"=X.intact$Chain.2[!is.na(X.intact$Chain.2)],"Category"=X.intact$Category[!is.na(X.intact$Chain.2)],"Main class"=X.intact$Main.class[!is.na(X.intact$Chain.2)],"Sub class"=X.intact$Sub.class[!is.na(X.intact$Chain.2)]))
  }
  if (sum(!is.na(X.intact$Chain.3))!=0){
    X.allchains.results<- rbind(X.allchains.results, data.frame("Lipid name"=X.intact$Lipid[!is.na(X.intact$Chain.3)],"Classifier"="Specific chain","Class"=X.intact$Chain.3[!is.na(X.intact$Chain.3)],"Category"=X.intact$Category[!is.na(X.intact$Chain.3)],"Main class"=X.intact$Main.class[!is.na(X.intact$Chain.3)],"Sub class"=X.intact$Sub.class[!is.na(X.intact$Chain.3)]))
  }
  if (sum(!is.na(X.intact$Chain.4))!=0){
    X.allchains.results<- rbind(X.allchains.results, data.frame("Lipid name"=X.intact$Lipid[!is.na(X.intact$Chain.4)],"Classifier"="Specific chain","Class"=X.intact$Chain.4[!is.na(X.intact$Chain.4)],"Category"=X.intact$Category[!is.na(X.intact$Chain.4)],"Main class"=X.intact$Main.class[!is.na(X.intact$Chain.4)],"Sub class"=X.intact$Sub.class[!is.na(X.intact$Chain.4)]))
  }
  #remove duplicates
  X.allchains.results<-X.allchains.results[!duplicated(X.allchains.results),]
  #add the "specific chains"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.allchains.results$Lipid.name,"Classifier"="Specific chain","Class"=X.allchains.results$Class))
  #add the "specific chains by category"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.allchains.results$Lipid.name,"Classifier"="Specific chains by category","Class"=paste0(X.allchains.results$Category," with the chain ",X.allchains.results$Class)))
  #add the "Specific chains by mainclass"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.allchains.results$Lipid.name,"Classifier"="Specific chains by mainclass","Class"=paste0(X.allchains.results$Main.class," with the chain ",X.allchains.results$Class)))
  #add the "Specific chains by subclass"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.allchains.results$Lipid.name,"Classifier"="Specific chains by subclass","Class"=paste0(X.allchains.results$Sub.class," with the chain ",X.allchains.results$Class)))

  #"Total chain carbon by category"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Total chain carbon by category","Class"=paste0(X.intact$Category," with a total number of chain carbon of ",X.intact$Total.Number.of.Carbon)))

  #"Total chain carbon by mainclass"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Total chain carbon by mainclass","Class"=paste0(X.intact$Main.class," with a total number of chain carbon of ",X.intact$Total.Number.of.Carbon)))

  #"Total chain carbon by subclass"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Total chain carbon by subclass","Class"=paste0(X.intact$Sub.class," with a total number of chain carbon of ",X.intact$Total.Number.of.Carbon)))

  #"Total number of DB by category"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Total number of DB by category","Class"=paste0(X.intact$Category," with a total number of unsaturation of ",X.intact$Double.Bonds)))

  #"Total number of DB by mainclass"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Total number of DB by mainclass","Class"=paste0(X.intact$Main.class," with a total number of unsaturation of ",X.intact$Double.Bonds)))

  #"Total number of DB by subclass"
  X.results<- rbind(X.results,data.frame("Lipid name"=X.intact$Lipid,"Classifier"="Total number of DB by subclass","Class"=paste0(X.intact$Sub.class," with a total number of unsaturation of ",X.intact$Double.Bonds)))

  #specific chains
  X.chain.results<- data.frame(matrix(ncol=3, nrow=0))
  colnames(X.chain.results) <- c("Lipid name","Classifier","Class")
  if (sum(X.chain$SCFA>0)!=0){
  X.chain.results<- rbind(X.chain.results,data.frame("Lipid name"= as.character(X.chain$Lipid[X.chain$SCFA>0]),"Classifier"="Chain(s) characteristics", "Class"="SCFA"))
  }
  if (sum(X.chain$MCFA>0)!=0){
  X.chain.results<- rbind(X.chain.results,data.frame("Lipid name"= as.character(X.chain$Lipid[X.chain$MCFA>0]),"Classifier"="Chain(s) characteristics", "Class"="MCFA"))
  }
  if (sum(X.chain$LCFA>0)!=0){
  X.chain.results<- rbind(X.chain.results,data.frame("Lipid name"= as.character(X.chain$Lipid[X.chain$LCFA>0]),"Classifier"="Chain(s) characteristics", "Class"="LCFA"))
  }
  if (sum(X.chain$VLCFA>0)!=0){
    X.chain.results<- rbind(X.chain.results,data.frame("Lipid name"= as.character(X.chain$Lipid[X.chain$VLCFA>0]),"Classifier"="Chain(s) characteristics", "Class"="VLCFA"))
  }
  if (sum(X.chain$Saturated>0)!=0){
    X.chain.results<- rbind(X.chain.results,data.frame("Lipid name"= as.character(X.chain$Lipid[X.chain$Saturated>0]),"Classifier"="Chain(s) characteristics", "Class"="Saturated"))
  }
  if (sum(X.chain$Monounsaturated>0)!=0){
    X.chain.results<- rbind(X.chain.results,data.frame("Lipid name"= as.character(X.chain$Lipid[X.chain$Monounsaturated>0]),"Classifier"="Chain(s) characteristics", "Class"="Monounsaturated"))
  }
  if (sum(X.chain$Diunsaturated>0)!=0){
    X.chain.results<- rbind(X.chain.results,data.frame("Lipid name"= as.character(X.chain$Lipid[X.chain$Diunsaturated>0]),"Classifier"="Chain(s) characteristics", "Class"="Diunsaturated"))
  }
  if (sum(X.chain$Polyunsaturated>0)!=0){
    X.chain.results<- rbind(X.chain.results,data.frame("Lipid name"= as.character(X.chain$Lipid[X.chain$Polyunsaturated>0]),"Classifier"="Chain(s) characteristics", "Class"="Polyunsaturated"))
  }

  X.results<-rbind(X.results,X.chain.results)
  X.results

}

#' lipid_network_maker()
#'
#' Enable the construction of enrichment networks
#'
#' @param X either a edge data.frame created with the function edge_maker() or Any lipid character object, preferably cleaned using the clean.lipid.list() function
#' @param Y any enrichment table created with Lipid Mini-On
#' @param pval unadjusted pvalue cutoff
#' @param adjpval BH adjusted pvalue cutoff
#' @param name This will be the name of the output of the function ("network by default")
#'
#' @details Create multiple enrichment network objects based on a query object and a Lipid-Mini-On enrichment analysis table result
#'
#' @return This function will return 3 object "name".edge (list of edges), "name".nodes_attributes(list of nodes, type of node, color), "name".edge_attribute (list of edges, colors)
#'
#' @examples lipid_network_maker(cleaned.queryExample,run_the_tests(lipid.miner(cleaned.queryExample,output.list = T), lipid.miner(cleaned.universeExample,output.list = T), test.type="EASE",general.select=c(T,T,T,T,T),subset.select=c(T,T,T),enrich=F,subset.by = "category"),adjpval=0.05)
#'
#' @author Geremy Clair
#' @export
#'

lipid_network_maker<-function(X,Y,pval=1,adjpval=1,name="network"){
   if(!exists("pval")){pval<-1}
  if(!exists("adjpval")){adjpval<-1}
  if(!exists("name")){name<-"network"}
  if(is.character(X)){X<-lipid_edge_maker(X)}

  #reduce the edge table to the two sides of the edge
  X.filtered<-X[,c(1,3)]
  #filter the result table of a test based on pvalue
  table.filtered<-Y[Y$fold.change>=1,]
  table.filtered<-table.filtered[table.filtered$Pvalue<=pval,]
  table.filtered<-table.filtered[table.filtered$Pvalue<=adjpval,]

  #perform the filtering
  X.filtered<- X.filtered[X.filtered$Class%in% table.filtered$Classifier,]
  X.filtered<- X.filtered[order(X.filtered$Class),]

  X.filtered$Color[1]<-1
  for (i in 2:nrow(X.filtered)){
    if(X.filtered$Class[i]==X.filtered$Class[i-1]){X.filtered$Color[i]<-X.filtered$Color[i-1]}
    else
    {X.filtered$Color[i]<-X.filtered$Color[i-1]+1}
  }
  pal<- colorRampPalette(c("#333399","#ff6633","#006666","#cc3366","#33cc33","#663399","#ffcc00","#009999","#ff3333"))

  X.filtered$Color<- pal(max(X.filtered$Color))[X.filtered$Color]

  class.attribute<-unique(X.filtered[,2:3])
  colnames(class.attribute)<-c("Node","Color")

  nodes.attributes<- rbind(unique(data.frame("Node"=X.filtered[,1],"Classifier"="Lipid")),unique(data.frame("Node"=X.filtered[,2],"Classifier"="Classifier")))
  nodes.attributes$Color<- class.attribute[match(nodes.attributes$Node,class.attribute$Node),2]

  nodes.attributes$Color[is.na(nodes.attributes$Color)]<-"#cccccc"

  assign(paste(name,".edges",sep=""),X.filtered[,1:2], envir = globalenv())
  assign(paste(name,".edges_attributes",sep=""),X.filtered, envir = globalenv())
  assign(paste(name,".nodes_attributes",sep=""),nodes.attributes, envir = globalenv())

}

#' lipid_visnetwork_plot()
#'
#' plots an interactive network of the lipids and enrichment terms
#'
#' @param X either a edge data.frame created with the function edge_maker() or Any lipid character object, preferably cleaned using the clean.lipid.list() function
#' @param Y any enrichment table created with Lipid Mini-On
#' @param pval unadjusted pvalue cutoff
#' @param adjpval BH adjusted pvalue cutoff
#'
#' @details plots an interactive network of the lipids and enrichment terms
#'
#' @return This function will plot an interactive network
#'
#' @examples lipid_visnetwork_plot(cleaned.queryExample,run_the_tests(lipid.miner(cleaned.queryExample,output.list = T), lipid.miner(cleaned.universeExample,output.list = T), test.type="EASE",general.select=c(T,T,T,T,T),subset.select=c(T,T,T),enrich=F,subset.by = "category"),adjpval=0.05)
#'
#' @author Geremy Clair
#' @export
#'

lipid_visnetwork_plot<-function(X,Y,pval=1,adjpval=1){
  if(!exists("pval")){pval<-1}
  if(!exists("adjpval")){adjpval<-1}
  if(is.character(X)){X<-lipid_edge_maker(X)}
  p<-pval
  ap<- adjpval
  lipid_network_maker(X,Y,pval=p,adjpval=ap,name="network")
  colnames(network.nodes_attributes)<-c("label","title","color.background")
  network.nodes_attributes2<-cbind(id=paste0("s",1:nrow(network.nodes_attributes)),network.nodes_attributes,shape=c("dot", "diamond")[as.numeric(network.nodes_attributes$title)],size=c(5, 30)[as.numeric(network.nodes_attributes$title)],borderWidth=0)
  network.nodes_attributes2$title<- network.nodes_attributes2$label
  network.edges_attributes2<-data.frame(from=network.nodes_attributes2$id[match(network.edges_attributes$Lipid.name,network.nodes_attributes2$label)],to=network.nodes_attributes2$id[match(network.edges_attributes$Class,network.nodes_attributes2$label)],color=network.edges_attributes$Color,width=1)
  visNetwork(network.nodes_attributes2, network.edges_attributes2, width="100%", height="1200px") %>% visOptions(highlightNearest = TRUE, selectedBy = "type.label",manipulation=T) %>% visEvents(stabilized = "function() { this.setOptions({nodes : {physics : false}})}")%>% visNodes(font = '18px arial #343434')
  }



