################################
# Rodin 0.1.0
# a Shiny app version of the tool is available at
# https://omicstools.pnnl.gov/shiny/lipid-mini-on/
# for help or troubleshouting contact geremy.clair@pnnl.gov

################################ STEP I : loading the data

#load the library
library("Rodin")

#first load the query data and the universe data
Query <- QueryExample
Universe <- UniverseExample

############################### STEP II : generate the enrichment analysis

cleaned.queryExample <- clean.lipid.list(Query)
cleaned.universeExample <- clean.lipid.list(Universe)

#the next step is to use the lipid miner function to create the objects 1- .intact, 2- .chain and 3- .allchains
lipid.miner(cleaned.queryExample, name="Query", TGcollapse.rm = TRUE)
lipid.miner(cleaned.universeExample, name="Universe", TGcollapse.rm = TRUE)

# run the Fisher (and create results objects)
intact_fisher_cat_result <- intact.fisher(Query.intact$Category, Universe.intact$Category)
intact_fisher_main_result <- intact.fisher(Query.intact$`Main class`, Universe.intact$`Main class`)
intact_fisher_sub_result <- intact.fisher(Query.intact$`Sub class`, Universe.intact$`Sub class`)
chain_fisher_result <- chain.fisher(Query.chain, Universe.chain)
allchains_fisher_result <- allchains.fisher(Query.allchains, Universe.allchains)

# Run the category specific total number of chain carbon enrichment
total_carbon_cat_result<-total.carbon.cat(Query.intact,Universe.intact, enrich=FALSE)
# Run the Main class specific total number of chain carbon enrichment
total_carbon_main_result<-total.carbon.main(Query.intact,Universe.intact, enrich=FALSE)
# Run the sub class specific total number of chain carbon enrichment
total_carbon_sub_result<-total.carbon.sub(Query.intact,Universe.intact, enrich=FALSE)

# Run the category specific number of chain DB enrichment
total_DB_cat_result<-total.DB.cat(Query.intact,Universe.intact, enrich=FALSE)
# Run the Main class specific number of chain DB enrichment
total_DB_main_result<-total.DB.main(Query.intact,Universe.intact, enrich=FALSE)
# Run the Main class specific number of chain DB enrichment
total_DB_sub_result<-total.DB.sub(Query.intact,Universe.intact, enrich=FALSE)

# Run the category specific allchains enrichment
allchains_cat_result<-allchains.cat(Query.intact,Universe.intact, enrich=FALSE)
# Run the Main class specific allchains enrichment
allchains_main_result<-allchains.main(Query.intact,Universe.intact, enrich=FALSE)
# Run the sub class specific allchains enrichment
allchains_sub_result<-allchains.sub(Query.intact,Universe.intact, enrich=FALSE)

# subclass enrichment by mainclass (fisher)
subclass_mainclass_result<- subclass.mainclass(Query.intact,Universe.intact,test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0)

#reformat the tables to add the type of test
intact_fisher_cat_result<- cbind(Type=(rep("Intact")),intact_fisher_cat_result)
intact_fisher_main_result<- cbind(Type=(rep("Main Class")),intact_fisher_main_result)
intact_fisher_sub_result<- cbind(Type=(rep("Sub Class")),intact_fisher_sub_result)
chain_fisher_result<- cbind(Type=(rep("Chain characteristics")),chain_fisher_result)
allchains_fisher_result<- cbind(Type=(rep("Specific chain")),allchains_fisher_result)
total_carbon_cat_result<- cbind(Type=(rep("TotalCarbon by Cat")),total_carbon_cat_result)
total_carbon_main_result<- cbind(Type=(rep("TotalCarbon by Main")),total_carbon_main_result)
total_carbon_sub_result<- cbind(Type=(rep("TotalCarbon by Sub")),total_carbon_sub_result)
total_DB_cat_result<- cbind(Type=(rep("total DB by Cat")),total_DB_cat_result)
total_DB_main_result<- cbind(Type=(rep("total DB by Main")),total_DB_main_result)
total_DB_sub_result<- cbind(Type=(rep("total DB by Sub")),total_DB_sub_result)
allchains_cat_result<- cbind(Type=(rep("Specific Chain by Cat")),allchains_cat_result)
allchains_main_result<- cbind(Type=(rep("Specific Chain by Main")),allchains_main_result)
allchains_sub_result<- cbind(Type=(rep("Specific Chain by Sub")),allchains_sub_result)
subclass_mainclass_result<-cbind(Type=(rep("Subclass by Mainclass")),subclass_mainclass_result)

# place all the tables in a globaloutput table
globaloutput<-rbind(intact_fisher_cat_result,intact_fisher_main_result,intact_fisher_sub_result,chain_fisher_result,total_carbon_cat_result,total_carbon_main_result,total_carbon_sub_result,total_DB_cat_result,total_DB_main_result,total_DB_sub_result,allchains_cat_result,allchains_main_result,allchains_sub_result,subclass_mainclass_result, allchains_fisher_result)

#filter this output with a pvalue cutoff
filteredGlobalOutput<-globaloutput[globaloutput$Pvalue<0.05&globaloutput$fold.change>1,]

#display the global table after filtering it
filteredGlobalOutput

############################# STEP III : generate the graphics#############

#category stack barplot
grid.arrange(intact.cat.stack(Universe.intact)+ggtitle("Universe"),intact.cat.stack(Query.intact)+ggtitle("Query"),ncol=2)

#Main class stack barplot
grid.arrange(intact.main.stack(Universe.intact)+ggtitle("Universe"),intact.main.stack(Query.intact)+ggtitle("Query"),ncol=2)

#Sub class stack barplot
grid.arrange(intact.sub.stack(Universe.intact)+ggtitle("Universe"),intact.sub.stack(Query.intact)+ggtitle("Query"),ncol=2)

# Chain Length stack barplot
grid.arrange(chain.length.stack(Universe.chain)+ggtitle("Universe"),chain.length.stack(Query.chain)+ggtitle("Query"),ncol=2)

#Unsaturation stack barplot
grid.arrange(chain.unsat.stack(Universe.chain)+ggtitle("Universe"),chain.unsat.stack(Query.chain)+ggtitle("Query"),ncol=2)

#Unique chains barplot (query and uninverse)
allchains.barplot(Query.allchains,Universe.allchains)+ggtitle("All chains")


#Barplot of the Unique chains with PC as a main class
mainclass.test<-"PC"
allchains.barplot(subsetmainclass(Query.allchains,mainclass = mainclass.test), subsetmainclass(Universe.allchains,mainclass = mainclass.test))+ggtitle(paste("All chains of the mainclass", mainclass.test))

#Barplot of the Unique chains with LPC( as a subclass
subclass.test<-"LPC("
allchains.barplot(subsetsubclass(Query.allchains,subclass = subclass.test), subsetsubclass(Universe.allchains,subclass = subclass.test))+ggtitle(paste("All chains of the mainclass", subclass.test))

#make Network
lipid_visnetwork_plot(cleaned.queryExample,filteredGlobalOutput,pval = 0.05)
