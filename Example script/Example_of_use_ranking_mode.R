################################
# Rodin 0.1.1
# a Shiny app version of the tool is available at
# https://omicstools.pnnl.gov/shiny/lipid-mini-on/
# for help or troubleshouting contact geremy.clair@pnnl.gov

################################ STEP I : loading the data

#load the library
library("Rodin")

#first load the rankingTable
RankingTable <- RankingTableExample

############################### STEP II : generate the enrichment analysis
#clean the lipid list to remove duplicates and empty
cleaned.RT <- clean.rankingTable(RankingTable)

#the next step is to use the lipid miner function to create the objects 1- .intact, 2- .chain and 3- .allchains
lipid.miner(cleaned.RT[,1],name="RT")

# Kolmogorov–Smirnov test based on the ranking and on the ranking values (pvalues, Fold-changes)
intact_KS_cat_result <- intact.KS(RT.intact, colname = "Category", cleaned.RT)
intact_KS_main_result <- intact.KS(RT.intact, colname = "Main class", cleaned.RT)
intact_KS_sub_result <- intact.KS(RT.intact, colname = "Sub class", cleaned.RT)
chain_KS_result <- chain.KS(RT.chain, cleaned.RT)
allchains_KS_result <- allchains.KS(RT.allchains, cleaned.RT)

# Run the category specific total number of chain carbon enrichment
total_carbon_cat_result<-total.carbon.cat.KS(RT.intact,cleaned.RT, enrich=FALSE)
# Run the Main class specific total number of chain carbon enrichment
total_carbon_main_result<-total.carbon.main.KS(RT.intact,cleaned.RT, enrich=FALSE)
# Run the sub class specific total number of chain carbon enrichment
total_carbon_sub_result<-total.carbon.sub.KS(RT.intact,cleaned.RT, enrich=FALSE)

# Run the category specific number of chain DB enrichment
total_DB_cat_result<-total.DB.cat.KS(RT.intact,cleaned.RT, enrich=FALSE)
# Run the Main class specific number of chain DB enrichment
total_DB_main_result<-total.DB.main.KS(RT.intact,cleaned.RT, enrich=FALSE)
# Run the Main class specific number of chain DB enrichment
total_DB_sub_result<-total.DB.sub.KS(RT.intact,cleaned.RT, enrich=FALSE)


# Run the category specific allchains enrichment
allchains_cat_result<-allchains.cat.KS(RT.intact,cleaned.RT, enrich=FALSE)
# Run the Main class specific allchains enrichment
allchains_main_result<-allchains.main.KS(RT.intact,cleaned.RT, enrich=FALSE)
# Run the sub class specific allchains enrichment
allchains_sub_result<-allchains.sub.KS(RT.intact,cleaned.RT, enrich=FALSE)

# subclass enrichment by mainclass (KS)
subclass_mainclass_result<- subclass.mainclass.KS(RT.intact,cleaned.RT, enrich=FALSE)

#reformat the tables to add the type of test
intact_KS_cat_result<- cbind(Type=(rep("Intact")),intact_KS_cat_result)
intact_KS_main_result<- cbind(Type=(rep("Main Class")),intact_KS_main_result)
intact_KS_sub_result<- cbind(Type=(rep("Sub Class")),intact_KS_sub_result)
chain_KS_result<- cbind(Type=(rep("Chain characteristics")),chain_KS_result)
allchains_KS_result<- cbind(Type=(rep("Specific chain")),allchains_KS_result)
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
globaloutput<-rbind(intact_KS_cat_result,intact_KS_main_result,intact_KS_sub_result,chain_KS_result,total_carbon_cat_result,total_carbon_main_result,total_carbon_sub_result,total_DB_cat_result,total_DB_main_result,total_DB_sub_result,allchains_cat_result,allchains_main_result,allchains_sub_result,subclass_mainclass_result, allchains_KS_result)

#filter this output with a pvalue cutoff
filteredGlobalOutput<-globaloutput[globaloutput$`p-value`<0.05,]

#display the global table after filtering it
filteredGlobalOutput

############################# STEP III : generate the graphics#############

#category stack barplot
intact.cat.stack(RT.intact)+ggtitle("Lipids")
#Main class stack barplot
intact.main.stack(RT.intact)+ggtitle("Lipids")

#Sub class stack barplot
intact.sub.stack(RT.intact)+ggtitle("Lipids")

# Chain Length stack barplot
chain.length.stack(Universe.chain)+ggtitle("Lipids")

#Unsaturation stack barplot
chain.unsat.stack(Universe.chain)+ggtitle("Lipids")

#Unique chains barplot (query and uninverse)
allchains.barplot(RT.allchains)+ggtitle("All chains")

#Barplot of the Unique chains with PC as a main class
mainclass.test<-"PC"
allchains.barplot(subsetmainclass(RT.allchains,mainclass = mainclass.test))+ggtitle(paste("All chains of the mainclass", mainclass.test))

#Barplot of the Unique chains with LPC( as a subclass
subclass.test<-"LPC("
allchains.barplot(subsetsubclass(RT.allchains,subclass = subclass.test))+ggtitle(paste("All chains of the mainclass", subclass.test))

#make Network
lipid_visnetwork_plot(cleaned.queryExample,filteredGlobalOutput,pval = 0.05)
