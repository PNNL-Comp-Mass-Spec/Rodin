#' run_the_tests()
#'
#' Runs the specified enrichment test(s) and returns a table of results
#'
#' @param Output.miner list output from lipid.miner() function run on the Query dataset
#' @param Y has to be either a second output.miuner (e.g. Universe.miner) list output from lipid.miner() function run on the Universe dataset OR a RankingTable if the test is KS
#' @param test.type character string specifying the type of enrichment test. Valid values are: "Fisher", "Binom", "Hyper", "EASE", or "KS".
#' @param general.select vector of logical values of length ?, where the vector specifies ? in the following order ?
#' @param subset.by character string specifying one of "category", "mainclass", or "subclass", which denote the ?
#' @param subset.select vector of logical values of length ?, where the vector specifies ? in the following order ?
#' @param enrich logical value specifying whether to ?
#' @param p is a regular pvalue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @example run_the_tests(lipid.miner(cleaned.queryExample,output.list = T), lipid.miner(cleaned.universeExample,output.list = T), test.type="Fisher",general.select=c(T,T,T,T,T),subset.select=c(T,T,T),enrich=F,subset.by = "category")
#' @example run_the_tests(lipid.miner(cleaned.RTexample[,1],output.list = T), cleaned.RTexample, test.type="KS",general.select=c(T,T,T,T,T),subset.select=c(T,T,T),enrich=F,subset.by = "category")
#'
#' @author Geremy Clair, Kelly Stratton
#' @export
#'
run_the_tests <- function(Output.miner, Y, test.type, general.select, subset.by, subset.select, enrich, p, q, order="ascending"){
  if(missing(order)){order="ascending"}

  #verify that Query.miner is of the right type
  if(!is.list(Y) & !names(Y)== c("intact","chain","allchains")){stop("Query.miner is of the wrong type")}


  #test if the type of Y and the type of test are compatible
  if(is.list(Y) & names(Y)== c("intact","chain","allchains")){
    Universe.miner<-Y
    if(test.type=="KS"){stop("Y is a an object generated with lipid.miner, 'KS' cannot be used as test")}
  }else{
    rankingTable<-Y
    if(test.type!="KS"){stop("Y seems to be a rankingTable, 'KS' is the only type of test working with such object")}}

  #Fisher
  if(test.type=="Fisher"){
    #intact cat
    if(general.select[1]==T){
      if (enrich==T){
        intact_cat_result <- intact.fisher.enrich(Output.miner$intact$Category, Universe.miner$intact$Category,p=p,q=q)
      }else{
        intact_cat_result <- intact.fisher(Output.miner$intact$Category, Universe.miner$intact$Category)
      }}
    #intact main
    if(general.select[2]==T){
      if (enrich==T){
        intact_main_result <- intact.fisher.enrich(Output.miner$intact$`Main class`, Universe.miner$intact$`Main class`,p=p,q=q)
      }else{
        intact_main_result <- intact.fisher(Output.miner$intact$`Main class`, Universe.miner$intact$`Main class`)
      }}
    #intact sub
    if(general.select[3]==T){
      if (enrich==T){
        intact_sub_result <- intact.fisher.enrich(Output.miner$intact$`Sub class`, Universe.miner$intact$`Sub class`,p=p,q=q)
      }else{
        intact_sub_result <- intact.fisher(Output.miner$intact$`Sub class`, Universe.miner$intact$`Sub class`)
      }}
    #chain
    if(general.select[4]==T){
      if (enrich==T){
        chain_result <- chain.fisher.enrich(Output.miner$chain, Universe.miner$chain,p=p,q=q)
      }else{
        chain_result <- chain.fisher(Output.miner$chain, Universe.miner$chain)
      }}
    #allchains
    if(general.select[5]==T){
      if (enrich==T){
        allchains_result <- allchains.fisher.enrich(Output.miner$allchains, Universe.miner$allchains,p=p,q=q)
      }else{
        allchains_result <- allchains.fisher(Output.miner$allchains, Universe.miner$allchains)
      }}
    #sub total carbon
    if(subset.select[1]==T){
      if (subset.by=="category"){
        sub_total_carbon_result<- total.carbon.cat(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="mainclass"){
        sub_total_carbon_result<- total.carbon.main(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="subclass"){
        sub_total_carbon_result<- total.carbon.sub(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
    }
    #sub total DB
    if(subset.select[2]==T){
      if (subset.by=="category"){
        sub_total_DB_result<- total.DB.cat(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="mainclass"){
        sub_total_DB_result<- total.DB.main(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="subclass"){
        sub_total_DB_result<- total.DB.sub(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
    }
    #sub allchains
    if(subset.select[2]==T){
      if (subset.by=="category"){
        sub_allchains_result<- allchains.cat(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="mainclass"){
        sub_allchains_result<- allchains.main(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="subclass"){
        sub_allchains_result<- allchains.sub(Output.miner$intact,Universe.miner$intact, enrich=enrich,p=p,q=q)
      }
    }
  }


    #Binom
    if(test.type=="Binom"){
      #intact cat
      if(general.select[1]==T){
        if (enrich==T){
          intact_cat_result <- intact.binom.enrich(Output.miner$intact$Category, Universe.miner$intact$Category,p=p,q=q)
        }else{
          intact_cat_result <- intact.binom(Output.miner$intact$Category, Universe.miner$intact$Category)
        }}
      #intact main
      if(general.select[2]==T){
        if (enrich==T){
          intact_main_result <- intact.binom.enrich(Output.miner$intact$`Main class`, Universe.miner$intact$`Main class`,p=p,q=q)
        }else{
          intact_main_result <- intact.binom(Output.miner$intact$`Main class`, Universe.miner$intact$`Main class`)
        }}
      #intact sub
      if(general.select[3]==T){
        if (enrich==T){
          intact_sub_result <- intact.binom.enrich(Output.miner$intact$`Sub class`, Universe.miner$intact$`Sub class`,p=p,q=q)
        }else{
          intact_sub_result <- intact.binom(Output.miner$intact$`Sub class`, Universe.miner$intact$`Sub class`)
        }}
      #chain
      if(general.select[4]==T){
        if (enrich==T){
          chain_result <- chain.binom.enrich(Output.miner$chain, Universe.miner$chain,p=p,q=q)
        }else{
          chain_result <- chain.binom(Output.miner$chain, Universe.miner$chain)
        }}
      #allchains
      if(general.select[5]==T){
        if (enrich==T){
          allchains_result <- allchains.binom.enrich(Output.miner$allchains, Universe.miner$allchains,p=p,q=q)
        }else{
          allchains_result <- allchains.binom(Output.miner$allchains, Universe.miner$allchains)
        }}
      #sub total carbon
      if(subset.select[1]==T){
        if (subset.by=="category"){
          sub_total_carbon_result<- total.carbon.cat(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="mainclass"){
          sub_total_carbon_result<- total.carbon.main(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="subclass"){
          sub_total_carbon_result<- total.carbon.sub(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
      }
      #sub total DB
      if(subset.select[2]==T){
        if (subset.by=="category"){
          sub_total_DB_result<- total.DB.cat(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="mainclass"){
          sub_total_DB_result<- total.DB.main(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="subclass"){
          sub_total_DB_result<- total.DB.sub(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
      }
      #sub allchains
      if(subset.select[2]==T){
        if (subset.by=="category"){
          sub_allchains_result<- allchains.cat(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="mainclass"){
          sub_allchains_result<- allchains.main(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="subclass"){
          sub_allchains_result<- allchains.sub(Output.miner$intact,Universe.miner$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
      }
    }


      #Hyper
      if(test.type=="Hyper"){
        #intact cat
        if(general.select[1]==T){
          if (enrich==T){
            intact_cat_result <- intact.hyper.enrich(Output.miner$intact$Category, Universe.miner$intact$Category,p=p,q=q)
          }else{
            intact_cat_result <- intact.hyper(Output.miner$intact$Category, Universe.miner$intact$Category)
          }}
        #intact main
        if(general.select[2]==T){
          if (enrich==T){
            intact_main_result <- intact.hyper.enrich(Output.miner$intact$`Main class`, Universe.miner$intact$`Main class`,p=p,q=q)
          }else{
            intact_main_result <- intact.hyper(Output.miner$intact$`Main class`, Universe.miner$intact$`Main class`)
          }}
        #intact sub
        if(general.select[3]==T){
          if (enrich==T){
            intact_sub_result <- intact.hyper.enrich(Output.miner$intact$`Sub class`, Universe.miner$intact$`Sub class`,p=p,q=q)
          }else{
            intact_sub_result <- intact.hyper(Output.miner$intact$`Sub class`, Universe.miner$intact$`Sub class`)
          }}
        #chain
        if(general.select[4]==T){
          if (enrich==T){
            chain_result <- chain.hyper.enrich(Output.miner$chain, Universe.miner$chain,p=p,q=q)
          }else{
            chain_result <- chain.hyper(Output.miner$chain, Universe.miner$chain)
          }}
        #allchains
        if(general.select[5]==T){
          if (enrich==T){
            allchains_result <- allchains.hyper.enrich(Output.miner$allchains, Universe.miner$allchains,p=p,q=q)
          }else{
            allchains_result <- allchains.hyper(Output.miner$allchains, Universe.miner$allchains)
          }}
        #sub total carbon
        if(subset.select[1]==T){
          if (subset.by=="category"){
            sub_total_carbon_result<- total.carbon.cat(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="mainclass"){
            sub_total_carbon_result<- total.carbon.main(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="subclass"){
            sub_total_carbon_result<- total.carbon.sub(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
        }
        #sub total DB
        if(subset.select[2]==T){
          if (subset.by=="category"){
            sub_total_DB_result<- total.DB.cat(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="mainclass"){
            sub_total_DB_result<- total.DB.main(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="subclass"){
            sub_total_DB_result<- total.DB.sub(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
        }
        #sub allchains
        if(subset.select[2]==T){
          if (subset.by=="category"){
            sub_allchains_result<- allchains.cat(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="mainclass"){
            sub_allchains_result<- allchains.main(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="subclass"){
            sub_allchains_result<- allchains.sub(Output.miner$intact,Universe.miner$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
        }
      }


        #EASE
        if(test.type=="EASE"){
          #intact cat
          if(general.select[1]==T){
            if (enrich==T){
              intact_cat_result <- intact.EASE.enrich(Output.miner$intact$Category, Universe.miner$intact$Category,p=p,q=q)
            }else{
              intact_cat_result <- intact.EASE(Output.miner$intact$Category, Universe.miner$intact$Category)
            }}
          #intact main
          if(general.select[2]==T){
            if (enrich==T){
              intact_main_result <- intact.EASE.enrich(Output.miner$intact$`Main class`, Universe.miner$intact$`Main class`,p=p,q=q)
            }else{
              intact_main_result <- intact.EASE(Output.miner$intact$`Main class`, Universe.miner$intact$`Main class`)
            }}
          #intact sub
          if(general.select[3]==T){
            if (enrich==T){
              intact_sub_result <- intact.EASE.enrich(Output.miner$intact$`Sub class`, Universe.miner$intact$`Sub class`,p=p,q=q)
            }else{
              intact_sub_result <- intact.EASE(Output.miner$intact$`Sub class`, Universe.miner$intact$`Sub class`)
            }}
          #chain
          if(general.select[4]==T){
            if (enrich==T){
              chain_result <- chain.EASE.enrich(Output.miner$chain, Universe.miner$chain,p=p,q=q)
            }else{
              chain_result <- chain.EASE(Output.miner$chain, Universe.miner$chain)
            }}
          #allchains
          if(general.select[5]==T){
            if (enrich==T){
              allchains_result <- allchains.EASE.enrich(Output.miner$allchains, Universe.miner$allchains,p=p,q=q)
            }else{
              allchains_result <- allchains.EASE(Output.miner$allchains, Universe.miner$allchains)
            }}
          #sub total carbon
          if(subset.select[1]==T){
            if (subset.by=="category"){
              sub_total_carbon_result<- total.carbon.cat(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="mainclass"){
              sub_total_carbon_result<- total.carbon.main(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="subclass"){
              sub_total_carbon_result<- total.carbon.sub(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
          }
          #sub total DB
          if(subset.select[2]==T){
            if (subset.by=="category"){
              sub_total_DB_result<- total.DB.cat(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="mainclass"){
              sub_total_DB_result<- total.DB.main(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="subclass"){
              sub_total_DB_result<- total.DB.sub(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
          }
          #sub allchains
          if(subset.select[2]==T){
            if (subset.by=="category"){
              sub_allchains_result<- allchains.cat(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="mainclass"){
              sub_allchains_result<- allchains.main(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="subclass"){
              sub_allchains_result<- allchains.sub(Output.miner$intact,Universe.miner$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
          }
        }

  #KS
  if(test.type=="KS"){
    #intact cat
    if(general.select[1]==T){
      if (enrich==T){
        intact_cat_result <- intact.KS.enrich(Output.miner$intact, colname= "Category", rankingTable , p=p , q=q , order=order)
      }else{
        intact_cat_result <- intact.KS(Output.miner$intact, colname= "Category", rankingTable , order=order)
      }}
    #intact main
    if(general.select[2]==T){
      if (enrich==T){
        intact_main_result <- intact.KS.enrich(Output.miner$intact, colname= "Main class", rankingTable , p=p , q=q , order=order)
      }else{
        intact_main_result <- intact.KS(Output.miner$intact, colname= "Main class", rankingTable , order=order)
      }}
    #intact sub
    if(general.select[3]==T){
      if (enrich==T){
        intact_sub_result <- intact.KS.enrich(Output.miner$intact, colname= "Sub class", rankingTable , p=p , q=q , order=order)
      }else{
        intact_sub_result <- intact.KS(Output.miner$intact, colname= "Sub class", rankingTable , order=order)
      }}
    #chain
    if(general.select[4]==T){
      if (enrich==T){
        chain_result <- chain.KS.enrich(Output.miner$chain, rankingTable , p=p , q=q , order=order)
      }else{
        chain_result <- chain.KS(Output.miner$chain, rankingTable, order=order)
      }}
    #allchains
    if(general.select[5]==T){
      if (enrich==T){
        allchains_result <- allchains.KS.enrich(Output.miner$allchains, rankingTable , p=p , q=q , order=order)
      }else{
        allchains_result <- allchains.KS(Output.miner$allchains, rankingTable , order=order)
      }}

    #sub total carbon
    if(subset.select[1]==T){
      if (subset.by=="category"){
        sub_total_carbon_result<- total.carbon.cat.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="mainclass"){
        sub_total_carbon_result<- total.carbon.main.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="subclass"){
        sub_total_carbon_result<- total.carbon.sub.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
    }
    #sub total DB
    if(subset.select[2]==T){
      if (subset.by=="category"){
        sub_total_DB_result<- total.DB.cat.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="mainclass"){
        sub_total_DB_result<- total.DB.main.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="subclass"){
        sub_total_DB_result<- total.DB.cat.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
    }
    #sub allchains
    if(subset.select[2]==T){
      if (subset.by=="category"){
        sub_allchains_result<- allchains.cat.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="mainclass"){
        sub_allchains_result<- allchains.main.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="subclass"){
        sub_allchains_result<- allchains.sub.KS(Output.miner$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
    }
  }

  #Output tables if the corresponding object is not NA
  #section titles object
  title.enrich<-NA

  #create a global output object (for a global download of the results)
  global.output<-data.frame(matrix(nrow=0,ncol=9))
  colnames(global.output)<-c("Test.performed","Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")

  options(warn=-1)
  if(intact_cat_result!=F){
    title.enrich<- paste("Category(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    intact_cat_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,intact_cat_result))
  }

  if(intact_main_result!=F){
    title.enrich<- paste("Main class(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    intact_main_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,intact_main_result))
  }

  if(intact_sub_result!=F){
    title.enrich<- paste("Sub class(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    intact_sub_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,intact_sub_result))
  }

  if(chain_result!=F){
    title.enrich<- paste("Chain(s) characteristics(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    chain_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,chain_result))
  }

  if(allchains_result!=F){
    title.enrich<- paste("Specific chain(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    allchains_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,allchains_result))
  }

  if(sub_total_carbon_result!=F){
    title.enrich<- paste("Total chain carbon by ",subset.by,"(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    sub_total_carbon_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich, sub_total_carbon_result))
  }

  if(sub_total_DB_result!=F){
    title.enrich<- paste("Total number of DB by ",subset.by,"(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    sub_total_DB_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich, sub_total_DB_result))
  }

  if(sub_allchains_result!=F){
    title.enrich<- paste("Specific chains by ",subset.by,"(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    sub_allchains_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich, sub_allchains_result))
  }

  # download the global output
  return(global.output)

}
