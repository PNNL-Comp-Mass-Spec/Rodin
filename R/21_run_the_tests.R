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
run_the_tests <- function(X, Y, test.type, general.select, subset.by, subset.select, enrich, p, q, order="ascending"){
  if(missing(order)){order="ascending"}

  #verify that Query.miner is of the right type
  if(test.type!="KS"){ 
    if(!is.list(Y) & !names(Y)== c("intact","chain","allchains")){stop("Query.miner is of the wrong type")}
    if(is.list(Y) & names(Y)== c("intact","chain","allchains")){
      Universe.miner<-Y
    } else{stop("Y seems to be a rankingTable, 'KS' is the only type of test working with such object")}
    
  } else if(test.type=="KS"){
    if(!is.null(names(Y))) {
        if(all(names(Y) %in% c("intact","chain","allchains"))){
      stop("Y is a an object generated with lipid.miner, 'KS' cannot be used as test")
        } else {
          rankingTable<-Y 
        }
      }else{
        stop("Improper ranking table")
      }
}
  #Fisher
  if(test.type=="Fisher"){
    #intact cat
    if("cat" %in% general.select){
      if (enrich==T){
        intact_cat_result <- intact.fisher.enrich(X$intact$Category, Y$intact$Category,p=p,q=q)
      }else{
        intact_cat_result <- intact.fisher(X$intact$Category, Y$intact$Category)
      }}
    #intact main
    if("main" %in% general.select){
      if (enrich==T){
        intact_main_result <- intact.fisher.enrich(X$intact$`Main class`, Y$intact$`Main class`,p=p,q=q)
      }else{
        intact_main_result <- intact.fisher(X$intact$`Main class`, Y$intact$`Main class`)
      }}
    #intact sub
    if("sub" %in% general.select){
      if (enrich==T){
        intact_sub_result <- intact.fisher.enrich(X$intact$`Sub class`, Y$intact$`Sub class`,p=p,q=q)
      }else{
        intact_sub_result <- intact.fisher(X$intact$`Sub class`, Y$intact$`Sub class`)
      }}
    #chain
    if("chains" %in% general.select){
      if (enrich==T){
        chain_result <- chain.fisher.enrich(X$chain, Y$chain,p=p,q=q)
      }else{
        chain_result <- chain.fisher(X$chain, Y$chain)
      }}
    #allchains
    if("length" %in% general.select){
      if (enrich==T){
        allchains_result <- allchains.fisher.enrich(X$allchains, Y$allchains,p=p,q=q)
      }else{
        allchains_result <- allchains.fisher(X$allchains, Y$allchains)
      }}
    #sub total carbon
    if("total_carbon" %in% subset.select){
      if (subset.by=="category"){
        sub_total_carbon_result<- total.carbon.cat(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="mainclass"){
        sub_total_carbon_result<- total.carbon.main(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="subclass"){
        sub_total_carbon_result<- total.carbon.sub(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="all"){
        sub_total_carbon_result<- total.carbon.cat(X$intact,Y$intact, enrich=enrich,p=p,q=q)
        sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.main(X$intact,Y$intact, enrich=enrich,p=p,q=q))
        sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.sub(X$intact,Y$intact, enrich=enrich,p=p,q=q))
      }
    }
    #sub total DB
    if("total_insaturation" %in% subset.select){
      if (subset.by=="category"){
        sub_total_DB_result<- total.DB.cat(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="mainclass"){
        sub_total_DB_result<- total.DB.main(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="subclass"){
        sub_total_DB_result<- total.DB.sub(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="all"){
        sub_total_DB_result<- total.DB.cat(X$intact,Y$intact, enrich=enrich,p=p,q=q)
        sub_total_DB_result<- rbind(sub_total_DB_result,total.DB.main(X$intact,Y$intact, enrich=enrich,p=p,q=q))
        sub_total_DB_result<- rbind(sub_total_DB_result,total.DB.sub(X$intact,Y$intact, enrich=enrich,p=p,q=q))
      }
    }
    #sub allchains
    if("specific_chains" %in% subset.select){
      if (subset.by=="category"){
        sub_allchains_result<- allchains.cat(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="mainclass"){
        sub_allchains_result<- allchains.main(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="subclass"){
        sub_allchains_result<- allchains.sub(X$intact,Y$intact, enrich=enrich,p=p,q=q)
      }
      if (subset.by=="all"){
        sub_allchains_result<- allchains.cat(X$intact,Y$intact, enrich=enrich,p=p,q=q)
        sub_allchains_result<- rbind(sub_allchains_result,allchains.main(X$intact,Y$intact, enrich=enrich,p=p,q=q))
        sub_allchains_result<- rbind(sub_allchains_result,allchains.sub(X$intact,Y$intact, enrich=enrich,p=p,q=q))
      }
    }
  }


    #Binom
    if(test.type=="Binom"){
      #intact cat
      if("cat" %in% general.select){
        if (enrich==T){
          intact_cat_result <- intact.binom.enrich(X$intact$Category, Y$intact$Category,p=p,q=q)
        }else{
          intact_cat_result <- intact.binom(X$intact$Category, Y$intact$Category)
        }}
      #intact main
      if("main" %in% general.select){
        if (enrich==T){
          intact_main_result <- intact.binom.enrich(X$intact$`Main class`, Y$intact$`Main class`,p=p,q=q)
        }else{
          intact_main_result <- intact.binom(X$intact$`Main class`, Y$intact$`Main class`)
        }}
      #intact sub
      if("sub" %in% general.select){
        if (enrich==T){
          intact_sub_result <- intact.binom.enrich(X$intact$`Sub class`, Y$intact$`Sub class`,p=p,q=q)
        }else{
          intact_sub_result <- intact.binom(X$intact$`Sub class`, Y$intact$`Sub class`)
        }}
      #chain
      if("chains" %in% general.select){
        if (enrich==T){
          chain_result <- chain.binom.enrich(X$chain, Y$chain,p=p,q=q)
        }else{
          chain_result <- chain.binom(X$chain, Y$chain)
        }}
      #allchains
      if("length" %in% general.select){
        if (enrich==T){
          allchains_result <- allchains.binom.enrich(X$allchains, Y$allchains,p=p,q=q)
        }else{
          allchains_result <- allchains.binom(X$allchains, Y$allchains)
        }}
      #sub total carbon
      if("total_carbon" %in% subset.select){
        if (subset.by=="category"){
          sub_total_carbon_result<- total.carbon.cat(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="mainclass"){
          sub_total_carbon_result<- total.carbon.main(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="subclass"){
          sub_total_carbon_result<- total.carbon.sub(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="all"){
          sub_total_carbon_result<- total.carbon.cat(X$intact,Y$intact, test= "Binom", enrich=enrich,p=p,q=q)
          sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.main(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q))
          sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.sub(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q))
        }
        
      }
      #sub total DB
      if("total_insaturation" %in% subset.select){
        if (subset.by=="category"){
          sub_total_DB_result<- total.DB.cat(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="mainclass"){
          sub_total_DB_result<- total.DB.main(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="subclass"){
          sub_total_DB_result<- total.DB.sub(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="all"){
          sub_total_DB_result<- total.DB.cat(X$intact,Y$intact, test= "Binom", enrich=enrich,p=p,q=q)
          sub_total_DB_result<- rbind(sub_total_DB_result,total.DB.main(X$intact,Y$intact, test= "Binom", enrich=enrich,p=p,q=q))
          sub_total_DB_result<- rbind(sub_total_DB_result,total.DB.sub(X$intact,Y$intact, test= "Binom", enrich=enrich,p=p,q=q))
        }
      }
      #sub allchains
      if("specific_chains" %in% subset.select){
        if (subset.by=="category"){
          sub_allchains_result<- allchains.cat(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="mainclass"){
          sub_allchains_result<- allchains.main(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="subclass"){
          sub_allchains_result<- allchains.sub(X$intact,Y$intact, test= "Binom",enrich=enrich,p=p,q=q)
        }
        if (subset.by=="all"){
          sub_allchains_result<- allchains.cat(X$intact,Y$intact,  test= "Binom",enrich=enrich,p=p,q=q)
          sub_allchains_result<- rbind(sub_allchains_result,allchains.main(X$intact,Y$intact,  test= "Binom",enrich=enrich,p=p,q=q))
          sub_allchains_result<- rbind(sub_allchains_result,allchains.sub(X$intact,Y$intact,  test= "Binom",enrich=enrich,p=p,q=q))
        }
      }
    }


      #Hyper
      if(test.type=="Hyper"){
        #intact cat
        if("cat" %in% general.select){
          if (enrich==T){
            intact_cat_result <- intact.hyper.enrich(X$intact$Category, Y$intact$Category,p=p,q=q)
          }else{
            intact_cat_result <- intact.hyper(X$intact$Category, Y$intact$Category)
          }}
        #intact main
        if("main" %in% general.select){
          if (enrich==T){
            intact_main_result <- intact.hyper.enrich(X$intact$`Main class`, Y$intact$`Main class`,p=p,q=q)
          }else{
            intact_main_result <- intact.hyper(X$intact$`Main class`, Y$intact$`Main class`)
          }}
        #intact sub
        if("sub" %in% general.select){
          if (enrich==T){
            intact_sub_result <- intact.hyper.enrich(X$intact$`Sub class`, Y$intact$`Sub class`,p=p,q=q)
          }else{
            intact_sub_result <- intact.hyper(X$intact$`Sub class`, Y$intact$`Sub class`)
          }}
        #chain
        if("chains" %in% general.select){
          if (enrich==T){
            chain_result <- chain.hyper.enrich(X$chain, Y$chain,p=p,q=q)
          }else{
            chain_result <- chain.hyper(X$chain, Y$chain)
          }}
        #allchains
        if("length" %in% general.select){
          if (enrich==T){
            allchains_result <- allchains.hyper.enrich(X$allchains, Y$allchains,p=p,q=q)
          }else{
            allchains_result <- allchains.hyper(X$allchains, Y$allchains)
          }}
        #sub total carbon
        if("total_carbon" %in% subset.select){
          if (subset.by=="category"){
            sub_total_carbon_result<- total.carbon.cat(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="mainclass"){
            sub_total_carbon_result<- total.carbon.main(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="subclass"){
            sub_total_carbon_result<- total.carbon.sub(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="all"){
            sub_total_carbon_result<- total.carbon.cat(X$intact,Y$intact, test= "Hyper", enrich=enrich,p=p,q=q)
            sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.main(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q))
            sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.sub(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q))
          }
        }
        #sub total DB
        if("total_insaturation" %in% subset.select){
          if (subset.by=="category"){
            sub_total_DB_result<- total.DB.cat(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="mainclass"){
            sub_total_DB_result<- total.DB.main(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="subclass"){
            sub_total_DB_result<- total.DB.sub(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="all"){
            sub_total_DB_result<- total.DB.cat(X$intact,Y$intact, test= "Hyper", enrich=enrich,p=p,q=q)
            sub_total_DB_result<- rbind(sub_total_DB_result,total.DB.main(X$intact,Y$intact, test= "Hyper", enrich=enrich,p=p,q=q))
            sub_total_DB_result<- rbind(sub_total_DB_result,total.DB.sub(X$intact,Y$intact, test= "Hyper", enrich=enrich,p=p,q=q))
          }
        }
        #sub allchains
        if("specific_chains" %in% subset.select){
          if (subset.by=="category"){
            sub_allchains_result<- allchains.cat(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="mainclass"){
            sub_allchains_result<- allchains.main(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="subclass"){
            sub_allchains_result<- allchains.sub(X$intact,Y$intact, test= "Hyper",enrich=enrich,p=p,q=q)
          }
          if (subset.by=="all"){
            sub_allchains_result<- allchains.cat(X$intact,Y$intact,  test= "Hyper",enrich=enrich,p=p,q=q)
            sub_allchains_result<- rbind(sub_allchains_result,allchains.main(X$intact,Y$intact,  test= "Hyper",enrich=enrich,p=p,q=q))
            sub_allchains_result<- rbind(sub_allchains_result,allchains.sub(X$intact,Y$intact,  test= "Hyper",enrich=enrich,p=p,q=q))
          }
        }
      }


        #EASE
        if(test.type=="EASE"){
          #intact cat
          if("cat" %in% general.select){
            if (enrich==T){
              intact_cat_result <- intact.EASE.enrich(X$intact$Category, Y$intact$Category,p=p,q=q)
            }else{
              intact_cat_result <- intact.EASE(X$intact$Category, Y$intact$Category)
            }}
          #intact main
          if("main" %in% general.select){
            if (enrich==T){
              intact_main_result <- intact.EASE.enrich(X$intact$`Main class`, Y$intact$`Main class`,p=p,q=q)
            }else{
              intact_main_result <- intact.EASE(X$intact$`Main class`, Y$intact$`Main class`)
            }}
          #intact sub
          if("sub" %in% general.select){
            if (enrich==T){
              intact_sub_result <- intact.EASE.enrich(X$intact$`Sub class`, Y$intact$`Sub class`,p=p,q=q)
            }else{
              intact_sub_result <- intact.EASE(X$intact$`Sub class`, Y$intact$`Sub class`)
            }}
          #chain
          if("chains" %in% general.select){
            if (enrich==T){
              chain_result <- chain.EASE.enrich(X$chain, Y$chain,p=p,q=q)
            }else{
              chain_result <- chain.EASE(X$chain, Y$chain)
            }}
          #allchains
          if("length" %in% general.select){
            if (enrich==T){
              allchains_result <- allchains.EASE.enrich(X$allchains, Y$allchains,p=p,q=q)
            }else{
              allchains_result <- allchains.EASE(X$allchains, Y$allchains)
            }}
          #sub total carbon
          if("total_carbon" %in% subset.select){
            if (subset.by=="category"){
              sub_total_carbon_result<- total.carbon.cat(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="mainclass"){
              sub_total_carbon_result<- total.carbon.main(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="subclass"){
              sub_total_carbon_result<- total.carbon.sub(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="all"){
              sub_total_carbon_result<- total.carbon.cat(X$intact,Y$intact, test= "EASE", enrich=enrich,p=p,q=q)
              sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.main(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q))
              sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.sub(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q))
            }
          }
          #sub total DB
          if("total_insaturation" %in% subset.select){
            if (subset.by=="category"){
              sub_total_DB_result<- total.DB.cat(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="mainclass"){
              sub_total_DB_result<- total.DB.main(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="subclass"){
              sub_total_DB_result<- total.DB.sub(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="all"){
              sub_total_DB_result<- total.DB.cat(X$intact,Y$intact, test= "EASE", enrich=enrich,p=p,q=q)
              sub_total_DB_result<- rbind(sub_total_DB_result,total.DB.main(X$intact,Y$intact, test= "EASE", enrich=enrich,p=p,q=q))
              sub_total_DB_result<- rbind(sub_total_DB_result,total.DB.sub(X$intact,Y$intact, test= "EASE", enrich=enrich,p=p,q=q))
            }
          }
          #sub allchains
          if("specific_chains" %in% subset.select){
            if (subset.by=="category"){
              sub_allchains_result<- allchains.cat(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="mainclass"){
              sub_allchains_result<- allchains.main(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="subclass"){
              sub_allchains_result<- allchains.sub(X$intact,Y$intact, test= "EASE",enrich=enrich,p=p,q=q)
            }
            if (subset.by=="all"){
              sub_allchains_result<- allchains.cat(X$intact,Y$intact,  test= "EASE",enrich=enrich,p=p,q=q)
              sub_allchains_result<- rbind(sub_allchains_result,allchains.main(X$intact,Y$intact,  test= "EASE",enrich=enrich,p=p,q=q))
              sub_allchains_result<- rbind(sub_allchains_result,allchains.sub(X$intact,Y$intact,  test= "EASE",enrich=enrich,p=p,q=q))
            }
          }
        }

  #KS
  if(test.type=="KS"){
    #intact cat
    if("cat" %in% general.select){
      if (enrich==T){
        intact_cat_result <- intact.KS.enrich(X$intact, colname= "Category", rankingTable , p=p , q=q , order=order)
      }else{

        intact_cat_result <- intact.KS(X$intact, colname= "Category", rankingTable , order=order)
      }}
    #intact main
    if("main" %in% general.select){
      if (enrich==T){
        intact_main_result <- intact.KS.enrich(X$intact, colname= "Main class", rankingTable , p=p , q=q , order=order)
      }else{
        intact_main_result <- intact.KS(X$intact, colname= "Main class", rankingTable , order=order)
      }}
    #intact sub
    if("sub" %in% general.select){
      if (enrich==T){
        intact_sub_result <- intact.KS.enrich(X$intact, colname= "Sub class", rankingTable , p=p , q=q , order=order)
      }else{
        intact_sub_result <- intact.KS(X$intact, colname= "Sub class", rankingTable , order=order)
      }}
    #chain
    if("chains" %in% general.select){
      if (enrich==T){
        chain_result <- chain.KS.enrich(X$chain, rankingTable , p=p , q=q , order=order)
      }else{
        chain_result <- chain.KS(X$chain, rankingTable, order=order)
      }}
    #allchains
    if("length" %in% general.select){
      if (enrich==T){
        allchains_result <- allchains.KS.enrich(X$allchains, rankingTable , p=p , q=q , order=order)
      }else{
        allchains_result <- allchains.KS(X$allchains, rankingTable , order=order)
      }}

    #sub total carbon
    if("total_carbon" %in% subset.select){
      if (subset.by=="category"){
        sub_total_carbon_result<- total.carbon.cat.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="mainclass"){
        sub_total_carbon_result<- total.carbon.main.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="subclass"){
        sub_total_carbon_result<- total.carbon.sub.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="all"){
        sub_total_carbon_result<- total.carbon.cat.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
        sub_total_carbon_result<- rbind(sub_total_carbon_result, total.carbon.main.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q))
        sub_total_carbon_result<- rbind(sub_total_carbon_result,total.carbon.sub.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q))
      }
    }
    #sub total DB
    if("total_insaturation" %in% subset.select){
      if (subset.by=="category"){
        sub_total_DB_result<- total.DB.cat.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="mainclass"){
        sub_total_DB_result<- total.DB.main.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="subclass"){
        sub_total_DB_result<- total.DB.cat.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="all"){
        sub_total_DB_result <- total.DB.cat.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
        sub_total_DB_result <- rbind(sub_total_DB_result, total.DB.main.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q))
        sub_total_DB_result <- rbind(sub_total_DB_result, total.DB.cat.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q))
      }
    }
    #sub allchains
    if("specific_chains" %in% subset.select){
      if (subset.by=="category"){
        sub_allchains_result<- allchains.cat.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="mainclass"){
        sub_allchains_result<- allchains.main.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="subclass"){
        sub_allchains_result<- allchains.sub.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
      }
      if (subset.by=="all"){
        sub_allchains_result <- allchains.cat.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q)
        sub_allchains_result <- rbind(sub_allchains_result, allchains.main.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q))
        sub_allchains_result <- rbind(sub_allchains_result, allchains.sub.KS(X$intact , rankingTable , order = order , enrich=enrich , p=p , q=q))
      }
    }
  }

  #Output tables if the corresponding object is not NA
  #section titles object
  title.enrich<-NA
  #create a global output object (for a global download of the results)
  if(test.type=="KS"){
    global.output<-data.frame(matrix(nrow=0,ncol=5))
    colnames(global.output)<-c("Test.performed","Classifier","Count.in.list", "p-value", "FDR.q-value")
  } else {
    global.output<-data.frame(matrix(nrow=0,ncol=9))
    colnames(global.output)<-c("Test.performed","Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")
  }
  
  options(warn=-1)
  if(exists("intact_cat_result")){
    title.enrich<- paste("Category(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    intact_cat_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,intact_cat_result))
  }

  if(exists("intact_main_result")){
    title.enrich<- paste("Main class(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    intact_main_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,intact_main_result))
  }

  if(exists("intact_sub_result")){
    title.enrich<- paste("Sub class(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    intact_sub_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,intact_sub_result))
  }

  if(exists("chain_result")){
    title.enrich<- paste("Chain(s) characteristics(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    chain_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,chain_result))
  }

  if(exists("allchains_result")){
    title.enrich<- paste("Specific chain(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    allchains_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich,allchains_result))
  }

  if(exists("sub_total_carbon_result")){
    title.enrich<- paste("Total chain carbon by ",subset.by,"(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    sub_total_carbon_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich, sub_total_carbon_result))
  }

  if(exists("sub_total_DB_result")){
    title.enrich<- paste("Total number of DB by ",subset.by,"(",test.type,")",sep="")
    #display title
    title.enrich
    #display this on the right panel
    sub_total_DB_result
    global.output<-rbind(global.output,cbind(Test.performed=title.enrich, sub_total_DB_result))
  }

  if(exists("sub_allchains_result")){
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
