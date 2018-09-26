#' QueryWithinUniverse()
#'
#' Verifies if the lipids of the a <Query> character vector are contained in a <Universe> character vector.
#'
#' @param X a character vector containing a lipid list
#' @param Y a character vector containing a lipid list
#'
#' @details Verifies if the lipids of the a <Query> character vector are contained in a <Universe> character vector.
#'
#' @return This function will return an error message if the all the elements of a <Query>character vector are not a subset of an <Universe>character
#'
#' @examples QueryWithinUniverse(cleaned.queryExample,cleaned.universeExample)
#'
#' @author Geremy Clair
#' @export
#'
QueryWithinUniverse<-function(X,Y){
  if(!is.character(X)){stop("X has to be a character vector")}
  if(!is.character(Y)){stop("Y has to be a character vector")}

  if(length(X[!X %in% Y])>0) {
    warning("some lipids are found in the query but not in the universe")
    X[!X %in% Y]
    }

  }
