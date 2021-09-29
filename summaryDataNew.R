#' @title Compute attribute level plot data
#' @description function to compute the distribution of the target with respect to an attribute.
#' @details This function will compute the weighted averages of the observed response, the model response, model coefficients, and model residuals.
#' @param obj The model object
#' @param dat data.table containing the attribute for visualization
#' @param selectedVar Visualization attribute name (string)
#' @param linearPredictor Flag for computing results in linear/link space(Depreciated)
#'
#' @return List containing the linear and link observed, fitted average, model, and resudual weighted averages for \code{selectedVar}
#'
#' @examples \dontrun{
#'     myData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
#'     model <- glm(admit ~gre+rank,data=myData,family = 'binomial')
#'     summaryDataNew(model,myData,"gre",linearPredictor = TRUE)
#' }
#'
#' @author Luke
#' @importFrom dplyr %>% group_by summarise
#' @importFrom stats weighted.mean
#' @importFrom forcats fct_explicit_na
#' @rawNamespace improt(data.table,except = c(first,melt))
#'
#' @export

summaryDataNew <- function(obj,dat,selectedVar,linearPredictor=TRUE){

    asNumericFactor <- function(x) {as.numeric(levels(x))[x]}

    weight <- NULL
    resp <- NULL
    pred <- NULL

    num_var <-setdiff(names(dat),names(dat[,which(sapply(dat,is.factor))]))
    datnew <- data.table(
        var = dat[[selectedVar]],
        resp = obj$y,
        weight = obj$prior.weights,
        pred = obj$fitted.values
    )

    if(any(is.na(datnew$var)) && is.factor(datnew$var)){
        datnew$var <- forcats::fct_explicit_na(data$var,'NaN')
    }

    uniqueVals <- length(unique(datnew$var))

    totalWeight <- sum(datnew$weight)

    if(is.numeric(datnew$var) && uniqueVals >250){
        oldValue <- datnew$var
        oldValue[which(datnew$ar <0)] <- NA
        specials <- c(991,994,995,999,9992,9999999.9)
        excludes <- c()
        uniqueOld <- unique(oldValue)
        for(i in rev(specials)){
            if (i %in% oldValue && i >= max(uniqueOld,na.rm=TRUE)){
                excludes <- c(excludes,i)
                uniqueOld <- uniqueOld[which(uniqueOld != max(uniqueOld,na.rm = TRUE))]
            }
            excludes <- rev(excludes)
        }
        rm(uniqueOld)
        oldValue[which(datnew$var %in% excludes)] <- NA
        oldMin <- min(oldValue, na.rm = TRUE)
        oldMax <- max(oldValue, na.rm = TRUE)
        newMin <- 1
        newMax <- 250
        oldRange <- (oldMax - oldMin)
        newRange <- (newMax - newMin)
        newValue <- (((oldValue-oldMin)*newRange)/oldRange)+newMin
        newValue <- round(newValue)
        finalValue <- (((newValue - newMin)* oldRange)/newRange)+oldMin
        finalValue[which(is.na(finalValue))] <- datnew$var[which(is.na(finalValue))]
        datnew$var <- finalValue
    }

    if(nrow(datnew) <5e6){
        varDF <- datnew %>% group_by(var,.drop=FALSE)

        EE <- as.data.table(varDF %>% summarise(average=sum(weight)/totalWeight))
        obs <- as.data.table(varDF %>% summarise(average = obj$family$linkfun(weighted.mean(resp,weight))))
        ca <- as.data.table(varDF %>% summarise(average = obj$family$linkfun(weighted.mean(pred,weight))))
    } else {
        res <- datnew[,.(EE=sum(weight)/totalWeight,obs=obj$family$linkfun(weighted.mean(resp,weight)),ca=obj$family$linkfun(weighted.mean(pred,weight))),by=var]
        if (is.factor(datnew$var)){
            res <- res[data.table::data.table("var" = levels(datnew$var)),on=.(var=var)][,EE := ifelse(is.na(EE),0,EE)]
            res[is.na(res)] <- NaN
        }

        data.table::setorder(res,var)

        EE <- data.table::data.table("var" = res[["var"]],"average" = res[["EE"]])
        obs <- data.table::data.table("var" = res[["var"]],"average" = res[["obs"]])
        ca <- data.table::data.table("var" = res[["var"]], "average" = res[["ca"]])
    }

    if(is.null(obj$offset)){
        offset <- 0
    } else {
        offset <- weighted.mean(obj$offset,obj$prior.weights)
    }

    if(selectedVar %in% num_var){
        if(selectedVar %in% names(obj$coefficients)){
            numLevs <- obs$var
            cm <-
                data.table(
                    x = numLevs,
                    y = obj$coefficients[which(names(obj$coefficients)==selectedVar)] * numLevs +obj$coefficients[1] +offset)
        } else {
            numLevs <- obs$var
            cm <-
                data.table(
                    x = numLevs,
                    y = 0 * numLevs + obj$coefficients[1] +offset)
            }
        } else {
            if(length(grep(paste0('^',selectedVar),names(obj$coefficients))) != 0){
                substring <- substr(selectedVar,1,nchar(selectedVar) - 1)
                levs <- levels(datnew$var)
                coefs <- c(levs[1],0)
                for (i in 2:length(levs)){
                    exp <- paste0('^',substring,'.',levs[i])
                    coefs <- rbind(coefs,c(levs[i],as.numeric(obj$coefficients[grep(exp,names(obj$coefficients))])))
                }
                coefs <- as.data.frame(coefs)
                if(all(InsureR::checkNumeric(coefs$V2))){
                    names(coefs) <- c("level","coefficient")
                    rownames(coefs) <-c()

                    cm <-
                        data.table(x = as.numeric(coefs$level),
                                   y = asNumericFactor(coefs$coefficient) + obj$coefficients[1] +offset)
                }
                } else{
                    temp <- levels(as.factor(datnew$var))
                    temp <- as.numeric(as.data.frame(temp)$temp)
                    cm <-
                        data.table(x = temp,
                                   y = 0 * temp + obj$coefficients[1] + offset)
                }
            }

            cu <- data.table(x = cm$x,
                             y = cm$y + obs$average - ca$average)
            names(obs) <- c("x","y")
            names(ca) <- c("x","y")
            names(cm) <- c("x","y")
            names(cu) <- c("x","y")
            names(EE) <- c("x","y")

            res <- list(
                linear = data.frame(
                    x = obs$x,
                    obs = obs$y,
                    ca = ca$y,
                    cm = cm$y,
                    cu = cu$y,
                    EE = EE$y
                ),
                link = data.frame(
                    x = obs$x,
                    obs = obj$family$linkinv(obs$y),
                    ca = obj$family$linkinv(ca$y),
                    cm = obj$family$linkinv(cm$y),
                    cu = obj$family$linkinv(cu$y),
                    EE = EE$y
                )
            )

            res$linear[res$linear == -Inf] <- NA
            res$link[res$link == -Inf] <- NA

            res$linear[is.na(res$linear)] <- NA
            res$link[is.na(res$link)] <- NA

            return(res)
}
