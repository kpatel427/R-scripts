# script to generate survival (Kaplan-Meier) plots using expression data
# Author: Khushbu Patel
#--------------------------------------------------------------------------------------------------#

getTimeStatus <- function(endpoint){

  if(endpoint == "overall survival") {
      time <- 'nti_surv_overall'
      status <- 'nti_event_overall_num'
    } else {
      time <- 'nti_surv_progrfree'
      status <- 'nti_event_progrfree_num'
    }
  
  result <- c(time,status)
  return(result)

}
survivalKM <- function(exprData, metadata, gene, endpoint, Risk) {

  result <- getTimeStatus(endpoint)
  time <- result[1]
  status <- result[2]
  
  #timeSurv <- metadata[,time]
  #statusSurv <- metadata[,status]
  
  # Step1: Prepare data ----------------------------------
  # condition for specific RISK group
  if (Risk != 'All') {
    metadata <- metadata[metadata$RISK == Risk,]
    exprData <- exprData[, which(colnames(exprData) %in% rownames(metadata))]
  } else {
    metadata <- metadata
    exprData <- exprData
  }
  
  # subset metadata
  metadata <- metadata[,c("RISK",time,status)]
  
  # subset expr data
  exprData <- exprData[gene,]
  
  
  # Step2: Normalizing expression data ----------------------------------
  # z-score > 0 = high expression
  # z-score < 0 = low expression
  exprData <- data.frame("gene.Z" = apply(X = exprData, MARGIN = 1, FUN = scale))
  

  # adding expr data to metadata
  metadata <- cbind(metadata, exprData)
  metadata <- metadata %>%
      gather(key = 'gene', value = 'zScore', -c("RISK",time,status))

  metadata$gene <- gsub('gene.Z.','',metadata$gene)
  metadata$gene <- as.factor(metadata$gene)
  # assigning exprStatus
  # 1 = high
  # 0 = low
  # 2 = mid expression
  # ignore mid expressions
  metadata$exprStatus <- ifelse(metadata$zScore > 0, 1,
                                ifelse(metadata$zScore < 0, 0, 2))


  # remove the ones with mid expressions
  metadata <- metadata[metadata$exprStatus != 2,]

  # Step3: Fitting survival curves ----------------------------------
  group = "exprStatus"

  #fit <- survfit(Surv(timeSurv, statusSurv) ~ exprStatus, data = metadata)
  fit <- do.call(survfit, list(Surv(eval(parse(text=time)), eval(parse(text=status))) ~ exprStatus, data = metadata))


  # Step4: Comparing survival times between groups ----------------------------------
  # conduct between-group significance tests using a log-rank test
  # between high-low expression groups

  diff <- survdiff(formula = Surv(eval(parse(text=time)), eval(parse(text=status))) ~ exprStatus, data = metadata)

  pval <- pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)

  adjpval <- p.adjust(pval, method = "fdr",  n = (fit$n[1] + fit$n[2]))

  plotTitle <- paste0(paste(gene, collapse = " | "), " P-val(Adj) :", format(pval, scientific=T, digits=3), "(", format(adjpval, scientific=T, digits=3), ")")



  # Step5: Generate plot ----------------------------------
  # change strata names
  low <- paste0("Low : n = ", fit$n[1])
  high <- paste0("High : n = ", fit$n[2])
  names(fit$strata) <- c(low,high)

  plot <- ggsurvplot(fit,
             data = metadata,
             pval = TRUE,
             pval.size = 5,
             conf.int = TRUE,
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             surv.median.line = "hv", # Specify median survival
             ggtheme = theme_Publication_scatter(base_size = 14), # Change ggplot2 theme
             risk.table.fontsize = 5,
             palette = c("#E7B800", "#2E9FDF"),
             title = plotTitle) + xlab('Survival Time')

  return(plot)
} # function ends here
