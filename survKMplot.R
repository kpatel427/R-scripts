# script to ge survival data for AML TARGET data (open access data)
# generates KM plots for EIF4E gene
# plot shows survival probability over days among various RISK groups

library(data.table)
library(tidyverse)
library(dplyr)
library(readxl)
library(survival)
#install.packages('survminer')
library(survminer)



# making survival (KM) plot
# OS
ggsurv_os <- ggsurvplot(
  fit = survfit(Surv(time = overall_survival_days, event = overall_survival_event) ~ RISK, data = EIF4E_subset), 
  title = 'Overal survival across various risk groups for EIF4E',
  xlab = 'Days', 
  ylab = 'Survival probability',
  legend.title = 'RISK groups',
  legend.labs = c('High','Low','Standard','Unknown'),
  #conf.int=TRUE,
  pval = TRUE,
  pval.size = 4,
  fun = 'pct', 
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.35,
  risk.table.y.text = FALSE, # show only legend bars on y.axis of risk table
  #linetype = "strata",
  palette = c("#E7B800","#2E9FDF", "#228B22","#D2691E"),
  legend = 'bottom',
  surv.median.line = "hv")


ggsurv_os$plot <- ggsurv_os$plot + theme(plot.title = element_text(hjust = 0.5))
ggsurv_os$table <- ggsurv_os$table + theme(plot.title = element_text(hjust = 0.5))

print(ggsurv_os)

# saving the plot
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file=paste0(Sys.Date(),'_OS_survival_EIF4E_TARGET_AML.pdf'), width = 8, height = 6)
print(ggsurv_os)
dev.off()


# EFS
ggsurv_efs <- ggsurvplot(
  fit = survfit(Surv(time = event_free_survival_days, event = event_free_survival_event) ~ RISK, data = EIF4E_subset), 
  title = 'Event free survival across various risk groups for EIF4E',
  xlab = 'Days', 
  ylab = 'Survival probability',
  legend.title = 'RISK groups',
  legend.labs = c('High','Low','Standard','Unknown'),
  #conf.int=TRUE,
  pval = TRUE,
  pval.size = 4,
  #pval.coord = c(1000, 50),
  fun = 'pct', 
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.35,
  risk.table.y.text = FALSE, # show only legend bars on y.axis of risk table
  #linetype = "strata",
  palette = c("#E7B800","#2E9FDF", "#228B22","#D2691E"),
  legend = 'bottom',
  surv.median.line = "hv")


ggsurv_efs$plot <- ggsurv_efs$plot + theme(plot.title = element_text(hjust = 0.5))
ggsurv_efs$table <- ggsurv_efs$table + theme(plot.title = element_text(hjust = 0.5))

print(ggsurv_efs)

# saving the plot
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file=paste0(Sys.Date(),'_EFS_survival_EIF4E_TARGET_AML.pdf'), width = 8, height = 6)
print(ggsurv_efs)
dev.off()
