setwd("R:/BSRBR/Analyses/lucy_bull/PhD work/Simulation study/Simulation_output_Feb22/")
#setwd("R:/BSRBR/Analyses/lucy_bull/Mark_CSF_output/")
require(rlist)

performance_files <- list.files(pattern = "performance_results_")
stability_files <- list.files(pattern = "stability_results_")


performance_list <- lapply(performance_files, function(file.name)
{
  df           <- read.csv(file.name)
  df$file.name <- file.name
  return(df)
})
performance_results_all <- list.rbind(performance_list)


stability_list <- lapply(stability_files, function(file.name)
{
  df           <- read.csv(file.name)
  df$file.name <- file.name
  return(df)
})

stability_results_all <- list.rbind(stability_list)

# local
# write.table(performance_results, file = "R:/BSRBR/Analyses/lucy_bull/PhD work/Simulation study/Functions/Cleaner_code/performance_results_all.csv", sep = ",", col.names = NA,
#            qmethod = "double")
# write.table(stability_results, file = "R:/BSRBR/Analyses/lucy_bull/PhD work/Simulation study/Functions/Cleaner_code/stability_results_all.csv", sep = ",", col.names = NA,
#            qmethod = "double")

write.table(performance_results_all, file = "./performance_results_all.csv", sep = ",", col.names = NA,
            qmethod = "double")
write.table(stability_results_all, file = "./stability_results_all.csv", sep = ",", col.names = NA,
            qmethod = "double")