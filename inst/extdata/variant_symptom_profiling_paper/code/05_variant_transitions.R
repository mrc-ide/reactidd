
#' First clear the environment of variables
rm(list=ls(all=TRUE))
# get root director of project
root.dir <- getwd()
# setwd(dir = "/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/")
outpath <- paste0(root.dir,"/output/")
figpath <-  paste0(root.dir,"/plots/")



source("E:/Group/functions/load_packages.R", local = T)
source("E:/Group/functions/full_join_multiple_datasets.R", local = T)
source("E:/Group/functions/wrangle_cleaner_functions.R", local = T)
source("E:/Group/functions/cats_and_covs.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/create_subfolder.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/forest_plot.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/save_styled_table.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/stability_selection.R", local = T)


#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr",
                  "tidyr", "pheatmap","ggbrace","Tabulizer","datapasta",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus","OverReact"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------

source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_data_prep.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "variant_transitions")


### DATA from www.covid19.sanger.ac.uk/lineages/raw
variant_data=read_csv("Proportion_in_England.csv")
variant_data$date <- as.Date(variant_data$date)


dataCreator <- function(startdate=as.Date("2020-09-05"),numweeks=10,by="-1 week",lineage = "B", rev=T){
  datelist = (seq(as.Date(startdate),by=by,length.out=numweeks))
  if(rev){
    datelist=rev(datelist)
    datelist <- datelist[1:(length(datelist)-1)]
  }else{
    datelist <- datelist[2:(length(datelist))]
  }
  df=data.frame(date=datelist,lineage = lineage, value = 100, upper=100,lower=100)
  return(df)
}

## add dummy data
variant_data <- rbind(dataCreator(startdate = as.Date("2020-09-05"),numweeks=36),
                      variant_data)
                      # (dataCreator(startdate = as.Date("2022-04-12"),numweeks=2,by = "1 week",lineage = "BA.2", rev=F)))

# wrangle variant data
p_variant_transition <- variant_data %>% 
  mutate(lineage=factor(case_when(lineage=="B.1.1.7"~"Alpha",
                                  lineage%in%c("AY.4.2","B.1.617.2")~"Delta",
                                  grepl("BA.1", lineage) ~ "BA.1",
                                  grepl("BA.2", lineage) ~ "BA.2",
                                  # lineage%in%c("BA.1","BA.2","BA.1.1")~"Omicron",
                                  T ~ "Wildtype/other"
  ),levels = c("Wildtype/other","Alpha","Delta","BA.1","BA.2"))) %>% 
  group_by(date,lineage) %>% 
  summarise(value=sum(value)) %>% 
  ggplot(aes(x=date, y=value, fill = lineage)) +
  geom_area(position = "fill", colour = "white", alpha=0.8) +
  # scale_fill_imperial(palette = "extended") +
  OverReact::scale_fill_imperial(palette = "default") +
  OverReact::theme_react() + 
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  theme_adjust +
  labs(x="Date", y="Proportion of typed \ncases in UK", fill = "") +
  theme(axis.text.x = element_text(angle=45, vjust =0.5, hjust =0.5))

  
  
  # define a function to add a date annotation bracket to the plot, for iteration
roundAnnotator <- function(p, r){
  start_date = round_dates$Fieldwork_start[[r]]
  end_date = round_dates$Fieldwork_end[[r]]
  mid_date =round_dates$Fieldwork_midpoint[[r]]
  
  p=p+geom_brace(aes(x=c(as.Date(start_date),as.Date(end_date)), y = c(1,1.05)), inherit.data = F,bending = 0.1)+
    geom_vline(xintercept = as.Date(start_date), linetype="dashed", size=0.05, col = "white") +
    geom_vline(xintercept = as.Date(end_date), linetype="dashed", size=0.05, col = "white") +
    annotate("rect",xmin=start_date,xmax=end_date,ymin=0,ymax=1, fill = "white", alpha=0.1) +
    annotate("text",x=as.Date(mid_date), y=1.07,  size=2,
             label = paste0("R",r), colour="black")
    
  p
  return(p)
}
  
# Loop over all rounds
# p=p_variant_transition
for (r in 1:19){
  p_variant_transition=roundAnnotator(p = p_variant_transition,r = r)
}
  
p_variant_transition

OverReact::saveREACTplot(p = p_variant_transition,figpath = figpath,filename = "variant_transition_plot",
                         width = 11,height = 3)

OverReact::saveREACTplot(p = p_variant_transition,figpath = figpath,filename = "variant_transition_plot_tall",
                         width = 11,height = 6)




