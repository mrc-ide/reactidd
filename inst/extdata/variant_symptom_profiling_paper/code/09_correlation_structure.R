
#' First clear the environment of variables
rm(list=ls(all=TRUE))
# get root director of project
root.dir <- getwd()
setwd(dir = "/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/")
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

# 
# 
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_bits_and_pieces.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_functions.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/symptom_prediction_children/code/00_bits_and_pieces.R", local = T)


#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr",
                  "tidyr", "pheatmap","scales","OverReact","ggstance",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus","ComplexHeatmap"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "correlation_structure")



# Create correlation matrices ---------------------------------------------

corlist <- list()
annolist <- list()

variants <- levels(dfRes$variant_inferred_detail)

for(i in 1:length(variants)){
  print(i)
  myvariant=variants[[i]]
  corlist[[i]] <- dfRes %>% filter(variant_inferred_detail==myvariant, estbinres==1) %>% 
    select(all_of(sympnames_type_df$symptom_code[1:26])) %>% 
    cor(use = "pairwise.complete.obs") %>% round(2)
  rownames(corlist[[i]])=colnames(corlist[[i]])=sympnames_type_df$symptom[1:26]
  annolist[[i]]=dfRes %>% filter(variant_inferred_detail==myvariant, estbinres==1) %>% 
    select(all_of(sympnames_type_df$symptom_code[1:26])) %>% 
    pivot_longer(cols =sympnames_type_df$symptom_code[1:26]) %>% 
    group_by(name) %>% 
    summarise(nobs=n(),
              n=sum(value,na.rm=T),
              perc=100*n/nobs) %>% 
    left_join(sympnames_type_df, by=c("name"="symptom_code"))
}
names(corlist)=names(annolist)=variants


### PLot heatmap
col_fun=circlize::colorRamp2(breaks = c(0,0.5,1),
                             colors = c("white",myCols[[2]],myCols[[1]]))



createHeatmap <- function(mymatrix,mymatrixname,percs){
  
  myheatmap=ComplexHeatmap::Heatmap(name = "Correlation",
                                    mymatrix,
                                    cluster_rows = F,cluster_columns = F,
                                    cell_fun = function(j,i,x,y,width, height, fill){
                                      grid.text(sprintf("%.2f", mymatrix[i,j]),
                                                x, y, gp=gpar(fontsize=4, col =
                                                                ifelse(mymatrix[i,j]==1,myCols[[1]],
                                                                       ifelse(mymatrix[i,j]>0.4,
                                                                       "white",
                                                                       "black"))))
                                    },
                                    row_title_gp = gpar(fontsize=7, fontface="bold"),
                                    column_title_gp = gpar(fontsize=11, fontface="bold"),
                                    row_split = sympnames_type_df$symptom_type[1:26],
                                    column_split = sympnames_type_df$symptom_type[1:26],
                                    row_names_gp = gpar(fontsize=7),
                                    column_names_gp = gpar(fontsize=7),
                                    top_annotation = columnAnnotation(
                                      `%`=anno_barplot(percs$perc, gp=gpar(fill=myCols[[1]]),
                                                        border=F,height = unit(1.7,"cm"),
                                                            ylim = c(0,50),
                                                          # axis_param = list(direction = "reverse")
                                                          )),
                                    # title = mymatrixname,
                                    # row_names_side = "left",
                                    # column_names_side = "top",
                                    column_title = mymatrixname,
                                    row_gap = unit(1.5,"mm"),
                                    column_gap = unit(1.5,"mm"),
                                    rect_gp = gpar(col = "white", lwd = 1),
                                    col = col_fun)
  
  return(myheatmap)
}

heatmaps=NULL
for(i in 1:5){
  heatmaps=heatmaps+createHeatmap(corlist[[i]],names(corlist)[[i]],annolist[[i]])
}

draw(heatmaps,ht_gap=unit(1,"cm"))


png(filename = paste0(figpath,"heatmaps_variants_5panel.png"),
    width = 26,height = 7,units = "in",res = 300)
draw(heatmaps,ht_gap=unit(1.2,"cm"))
dev.off()


# Heatmap of all positives ------------------------------------------------

cormat <- dfRes %>% filter(estbinres==1) %>% 
  select(all_of(sympnames_type_df$symptom_code[1:26])) %>% 
  cor(use = "pairwise.complete.obs") %>% round(2)
percs=dfRes %>% filter(estbinres==1) %>% 
  select(all_of(sympnames_type_df$symptom_code[1:26])) %>% 
  pivot_longer(cols =sympnames_type_df$symptom_code[1:26]) %>% 
  group_by(name) %>% 
  summarise(nobs=n(),
            n=sum(value,na.rm=T),
            perc=100*n/nobs) %>% 
  left_join(sympnames_type_df, by=c("name"="symptom_code"))

rownames(cormat)=colnames(cormat)=sympnames_type_df$symptom[1:26]
heatmap_main=createHeatmap(mymatrix = cormat,mymatrixname = "All positives",percs = percs)

png(filename = paste0(figpath,"heatmaps_all_positives.png"),
    width = 7,height = 7,units = "in",res = 300)
heatmap_main
dev.off()



# Correlation between correlation matrices --------------------------------
i=1
j=2
cor_cors=matrix(data = NA,nrow = 5,ncol=5,dimnames = list(variants,variants))
cor_ps=matrix(data = NA,nrow = 5,ncol=5,dimnames = list(variants,variants))

for(i in 1:5){
  for(j in 1:5){
  speartest=cor.test(x = corlist[[i]][upper.tri(corlist[[i]], diag=F)],
                y = corlist[[j]][upper.tri(corlist[[j]], diag=F)],conf.level = 0.95,
                method = "pearson", null.value=1)
  se=speartest$estimate -speartest$conf.int[1]
  zscore=(1-speartest$estimate)/se
  
  cor_cors[i,j] <- speartest$estimate
  cor_ps[i,j] <- pnorm(q = zscore,lower.tail = F)
  print(speartest)
  }
}




printmat=matrix(data = NA,nrow = 5,ncol=5,dimnames = list(variants,variants))

for(i in 1:5){
  for(j in 1:5){
    printmat[i,j]=paste0("r=",round(cor_cors[i,j],2),"\n\np=", 
                         # scales::scientific(
                           round(cor_ps[i,j],6)
                         # ,digits=2
                         # )
    )
  }
}
printmat


printmat_out=as.data.frame(cbind(Variant=variants,printmat))

savePrettyExcelWorkbook(listOfTables = list(correlations=as.data.frame(printmat_out)),
                        workbookName = "correlation_structure_tests",
                        outpath = outpath)


# Plot heatmap ------------------------------------------------------------



myheatmap_corcors=ComplexHeatmap::Heatmap(name = "Correlation\nbetween\ncorrelation\nmatrices",
                                          cor_cors,
                                  cluster_rows = F,cluster_columns = F,
                                  cell_fun = function(j,i,x,y,width, height, fill){
                                    grid.text(printmat[i,j],
                                              x, y, gp=gpar(fontsize=7, col =
                                                              ifelse(i==j,NA,
                                                                     ifelse(cor_cors[i,j]>0.75,
                                                                            "white",
                                                                            "black"))))
                                  },
                                  row_title_gp = gpar(fontsize=9, fontface="bold"),
                                  column_title_gp = gpar(fontsize=11, fontface="bold"),
                                  # row_split = sympnames_type_df$symptom_type[1:26],
                                  # column_split = sympnames_type_df$symptom_type[1:26],
                                  row_names_gp = gpar(fontsize=9),
                                  column_names_gp = gpar(fontsize=9),
                                  
                                  # column_title = mymatrixname,
                                  row_gap = unit(1.5,"mm"),
                                  column_gap = unit(1.5,"mm"),
                                  rect_gp = gpar(col = "white", lwd = 1),
                                  col = circlize::colorRamp2(breaks = c(0.6,0.8,1),
                                                             colors = c("white",myCols[[2]],myCols[[1]])))
myheatmap_corcors


png(filename = paste0(figpath,"heatmaps_correlation_between_matrices.png"),
    width = 6,height = 5,units = "in",res = 300)
myheatmap_corcors
dev.off()



# 
# # Modelling sense of smell/taste asfunction of round + taste --------------
# 
# sympnames_type_df
# class(dfRes$variant_inferred_detail)
# mod_dat=dfRes %>% filter(estbinres==1, variant_inferred_detail %in% variants[c(3,5)]) %>% 
#   mutate(omicron=case_when(variant_inferred_detail%in%variants[4:5] ~ 1,
#                            T~ 0))
# f=as.formula("symptnowaw_2~symptnowaw_1*omicron")
# mod <- glm(formula = f,data = mod_dat,family = "binomial")
# 
# jtools::summ(mod)
# table(mod_dat$omicron)







# Create coccurrence matrices ---------------------------------------------

corlist <- list()
annolist <- list()

variants <- levels(dfRes$variant_inferred_detail)

for(i in 1:length(variants)){
  print(i)
  myvariant=variants[[i]]
  df_for_cp <-dfRes %>% 
    filter(variant_inferred_detail==myvariant) %>% 
    filter(estbinres==1, symptomatic==1) %>% 
    select(sympnames_type_df$symptom_code[1:26]) %>% 
    as.matrix()
  df_for_cp[is.na(df_for_cp)] <- 0
  cooccurmat <- df_for_cp %>% crossprod()
  cooccurmat=cooccurmat/nrow(df_for_cp)
  rownames(cooccurmat)=colnames(cooccurmat)=sympnames_type_df$symptom[1:26]
  cooccurmat <- 100*round(cooccurmat,2)
  marginal_counts=diag(cooccurmat)
  
  
  corlist[[i]] <- cooccurmat
  annolist[[i]]=marginal_counts
}
names(corlist)=names(annolist)=variants




createHeatmapCoccur <- function(cooccurmat,mymatrixname,marginal_counts){
  
  
  myheatmap_cooccur_first=ComplexHeatmap::Heatmap(name = "% symptom \noccurrence or\nco-occurrence",
                                                  cooccurmat,
                                                  cluster_rows = F,cluster_columns = F,
                                                  cell_fun = function(j,i,x,y,width, height, fill){
                                                    grid.text(cooccurmat[i,j],
                                                              x, y, gp=gpar(fontsize=6, col =
                                                                              ifelse(cooccurmat[i,j]>40,
                                                                                     "white",
                                                                                     "black")))
                                                  },
                                                  row_title_gp = gpar(fontsize=7, fontface="bold"),
                                                  column_title_gp = gpar(fontsize=7, fontface="bold"),
                                                  row_split = sympnames_type_df$symptom_type[1:26],
                                                  column_split = sympnames_type_df$symptom_type[1:26],
                                                  row_names_gp = gpar(fontsize=9),
                                                  column_names_gp = gpar(fontsize=9),
                                                  top_annotation = columnAnnotation(
                                                    `%`=anno_barplot(marginal_counts, gp=gpar(fill=myCols[[1]]),
                                                                     border=F,height = unit(1.7,"cm"),
                                                                     ylim = c(0,70),
                                                                     # axis_param = list(direction = "reverse")
                                                    )),
                                                  column_title = mymatrixname,
                                                  row_gap = unit(1.5,"mm"),
                                                  column_gap = unit(1.5,"mm"),
                                                  rect_gp = gpar(col = "white", lwd = 1),
                                                  col = circlize::colorRamp2(breaks = c(00,40,80),
                                                                             colors = c("white",
                                                                                        myCols[[2]],
                                                                                        myCols[[1]])))
  
  return(myheatmap_cooccur_first)
}

heatmaps_coccur=NULL
for(i in 1:5){
  heatmaps_coccur=heatmaps_coccur+createHeatmapCoccur(cooccurmat = corlist[[i]],
                                                      mymatrixname  = names(corlist)[[i]],
                                                      marginal_counts = annolist[[i]])
}

draw(heatmaps_coccur,ht_gap=unit(1,"cm"))


png(filename = paste0(figpath,"heatmaps_variant_cooccurrence_5panel.png"),
    width = 29,height = 7,units = "in",res = 300)
draw(heatmaps_coccur,ht_gap=unit(1.2,"cm"))
dev.off()


# Heatmap of all positives ------------------------------------------------

cormat <- dfRes %>% filter(estbinres==1) %>% 
  select(all_of(sympnames_type_df$symptom_code[1:26])) %>% 
  cor(use = "pairwise.complete.obs") %>% round(2)
percs=dfRes %>% filter(estbinres==1) %>% 
  select(all_of(sympnames_type_df$symptom_code[1:26])) %>% 
  pivot_longer(cols =sympnames_type_df$symptom_code[1:26]) %>% 
  group_by(name) %>% 
  summarise(nobs=n(),
            n=sum(value,na.rm=T),
            perc=100*n/nobs) %>% 
  left_join(sympnames_type_df, by=c("name"="symptom_code"))

rownames(cormat)=colnames(cormat)=sympnames_type_df$symptom[1:26]
heatmap_main=createHeatmap(mymatrix = cormat,mymatrixname = "All positives",percs = percs)

png(filename = paste0(figpath,"heatmaps_all_positives.png"),
    width = 7,height = 7,units = "in",res = 300)
heatmap_main
dev.off()










# First symptom cooccurrence ----------------------------------------------


# add first symptoms to covnamelist
cov_name_list[sympnames_type_df$first_symptom_code[1:26]] <- sympnames_type_df$symptom[1:26]

# Sort out symptoms
dfRes <- dfRes %>% 
  mutate_at(all_of(sympnames_type_df$first_symptom_code[1:26]),binaryCleaner_1_0) 


### Replace NA with 0 in symptoms
dfRes[sympnames_type_df$first_symptom_code[1:26]][is.na(dfRes[sympnames_type_df$first_symptom_code[1:26]])] <- 0
dfRes$variant_inferred %>% table()
df_for_cp <-dfRes %>% 
  filter(round %in% c(17:19)) %>% 
  filter(estbinres==1, symptomatic==1) %>% 
  select(sympnames_type_df$first_symptom_code[1:26]) %>% 
  as.matrix()
df_for_cp[is.na(df_for_cp)] <- 0
cooccurmat <- df_for_cp %>% crossprod()
cooccurmat=cooccurmat/nrow(df_for_cp)
rownames(cooccurmat)=colnames(cooccurmat)=sympnames_type_df$symptom[1:26]
# diag(cooccurmat_L1) <- 0
cooccurmat <- 100*round(cooccurmat,3)
marginal_counts=diag(cooccurmat)


myheatmap_cooccur_first=ComplexHeatmap::Heatmap(name = "% co-occurrence \nas first symptom",
                                          cooccurmat,
                                          cluster_rows = F,cluster_columns = F,
                                          cell_fun = function(j,i,x,y,width, height, fill){
                                            grid.text(cooccurmat[i,j],
                                                      x, y, gp=gpar(fontsize=7, col =
                                                                      ifelse(cooccurmat[i,j]>10,
                                                                                    "white",
                                                                                    "black")))
                                          },
                                          row_title_gp = gpar(fontsize=7, fontface="bold"),
                                          column_title_gp = gpar(fontsize=7, fontface="bold"),
                                          row_split = sympnames_type_df$symptom_type[1:26],
                                          column_split = sympnames_type_df$symptom_type[1:26],
                                          row_names_gp = gpar(fontsize=9),
                                          column_names_gp = gpar(fontsize=9),
                                          top_annotation = columnAnnotation(
                                            `%`=anno_barplot(marginal_counts, gp=gpar(fill=myCols[[1]]),
                                                             border=F,height = unit(1.7,"cm"),
                                                             ylim = c(0,50),
                                                             # axis_param = list(direction = "reverse")
                                            )),
                                          # column_title = mymatrixname,
                                          row_gap = unit(1.5,"mm"),
                                          column_gap = unit(1.5,"mm"),
                                          rect_gp = gpar(col = "white", lwd = 1),
                                          col = circlize::colorRamp2(breaks = c(00,15,40),
                                                                     colors = c("white",
                                                                                myCols[[2]],
                                                                                myCols[[1]])))
myheatmap_cooccur_first


png(filename = paste0(figpath,"heatmap_cooccurrence_first_symptom.png"),
    width = 10,height = 9,units = "in",res = 300)
myheatmap_cooccur_first
dev.off()



# Relationship between Ct Values and time since symtom onset --------------
dfRes$variant_inferred_detail %>% unique()

table(dfRes$days_since_symptom_onset)
ct_plot=dfRes %>% filter(estbinres==1, !is.na(variant_inferred_detail), 
                         variant_inferred_detail %in% c("BA.1 (Omicron)", "BA.2 (Omicron)")) %>% 
  mutate(days_since_symptom_onset=case_when(days_since_symptom_onset==12 ~ NA_real_,
                           T ~days_since_symptom_onset)) %>% 
  ggplot(aes(x=days_since_symptom_onset,y=ct1)) +
  geom_jitter(alpha=0.1, size=0.5) +
  OverReact::theme_react(strip_text_size = 8) +
  facet_wrap(.~variant_inferred_detail) +
  geom_smooth(method="lm", col = "firebrick2", size=0.5) +
  ggpubr::stat_regline_equation(size=2,label.y = 43) +
  ggpubr::stat_cor(size=2,label.y=47) +
  scale_x_continuous(breaks = 0:11) +
  labs(x="Days since symptom onset", y="N-gene Ct value") 


### PLot day distributions
# get means
means_days=dfRes %>% filter(estbinres==1,  variant_inferred_detail %in% c("BA.1 (Omicron)", "BA.2 (Omicron)")) %>% 
  filter(round %in% c(2:7,8:10,13:15,16,17,18,19)) %>% 
  # mutate(days_since_symptom_onset=case_when(days_since_symptom_onset==12 ~ NA_real_,
  #                                           T ~days_since_symptom_onset)) %>% 
  group_by(variant_inferred_detail) %>% 
  summarise(mean_days=mean(days_since_symptom_onset,na.rm=T))
means_days

ct_plot_daydist=dfRes %>% filter(estbinres==1, !is.na(variant_inferred_detail), 
                         variant_inferred_detail %in% c("BA.1 (Omicron)", "BA.2 (Omicron)")) %>% 
  mutate(days_since_symptom_onset=case_when(days_since_symptom_onset==12 ~ NA_real_,
                                            T ~days_since_symptom_onset)) %>% 
  ggplot(aes(x=days_since_symptom_onset,fill=variant_inferred_detail)) +
  geom_density(outline.type = "full",col="white", alpha=0.5) +
  theme_void()+
  facet_wrap(.~variant_inferred_detail) +
  scale_x_continuous(breaks = 0:11) +
  scale_fill_manual(values=c( "#EC7300","#960078")) +
  labs(x="Days since symptom onset")  +
  geom_vline(data = data.frame(variant_inferred_detail="BA.1 (Omicron)", xint=means_days$mean_days[[1]]),
             aes(xintercept=xint),size=.1,
            linetype="dashed", col="white")+
  geom_text(data = data.frame(variant_inferred_detail="BA.1 (Omicron)", xint=means_days$mean_days[[1]]),
            aes(x=xint+1.05,y=0.01),label=paste0("Mean: ", round(means_days$mean_days[[1]],2)," days"),
            size=2,col="white") +
  geom_vline(data = data.frame(variant_inferred_detail="BA.2 (Omicron)", xint=means_days$mean_days[[2]]),
             aes(xintercept=xint),size=.1,
             linetype="dashed", col="white")+
  geom_text(data = data.frame(variant_inferred_detail="BA.2 (Omicron)", xint=means_days$mean_days[[2]]),
            aes(x=xint+1.05,y=0.01),label=paste0("Mean: ", round(means_days$mean_days[[2]],2)," days"), 
            size=2,col="white") +
  theme(strip.text = element_blank(),
        legend.position = "none")
ct_plot_daydist

ct_plot_ctdist=dfRes %>% filter(estbinres==1, !is.na(variant_inferred_detail), 
                                 variant_inferred_detail %in% c("BA.1 (Omicron)", "BA.2 (Omicron)")) %>% 
  mutate(days_since_symptom_onset=case_when(days_since_symptom_onset==12 ~ NA_real_,
                                            T ~days_since_symptom_onset)) %>% 
  ggplot(aes(x=ct1,fill=variant_inferred_detail)) +
  geom_density(outline.type = "full", alpha=0.5, col="white") +
  theme_void()+
  # facet_wrap(.~variant_inferred_detail) +
  coord_flip() +
  scale_fill_manual(values=c( "#EC7300","#960078")) +
  scale_x_continuous(breaks = 0:11) +
  labs(x="Days since symptom onset",fill="")  +
  theme(strip.text = element_blank(),
        legend.position = "none")
ct_plot_ctdist
ct_plot_daydist



layout="
AAAAAAAAAB
CCCCCCCCCD
CCCCCCCCCD
CCCCCCCCCD
CCCCCCCCCD
"



plots <- list(ct_plot_daydist, plot_spacer(),ct_plot,ct_plot_ctdist)
myplot=patchwork::wrap_plots(plots, design = layout)

OverReact::saveREACTplot(p = myplot,figpath = figpath,
                         filename = "ct_vs_days_since_sx_onset_scatter",
                         width = 8,height = 5)

