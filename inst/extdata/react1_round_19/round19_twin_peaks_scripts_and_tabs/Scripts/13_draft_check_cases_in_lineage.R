n_negatives=rep(NA,19)
names(n_negatives)=1:19
for (i in 8:19){
lineage=readRDS(paste0("E:/Group/saved_objects/rep",i,"_lineage.rds"))
print(table(lineage$estbinres))
n_negatives[i]=sum(lineage$estbinres==0)
}
