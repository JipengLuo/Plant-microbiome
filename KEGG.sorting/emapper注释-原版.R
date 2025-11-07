####################################构建OrgDb
#install.packages('AnnotationForge')
library(dplyr)
library(stringr)
library(jsonlite)
library(AnnotationForge)
#顺手设置一下options
options(stringsAsFactors = F)

#读取生成的annotations文件
emapper <- read.table("eggnog.emapper.annotations", header = T, sep = "\t",quote = "")
emapper[emapper==""] <- NA

#提取GO信息;提取query列、Preferred name、GOs列，其中%>%为管道符
# library(dplyr)
gene_info <- emapper %>% dplyr::select(GID=query, GENENAME=Preferred_name) %>% na.omit()
gos <- emapper %>% dplyr::select(query, GOs) %>% na.omit()

#构建一个空的【gene2go】数据框，为后面填充数据用
gene2go = data.frame(GID = character(),
                     GO = character(),
                     EVIDENCE = character())

# library(stringr)
gos_list <- function(x){
  the_gos <- str_split(x[2], ",", simplify = FALSE)[[1]]
  df_temp <- data.frame(GID = rep(x[1], length(the_gos)),
                        GO = the_gos,
                        EVIDENCE = rep("IEA", length(the_gos)))
  return(df_temp)
}

gene2gol <- apply(as.matrix(gos), 1, gos_list)
gene2gol_df <- do.call(rbind.data.frame, gene2gol)
gene2go <- gene2gol_df
gene2go$GO[gene2go$GO=="-"] <- NA
gene2go <- na.omit(gene2go)


####提取KEGG信息
gene2ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko)
gene2ko$Ko[gene2ko$Ko=="-"]<-NA
gene2ko<-na.omit(gene2ko)
gene2kol <- apply(as.matrix(gene2ko), 1, gos_list)
gene2kol_df <- do.call(rbind.data.frame, gene2kol)
gene2ko <- gene2kol_df[,1:2]
colnames(gene2ko) <- c("GID","Ko")
gene2ko$Ko <- gsub("ko:", "", gene2ko$Ko)


# 下载KO的json文件，直接运行即可
if(T){
  library(jsonlite)
  # 下面的json = "ko00001.json"，如果你下载到其他地方，记得加上路径
  update_kegg <- function(json = "ko00001.json") {
    pathway2name <- tibble(Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())
    kegg <- fromJSON(json)
    for (a in seq_along(kegg[["children"]][["children"]])) {
      A <- kegg[["children"]][["name"]][[a]]
      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          kos <- str_match(kos_info, "K[0-9]*")[,1]
          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))}}}
    save(pathway2name, ko2pathway, file = "kegg_info.RData")}
}

# 调用函数后在本地创建kegg_info.RData文件，以后只需要载入 "kegg_info.RData"即可
update_kegg()
# 载入kegg_info.RData文件
load(file = "kegg_info.RData")

# 提取 pathway 信息
gene2pathway <- gene2ko %>% 
  left_join(ko2pathway, by = "Ko") %>% 
  dplyr::select(GID, Pathway) %>% na.omit()


# 烟草
tax_id = '4097'
genus = 'Nicotiana'
species = 'tabacum'

# 可以去重，以防万一，或者等下面报错提示去重再去跑代码也行。
# 原理是，makeOrgPackage不允许有重复的行，因此需要删除，
# 除了gene2pathway可能有重复的行，其他几乎不可能有重复的情况。

if(T){
  gene2go <- unique(gene2go)
  gene2go <- gene2go[!duplicated(gene2go),]
  gene2ko <- gene2ko[!duplicated(gene2ko),]
  gene2pathway <- gene2pathway[!duplicated(gene2pathway),]
  gene_info <- gene_info[!duplicated(gene_info),]
}

gene2go <- unique(gene2go)
gene2go <- gene2go[!duplicated(gene2go),]

########## 构建OrgDb ##########
# library(AnnotationForge)
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               pathway=gene2pathway,
               version="1.38.0",  #版本，使用?makeOrgPackage，拉到最下面查看
               maintainer = "luojp <luojp900510@gmail.com>",  #修改为你的名字和邮箱
               author = "luojp <luojp900510@gmail.com>",  #修改为你的名字和邮箱
               outputDir = "./",  #输出文件位置
               tax_id=tax_id,  #你在NCBI上查并定义的id
               genus=genus,  #你在NCBI上查并定义的属名
               species=species,  #你在NCBI上查并定义的种名
               goTable="go")
# 构建完成，生成一个 eg.db 的文件夹

########## 导入建OrgDb ##########
install.packages("./org.Ntabacum.eg.db", repos=NULL, type="source")
library(org.Ntabacum.eg.db)
# 查看所有列的信息
columns(org.Ntabacum.eg.db)
# 查看所有基因
keys(org.Ntabacum.eg.db)
# 查看特定基因的信息
# library(dplyr)
select(org.Ntabacum.eg.db, keys = "Nitab4.5_0000621g0130.1", 
       columns = c("GO"))

########## 导出背景文件 ##########

# 因为我们在后面进行kegg分析的时候需要背景文件
# 因此为了方便即可导出pathway2name和pathway2gene

# 做一些常规格式化
pathway2name$Name <- gsub(" \\[BR:ko[0-9]{5}\\]", "", pathway2name$Name)
pathway2name<- na.omit(pathway2name)
pathway2gene <-gene2pathway[, c("Pathway","GID")]
# 输出
write.table(pathway2name, file = "./pathway2name", sep = "\t", quote = F, row.names = F)
write.table(pathway2gene, file = "./pathway2gene", sep = "\t", quote = F, row.names = F)


########## clusterProfiler 富集分析 ##########
# KEGG（不需要OrgDb）
# library(clusterProfiler)
# 每次只需导入下面两个文件即可（同一物种）
pathway2gene <- read.table("./pathway2gene",header = T,sep = "\t")
pathway2name <- read.table("./pathway2name",header = T,sep = "\t")

# KEGG pathway 富集
ekp_up <- enricher(gene_up, 
                   TERM2GENE = pathway2gene, 
                   TERM2NAME = pathway2name, 
                   pvalueCutoff = 1,  # 表示全部保留，可以设为0.05作为阈值
                   qvalueCutoff = 1, # 表示全部保留，可以设为0.05作为阈值
                   pAdjustMethod = "BH",
                   minGSSize = 1)
ekp_down <- enricher(gene_down, 
                     TERM2GENE = pathway2gene, 
                     TERM2NAME = pathway2name, 
                     pvalueCutoff = 1,  # 表示全部保留，可以设为0.05作为阈值
                     qvalueCutoff = 1, # 表示全部保留，可以设为0.05作为阈值
                     pAdjustMethod = "BH",
                     minGSSize = 1)
p1 <- dotplot(ekp_up)
p2 <- dotplot(ekp_down)
p1+p2


# GO富集（需要OrgDb）
ego_up <- enrichGO(gene=gene_up,
                   OrgDb=org.Ntabacum.eg.db,
                   keyType="GID",
                   ont="ALL",   #CC/BP/MF可选
                   qvalueCutoff = 1,
                   pvalueCutoff = 1)
ego_down <- enrichGO(gene=gene_down,
                     OrgDb=org.Ntabacum.eg.db,
                     keyType="GID",
                     ont="ALL",   #CC/BP/MF可选
                     qvalueCutoff = 1,
                     pvalueCutoff = 1)
p3 <- dotplot(ego_up)
p4 <- dotplot(ego_down)
p3+p4


# 无论是ekp还是ego都可以选择将数据导出，然后可以自己使用ggplot画
ego_results <- as.data.frame(ego)
ekp_result <- as.data.frame(ekp)
write.table(XXX_results, file = "XXX.txt", quote = F)






