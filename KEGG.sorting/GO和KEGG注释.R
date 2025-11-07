library(stringr)
library(dplyr)
library(clusterProfiler)
egg <- read.table("eggnog.emapper.annotations", sep="\t", header=T, quote = "")
colnames(egg)

koterms <- egg %>%
           dplyr::select(GID = query, KO=KEGG_ko) %>% na.omit() %>% 
           filter(str_detect(KO, "ko")) #根据egg第一列和KEGG_ko列的标题提取基因的KEGG注释。
head(koterms)

load("kegg_info.RData")  #读取kegg_info.RData文件
head(ko2pathway)
head(pathway2name)

# 获得gene2pathway
library(stringr)
# 把ko2pathway的列名改为KO和Pathway，与koterms一致
colnames(ko2pathway)=c('Ko','Pathway')
head(ko2pathway)
#把koterms的KO值的前缀ko:去掉，与ko2pathway格式一致
koterms$KO <- str_replace_all(koterms$KO, "ko:", "")
head(koterms)
colnames(koterms)=c('GID', 'Ko')
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

#合并koterms和ko2pathway到gene2pathway，将基因与pathway的对应关系整理出来
gene2pathway <- koterms %>% left_join(ko2pathway, by = "Ko") %>% 
  dplyr::select(GID, Pathway) %>% na.omit() 




