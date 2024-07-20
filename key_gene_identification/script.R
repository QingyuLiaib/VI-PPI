library(readr)
library(stringr)
library(igraph)
setwd('')

pneumonia_ctd <- read_csv("./pneumonia_CTD.csv")
virus_human <- read_delim("./virus_human_HVIDB.csv",
                          delim = "\t", escape_double = FALSE,
                          trim_ws = TRUE)
a = str_split(virus_human$Human_GeneName,';')
re = c()
for(i in 1:length(a)){
  re = c(re,a[[i]][1])
}
virus_human$genename = re
###############################################################################
infla_net <- readRDS("./infla_net.rds")
infla_net.igraph = graph_from_data_frame(infla_net,directed = F)
IAV = virus_human[virus_human$Short=='H1N1',]
result = data.frame()
for(i in 1:length(unique(IAV$Organism_virus))){
  name = unique(IAV$Organism_virus)[i]
  tmp = IAV[IAV$Organism_virus%in%name,]
  result = rbind(result,data.frame(names = name,genes = unique(tmp$genename)))
}

sta = data.frame(table(result$genes))
IAVgene = sta[sta$Freq>3,]$Var1
############################################################
dis_gene = pneumonia_ctd[pneumonia_ctd$`Direct Evidence`=='marker/mechanism',]$`Gene Symbol`
dis_gene = dis_gene[!is.na(dis_gene)]

list1 = intersect(IAVgene,c(infla_net$node1,infla_net$node2))
list2 = intersect(dis_gene, c(infla_net$node1,infla_net$node2))

all  = unique(c(infla_net$node1,infla_net$node2))

betw=list()
bet_cal = function(list1,list2,net){
  genes = unique(c(list1,list2))
  infla_net.igraph = graph_from_data_frame(net,directed = F)
  neigh = c()
  for(i in 1:length(list1)){
    tmp = all_shortest_paths(graph = infla_net.igraph,from = list1[i],to = list2)
    neigh = c(neigh,unique(names(unlist(tmp$res))))
  }
  neigh = unique(neigh)
  betw = data.frame(row.names = neigh)
  betw$values = 0
  for(j in 1:length(list1)){
    print(j)
    for(k in 1:length(list2)){
      tmp = all_shortest_paths(graph = infla_net.igraph,from = list1[j],to = list2[k])
      count = table(names(unlist(tmp$res)))
      count = count[!names(count) %in% c(list1[j],list2[k])]/length(tmp$res)
      betw[names(count),]=betw[names(count),]+count
    }
  }
  betw$ratio = betw$values/(length(list1)*length(list2))
  betw = betw[order(betw$values,decreasing = T),]
  return(betw)
}

sample_nodelist = function(source_net_node,background){
  #source_net: data.frame(n*2)_colnames = c(node1,node2)
  #background: data.frame(n*2)_colnames = c(node1,node2)
  node_list = c()
  source_net_node = source_net_node
  background_degree=as.data.frame(table(c(background$node1,background$node2)))
  source_net_degree_count = background_degree[background_degree$Var1%in%source_net_node,]
  if(nrow(source_net_degree_count)==0){
    return(0)
  }
  source_net_degree_count = as.data.frame(table(source_net_degree_count$Freq))
  for(i in 1:nrow(source_net_degree_count)){
    deg = as.numeric(as.character(source_net_degree_count[i,1]))
    count = as.numeric(as.character(source_net_degree_count[i,2]))
    tmp = background_degree[background_degree[,2]==deg,]

    node_sample = tmp[sample(1:nrow(tmp),count),1]
    node_sample = as.character(node_sample)
    node_list = c(node_list,node_sample)
  }
  return(node_list)
}

for(j in 1:length(all)){
  print(j/length(all)*100)
  tmp1 = all_shortest_paths(graph = infla_net.igraph,from =all[j],to = list1)
  neig = unique(names(unlist(tmp1$res)))
  ta = data.frame(row.names = neig)
  ta$value=0
  for ( i in 1:length(list1)){
    tmp = all_shortest_paths(graph = infla_net.igraph,from =list1[i],to = all[j])
    count = table(names(unlist(tmp$res)))
    count = count[!names(count) %in% c(list1[i],all[j])]/length(tmp$res)
    ta[names(count),] = ta[names(count),]+count
  }
  betw[[all[j]]] =ta
}
saveRDS(betw,'./betw.rds')


nodes = unique(c(infla_net$node1,infla_net$node2))
con = data.frame(row.names = nodes)
for ( t in 1:10000){
  print(t)
  sam = sample_nodelist(list2,infla_net)
  alldata = data.frame(row.names = nodes)
  alldata$value = 0
  for(i in 1:length(sam)){
    tmp = betw[[sam[i]]]
    alldata[rownames(tmp),]=alldata[rownames(tmp),]+tmp$value
  }
  con = cbind(con,alldata)
}
saveRDS(con,'./IAV_con_sample.rds')

IAV_betw = bet_cal(list1,list2,infla_net)
keygene =data.frame()
for (i in 1: nrow(IAV_betw)){
  gene = rownames(IAV_betw)[i]
  distri = con[gene,]
  q3 = as.numeric(quantile(as.numeric(distri),probs = c(0.75)))
  q1 = as.numeric(quantile(as.numeric(distri),probs = c(0.25)))
  maxbet = q3+3*(q3-q1)
  print(maxbet)
  if(IAV_betw[i,'values']>maxbet){
    print(gene)
    keygene = rbind(keygene,
                    data.frame(gene = gene,betw = IAV_betw[i,'values']))
  }
}
write_csv(keygene,'./IAVKEYgene.csv')
