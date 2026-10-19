args<-commandArgs(trailingOnly=TRUE);root<-args[1];out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(microeco));suppressPackageStartupMessages(library(igraph))
stopifnot(as.character(packageVersion('microeco'))=='2.0.0')
nodes<-letters[1:12]
edges<-rbind(c('a','b'),c('b','c'),c('c','a'),c('c','d'),c('e','f'),c('f','g'),c('g','h'),c('i','j'),c('j','k'),c('k','i'),c('a','e'),c('a','i'),c('b','f'),c('c','g'),c('d','h'),c('d','i'),c('e','j'))
g<-graph_from_data_frame(as.data.frame(edges),directed=FALSE,vertices=nodes);V(g)$module<-c(rep('M1',4),rep('M2',4),rep('M3',3),'M4');E(g)$weight<-1
# Existing package role code on fixed modules: no community optimization, no plots.
otu<-matrix(1,12,3,dimnames=list(nodes,c('s1','s2','s3')))
dataset<-microtable$new(otu_table=as.data.frame(otu),auto_tidy=FALSE)
t1<-trans_network$new(dataset=dataset,cor_method=NULL);t1$res_network<-g
private<-t1$.__enclos_env__$private
wmd<-private$within_module_degree(g);z<-private$zscore(wmd)
for(mod in unique(wmd$module)){idx<-wmd$module==mod;d<-wmd$mod_links[idx];expected<-(d-mean(d))/sqrt(sum((d-mean(d))^2)/(length(d)-1));stopifnot(isTRUE(all.equal(z$z[idx],expected,tolerance=1e-15)))}
amc<-private$among_module_connectivity(g);td<-data.frame(taxa=nodes,total_links=degree(g));pi<-private$participation_coeffiecient(amc,td)
for(n in nodes){nb<-as_ids(neighbors(g,n));exp<-if(!length(nb))NaN else 1-sum((table(V(g)$module[match(nb,nodes)])/length(nb))^2);stopifnot(isTRUE(all.equal(unname(pi[n,'p']),unname(exp),tolerance=1e-15)))}
roles<-private$assign_module_roles(data.frame(z=c(2.5,2.5,2.50001,2.50001,NA),p=c(.62,.62001,.62,.62001,.7),row.names=c('boundary','connector','module','network','invalid')))
stopifnot(identical(roles$taxa_roles,c('Peripheral nodes','Connectors','Module hubs','Network hubs')),!('invalid'%in%rownames(roles)))
t1$get_node_table(node_roles=TRUE);tab<-t1$res_node_table
stopifnot(setequal(tab$name,nodes),all(is.na(tab$z[tab$name%in%c('i','j','k','l')])),all(is.na(tab$taxa_roles[tab$name%in%c('i','j','k','l')])))
write.table(tab,file.path(out,'ZIPI_ACTUAL_NODE_TABLE.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(cbind(name=rownames(roles),roles),file.path(out,'ZIPI_ROLE_BOUNDARIES.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS Zi/Pi actual microeco2.0.0 fixed-module formulas; sample SD; <=2.5 and <=0.62; NA preserved; no optimization/plot\n')
