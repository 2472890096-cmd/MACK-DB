args<-commandArgs(trailingOnly=TRUE);root<-args[1];out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
core<-file.path(root,'pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R')
# Execute every original program expression except its CLI acquisition; supply exactly the production args.
e<-new.env(parent=globalenv());parsed<-parse(core)
fn<-function(x) is.call(x)&&identical(x[[1]],as.name('<-'))&&is.symbol(x[[2]])
for(x in parsed)if(fn(x)&&identical(x[[2]],as.name('parse_args')))eval(x,e)
mat<-rbind(c(1,0,0,0,0,0,0,0),c(1,1,0,0,0,0,0,0),c(2,0,3,0,4,0,1,0),c(1,2,3,4,5,6,0,0),c(2,2,2,2,2,2,2,0),rep(2,8),rep(0,8),c(0,0,1,3,0,0,0,0),c(0,2,0,3,0,1,0,1),c(1,2,3,2,1,3,2,1))
rownames(mat)<-0:(nrow(mat)-1);colnames(mat)<-paste0('S',1:8)
a<-file.path(out,'toy_abundance.matrix');m<-file.path(out,'toy_species.matrix')
write.table(cbind(OTU_id=rownames(mat),mat),a,sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(id=0:9,name=paste0('Synthetic_taxon_',0:9)),m,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
e$opts<-e$parse_args(c('-a',a,'-s',m,'-o','PRJTOYEPI','-d',out,'--mode','prevalence','--prevalence-min','0.25','--prevalence-max','0.75','--network-quantile','0.1'))
for(x in parsed)if(!(fn(x)&&identical(x[[2]],as.name('opts'))))eval(x,e)
near<-function(a,b)stopifnot(identical(is.na(a),is.na(b)),max(abs(a-b),na.rm=TRUE)<=1e-15)
# Independent scalar pairwise Bray-Curtis and discrete-ECDF oracle, without vegan or the production helper.
normalize<-function(x){for(j in seq_len(ncol(x))){s<-sum(x[,j]);x[,j]<-if(s>0)x[,j]/s else 0};x}
bc<-function(a,b){d<-sum(a+b);if(d>0)sum(abs(a-b))/d else 0}
ec<-function(x,p){x<-x[is.finite(x)];if(!length(x))return(NA_real_);u<-sort(unique(x));ok<-u[vapply(u,function(t)mean(x<=t)<=p,logical(1))];if(length(ok))max(ok) else min(u)}
d1<-d2<-q<-rep(NA_real_,nrow(mat));thresholds<-rep(NA_real_,nrow(mat))
for(i in seq_len(nrow(mat))){pres<-mat[i,]>0;if(mean(pres)<.25||mean(pres)>.75)next
 z<-normalize(mat[-i,,drop=FALSE]);b<-matrix(0,8,8)
 for(j in 1:7)for(k in (j+1):8)b[j,k]<-b[k,j]<-bc(z[,j],z[,k])
 d1[i]<-mean(b[pres,!pres]);d2[i]<-bc(rowMeans(z[,pres,drop=FALSE]),rowMeans(z[,!pres,drop=FALSE]))
 th<-ec(b[upper.tri(b)],.1);thresholds[i]<-th;A<-b<=th;diag(A)<-FALSE;degree<-rowSums(A);twom<-sum(A)
 q[i]<-if(twom==0)0 else sum((A-outer(degree,degree)/twom)*outer(pres,pres,'=='))/twom
}
near(e$D1,d1);near(e$D2,d2);near(e$Q,q)
stopifnot(all(is.finite(e$D1[c(2,3,4,8,9)])),all(is.na(e$D1[c(1,5,6,7,10)])))
for(x in list(c(.1,.2,.2,.2,.2,.2,.3,.4,.5,.6),c(rep(.2,4),rep(.8,6)),c(0,0,0,0,1,1),c(NA,Inf,.3,.7)))for(p in c(0,.1,.3,.5,1))near(e$ecdf_threshold_fast(x,p),ec(x,p))
stopifnot(is.na(e$ecdf_threshold_fast(numeric(),.1)),e$ecdf_threshold_fast(c(rep(.2,4),rep(.8,6)),.1)==.2)
stopifnot(!any(e$select_outliers(rep(0,8))),!any(e$select_outliers(c(NA,1))),identical(e$select_outliers(c(rep(0,9),10)),c(rep(FALSE,9),TRUE)))
stopifnot(e$normalize_columns(matrix(c(0,0,1,1),2))[1,1]==0,e$bray_two(c(0,0),c(0,0))==0)
write.table(data.frame(taxon=0:9,prevalence=rowMeans(mat>0),D1=e$D1,D2=e$D2,Q=e$Q,threshold=thresholds),file.path(out,'EPI_ACTUAL_VS_SCALAR_ORACLE.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PASS EPI v3: original whole-program expressions; D1/D2/Q scalar oracle <=1e-15; prevalence boundaries; ECDF ties/fallback; strict > and legal zero\n')
