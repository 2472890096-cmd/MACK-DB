#!/usr/bin/env Rscript
# Scores existing predictions only; evaluate unchanged scorer expressions in bounded chunks.
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L)
project<-args[1]; data_dir<-args[2]; result_dir<-args[3]; out<-args[4]; total_cases<-as.integer(args[5])
stopifnot(dir.exists(out))
gate<-file.path(out,'LOWMEM_REPRODUCTION.json')
stopifnot(file.exists(gate),any(grepl('"status": "PASS"',readLines(gate),fixed=TRUE)))
script_path<-sub('^--file=', '', commandArgs()[grepl('^--file=',commandArgs())][1])
# Rscript encodes spaces in its own --file argument; trailing scientific arguments are unchanged.
script_path<-gsub('~+~',' ',script_path,fixed=TRUE)
source_path<-normalizePath(file.path(dirname(script_path),'../calculate_structural_keystoneness_rk4_fast.R'),mustWork=TRUE)
source_text<-readLines(source_path,warn=FALSE)
parsed<-parse(text=source_text)
helpers<-c('read_matrix_csv','safe_norm_cols','safe_norm_rows','detect_index_base_and_adjust')
for(ex in parsed){
 if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&is.symbol(ex[[2]])&&as.character(ex[[2]])%in%helpers)eval(ex,envir=.GlobalEnv)
}
stopifnot(all(vapply(helpers,exists,logical(1))))
extract<-function(begin,end){
 i<-which(grepl(begin,source_text,fixed=TRUE));j<-which(grepl(end,source_text,fixed=TRUE));stopifnot(length(i)==1L,length(j)==1L,j>i)
 source_text[i:(j-1L)]
}
null_lines<-extract('    # Predicted null composition','    # Reference null composition')
bc_lines<-extract('    # Bray-Curtis predicted','    # Bray-Curtis reference')
summary_lines<-extract('      summary_df <- valid_df %>%','      write.csv(summary_df')
scalar_lines<-c(source_text[which(grepl('    original_abundance <- Ptrain[cbind(sp_idx, s_idx)]',source_text,fixed=TRUE))],source_text[which(grepl('    Ks_pred <- BC_pred * (1 - original_abundance)',source_text,fixed=TRUE))],source_text[which(grepl('    valid <- is.finite(Ks_pred)',source_text,fixed=TRUE))])
stopifnot(length(scalar_lines)==3L)
writeLines(c('# Unmodified expression excerpts; source:',paste0('# ',source_path),null_lines,bc_lines,scalar_lines,summary_lines),file.path(out,'EXACT_R_SOURCE_EXPRESSIONS.R.txt'))
null_expr<-parse(text=null_lines);bc_expr<-parse(text=bc_lines);scalar_expr<-parse(text=scalar_lines);sum_expr<-parse(text=summary_lines)
Ptrain<-safe_norm_cols(read_matrix_csv(file.path(data_dir,'Ptrain.csv')))$mat
qtrn<-safe_norm_rows(read_matrix_csv(file.path(result_dir,'qtrn.csv')))$mat
species_id<-scan(file.path(data_dir,'species_id.csv'),what=integer(),quiet=TRUE,sep=',')
sample_id<-scan(file.path(data_dir,'sample_id.csv'),what=integer(),quiet=TRUE,sep=',')
sp_all<-detect_index_base_and_adjust(species_id,nrow(Ptrain),'species_id')$idx
sa_all<-detect_index_base_and_adjust(sample_id,nrow(qtrn),'sample_id')$idx
stopifnot(length(sp_all)==total_cases,length(sa_all)==total_cases,ncol(qtrn)==nrow(Ptrain),nrow(qtrn)==ncol(Ptrain))
ks_all<-bc_all<-rep(NA_real_,total_cases);qn_valid_all<-qt_valid_all<-rep(FALSE,total_cases)
qn_sum_all<-rep(NA_real_,total_cases)
con<-file(file.path(result_dir,'qtst.csv'),open='r')
count<-0L;bad_input<-negative_input<-0;started<-proc.time();chunk_rows<-256L
name_map<-read.delim(file.path(out,'taxon_name_map.tsv'),check.names=FALSE)
focal_ids<-name_map$species_id[name_map$name=='Vescimonas_sp934591815']
probe_cases<-unique(c(1L,total_cases,as.integer(seq(1,total_cases,length.out=10)),which(sp_all%in%focal_ids)[c(1,2)]))
probe_cases<-probe_cases[is.finite(probe_cases)&probe_cases>=1&probe_cases<=total_cases]
probe_q<-list();probe_info<-list()
repeat{
 lines<-readLines(con,n=chunk_rows,warn=FALSE)
 if(!length(lines))break
 qtst_raw<-as.matrix(data.table::fread(text=paste(lines,collapse='\n'),header=FALSE,data.table=FALSE,showProgress=FALSE))
 storage.mode(qtst_raw)<-'numeric'
 stopifnot(ncol(qtst_raw)==nrow(Ptrain),nrow(qtst_raw)==length(lines),count+nrow(qtst_raw)<=total_cases)
 bad_input<-bad_input+sum(!is.finite(qtst_raw));negative_input<-negative_input+sum(qtst_raw<0,na.rm=TRUE)
 qtst<-safe_norm_rows(qtst_raw)$mat
 N_exp<-nrow(qtst);global_rows<-count+seq_len(N_exp)
 sp_idx<-sp_all[global_rows];s_idx<-sa_all[global_rows];idx_mat<-cbind(seq_len(N_exp),sp_idx)
 eval(null_expr);qtst_valid<-rowSums(qtst)>1e-12;eval(bc_expr);eval(scalar_expr)
 ks_all[global_rows]<-Ks_pred;bc_all[global_rows]<-BC_pred
 qn_valid_all[global_rows]<-q_null_valid;qt_valid_all[global_rows]<-qtst_valid;qn_sum_all[global_rows]<-q_null_rs
 for(i in which(global_rows%in%probe_cases)){
  probe_q[[length(probe_q)+1L]]<-qtst_raw[i,]
  probe_info[[length(probe_info)+1L]]<-data.frame(case_id=global_rows[i],species_id=sp_idx[i],sample_id=s_idx[i],original_abundance=original_abundance[i],BC_pred=BC_pred[i],Ks_pred=Ks_pred[i],valid=valid[i])
 }
 count<-count+N_exp
 if(count%%(256L*200L)==0L)cat(project,'R score rows',count,'/',total_cases,'\n')
}
close(con);stopifnot(count==total_cases)
detail<-data.frame(species=sp_all,sample=sa_all,BC_pred=bc_all,Ks_pred=ks_all,valid=is.finite(ks_all))
valid_df<-detail %>% filter(valid,is.finite(Ks_pred))
stopifnot(nrow(valid_df)>0)
eval(sum_expr)
write_precise<-function(x,path){
 for(k in names(x))if(is.double(x[[k]]))x[[k]]<-format(x[[k]],digits=17,scientific=TRUE,trim=TRUE)
 write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
}
write_precise(summary_df,file.path(out,'R_taxon_summary.tsv'))
saveRDS(list(species=sp_all,sample=sa_all,BC_pred=bc_all,Ks_pred=ks_all,q_null_valid=qn_valid_all,qtst_valid=qt_valid_all,q_null_sum_before_renormalization=qn_sum_all),file.path(out,'R_case_scores.rds'),compress=FALSE)
# Small exact-source spot-check payload for an independent scalar implementation.
write.table(do.call(rbind,probe_q),file.path(out,'R_probe_qtst.csv'),sep=',',quote=FALSE,row.names=FALSE,col.names=FALSE)
write_precise(do.call(rbind,probe_info),file.path(out,'R_probe_scores.tsv'))
write.table(data.frame(field=c('project','status','rows','valid','invalid_q_null','invalid_qtst','nonfinite_input_cells','negative_input_cells','taxa_with_valid_scores','elapsed_seconds','R_version','dplyr_version','data_table_version'),value=c(project,'PASS',count,sum(is.finite(ks_all)),sum(!qn_valid_all),sum(!qt_valid_all),bad_input,negative_input,nrow(summary_df),(proc.time()-started)[['elapsed']],R.version.string,as.character(packageVersion('dplyr')),as.character(packageVersion('data.table')))),file.path(out,'R_SCORING_STATUS.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat(project,'R SCORING PASS',count,'rows',nrow(summary_df),'taxa\n')
