# Test with future package framework
library(tidyverse)
library(readr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(valr)
library(future)
library(future.apply)
options(scipen = 999999999)
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set
#-----------------------------------------
#Utils. Fn
data_in_fn<-function(file){
  dat.out <- get(load(file))
  tmp_obj <- names(mget(load(file)))
  rm(list = tmp_obj)
  rm(tmp_obj)
  return(dat.out)
}
#-----------------------------------------
#Trivial example to illustrate the pattern for future package parallel computing

test_dat_file<-"./data/cluster/raw/chr1_spec_res.Rda"

chr_spec_res<-data_in_fn(test_dat_file)

tmp_set_l<-chr_spec_res$cl_member
# Set-up parallelisation (method and number of workers)
plan(multisession,workers=3)
# Run parallel computation, specifying the required packages
parallel_test_l<-future_lapply(1:100,future.seed = NULL,future.packages = c("dplyr"),function(x){
  check<-tmp_set_l[[sample(1:length(tmp_set_l),1)]]
  return(tibble(range=diff(range(as.numeric(check)))))
})
# Return to sequential computation to close the parallelisation
plan(sequential)
#-----------------------------------------
#Trial closer to current application with smallest chromosome
feature_file<-"./data/feature/feature_wrapped.Rda"
feature_Grange<-data_in_fn(feature_file)

hg19_coord <- read_delim("./data/annotation/hg19.genome", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")

fn_repo<-"./data/annotation/fn_BED/"
chr_set<-list.files(fn_repo)

chromo<-"chr22"

chr_feature_Grange<-feature_Grange[seqnames(feature_Grange)==chromo]
#Generate random CAGE coordinate
fn_folder<-paste0(fn_repo,chromo,"/")
fn_file<-grep('BED$',list.files(fn_folder),value = T)
fn_bed_l<-lapply(fn_file,function(f){
  read_bed(paste0(fn_folder,f),n_fields = 3)
})
names(fn_bed_l)<-fn_file

#Generate the random peak coordinates
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)  <- chromo
peakAnno <- annotatePeak(chr_feature_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)

rn_annotation<-sample(peakAnno@annoStat$Feature,size = length(chr_feature_Grange),prob = peakAnno@annoStat$Frequency/100,replace = T)
#check number of peaks from that category
n5<-length(grep("5'",as.character(rn_annotation)))
n3<-length(grep("3'",as.character(rn_annotation)))
nexon<-length(grep("Exon",as.character(rn_annotation)))
nintron<-length(grep("Intron",as.character(rn_annotation)))
n1kb<-length(grep("1kb",as.character(rn_annotation)))
n2kb<-length(grep("2kb",as.character(rn_annotation)))
n3kb<-length(grep("3kb",as.character(rn_annotation)))
ndown<-length(grep("Down",as.character(rn_annotation)))
ninter<-length(grep("Inter",as.character(rn_annotation)))
n_vec<-c(n3,n5,ndown,nexon,ninter,nintron,n1kb,n2kb,n3kb)
names(n_vec)<-fn_file
rm(ChIPseekerEnv)
rm(n5,n3,nexon,nintron,n1kb,n2kb,n3kb,ndown,ninter)
n_vec<-n_vec[n_vec>0]
#----------------------------------------------------
# Run parallel computation in Global environment
tmp_cage_tbl<-chr_feature_Grange %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)


rn_fn_coord_l<-vector('list',length(n_vec))
names(rn_fn_coord_l)<-names(n_vec)
for(f in names(n_vec)){
  message(f)
  tmp_n<-n_vec[f]
  
  if(f==fn_file[5]){
    plan(multisession,workers=5)
    
    rn_fn_coord_l[[f]]<-future_lapply(1:10,future.seed = NULL,future.packages = c("dplyr","valr"),function(x){
      rn_pol<-valr::bed_shuffle(x = tmp_cage_tbl%>%sample_n(tmp_n),genome = hg19_coord,excl = fn_bed_l[[f]],within = T,max_tries=1e6)
      return(rn_pol)
    })
    
    plan(sequential)

  }
  if(f!=fn_file[5]){
    plan(multisession,workers=5)
    rn_fn_coord_l[[f]]<-future_lapply(1:10,future.seed = NULL,future.packages = c("dplyr","valr"),function(x){
      rn_pol<-try(valr::bed_shuffle(x = tmp_cage_tbl%>%sample_n(tmp_n),genome = hg19_coord,incl = fn_bed_l[[f]],within=T,max_tries=1e6),silent=T)
      return(rn_pol)
    })
    plan(sequential)
    rn_fn_coord_l[[f]]<-rn_fn_coord_l[[f]][!(unlist(lapply(rn_fn_coord_l[[f]],function(x)any(class(x) %in% "try-error"))))]
  }
}

#Assemble these blocks into 10000 rnadom combo
plan(multisession,workers=5)

rn_peak_coord_tbl_l<-future_lapply(1:10000,future.packages = c("dplyr"),future.seed = NULL,function(x){
  return(do.call(bind_rows,lapply(rn_fn_coord_l,function(f)f[[sample(1:length(f),1)]])))
})
plan(sequential)

#-----------------------------------------
# Run parallel computation within function

produce_rn_coord_tbl_fn<-function(tmp_cage_tbl,n_vec,hg19_coord,fn_bed_l,nboot){
  fn_env <- environment()
  
  rn_fn_coord_l<-vector('list',length(n_vec))
  names(rn_fn_coord_l)<-names(n_vec)
  for(f in names(n_vec)){
    message(f)
    tmp_n<-n_vec[f]
    
    if(f==fn_file[5]){
      plan(multisession,workers=5)
      
      rn_fn_coord_l[[f]]<-future_lapply(1:10,future.seed = NULL,future.packages = c("dplyr","valr"),future.envir = fn_env,function(x){
        rn_pol<-valr::bed_shuffle(x = tmp_cage_tbl%>%sample_n(tmp_n),genome = hg19_coord,excl = fn_bed_l[[f]],within = T,max_tries=1e6)
        return(rn_pol)
      })
      
      plan(sequential)
      
    }
    if(f!=fn_file[5]){
      plan(multisession,workers=5)
      rn_fn_coord_l[[f]]<-future_lapply(1:10,future.seed = NULL,future.packages = c("dplyr","valr"),future.envir = fn_env,function(x){
        rn_pol<-try(valr::bed_shuffle(x = tmp_cage_tbl%>%sample_n(tmp_n),genome = hg19_coord,incl = fn_bed_l[[f]],within=T,max_tries=1e6),silent=T)
        return(rn_pol)
      })
      plan(sequential)
      rn_fn_coord_l[[f]]<-rn_fn_coord_l[[f]][!(unlist(lapply(rn_fn_coord_l[[f]],function(x)any(class(x) %in% "try-error"))))]
    }
  }
  
  #Assemble these blocks into 10000 rnadom combo
  plan(multisession,workers=5)
  
  rn_peak_coord_tbl_l<-future_lapply(1:nboot,future.packages = c("dplyr"),future.seed = NULL,future.envir = fn_env,function(x){
    return(do.call(bind_rows,lapply(rn_fn_coord_l,function(f)f[[sample(1:length(f),1)]])))
  })
  plan(sequential)
  
  return(rn_peak_coord_tbl_l)
  
}

rn_peak_coord_tbl_l<-produce_rn_coord_tbl_fn(tmp_cage_tbl,n_vec,hg19_coord,fn_bed_l,5)
