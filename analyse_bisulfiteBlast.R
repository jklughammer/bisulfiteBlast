#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
theme_set(theme_bw())
options(stringsAsFActors=FALSE)

args=commandArgs(trailingOnly = TRUE)
readno=args[2]
db=args[3] 
setwd(args[1])

#Annotation table containing all sample annotations in your project. Must contain columsn "Abbreviateion" for species abbreviation/name used in RefFreeDMA, 
#"English" for English common name, and "scientific_name" for the scientific name of the species. 
annot=fread("<your path>/meta/Patholist_selection.tsv")
annot_red=annot[,list(common_name_query=English[1],scientific_name_query=scientific_name[1]),by="Abbreviation"]


# blast db taxa annotation needs to be created first (see createDB)
taxa=fread(paste0("<your path>/blastDB/",db,"_extracted_taxa.tsv"),sep="\t") 
setnames(taxa,names(taxa),c("id","taxid","length","common_name_query","species_complete","kingdom"))

#read stats files for all species
stats_files=system(paste0("ls */*",db,"*/*unmapped-blast.txt"),intern=TRUE)

#set colors
diff_cols=c("#e41a1c","#1D81F3","#4daf4a","#984ea3","#ff7f00","#F5DC3F","#964B00","pink","#999999")


#now loop through all stats of all species
all_stats=data.table()
for(file in stats_files){
  sample=unlist(strsplit(basename(file),"_unmapped-blast.txt"))[1]
  dir=dirname(file)
  rrbs_species=dirname(dirname(file))
  if(file.info(file)$size == 0){ 
    print(paste0("NO HITS: ", sample))
    next }
  print(sample)
  stats=fread(file,sep="\t",drop=c(13))
  setnames(stats,names(stats),c("query","target","species_complete","common_name","qlen","slen","sstart", "send",  "pident", "length", "evalue", "bitscore"))
  stats[,Abbreviation:=sample,]
  stats[,dir:=dir,]
  stats[,rrbs_species:=rrbs_species,]
  all_stats=rbindlist(list(all_stats,stats))  
  
}

all_stats=merge(all_stats,annot_red,by="Abbreviation",all.x=TRUE)
all_stats[,start_round:=round(sstart/100)*100,]
all_stats[,species:=unlist(lapply(species_complete,function(x){paste0(unlist(strsplit(x," "))[1:2],collapse=" ")})),]
topSel=8

taxa_sub=taxa[species_complete%in%unique(all_stats$species_complete)]
taxa_sub[,species:=unlist(lapply(species_complete,function(x){paste0(unlist(strsplit(x," "))[1:2],collapse=" ")})),]
taxa_sub[,amb_count:=length(unique(kingdom)),by="species"]
taxa_sub[amb_count>1,species:=species_complete,]
seqs_per_taxon=taxa_sub[,list(sequences=.N,bases=sum(length)),by=c("species","kingdom")]


all_species=unique(all_stats$rrbs_species)
pl=list()
suffix="filt"
for (sel_species in all_species){
  print(sel_species)
  sel_stats=all_stats[rrbs_species==sel_species]
  
  out_dir=sel_stats$dir[1]
  title=paste0(sel_stats$scientific_name_query[1]," (",sel_stats$common_name_query[1],")")
  species_count=merge(sel_stats[,list(N=.N,detected_seqs=length(unique(target))),by="species"],seqs_per_taxon,by="species",all.x=TRUE)[order(N,decreasing=TRUE)]
  species_count[,plotSpecies:=paste0(species," ",detected_seqs,"/",sequences," (",round(bases/1000,1),"kb)"),]
  sel_stats=merge(species_count,sel_stats,by="species",all=TRUE)
  sel_stats[,plotSpecies:=ifelse(species%in%species_count[1:topSel]$species,plotSpecies,"other")]
  sel_stats[,plotSpecies:=factor(plotSpecies,levels=c(unique(grep("^other$",plotSpecies[order(N,decreasing=TRUE)],invert=TRUE,value=TRUE)),"other"))]
  height=0.2*length(unique(sel_stats$Abbreviation))+2
  pdf(paste0(out_dir,"/",sel_species,"_contamination_species_",db,"_",suffix,".pdf"),height=height,width=12)
  ggpl1=ggplot(sel_stats,aes(x=Abbreviation,fill=plotSpecies))+geom_bar()+coord_flip()+scale_fill_manual(values=c(diff_cols))+ylim(0,as.numeric(readno))+ggtitle(label=title)+ylab("")
  print(ggpl1)
  dev.off()
  pl[[paste0(sel_species,"_species")]]=ggpl1+theme(plot.margin=unit(c(3.8-height/2,1,3.8-height/2,1),"inches") )
  
  target_count=sel_stats[,.N,by=c("target","sstart","species")][order(N,decreasing=TRUE)]
  target_count[,rank:=rank(-N,ties.method="random"),by="species"]
  target_count[,top_target:=ifelse(rank<=topSel,target,"other"),]
  target_count[,top_rank:=ifelse(rank<=topSel,rank,"other"),]
  target_count_red=target_count[,list(count=sum(N),target_len=max(sstart)),by=c("top_target","top_rank","species")]
  target_count_red[,sum:=sum(count),by="species"]
  
  target_count_red[,annot:=ifelse(top_rank=="1",paste0(top_target,round(target_len/1000,1),"kb"),NA),]
  target_count_red[,species:=factor(species,levels=unique(species[order(sum)])),]
  target_count_red=target_count_red[order(top_rank)]
  
  pdf(paste0(out_dir,"/",sel_species,"_contamination_target-distrib_",db,"_",suffix,".pdf"),height=4,width=8)
  ggpl2=ggplot(target_count_red[species%in%species_count[1:topSel]$species],aes(x=species,y=count,fill=as.factor(top_rank)))+geom_bar(stat="identity",position="fill")+geom_text(aes(label=annot),y=0.005,hjust=0)+coord_flip()+scale_fill_manual(values=diff_cols,guide = guide_legend(title = "target"))+ylab("target frequency")+ggtitle(label=title)
  print(ggpl2)
  dev.off()
  pl[[paste0(sel_species,"_distrib")]]=ggpl2+theme(plot.margin=unit(c(1.5,1,1.5,1),"inches") )
  write.table(target_count_red,paste0(out_dir,"/",sel_species,"_contamination_target_",suffix,".tsv"),sep="\t",quote=FALSE,row.names=FALSE)
}

pdf(paste0("../results_analysis/basicStats/bisulfiteBlast_contamination_",db,"_",suffix,".pdf"),height=7,width=12)
print(pl)
dev.off()
