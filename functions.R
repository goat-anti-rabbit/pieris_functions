# Functions

#########################
### THE SCAFFOLD FUNCTION
#########################


my_scaffold_fun <- function(subject,max_overlap=10,min_perc_id_second_hitOR=95,min_length_second_hitOR=50,min_perc_id_second_hitAND=90,min_length_second_hitAND=30) # the two limits are OR, the other two AND
{
    #print(subject)
    turn_to_coverage_matrix <- function(x){X<-rep(0,x[1]);X[x[2]:x[3]]<-1;return(X)}
    b <- B[B$sseq==subject,]
    b <- b[order(b$bitscore,decreasing=T),]

    if (length(unique(b$qseqid))>1)
    {
        # Now this ordering is a bit weird, but it is important
        # In stead of just ordering per bit score, we order per cumulative bitscore per contig
        # Next, we order per bitscore within contig
        b$qseqid <- factor(b$qseqid)
        order1   <- names(sort(by(b$bitscore,b$qseqid,sum),decreasing=T))
        b        <- b[unlist(sapply(order1, function(x) which (b$qseqid==x))),]
        
        # Now before we move on we check whether it is worth moving on. 
        # We will execute the function only if there are potentially interesting hits beyond the best scoring one:
        firsthit   <- as.character(b$qseqid[1])
        bs         <- b[b$qseqid!=firsthit,]
        bs         <- bs[bs$pident>min_perc_id_second_hitOR | bs$length>min_length_second_hitOR,]
        bs         <- bs[bs$pident>min_perc_id_second_hitAND & bs$length>min_length_second_hitAND,]
        b          <- rbind.data.frame(b[b$qseqid==firsthit,],bs)
        if (length(unique(b$qseqid))>1)
        {
            bb         <- cbind.data.frame(b$slen,b$sstart,b$send)
            
            covmat     <- t(apply(bb,1,turn_to_coverage_matrix))
            outmat     <- rbind(0,covmat[1,])
            whichrows  <- c(1)
            for(i in 2:nrow(covmat))
            {
                testmat <- rbind(outmat,covmat[i,])
                # If the overlap is smaller than 10 amino acids in total
                if(sum(colSums(testmat)==2)<=max_overlap)
                {
                    outmat<-testmat
                    whichrows<-c(whichrows, i)
                }
            }   
            b <- b[whichrows,]
        }
    }
    return(b)
}

##############################
### THE SCAFFOLD PLOT FUNCTION
##############################

plotblastmap<-function(subject,x)
{
    b        <- x[x$sseqid==subject,]   
    b        <- b[order(b$bitscore,decreasing=T),]
    b$qseqid <- factor(b$qseqid)
    order1   <- names(sort(by(b$bitscore,b$qseqid,sum),decreasing=T))
    b        <- b[unlist(sapply(order1, function(x) which(b$qseqid==x))),]
    b$qseqid <- factor(b$qseqid,levels=unique(b$qseqid))
          
    
    xmax   <- b$slen[1]
    height <- max(3,(length(unique(b$qseqid)) +1))
    plot(1,1,xlim=c(0,xmax),ylim=c(0,height),col="white",main=subject,yaxt="n",ylab="")
    for(line in 1:nrow(b))
    {
        ycoord <- as.numeric(b$qseqid)[line]
        B1     <- round(b$sstart[line]-((min(b$qstart[line],b$qend[line])/3)),0)
        B2     <- round(B1 + (b$qlen[line]/3),0)
        B      <- B1:B2
        points(B,rep(ycoord,length(B)),col="grey",pch="*")
    }
    
    for(line in 1:nrow(b))
    {
        ycoord <- as.numeric(b$qseqid)[line]
        P      <- b$sstart[line]:b$send[line]
        points(P,rep(ycoord,length(P)),col="red",pch="*")
        text(mean(P),ycoord-(0.5/height),labels=round(b$pident[line],1))
    }
}   


get_blast             <- function(x)
{unlist(strsplit(unlist(strsplit(x[1],split="Full="))[2],split=";"))[1]}

my_GO_from_trinotate2 <- function(x,id,blastp=T,blastx=T,pfam=T)
{
  readerfunction <- function(x)
  {
    as.vector(sapply(unlist(sapply(x,function(x){return(unlist(strsplit(as.character(x),split="`")))})),function(x){return(unlist(strsplit(x,split="\\^"))[1])}))
  }
  outputlist <- vector("list",length(id))
  names(outputlist) <- id
  for (name in id)
  {
    if(blastp==T){goblastp <- readerfunction(paste(x$gene_ontology_BLASTP[x$transcript_id==name],collapse="`"))}else{goblast<-NULL}
    if(blastx==T){goblastx <- readerfunction(paste(x$gene_ontology_BLASTX[x$transcript_id==name],collapse="`"))}else{goblast<-NULL}
    if(pfam==T)  {gopfam  <- readerfunction(paste(x$gene_ontology_pfam[x$transcript_id==name],collapse="`"))}else{gopfam<-NULL}
    output <- unique(c(goblastp,goblastx,gopfam))
    output <- output[grepl("GO",output)]
    if(length(output)>0)
    {output<-paste(output,collapse=";");outputlist[[name]] <- output}else
    {outputlist[[name]] <- "unknown"}
  }
  return(outputlist)
}

# A small function that selects within a txi object only those features above a coverage cutoff in a certain number of samples
my_coveragefilter <- function(txi_object, mincounts = 2, minsamples = 3) 
{
  mappings_in_all_samples <- NULL
  for(i in 0:25000)
  {
    mappings_in_all_samples<-c(mappings_in_all_samples,sum(rowSums(txi_object$counts > i) == ncol(txi_object$counts)))
  }
  plot(mappings_in_all_samples,type="l",log="x",xlab="number of mappings",ylab="number of contigs",main="Number of contigs with at least x mappings in all samples")
  
  
  tokeep <- !apply(txi_object$counts,1,function(x){sum(x<mincounts)})>= minsamples
  txi_object$abundance 	  <- txi_object$abundance [tokeep,]
  txi_object$counts 		  <- txi_object$counts    [tokeep,]
  txi_object$variance 	  <- txi_object$variance  [tokeep,]
  txi_object$length 		  <- txi_object$length    [tokeep,]
  return(txi_object)
}  
# function that puts everything in a simple object that I can easily work with:

my_deseqfun      <- function(dds,contrast,meta,trinot,remove_unknown=T)
{
  # Some functions only to be defined within this one (prevents clutter in global environment)

  
  # Return transfos and data
  vst.s   <- vst(dds, blind=FALSE) # this was called vst_something_supervised
  vst.u   <- vst(dds, blind=TRUE)  # this was called vst_something_unsupervised
  dat.s   <- assay(varianceStabilizingTransformation(dds, blind=FALSE)) # this was called DAT_vst_something_supervised
  dat.u   <- assay(varianceStabilizingTransformation(dds, blind=TRUE))  # this was called DAT_vst_something_unsupervised
  
  if(any(contrast=="sex")) 
  {
    col_names <- paste(meta$sample, meta$sex, toupper(meta$dispersotype),sep="_")
  }else{
    col_names <- paste(meta$sample, meta$month, toupper(meta$dispersotype),sep="_")
  }
  
  colnames(dat.s) <- col_names
  colnames(dat.u) <- col_names
  
  # Run binomial wald test:
  res     <- DESeq2::results(dds,contrast=contrast)
  res.o   <- res[order(res$pvalue),] # this was called resOrdered_something 
  res.o   <- cbind.data.frame("id"=rownames(res.o),res.o[,c(1,2,5,6)])
  # Now merge DE info with trinotate info
  # Careful: doubles are created because sometimes multiple protein predictions per transcript!
  res.o           <- merge(res.o,trinot,by.x="id",by.y="transcript_id",all.x=T,all.y=F)
  colnames(res.o) <- c("id","mean","log2FC","pval","padj","locus","sprot_Top_BLASTX_hit","RNAMMER","prot_id","prot_coords","sprot_Top_BLASTP_hit","Pfam","SignalP","TmHMM","eggnog","Kegg","gene_ontology_blast","gene_ontology_pfam","transcript","peptide")
  res.o_blastx    <- unlist(lapply(res.o$sprot_Top_BLASTX_hit,get_blast))
  res.o_blastp    <- unlist(lapply(res.o$sprot_Top_BLASTP_hit,get_blast))
  res.o           <- cbind.data.frame(res.o,"blastx"=res.o_blastx,"blastp"=res.o_blastp)
  res.o           <- res.o[order(res.o$pval),] # this was called resOrdered_something 
  
  # Generate output for GO_MWU:
  transcript_ids    <- unique(res.o$id)
  annot             <- unlist(my_GO_from_trinotate2(trinot,transcript_ids))
  annot            	<- cbind.data.frame(transcript_ids,annot)
  colnames(annot)  	<- c("id","go")
  
  ### IMPORTANT
  ### This is new since oct 2020
  ### Remove the unknowns!
  ### They turn your inferences overconfident!
  if(remove_unknown){annot <- annot[annot$go != "unknown",]}
  
  pvals				<-cbind.data.frame(rownames(res),-log10(res$pvalue))
  colnames(pvals)		  <-c("id","p")
  pvals <- pvals[!is.na(pvals$p),] 
  pvals <- pvals[pvals$id %in% annot$id,]
  
  FC           	<-cbind.data.frame(rownames(res),res$log2FoldChange)
  colnames(FC) 	<-c("id","foldchange")
  FC <- FC[!is.na(FC$foldchange),] 
  FC <- FC[FC$id %in% annot$id,]
  
  FC_abs <- FC
  FC_abs$foldchange	  <- abs(FC_abs$foldchange)
  colnames(FC_abs)  	<-c("id","foldchange_abs")
  
  MWU <- list("annot"=annot,"pvals"=pvals,"FC"=FC,"absFC"=FC_abs)
  
  return(list("dds"=dds,"res"=res,"res.o"=res.o,"vst.s"=vst.s,"vst.u"=vst.u,"dat.s"=dat.s,"dat.u"=dat.u,"meta"=meta,"contrast"=contrast,"MWU"=MWU))
}



# This is a function I can run on the output of my_deseqfun
my_deseq_summary <- function(x,padj_cutoff=0.05,log2FC_cutoff=1)
{
  
  xx <- x$res.o[!duplicated(x$res.o$locus),]
  # DEG stats:
  cat("number of tests:\n")
  ntests <- nrow(xx) - sum(is.na(xx$padj))
  ups    <- sum(xx$padj < padj_cutoff & xx$log2FC > log2FC_cutoff ,na.rm=T) 
  downs  <- sum(xx$padj < padj_cutoff & xx$log2FC < -log2FC_cutoff ,na.rm=T) 
  print(ntests)
  
  cat(paste("number of p.adj < ", padj_cutoff, "and abs (log2FC) >", log2FC_cutoff, ":\n", sep=" "))
  
  
  
  cat (paste ("up in", x$contrast[2], ":",sep=" "))
  print(ups)
  cat ("(",round(ups/ntests*100,2),"%)")
  cat("\n")
  cat (paste ("up in", x$contrast[3], ":",sep=" "))
  print(downs) 
  cat ("(",round(downs/ntests*100,2),"%)")
  
  # Annot stats:
  
  # And 4 informative plots:
  par(mfrow=c(2,2))
  
  # Shrinkage plot:
  plotDispEsts(x$dds)
  
  # density of transformed values
  plot(density(x$dat.u[,1],bw=0.25),main="density of transformed values",xlab="transformed values",ylim=c(0,0.25))
  for (i in 2:ncol(x$dat.u)){lines(density(x$dat.u[,i],bw=0.25))}
  
  fclim <- max(abs(xx$log2FC))
  
  
  ### Log 2 fold change direction:
  FC_lab <- paste("up in", x$contrast[3], "<- log2FC ->" , "up in", x$contrast[2],sep=" ")
  
  # MA plot
  plot(log2(xx$mean),xx$log2FC,pch=20,ylim=c(-fclim,fclim),col=scales::alpha(c(xx$padj<padj_cutoff)+1,0.5),xlab="log 2 mean counts",ylab=FC_lab,main = "MA plot")
  
  # Volcano plot:
  plot(-log10(xx$pval)~xx$log2FC,xlim=c(-fclim,fclim),pch=20,xlab=FC_lab,ylab="-log10 p-value",col=scales::alpha("black",0.5), main = "volcano plot")
  abline(h=-log10(padj_cutoff),col="red")
}


my_deseq_intersect <- function(x1,x2,padj_cutoff=0.05,log2FC_cutoff=0)
{  
  xx1 <- x1$res.o[!duplicated(x1$res.o$locus),]
  xx2 <- x2$res.o[!duplicated(x2$res.o$locus),]
  
  ID_all <- unique(c(xx1$id,xx2$id))
  
  ups1    <- as.vector(na.omit(xx1$id[xx1$padj < padj_cutoff & xx1$log2FC > log2FC_cutoff])) 
  downs1  <- as.vector(na.omit(xx1$id[xx1$padj < padj_cutoff & xx1$log2FC < -log2FC_cutoff])) 
  
  ups2    <- as.vector(na.omit(xx2$id[xx2$padj < padj_cutoff & xx2$log2FC > log2FC_cutoff])) 
  downs2  <- as.vector(na.omit(xx2$id[xx2$padj < padj_cutoff & xx2$log2FC < -log2FC_cutoff]))     
  
  TAB   <- table(ID_all %in% c(ups1,downs1),ID_all %in% c(ups2,downs2))
  CHISQ <- chisq.test(TAB)
  cat("Total passing criteria: ")
  cat(TAB[2,2])
  cat("\nExpected:             ")
  cat(round(CHISQ$expected[2,2]))
  cat("\n");
  
  
  upup     <- intersect(ups1,ups2)
  cat (length(upup)); cat (" up in "); cat(x1$contrast[2]); cat(" and up in "); cat(x2$contrast[2]); cat("\n");
  cat ("Expected: "); cat (round(chisq.test(table(ID_all %in% ups1 ,ID_all %in% ups2))$expected[2,2])); cat("\n");
  updown   <- intersect(ups1,downs2)
  cat (length(updown));cat (" up in "); cat(x1$contrast[2]); cat(" and up in "); cat(x2$contrast[3]); cat("\n");
  cat ("Expected: "); cat (round(chisq.test(table(ID_all %in% ups1 ,ID_all %in% downs2))$expected[2,2])); cat("\n");
  downup   <- intersect(downs1,ups2)
  cat (length(downup));cat (" up in "); cat(x1$contrast[3]); cat(" and up in "); cat(x2$contrast[2]); cat("\n");
  cat ("Expected: "); cat (round(chisq.test(table(ID_all %in% downs1 ,ID_all %in% ups2))$expected[2,2])); cat("\n");
  downdown <- intersect(downs1,downs2)
  cat (length(downdown));cat (" up in "); cat(x1$contrast[3]); cat(" and up in "); cat(x2$contrast[3]); cat("\n");
  cat ("Expected: "); cat (round(chisq.test(table(ID_all %in% downs1 ,ID_all %in% downs2))$expected[2,2])); cat("\n");
  
  return(list("upup"=upup,"updown"=updown,"downup"=downup,"downdown"=downdown,"chisq"=CHISQ))
}  


### This is an interesting object/function/object to quickly recover all info you need:
### A bit of a cumbersome way to transform this data structure (there must be easier, but it works)
### This function takes a pair ("id","go") where go can be a single term, or a semicolon separated one
### It will return a dataframe with "id" repeated at each line, and separate GO terms per line, but only in requested ontology
### It will also add all ancestor terms!
my_complete_go <- function(x,ontology="BP")
{		
  terms    <- unlist(strsplit(as.character(x[2]),split=";"))
  terms    <- terms[!terms %in% c("unknown","obsolete")]
  allterms <- terms
  for(term in terms)
  { 
    #print(term)
    TERM <- GOTERM[[term]]
    if(!is.null(TERM))
    {
      ONT     <- TERM@Ontology
      if(ONT=="BP" & ontology == "BP")
      {
        allterms <-  c(allterms,(get(term,GOBPANCESTOR)))
      }
      if(ONT=="CC" & ontology == "CC")
      {
        allterms <-  c(allterms,(get(term,GOCCANCESTOR)))
      }
      
      if(ONT=="MF" & ontology == "MF")
      {
        allterms <-  c(allterms,(get(term,GOMFANCESTOR)))
      }
    }
  }
  allterms <- unique(allterms)
  allterms <- allterms[allterms!="all"]
  id <- rep(x[1],length(allterms))
  #cat(".")
  return(cbind.data.frame("id"=id,"go"=allterms))		
}		



# Thi sis a bit of a wacky function because the function is wacky as it is
# It is also extremely slow for no good reason...
my_GO_MWU <- function(dds_object,statistic="pvals",ontology="BP",directory,largest=0.1, smallest=5, clusterCutHeight=0.25, alternative="g",clean=T)
{
  system(paste("git clone https://github.com/z0on/GO_MWU", directory,sep=" "))
  setwd(directory)
  
  write.table(x=dds_object$MWU$annot,file="annotations.go",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(x=dds_object$MWU[[statistic]], file="X",sep=",",quote=F,row.names=F,col.names=T)
  source("gomwu.functions.R")
  gomwuStats("X","go.obo", "annotations.go",ontology,perlPath="perl", largest=largest, smallest=smallest, clusterCutHeight=clusterCutHeight, Alternative=alternative)
  
  file <- paste("MWU",ontology,"X",sep="_")
  dds_object$go_results[[ontology]][[statistic]] <- read.table(file,header=T)
  
  setwd("..")
  if(clean) {system(paste("rm -rf", directory,sep=" "))}
  
  return(dds_object)
  
}

### Function to add long format frame with GO results and DEG results
### Careful. This function takes ALL_x objects directly from the environment.
### Not very pretty, but it works.



my_add_long_format <- function(dds,statistic="pvals",ontology="BP")
{
  
  # In a first step, take the results, and edit a little bit:
  # So, GO_MWU clusters terms that are overlapping
  # And gives them as combination terms, eg GO:0009394;GO:0009262;GO:0019692 deoxyribonucleotide metabolic process
  # The problem is that the order of the terms doesn't stick with the description
  go_results <- cbind.data.frame(dds$go_results[[ontology]][[statistic]],"clustered"=NA)
  ONT<-unlist(Term(GOTERM))
  
  for(i in 1:nrow(go_results))
  {
    x <- go_results[i,]
    if(grepl(";",x$term))
    {
      #cat("*")
      x$clustered <- x$term
      term        <- names(which(x$name==ONT))
      if(identical(term,character(0)))
      {
        cat("WARNING: ");cat(x$term);cat(" term not found\n")
        cat(i);cat("\n")
        term <- unlist(strsplit(x$term,";"))[1]  
      }  
      x$term      <- term
      go_results[i,] <- x   
      rm(x)
    }#else{cat(".")}  
  }  
  # Now merge with long format ontology information:
  if(ontology=="BP"){M <- merge (ALL_BP,go_results[,c(1,2,3,5,7,8)],by.x="go",by.y="term",all=T)}
  if(ontology=="MF"){M <- merge (ALL_MF,go_results[,c(1,2,3,5,7,8)],by.x="go",by.y="term",all=T)}
  if(ontology=="CC"){M <- merge (ALL_CC,go_results[,c(1,2,3,5,7,8)],by.x="go",by.y="term",all=T)}
  M <- merge (M,as.data.frame(dds$res[,c(1,2,5,6)]),by.x="id",by.y="row.names",all=T)
  M <- merge (M,dds$dat.u,by.x="id",by.y="row.names",all=T)
  colnames(M) <- c(c("id","GO", "dran_GO", "pval_GO","level_GO","padj_GO", "clustered","basemean", "log2FC", "pval_DEG", "padj_DEG"),colnames(M)[12:ncol(M)])
  M$dran_GO  <- as.numeric(M$dran_GO)    
  M$pval_GO  <- as.numeric(M$pval_GO)
  M$level_GO <- as.numeric(M$level_GO)
  M$padj_GO  <- as.numeric(M$padj_GO)
  M$basemean <- as.numeric(M$basemean)
  M$log2FC   <- as.numeric(M$log2FC)
  M$pval_DEG <- as.numeric(M$pval_DEG)
  M$padj_DEG <- as.numeric(M$padj_DEG)
  for(i in 12:ncol(M))
  {
    M[i,] <- as.numeric(M[i,])
  }  

  M <- unique(M)
  M[is.na(M$GO),"GO"] <- "unknown"
  M <- M[!is.na(M$basemean),]
  M <- M[order(M$pval_GO,decreasing=F),]
  
  go_results <- go_results[order(go_results$pval),]
  return(list("M"=M,"go_results"=go_results))
}



my.little.heatmap <- function(x, GOTERM, pval_cutoff=1, padj_cutoff=1, title=NA, col=cm.colors(256),indices=c(1,3,6,2,4,5))
{
  if(!grepl("GO:",GOTERM)){GOTERM<-as.character(x[x$name==GOTERM,"term"])}
  #print(GOTERM)
  if(is.na(title)){title<-paste(GOTERM,Term(GOTERM),sep=" ")}
  #print(title)
  
  x  <- x[x$GO==GOTERM,]
  x$pval_DEG[is.na(x$pval_DEG)] <- 1
  x$padj_DEG[is.na(x$padj_DEG)] <- 1
  x <- x[x$pval_DEG < pval_cutoff & x$padj_DEG < padj_cutoff,]
  
  x  <- x[!is.na(x$basemean),]
  # Scale per gene. Note the double transformation
  HM_input <- t(scale(t(x[12:ncol(x)])))[,indices]
  rownames(HM_input) <- x$id
  foo<-heatmap(HM_input,Colv=NA,scale="none",col = col,main=title)
  x <- x[order(x$pval_DEG),]
  return(x)
}


# function to quickly annotate contigs based on their LOC identifier:
my_annot<-function(x)
{
  annots <- read.delim("~/pieris/loci_annots.txt",sep="\t",header=F)
  A<-NULL
  for(i in x)
  {
    annot  <- annots[annots[,1]==i,2][1]
    A <- c(A,annot)
    
  }
  return(A)
}




