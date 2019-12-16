# @author Dinesh Barupal dinkumar@ucdavis.edu
# @version 1.0 July 2016
# @use call "getMetaMapp" function via a Ajax call from JS.
pacman::p_load(RColorBrewer, officer, rvg, ggplot2, dplyr, ggrepel, extrafont, dynamicTreeCut, plotly, htmlwidgets)
getChemSimNet <- function (cids, cutoff=0.7) {
  #cids <- c(1:5)
  #cutoff <- 0.7
  if (length(cids)>200) {
    hex <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f")
    bin <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
    cidlist <- split(cids, ceiling(seq_along(cids)/200))
    
    
    
    
    
    
    subkeys <- do.call(rbind,lapply(cidlist,function(x) { read.csv(paste(c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',paste(x,collapse=","),'/property/Fingerprint2D/csv'),collapse=""))}))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    m <- t(sapply(subkeys[,2],function(x) { as.integer(strsplit(paste(sapply(strsplit(paste(RCurl::base64Decode(x,"raw")[5:115],collapse=""),"")[[1]],function(x) {bin[which(hex==x)]}),collapse=""),"")[[1]][1:881] ) }))
    mat <- m%*%t(m)
    len <- length(m[,1])
    s <- mat.or.vec(len,len)
    for (i in 1:len) {
      for (j in 1:len){
        s[i,j] <- mat[i,j]/(mat[i,i]+mat[j,j]-mat[i,j])
      }
    }
    diag(s) <- 0
    dfmax <- cbind(cids,"tmsim",cids[sapply(1:length(cids), function (k) {which.max(s[k,]) })])
    s[lower.tri(s)]<-0
    chemsimdf <- do.call(rbind,sapply(1:length(cids), function (k) { if(length(which(s[k,]>cutoff))>0) {cbind(cids[k],"tmsim",cids[which(s[k,]>cutoff)])}} ))
    chemsimdf <- rbind(chemsimdf, cbind(cids,"tmsim",""), dfmax )
    
    ## Duplicated edges removal
    chemsimdf <- chemsimdf[-which(chemsimdf[,3]==""),]
    chemsimdf <- rbind(chemsimdf,cbind(chemsimdf[,3], chemsimdf[,2], chemsimdf[,1]))
    chemsimdf <-  chemsimdf[!duplicated( chemsimdf),]
    chemsimdf <- chemsimdf[order(chemsimdf[,3],decreasing = F),]
    chemsimdf <- chemsimdf[order(chemsimdf[,1],decreasing = F),]
    
    pmids_a <- cids
    
    for (i in 1:length(pmids_a)) {
      sind <- c((max(which(chemsimdf[,1]==pmids_a[i])) +1) :nrow(chemsimdf))
      chemsimdf[,3][which(chemsimdf[,3][sind]==pmids_a[i]) + (sind[1]-1) ] <- "XX"
    }
    chemsimdf <- chemsimdf[-which(chemsimdf[,3]=="XX"),]
    
    
    write.table(chemsimdf,file=paste(c("chemsim_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
    return(chemsimdf)
  } else{
    
    hex <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f")
    bin <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
    subkeys <- read.csv(paste(c('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',paste(cids,collapse=","),'/property/Fingerprint2D/csv'),collapse=""))
    m <- t(sapply(subkeys[,2],function(x) { as.integer(strsplit(paste(sapply(strsplit(paste(RCurl::base64Decode(x,"raw")[5:115],collapse=""),"")[[1]],function(x) {bin[which(hex==x)]}),collapse=""),"")[[1]][1:881] ) }))
    mat <- m%*%t(m)
    len <- length(m[,1])
    s <- mat.or.vec(len,len)
    for (i in 1:len) {
      for (j in 1:len){
        s[i,j] <- mat[i,j]/(mat[i,i]+mat[j,j]-mat[i,j])
      }
    }
    diag(s) <- 0
    dfmax <- cbind(cids,"tmsim",cids[sapply(1:length(cids), function (k) {which.max(s[k,]) })])
    s[lower.tri(s)]<-0
    chemsimdf <- do.call(rbind,sapply(1:length(cids), function (k) { if(length(which(s[k,]>cutoff))>0) {cbind(cids[k],"tmsim",cids[which(s[k,]>cutoff)])}} ))
    chemsimdf <- rbind(chemsimdf, cbind(cids,"tmsim",""), dfmax )
    
    ## Duplicated edges removal
    chemsimdf <- chemsimdf[-which(chemsimdf[,3]==""),]
    chemsimdf <- rbind(chemsimdf,cbind(chemsimdf[,3], chemsimdf[,2], chemsimdf[,1]))
    chemsimdf <-  chemsimdf[!duplicated( chemsimdf),]
    chemsimdf <- chemsimdf[order(chemsimdf[,3],decreasing = F),]
    chemsimdf <- chemsimdf[order(chemsimdf[,1],decreasing = F),]
    
    pmids_a <- cids
    
    for (i in 1:length(pmids_a)) {
      sind <- c((max(which(chemsimdf[,1]==pmids_a[i])) +1) :nrow(chemsimdf))
      chemsimdf[,3][which(chemsimdf[,3][sind]==pmids_a[i]) + (sind[1]-1) ] <- "XX"
    }
    chemsimdf <- chemsimdf[-which(chemsimdf[,3]=="XX"),]
    
    write.table(chemsimdf,file=paste(c("chemsim_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
    return(chemsimdf)
  }
}
## Tanimoto score calculation was adopted from http://data2quest.blogspot.com/2013/10/fast-tanimoto-similarity-calculation.html

krp <- read.table("https://raw.githubusercontent.com/barupal/metamapp/master/KRPlinks.txt", sep="\t")
# cids = df1[,1][is.na(df1[,1])==FALSE];keggids= df1[,2][is.na(df1[,1])==FALSE]
getKEGGRpairs <- function (cids, keggids, cutoff=0.7) {
  krp.1 <- match(krp[,1],keggids)
  krp.2 <- match(krp[,2],keggids)
  krp.cbind <- cbind (krp.1,krp.2)
  krp.net <- subset(krp.cbind, krp.1!="NA" & krp.2!="NA")
  cid.krp.2 <- cids[krp.net[,2]]
  cid.krp.1 <- cids[krp.net[,1]]
  krp.cid.net <- cbind(cid.krp.1,"krp",cid.krp.2)
  chemsim <- getChemSimNet(cids,cutoff)
  krp.cid.net <- rbind(krp.cid.net,chemsim)
  write.table(krp.cid.net,file=paste(c("chemsim_krp_",gsub("[.]","",as.character(cutoff)),".sif"),collapse=""), quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)  ## To write the cytoscape network file as an output
}


getNodeAttributes <- function (stat_file,samplecount) {
  names <- strsplit(strsplit(stat_file,"\n")[[1]][1],"\t")[[1]][-1]
  df1 <- do.call(rbind,lapply(strsplit(stat_file,"\n")[[1]][-1],function(x){as.integer(strsplit(x,"\t")[[1]])}))
  pvals <- sapply(1:length(df1[,1]),function(i){ t.test(df1[i,2:(samplecount+1)],df1[i,(samplecount+2):(samplecount+1+samplecount)],alternative = c("two.sided", "less", "greater"), var.equal=TRUE)$p.value  })
  cids <- df1[,1]
  statres <- do.call(rbind,lapply(1:length(pvals),function(z){ if(pvals[z]<0.05) {fc <- median(df1[z,2:(samplecount+1)])/median(df1[z,(samplecount+2):(samplecount+1+samplecount)]) ; if (fc>1.0){return( c(cids[z],pvals[z],round(fc,2),"Up")   )}else{  return( c(cids[z],pvals[z],round(1/fc,2),"Down")) }     } else {return(c(cids[z],pvals[z],1,"No Change"))} }))
  colnames(statres) <- c("CID","p.value","fold-change","direction")
  write.table(statres,file=paste(c("ttest_attr_",gsub("[.]| ","",names[1]),"_Vs_",gsub("[.]| ","",names[1+samplecount]),".txt"),collapse=""), quote=FALSE,sep="\t",col.names=T,row.names=FALSE)  ## To write the cytoscape network file as an output
}

getMetaMapp <- function (rsess, pvalind=1, cutoff=0.7) {
  statres <- paste(rsess,"R/.val/csv",sep="")
  df <- read.csv(statres)
  getKEGGRpairs(df$PubChem[is.na(df$PubChem)==FALSE], as.character(df$KEGG[is.na(df$PubChem)==FALSE]), cutoff)
  fcdf <- df[which(is.na(df$PubChem)==FALSE),grep("^Mean",colnames(df))]
  fcdf$meanratio <- fcdf$Mean.of.t1...None / fcdf$Mean.of.t2...Test.Compound
  pvaldf <- df[which(is.na(df$PubChem)==FALSE),grep("p_value",colnames(df))]
  exportdf <- data.frame(Pubchem=df$PubChem[is.na(df$PubChem)==FALSE])
  exportdf$KEGG <- as.character(df$KEGG[is.na(df$PubChem)==FALSE])
  exportdf$fold.change <- rep(1,length(exportdf[,1]))
  exportdf$pvaldirection <- rep("No Change",length(exportdf[,1]))
  exportdf <- cbind(exportdf,df[which(is.na(df$PubChem)==FALSE),grep("p_value",colnames(df))])
  exportdf$CpdName <- as.character(df$BinBase_name[is.na(df$PubChem)==FALSE])
  sigind <- which(exportdf[,grep("p_value",colnames(exportdf))[1]]<0.05)
  for( x in sigind)  {
    if(fcdf$meanratio[x]<1) {
      exportdf$fold.change[x] <- round(1/fcdf$meanratio[x],1)
      exportdf$pvaldirection[x] <- "Down"
    } else {
      exportdf$fold.change[x] <- round(fcdf$meanratio[x],1)
      exportdf$pvaldirection[x] <- "Up"
    }
  }
  #exportdf$meshanno <- sapply(exportdf$Pubchem, function (x) { paste(fromJSON( paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pccompound&db=mesh&id=",x,"&retmode=json",sep="") )$linksets[[1]]$linksetdbs[[1]]$links,collapse=",")  }   )
  write.table( exportdf, file=paste("node_attributes_chemsim_krp_",gsub("[.]","",as.character(cutoff)) ,".tsv", sep="" ), col.names = T, row.names = F, quote = F, sep = "\t" )
}

runMetaMapp <- function(stat_file, cutoff=0.7) {
  cfile <- strsplit(stat_file,"\n")[[1]]
  df1 <- do.call(rbind, lapply(cfile, function (x) { strsplit(x,"\t")[[1]]  } ))
  colnames(df1) <- df1[1,]
  df1 <- df1[-1,]
  getKEGGRpairs(df1[,1][is.na(df1[,1])==FALSE], df1[,2][is.na(df1[,1])==FALSE], cutoff)
  #exportdf <- data.frame(Pubchem_ID=df1[,1][is.na(df1[,1])==FALSE], KEGG_ID=df1[,2][is.na(df1[,1])==FALSE], CompoundName=df1[,3][is.na(df1[,1])==FALSE])
  exportdf <- as.data.frame(setNames(replicate(  length( grep("foldchange",colnames(df1) ) ), rep("No Change",length(df1[,1])), simplify = F), paste(colnames(df1)[grep("foldchange",colnames(df1))],"_direction",sep="") ), stringsAsFactors=FALSE)
  for (k in grep("pvalue",colnames(df1)  )) {
    df1[,k] <- as.numeric(df1[,k])
    sigind <- which(df1[,k]<0.05)
    df1[which(1:length(df1[,1])%in%sigind==FALSE),(k+1)] <- 1.0  ## convert all the non-significant fold changes to 1.00.
    for( x in sigind)  {
      if(df1[x,(k+1)]<1) {
        df1[x,(k+1)] <- round(1/as.numeric(df1[x,(k+1)] ),1)
        exportdf[x,paste(colnames(df1)[k+1],"_direction",sep="")] <- "Down"
      } else {
        df1[x,(k+1)] <- round(as.numeric(df1[x,(k+1)] ),1)
        exportdf[x,paste(colnames(df1)[k+1],"_direction",sep="")] <- "Up"
      }
    }
  }
  exportdf <- cbind(df1, exportdf)
  write.table( exportdf, file=paste("node_attributes_chemsim_krp_",gsub("[.]","",as.character(cutoff)) ,".tsv", sep="" ), col.names = T, row.names = F, quote = F, sep = "\t" )
}


chemrich.getChemicalClass = function(){
  if (length(which(colnames(chemrich.input.file) == "ChemicalClass")) != 
      0) {
    chemrich.input.file$ChemRICHClusters <- chemrich.input.file$ChemicalClass
    return(chemrich.input.file)
  }
  smi.all.fas <- as.character(sapply(chemrich.input.file$SMILES, 
                                     makeSmiles.clean))
  falabelvec <- sapply(smi.all.fas, function(x) {
    elecount <- table(strsplit(gsub("[0-9]|[)]|[(]|=", 
                                    "", x), "")[[1]])
    falabel <- ""
    if (length(table(c("c", "o") %in% tolower(names(elecount)))) == 
        1) {
      if (length(grep("n", x, ignore.case = T)) == 
          0) {
        if (elecount["C"] > 7 & length(grep("CCCC", 
                                            x)) == 1 & length(grep("C2", x)) != 1) {
          if (elecount["O"] == 2) {
            dlen <- length(strsplit(x, "=")[[1]]) - 
              2
            falabel <- paste(c("FA", elecount["C"], 
                               dlen), collapse = "_")
          }
          if (elecount["O"] >= 3) {
            if (length(grep("C1", x)) == 1) {
              if (length(strsplit(x, "C1")[[1]]) == 
                  3) {
                dlen <- length(strsplit(x, "=")[[1]]) - 
                  2
              }
              else {
                dlen <- length(strsplit(x, "=")[[1]]) - 
                  2
                falabel <- paste(c("Epoxy FA", 
                                   elecount["C"]), collapse = "_")
              }
            }
            else {
              if (length(strsplit(x, "=O|CO|OC")[[1]]) - 
                  2 == 0) {
                dlen <- length(strsplit(x, "=")[[1]]) - 
                  2
                falabel <- paste(c("OH-FA", elecount["C"], 
                                   dlen, (elecount["O"] - 2)), collapse = "_")
              }
              else {
                if (length(strsplit(x, "OC|CO")[[1]]) < 
                    3) {
                  dlen <- length(strsplit(x, "=")[[1]]) - 
                    2
                  falabel <- paste(c("O=FA", elecount["C"], 
                                     dlen), collapse = "_")
                }
              }
            }
          }
        }
      }
    }
    falabel
  })
  falabelvec[which(falabelvec == "OH-FA_20_3_2")] <- "DiHETrE"
  falabelvec[which(falabelvec == "OH-FA_20_4_2")] <- "DiHETE"
  falabelvec[which(falabelvec == "O=FA_18_3")] <- "oxo-ODE"
  falabelvec[which(falabelvec == "O=FA_20_5")] <- "oxo-ETE"
  falabelvec[which(falabelvec == "OH-FA_18_1_2")] <- "DiHOME"
  falabelvec[which(falabelvec == "OH-FA_18_1_3")] <- "TriHOME"
  falabelvec[which(falabelvec == "OH-FA_18_2_1")] <- "HODE"
  falabelvec[which(falabelvec == "OH-FA_18_2_2")] <- "DiHODE"
  falabelvec[which(falabelvec == "OH-FA_18_3_1")] <- "HOTrE"
  falabelvec[which(falabelvec == "OH-FA_20_3_1")] <- "HETrE"
  falabelvec[which(falabelvec == "OH-FA_20_4_1")] <- "HETE"
  falabelvec[which(falabelvec == "OH-FA_20_5_1")] <- "HEPE"
  falabelvec[which(falabelvec == "OH-FA_22_5_2")] <- "DiHDPE"
  falabelvec[which(falabelvec == "Epoxy FA_22")] <- "EpDPE"
  falabelvec[which(falabelvec == "Epoxy FA_18")] <- "EpETrE"
  falabelvec[which(falabelvec == "Epoxy FA_20")] <- "EpODE"
  falabelvec[grep("^FA_[0-9]{1,2}_0$", falabelvec)] <- "Saturated FA"
  falabelvec[grep("^FA_[0-9]{1,2}_[1-9]$", falabelvec)] <- "UnSaturated FA"
  cat("Computing sub-structure fingerprint\n")
  fps <- t(sapply(1:nrow(chemrich.input.file), function(x) {
    xy <- 0
    xy <- tryCatch({as.character(rcdk::get.fingerprint(rcdk::parse.smiles(chemrich.input.file$SMILES[x])[[1]], 
                                                       type = "pubchem"))
    }, error = function(er){
      0
    })
    xy
  }))
  cid.mesh.df <- data.frame(CID = df.mega.mesh$CID[df.mega.mesh$CID %in% 
                                                     chemrich.input.file$PubChemID], MESHTREE = df.mega.mesh$MESSTREE[df.mega.mesh$CID %in% 
                                                                                                                        chemrich.input.file$PubChemID], stringsAsFactors = F)
  
  
  print(sum(apply(fps, 2, nchar)==1))
  
  chemrich.input.file = chemrich.input.file[!apply(fps, 2, nchar)==1,]
  falabelvec = falabelvec[!apply(fps, 2, nchar)==1]
  
  
  
  fps <- t(sapply(1:nrow(chemrich.input.file), function(x) {
    xy <- 0
    xy <- tryCatch({as.character(rcdk::get.fingerprint(rcdk::parse.smiles(chemrich.input.file$SMILES[x])[[1]], 
                                                       type = "pubchem"))
    }, error = function(er){
      0
    })
    xy
  }))
  cid.mesh.df <- data.frame(CID = df.mega.mesh$CID[df.mega.mesh$CID %in% 
                                                     chemrich.input.file$PubChemID], MESHTREE = df.mega.mesh$MESSTREE[df.mega.mesh$CID %in% 
                                                                                                                        chemrich.input.file$PubChemID], stringsAsFactors = F)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  inmesh.vec <- rep("No", nrow(chemrich.input.file))
  inmesh.vec[chemrich.input.file$PubChemID %in% cid.mesh.df$CID] <- "Yes"
  df1.bitmat <- do.call(rbind, lapply(fps, function(x) as.integer(strsplit(x,"")[[1]][1:881])))
  df1.bitmat.location <- lapply(1:nrow(df1.bitmat), function(x) {
    which(df1.bitmat[x, ] == 1)
  })
  only_b <- sapply(1:length(bitloclist), function(x) {
    length(bitloclist[[x]])
  })
  bitmeans <- sapply(1:length(bitloclist), function(x) {
    median(bitloclist[[x]])
  })
  fpsmeans <- sapply(df1.bitmat.location, function(x) {
    median(x)
  })
  cat("Obtaining MESH class annotation\n")
  directlabels <- sapply(tolower(chemrich.input.file$CompoundName), 
                         function(x) {
                           clabel = "Not Found"
                           findind <- which(df.mega.mesh$CompoundName == x)
                           if (length(findind) > 0) {
                             classvec <- as.character(df.mega.mesh$MESSTREE[findind])
                             classvec <- strsplit(classvec[1], ";")[[1]]
                             if (length(grep("^D01[.]", classvec)) > 
                                 0) {
                               classvec <- classvec[-grep("^D01[.]|^D03[.]", 
                                                          classvec)]
                             }
                             clabel <- names(which.max(sapply(classvec, nchar)))
                           }
                           clabel
                         })
  
  
  
  
  
  labelvec <- sapply(1:nrow(chemrich.input.file), function(i) {
    clabel <- "Not Found"
    if (falabelvec[i] == "" & inmesh.vec[i] == "No" & 
        directlabels[i] == "Not Found") {
      meanindex <- which(bitmeans < (fpsmeans[i] + 5) & 
                           bitmeans > (fpsmeans[i] - 5))
      bitloclist.sb <- bitloclist[meanindex]
      only_b.sb <- only_b[meanindex]
      
      if(!length(bitloclist.sb) == 0){
        overlapvec <- sapply(1:length(bitloclist.sb), function(x) {
          length(which(bitloclist.sb[[x]] %in% df1.bitmat.location[[i]] == 
                         TRUE))
        })
        tmvec <- overlapvec/((length(df1.bitmat.location[[i]]) + 
                                only_b.sb) - overlapvec)
        if (length(which(tmvec > 0.9)) > 0) {
          if (length(which(tmvec > 0.98)) > 0) {
            cidindex <- meanindex[which(tmvec > 0.98)]
            if (length(cidindex) == 1) {
              clabel <- df.mega.mesh$MESSTREE[cidindex]
            }
            else {
              clabel.table <- sort(table(unlist(sapply(unique(df.mega.mesh[cidindex[order(tmvec[which(tmvec > 
                                                                                                        0.98)])], "MeSHUID"])[1:10], function(x) {
                                                                                                          if (!is.na(x)) {
                                                                                                            strsplit(df.mega.mesh$MESSTREE[which(df.mega.mesh$MeSHUID == 
                                                                                                                                                   x)][1], ";")
                                                                                                          }
                                                                                                        }))), decreasing = T)
              clabel.table <- which(clabel.table == clabel.table[1])
              clabel.table.sort <- sort(sapply(names(clabel.table), 
                                               nchar), decreasing = T)
              clabel.table.sort.max <- which.max(clabel.table.sort)
              if (length(clabel.table.sort.max == 1)) {
                clabel <- names(clabel.table.sort.max)
              }
              else {
                clabel <- sort(names(clabel.table.sort.max))[1]
              }
            }
          }
          else {
            cidindex <- meanindex[which(tmvec > 0.9)]
            if (length(cidindex) == 1) {
              clabel <- df.mega.mesh$MESSTREE[cidindex]
            }
            else {
              clabel.table <- sort(table(unlist(sapply(unique(df.mega.mesh[cidindex[order(tmvec[which(tmvec > 
                                                                                                        0.9)])], "MeSHUID"])[1:10], function(x) {
                                                                                                          if (!is.na(x)) {
                                                                                                            strsplit(df.mega.mesh$MESSTREE[which(df.mega.mesh$MeSHUID == 
                                                                                                                                                   x)][1], ";")
                                                                                                          }
                                                                                                        }))), decreasing = T)
              clabel.table <- which(clabel.table == clabel.table[1])
              clabel.table.sort <- sort(sapply(names(clabel.table), 
                                               nchar), decreasing = T)
              clabel.table.sort.max <- which.max(clabel.table.sort)
              if (length(clabel.table.sort.max == 1)) {
                clabel <- names(clabel.table.sort.max)
              }
              else {
                clabel <- sort(names(clabel.table.sort.max))[1]
              }
            }
          }
        }
      }
      
      
    }
    clabel
  })
  
  
  
  
  finalMesh.df <- cid.mesh.df
  finalMesh.df <- finalMesh.df[which(finalMesh.df$CID %in% 
                                       chemrich.input.file$PubChemID[which(falabelvec != "")] == 
                                       FALSE), ]
  finalMesh.df <- finalMesh.df[which(finalMesh.df$CID %in% 
                                       chemrich.input.file$PubChemID[which(directlabels != "Not Found")] == 
                                       FALSE), ]
  finalMesh.df <- rbind(finalMesh.df, data.frame(CID = chemrich.input.file$PubChemID, 
                                                 MESHTREE = labelvec))
  finalMesh.df$NewMesh <- finalMesh.df$MESHTREE
  finalMesh.df <- finalMesh.df[which(finalMesh.df$NewMesh != 
                                       "Not Found"), ]
  finalMesh.df <- finalMesh.df[which(finalMesh.df$NewMesh != 
                                       ""), ]
  finalMesh.df <- finalMesh.df[!duplicated(finalMesh.df), ]
  cat("Computing Chemical Similarity\n")
  m <- df1.bitmat
  mat <- m %*% t(m)
  len <- length(m[, 1])
  s <- mat.or.vec(len, len)
  for (i in 1:len) {
    for (j in 1:len) {
      s[i, j] <- mat[i, j]/(mat[i, i] + mat[j, j] - mat[i, 
                                                        j])
    }
  }
  diag(s) <- 0
  hc <- hclust(as.dist(1 - s), method = "ward.D2")
  clust1 <- cutreeDynamic(hc, distM = as.matrix(1 - s), deepSplit = 4, 
                          minClusterSize = 3)
  chemrich.input.file$ClusterNumber <- clust1
  chemrich.input.file$xlogp <- as.numeric(sapply(chemrich.input.file$SMILES, 
                                                 function(x) {
                                                   rcdk::get.xlogp(rcdk::parse.smiles(x)[[1]])
                                                 }))
  finalterm.df <- data.frame(CID = chemrich.input.file$PubChemID, 
                             Clabel = falabelvec, stringsAsFactors = F)
  directlabindex <- as.integer(which(directlabels != "Not Found"))[which(as.integer(which(directlabels != 
                                                                                            "Not Found")) %in% which(finalterm.df$Clabel == 
                                                                                                                       "") == TRUE)]
  finalterm.df$Clabel[directlabindex] <- as.character(directlabels[directlabindex])
  for (i in 1:nrow(finalterm.df)) {
    if (finalterm.df$Clabel[i] == "" & length(which(finalMesh.df$CID == 
                                                    chemrich.input.file$PubChemID[i])) > 0) {
      finalterm.df$Clabel[i] <- names(which.max(sapply(unlist(strsplit(finalMesh.df$NewMesh[which(finalMesh.df$CID == 
                                                                                                    chemrich.input.file$PubChemID[i])], ";")), 
                                                       nchar)))
    }
  }
  cat("Detecting New Compound Classes\n")
  finalterm.df.2 <- finalterm.df
  newClustVec <- names(which(table(chemrich.input.file$ClusterNumber[which(finalterm.df.2$Clabel == 
                                                                             "")]) > 3))
  clustMeanvec <- sapply(newClustVec, function(x) {
    mean(s[which(chemrich.input.file$ClusterNumber == x), 
           which(chemrich.input.file$ClusterNumber == x)])
  })
  newClustVec <- newClustVec[which(clustMeanvec > 0.7)]
  if (length(newClustVec) > 0) {
    for (i in which(finalterm.df.2$Clabel == "")) {
      if (chemrich.input.file$ClusterNumber[i] %in% newClustVec) {
        finalterm.df$Clabel[i] <- paste0("NewCluster_", 
                                         chemrich.input.file$ClusterNumber[i])
      }
    }
  }
  for (i in which(finalterm.df$Clabel == "")) {
    if (max(s[i, ]) > 0.75) {
      simorder <- order(s[i, ], decreasing = T)[which(s[i, 
                                                        ][order(s[i, ], decreasing = T)] > 0.75)]
      simorder.class <- sapply(simorder, function(x) {
        finalterm.df$Clabel[x]
      })
      simorder.class <- simorder.class[!is.na(simorder.class)]
      if (length(simorder.class) > 0) {
        if (simorder.class[1] != "") {
          finalterm.df$Clabel[i] <- simorder.class[which(simorder.class != 
                                                           "")][1]
        }
        else if (length(simorder.class) > 1) {
          finalterm.df$Clabel[i] <- simorder.class[which(simorder.class != 
                                                           "")][1]
        }
      }
    }
  }
  cat("Detecting non-overlapping class definition\n")
  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, 
                                          function(x) {
                                            length(which(finalterm.df$Clabel == x))
                                          }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, 
                                           function(x) {
                                             length(grep(x, finalterm.df$Clabel))
                                           }))
  exclusionVec <- c("D02", "D03.383", "D03.633.100", 
                    "D03.633.300", "D03.633.400", "D03.633", 
                    "D03.605", "D02.241.081")
  exclusionVec <- c(exclusionVec, unique(falabelvec)[-1])
  for (i in which(finalterm.df$gCount < 3)) {
    qpat <- gsub("[.][0-9]{2,3}$", "", finalterm.df$Clabel[i])
    if (length(grep(qpat, finalterm.df$Clabel)) > 2 & !qpat %in% 
        exclusionVec) {
      finalterm.df$Clabel[i] <- qpat
    }
  }
  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, 
                                          function(x) {
                                            length(which(finalterm.df$Clabel == x))
                                          }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, 
                                           function(x) {
                                             length(grep(x, finalterm.df$Clabel))
                                           }))
  for (i in which(finalterm.df$gCount < 3)) {
    if (max(s[i, ]) > 0.85) {
      simorder <- order(s[i, ], decreasing = T)[which(s[i, 
                                                        ][order(s[i, ], decreasing = T)] > 0.85)]
      simorder.class <- sapply(simorder, function(x) {
        finalterm.df$Clabel[x]
      })
      simorder.class <- simorder.class[!is.na(simorder.class)]
      if (simorder.class[1] != "") {
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class != 
                                                         "")][1]
      }
      else if (length(simorder.class) > 1) {
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class != 
                                                         "")][1]
      }
    }
  }
  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, 
                                          function(x) {
                                            length(which(finalterm.df$Clabel == x))
                                          }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, 
                                           function(x) {
                                             length(grep(x, finalterm.df$Clabel))
                                           }))
  for (i in which(finalterm.df$gCount < 3)) {
    qpat <- gsub("[.][0-9]{2,3}$", "", finalterm.df$Clabel[i])
    if (length(grep(qpat, finalterm.df$Clabel)) > 2 & !qpat %in% 
        exclusionVec) {
      finalterm.df$Clabel[i] <- qpat
    }
  }
  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, 
                                          function(x) {
                                            length(which(finalterm.df$Clabel == x))
                                          }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, 
                                           function(x) {
                                             length(grep(x, finalterm.df$Clabel))
                                           }))
  for (i in which(finalterm.df$gCount < 3)) {
    if (max(s[i, ]) > 0.85) {
      simorder <- order(s[i, ], decreasing = T)[which(s[i, 
                                                        ][order(s[i, ], decreasing = T)] > 0.85)]
      simorder.class <- sapply(simorder, function(x) {
        finalterm.df$Clabel[x]
      })
      simorder.class <- simorder.class[!is.na(simorder.class)]
      if (simorder.class[1] != "") {
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class != 
                                                         "")][1]
      }
      else if (length(simorder.class) > 1) {
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class != 
                                                         "")][1]
      }
    }
  }
  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, 
                                          function(x) {
                                            length(which(finalterm.df$Clabel == x))
                                          }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, 
                                           function(x) {
                                             length(grep(x, finalterm.df$Clabel))
                                           }))
  finalterm.df$Clabel[which(finalterm.df$Count < 3)] <- finalterm.df.2$Clabel[which(finalterm.df$Count < 
                                                                                      3)]
  finallabelvec <- finalterm.df$Clabel
  HasSaturatedFats <- names(which(table(finallabelvec[grep("D10|D09.400.410", 
                                                           finallabelvec)[which(sapply(grep("D10|D09.400.410", 
                                                                                            finallabelvec), function(x) {
                                                                                              length(grep("C=C", chemrich.input.file$SMILES[x]))
                                                                                            }) == 0)]]) > 2))
  for (i in 1:nrow(finalterm.df)) {
    if (finallabelvec[i] %in% HasSaturatedFats) {
      if (length(grep("C=C", chemrich.input.file$SMILES[i])) == 
          0) {
        finallabelvec[i] <- paste0("Saturated_", 
                                   getCNames(finallabelvec[i]))
      }
      else {
        finallabelvec[i] <- paste0("Unsaturated_", 
                                   getCNames(finallabelvec[i]))
      }
    }
  }
  clusterids <- sapply(as.character(finallabelvec), getCNames)
  cat("ChemRICH Class Estimation Finished. ChemRICHClusters column has been added to chemrich.input.file data.frame\n")
  chemrich.input.file$ChemRICHClusters <- clusterids
  chemrich.input.file
  
  
}
