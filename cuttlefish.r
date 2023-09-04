

args <- commandArgs(trailingOnly = TRUE)
dat1 = args[1]
dat2 = args[2]
name1 = args[3]
name2 = args[4]
pos = args[5]
end = args[6]
range = as.numeric(args[7])
dir = args[8]

print(args)

sim_function = function(dat1,dat2,name1,name2,pos,end,range,dir){
  
  #dat1 must be nv; dat2 must be other software data
  dat1 <- read.table(gzfile(dat1))
  dat2 <- read.table(gzfile(dat2))
  
  if (name1 == 'nv'){
    print('1 = nv')
    dat1_clean = parse_nv(dat1,1)
  } else if (name1 == 'snif'){
    print('1 = snif')
    dat1_clean = parse_snif(dat1,1)
  } else if (name1 == 'cuteSV'){
    print('1 = cuteSV')
    dat1_clean = parse_cuteSV(dat1,1)
  }
  
  if (name2 == 'nv'){
    dat2_clean = parse_nv(dat2,2)
    print('2 = nv')
  } else if (name2 == 'snif'){
    dat2_clean = parse_snif(dat2,2)
    print('2 = snif')
  } else if (name2 == 'cuteSV'){
    dat2_clean = parse_cuteSV(dat2,2)
    print('2 = cuteSV')
  }
  
  
  print('Parsing done')
  
  if (pos == TRUE){
    mylist <- pos_true(dat1,dat2,dat1_clean,dat2_clean,range)
  }
  
  within_pos_range <- mylist[[5]]
  
  
  
  
  if (end == TRUE){
    mylist <- end_true(dat1,dat2,dat1_clean,dat2_clean,range)
  }
  
  sim_dat1 <- mylist[[1]]
  unique_dat1 <- mylist[[2]]
  sim_dat2 <- mylist[[3]]
  unique_dat2 <- mylist[[4]]
  
  within_end_range <- mylist[[5]]
  
  
  
  
  if (pos == TRUE & end == TRUE){
    
    pos_dat1 = within_pos_range[which(within_pos_range$Dataset==1),]
    pos_dat2 = within_pos_range[which(within_pos_range$Dataset==2),]
    
    
    end_dat1 = within_end_range[which(within_end_range$Dataset==1),]
    end_dat2 = within_end_range[which(within_end_range$Dataset==2),]
    
    
    if (nrow(pos_dat1)>nrow(end_dat1)){
      index_sim_dat1 = match(end_dat1$Index,pos_dat1$Index)
      new_index_dat1 <- pos_dat1[index_sim_dat1,'Index']
      print(nrow(new_index_dat1))
    } else{
      index_sim_dat1 = match(pos_dat1$Index,end_dat1$Index)
      new_index_dat1 <- end_dat1[index_sim_dat1,'Index']
    } 
    
    
    if (nrow(pos_dat2)>nrow(end_dat2)){
      index_sim_dat2 = match(end_dat2$Index,pos_dat2$Index)
      new_index_dat2 <- pos_dat2[index_sim_dat2,'Index']
    } else{
      index_sim_dat2 = match(pos_dat2$Index,end_dat2$Index)
      new_index_dat2 <- end_dat2[index_sim_dat2,'Index']
    } 
    
    sim_dat1 <- dat1[new_index_dat1,]
    unique_dat1 <- dat1[-new_index_dat1,]
    sim_dat2 <- dat2[new_index_dat2,]
    unique_dat2 <- dat2[-new_index_dat2,]
    
    print('Position and end done')
    
  }
  
  write.table(sim_dat1, file="sim_dat1.csv", quote=F)
  write.table(unique_dat1, file="unique_dat1.csv", quote=F)
  write.table(sim_dat2, file="sim_dat2.csv", quote=F)
  write.table(unique_dat2, file="unique_dat2.csv", quote=F)
  
  print('Type done')
  
  print('End')
}


pos_true = function(dat1,dat2,dat1_clean,dat2_clean,range){
  dat1_position <- data.frame('Chromosome' = dat1_clean$Chromosome,'Type' = dat1_clean$Type, 'Index' = dat1_clean$Index, 'Position' = dat1_clean$Position,'Dataset' = dat1_clean$Dataset)
  dat2_position <- data.frame('Chromosome' = dat2_clean$Chromosome,'Type' = dat2_clean$Type,'Index' = dat2_clean$Index, 'Position' = dat2_clean$Position,'Dataset' = dat2_clean$Dataset)
  
  all_position = merge(dat1_position,dat2_position, all=TRUE)
  all_position <- all_position[order(all_position$Position),]
  all_position$Choromosome <- as.character(all_position$Chromosome)
  all_position$Position <- as.numeric(all_position$Position)
  
  within_pos_range = data.frame()
  last = nrow(all_position)-1
  
  for (i in 1:1242){
    for (j in i+1:1243){
      if (j == nrow(all_position))
        break
      
      if (all_position$Position[j]-all_position$Position[i]<=range){
        
        if (all_position$Chromosome[j] == all_position$Chromosome[i]){
          
          if (all_position$Dataset[j] !=all_position$Dataset[i]){
            
            within_pos_range <- merge(within_pos_range,merge(all_position[i,],all_position[j,],all=TRUE),all=TRUE)
            
          } 
        } else break
      } 
    }
  }
  
  
  
  if (sum(duplicated(within_pos_range))>=1){
    within_pos_range = within_pos_range[-which(duplicated(within_pos_range)==TRUE),]
  }
  
  index_dat1 <- within_pos_range[which(within_pos_range$Dataset==1),'Index']
  
  index_dat2 <- within_pos_range[which(within_pos_range$Dataset==2),'Index']
  
  
  sim_dat1 <- dat1[index_dat1,]
  unique_dat1 <- dat1[-index_dat1,]
  
  
  sim_dat2 <- dat2[index_dat2,]
  unique_dat2 <- dat2[-index_dat2,]
  
  print('Position done')
  
  return(list(sim_dat1,unique_dat1,sim_dat2,unique_dat2,within_pos_range))
}


end_true = function(dat1,dat2,dat1_clean,dat2_clean,range){
  dat1_end <- data.frame('Chromosome' = dat1_clean$Chromosome,'Type' = dat1_clean$Type,'Index' = dat1_clean$Index, 'End' = dat1_clean$End,'Dataset' = dat1_clean$Dataset)
  dat2_end <- data.frame('Chromosome' = dat2_clean$Chromosome,'Type' = dat2_clean$Type, 'Index' = dat2_clean$Index, 'End' = dat2_clean$End,'Dataset' = dat2_clean$Dataset)
  all_end = merge(dat1_end,dat2_end, all=TRUE)
  
  all_end <- all_end[order(all_end$End),]
  within_end_range = data.frame()
  
  
  for (i in 1:1242){
    for (j in i+1:1243){
      if (j == nrow(all_end))
        break
      if ((all_end$End[j]-all_end$End[i]<=range) & (all_end$Chromosome[j] == all_end$Chromosome[i])){
        if (all_end$Dataset[j] !=all_end$Dataset[i]){
          within_end_range <- merge(within_end_range,merge(all_end[i,],all_end[j,],all=TRUE),all=TRUE)
        } 
      } else break
    } 
  }
  
  if (sum(duplicated(within_end_range))>=1){
    within_end_range = within_end_range[-which(duplicated(within_end_range)==TRUE),]
  }
  
  index_dat1 <- within_end_range[which(within_end_range$Dataset==1),'Index']
  
  index_dat2 <- within_end_range[which(within_end_range$Dataset==2),'Index']
  
  
  sim_dat1 <- dat1[index_dat1,]
  unique_dat1 <- dat1[-index_dat1,]
  
  
  sim_dat2 <- dat2[index_dat2,]
  unique_dat2 <- dat2[-index_dat2,]
  
  print('End done')
  
  return(list(sim_dat1,unique_dat1,sim_dat2,unique_dat2,within_end_range))
}


parse_nv = function(dat1,dataset){
  dat1_clean <- dat1[,c('V1','V2','V5','V8')]
  names(dat1_clean) <- c('Chromosome','Position','Type','End')
  dat1_clean['End'] <- sub('.*END=','',dat1_clean[,'End'])
  dat1_clean['End'] <- sub(';.*','',dat1_clean[,'End'])
  dat1_clean$Type <- sub('.*<','',dat1_clean$Type)
  dat1_clean$Type <- sub('>.*','',dat1_clean$Type)
  dat1_clean[['End']] <- as.integer(dat1_clean[['End']])
  dat1_clean <-  data.frame(dat1_clean,c(rep(dataset,each=nrow(dat1_clean))),c(seq(from=1,to=nrow(dat1_clean))))
  names(dat1_clean)[5] <- 'Dataset'
  names(dat1_clean)[6] <- 'Index'
  return(dat1_clean)
}

parse_snif = function(dat2,dataset){
  dat2_clean <- dat2[,c('V1','V2','V8')]
  dat2_clean['Type'] <- sub('.*SVTYPE=','',dat2_clean[,'V8'])
  dat2_clean['Type'] <- sub(';.*','',dat2_clean[,'Type'])
  dat2_clean['End'] <- sub('.*END=','',dat2_clean[,'V8'])
  dat2_clean['End'] <- sub(';.*','',dat2_clean[,'End'])
  dat2_clean <-  subset(dat2_clean,select = -c(V8))
  names(dat2_clean)[1:2] <- c('Chromosome','Position') 
  dat2_clean[['End']] <- as.integer(dat2_clean[['End']])
  dat2_clean <-  data.frame(dat2_clean,c(rep(dataset,each=nrow(dat2_clean))),c(seq(from=1,to=nrow(dat2_clean))))
  names(dat2_clean)[5] <- 'Dataset'
  names(dat2_clean)[6] <- 'Index'
  if (sum(is.na(dat2_clean))!=0){
    row_index_na = which(is.na(dat2_clean),arr.ind=TRUE)[,1]
    dat2_clean[row_index_na,'End']<- dat2_clean[row_index_na,'Position']
  }
  return(dat2_clean)
}

parse_cuteSV = function(dat2,dataset){
  dat2_clean <- dat2[,c('V1','V2','V3','V8')]
  dat2_clean['Type'] <- sub('.*cuteSV.','',dat2_clean[,'V3'])
  dat2_clean['Type'] <- gsub("\\..*","",dat2_clean$Type)
  dat2_clean['End'] <- sub('.*END=','',dat2_clean[,'V8'])
  dat2_clean['End'] <- sub(';.*','',dat2_clean[,'End'])
  dat2_clean <-  subset(dat2_clean,select = -c(V8))
  dat2_clean <-  subset(dat2_clean,select = -c(V3))
  names(dat2_clean)[1:2] <- c('Chromosome','Position') 
  dat2_clean[['End']] <- as.integer(dat2_clean[['End']])
  dat2_clean <-  data.frame(dat2_clean,c(rep(dataset,each=nrow(dat2_clean))),c(seq(from=1,to=nrow(dat2_clean))))
  names(dat2_clean)[5] <- 'Dataset'
  names(dat2_clean)[6] <- 'Index'
  if (sum(is.na(dat2_clean))!=0){
    row_index_na = which(is.na(dat2_clean),arr.ind=TRUE)[,1]
    dat2_clean[row_index_na,'End']<- dat2_clean[row_index_na,'Position']
  }
  return(dat2_clean)
}

sim_function(dat1,dat2,name1,name2,pos,end,range,dir)
