#Step 1: search blast2GO output for transcripts containing certain annotation features (seqDesc)
#Step 2: extract the gene IDs from each species and determine their presence in MCL groups, singletons list, or neither
#step 3: extract expression information from edgeR outputs if available
#step 3: combine lines for Manch and white 1-to-1 ortholog members and output data
#step 4: collapse table to genes with at least some expression information

#Excell is then used to change names, transpose data, etc before feeding into the script: 
    #"~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/Loop_boxplot_CPM_by_trt_for_query.R"

#no required packages
#source script before running: 
  #source("~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/input_for_boxplots_from_annot_query.R")
#query should be entered in quotes and field corresponds to the column in the blas2GO output.
#So searching for "jasmon" in the Description field would be executed by typing: 
    #boxplots_from_query("jasmon","Description",output_folder) 
      #where output_folder should name a path of a folder that has already been created and excludes the final / 
        #e.g. output_folder<-"~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/automated_boxplot_outputs/jasmon"
#could further generalize the script to include entering paths of input files (such as orthoMCL and edgeR outputs) as inputs to the fuctions instead of specifiying them inside the function

boxplots_from_query<- function(query,field,output_folder){

  ####Step 1: search blast2GO output for transcripts containing certain annotation features (seqDesc)####  
  manch_b2g<-read.delim("~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/EBSeq/Manch_blast2go_table.txt")
  white_b2g<-read.delim("~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/EBSeq/white_blast2go_table.txt")
  
  ###old before adding the function to search any field
  #add column with logical flag indicating if the query string was found in the sequence description
  #manch_b2g$match<-grepl(pattern=query, manch_b2g$Description)
  #subset the dataframe to only those columns where the query was found and keep only column 1 containing manch gene IDs
  #keep_mgID<-subset(manch_b2g, match==TRUE, select = "SeqName")
  #rename column
  #colnames(keep_mgID)<-"Manch_ID"
  #repeat for white
  #white_b2g$match<-grepl(pattern=query, white_b2g$Description)
  #keep_wgID<-subset(white_b2g, match==TRUE, select = "SeqName")
  #colnames(keep_wgID)<-"white_ID"
  
  manch_b2g$match<-grepl(pattern=query, manch_b2g[,field])
  #subset the dataframe to only those columns where the query was found and keep only column 1 containing manch gene IDs
  keep_mgID<-subset(manch_b2g, match==TRUE, select = "SeqName")
  #rename column
  colnames(keep_mgID)<-"Manch_ID"
  
  white_b2g$match<-grepl(pattern=query, white_b2g[,field])
  #subset the dataframe to only those columns where the query was found and keep only column 1 containing manch gene IDs
  keep_wgID<-subset(white_b2g, match==TRUE, select = "SeqName")
  #rename column
  colnames(keep_wgID)<-"white_ID"
  
  
  ####Step 2: determine the presence of IDs in MCL groups, singletons list, or neither####  
  #import list of singletons
  singletonsMW <- read.table("~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/singletonsMW.txt", quote="\"", comment.char="")
  #trim species and isoform path information from gene identifiers
  #W|c112197_g1_i1:1-351(+) to c115197_g1_i1
  
  singletonsMW<- as.data.frame(gsub(".\\|c", "c", singletonsMW[,1]))
  #then the back ")" is also a metacharacter andany  character is denoted as "." with any number of items denoted as "*" so any number of any characters is ".*" 
  singletonsMW<- as.data.frame(gsub(":.*\\)", "", singletonsMW[,1]))
  #split by species. white is rows 1:10199, Manch is 10200:19156
  singletons_W<-as.data.frame(singletonsMW[1:10199,]) 
  colnames(singletons_W)<-"G"
  singletons_M<-as.data.frame(singletonsMW[10200:19156,]) 
  colnames(singletons_M)<-"G"
  
  #Add flag if gene ID appears in singletons list. example a$flag[a$id %in% temp] <- 1
  keep_mgID$singleton<-0
  keep_mgID$singleton[keep_mgID$Manch_ID %in% singletons_M$G]<-1
  keep_wgID$singleton<-0
  keep_wgID$singleton[keep_wgID$white_ID %in% singletons_W$G]<-1
  
  #import orthoMCL groups
  groupsMW <- read.csv("~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/groupsMW.csv")
  #subset groupsMW to 1-to-1 orthologs
  orthologs<-as.data.frame(groupsMW[(groupsMW$Manch_count=="1" & groupsMW$White_count=="1"),c(1,5,6)])
    #modify orthologs table to put species information in column names and remove from column contents
  colnames(orthologs)<- c("MCL_Group","Manch_ID", "white_ID")
  #trim species and isoform information from gene identifiers in both columns
  #W|c112197_g1_i1:1-351(+) to c115197_g1_i1
  #first the front "W|"  | is a metacharacter so must be preceded by \\
  orthologs<- as.data.frame(apply(orthologs,2,function(y) gsub(".\\|c", "c", y)))
  #then the back ")" is also a metacharacter andany  character is denoted as "." with any number of items denoted as "*" so any number of any characters is ".*" 
  orthologs<- as.data.frame(apply(orthologs,2,function(y) gsub(":.*\\)", "", y)))
  
  #Add flag if gene ID appears in ortholog list. example a$flag[a$id %in% temp] <- 1
  keep_mgID$ortholog<-0
  keep_mgID$ortholog[keep_mgID$Manch_ID %in% orthologs$Manch_ID]<-1
  keep_wgID$ortholog<-0
  keep_wgID$ortholog[keep_wgID$white_ID %in% orthologs$white_ID]<-1
  
  #add MCL_group information for orthologs
  keep_mgID <- merge(keep_mgID, orthologs, by.x = "Manch_ID", by.y = "Manch_ID", all.x=TRUE)
  keep_wgID <- merge(keep_wgID, orthologs, by.x = "white_ID", by.y = "white_ID", all.x=TRUE)
  
  
  #create subset of orthoMCL groups with just ortholog-paralog groups
  ortho_para= groupsMW[(groupsMW$Manch_count != 0 & groupsMW$White_count != 0 & groupsMW$Total_count > 2),]
  
  #NOTE: may also be useful to identify the composition of the ortho-para groups e.g. 2 manch 3 white. can get this information from "Manch_Count" and "White_Count" columns of ortho_para df 
  
  #combine all columns with gene identifiers into a single column
  member_1df=as.data.frame(ortho_para[,(c(1,5))])
  colnames(member_1df)[2]<-"G"
  member_2df=as.data.frame(ortho_para[,(c(1,6))])
  colnames(member_2df)[2]<-"G"
  member_3df=as.data.frame(ortho_para[,(c(1,7))])
  colnames(member_3df)[2]<-"G"
  member_4df=as.data.frame(ortho_para[,(c(1,8))])
  colnames(member_4df)[2]<-"G"
  member_5df=as.data.frame(ortho_para[,(c(1,9))])
  colnames(member_5df)[2]<-"G"
  member_6df=as.data.frame(ortho_para[,(c(1,10))])
  colnames(member_6df)[2]<-"G"
  member_7df=as.data.frame(ortho_para[,(c(1,11))])
  colnames(member_7df)[2]<-"G"
  member_8df=as.data.frame(ortho_para[,(c(1,12))])
  colnames(member_8df)[2]<-"G"
  member_9df=as.data.frame(ortho_para[,(c(1,13))])
  colnames(member_9df)[2]<-"G"
  member_10df=as.data.frame(ortho_para[,(c(1,14))])
  colnames(member_10df)[2]<-"G"
  member_11df=as.data.frame(ortho_para[,(c(1,15))])
  colnames(member_11df)[2]<-"G"
  member_12df=as.data.frame(ortho_para[,(c(1,16))])
  colnames(member_12df)[2]<-"G"
  member_13df=as.data.frame(ortho_para[,(c(1,17))])
  colnames(member_13df)[2]<-"G"
  member_14df=as.data.frame(ortho_para[,(c(1,18))])
  colnames(member_14df)[2]<-"G"
  member_15df=as.data.frame(ortho_para[,(c(1,19))])
  colnames(member_15df)[2]<-"G"
  member_16df=as.data.frame(ortho_para[,(c(1,20))])
  colnames(member_16df)[2]<-"G"
  member_17df=as.data.frame(ortho_para[,(c(1,21))])
  colnames(member_17df)[2]<-"G"
  member_18df=as.data.frame(ortho_para[,(c(1,22))])
  colnames(member_18df)[2]<-"G"
  member_19df=as.data.frame(ortho_para[,(c(1,23))])
  colnames(member_19df)[2]<-"G"
  member_20df=as.data.frame(ortho_para[,(c(1,24))])
  colnames(member_20df)[2]<-"G"
  member_21df=as.data.frame(ortho_para[,(c(1,25))])
  colnames(member_21df)[2]<-"G"
  member_22df=as.data.frame(ortho_para[,(c(1,26))])
  colnames(member_22df)[2]<-"G"
  member_23df=as.data.frame(ortho_para[,(c(1,27))])
  colnames(member_23df)[2]<-"G"
  member_24df=as.data.frame(ortho_para[,(c(1,28))])
  colnames(member_24df)[2]<-"G"
  member_25df=as.data.frame(ortho_para[,(c(1,29))])
  colnames(member_25df)[2]<-"G"
  member_26df=as.data.frame(ortho_para[,(c(1,30))])
  colnames(member_26df)[2]<-"G"
  member_27df=as.data.frame(ortho_para[,(c(1,31))])
  colnames(member_27df)[2]<-"G"
  member_28df=as.data.frame(ortho_para[,(c(1,32))])
  colnames(member_28df)[2]<-"G"
  member_29df=as.data.frame(ortho_para[,(c(1,33))])
  colnames(member_29df)[2]<-"G"
  member_30df=as.data.frame(ortho_para[,(c(1,34))])
  colnames(member_30df)[2]<-"G"
  member_31df=as.data.frame(ortho_para[,(c(1,35))])
  colnames(member_31df)[2]<-"G"
  member_32df=as.data.frame(ortho_para[,(c(1,36))])
  colnames(member_32df)[2]<-"G"
  member_33df=as.data.frame(ortho_para[,(c(1,37))])
  colnames(member_33df)[2]<-"G"
  member_34df=as.data.frame(ortho_para[,(c(1,38))])
  colnames(member_34df)[2]<-"G"
  member_35df=as.data.frame(ortho_para[,(c(1,39))])
  colnames(member_35df)[2]<-"G"
  member_36df=as.data.frame(ortho_para[,(c(1,40))])
  colnames(member_36df)[2]<-"G"
  member_37df=as.data.frame(ortho_para[,(c(1,41))])
  colnames(member_37df)[2]<-"G"
  member_38df=as.data.frame(ortho_para[,(c(1,42))])
  colnames(member_38df)[2]<-"G"
  member_39df=as.data.frame(ortho_para[,(c(1,43))])
  colnames(member_39df)[2]<-"G"
  member_40df=as.data.frame(ortho_para[,(c(1,44))])
  colnames(member_40df)[2]<-"G"
  member_41df=as.data.frame(ortho_para[,(c(1,45))])
  colnames(member_41df)[2]<-"G"
  member_42df=as.data.frame(ortho_para[,(c(1,46))])
  colnames(member_42df)[2]<-"G"
  member_43df=as.data.frame(ortho_para[,(c(1,47))])
  colnames(member_43df)[2]<-"G"
  member_44df=as.data.frame(ortho_para[,(c(1,48))])
  colnames(member_44df)[2]<-"G"
  member_45df=as.data.frame(ortho_para[,(c(1,49))])
  colnames(member_45df)[2]<-"G"
  member_46df=as.data.frame(ortho_para[,(c(1,50))])
  colnames(member_46df)[2]<-"G"
  member_47df=as.data.frame(ortho_para[,(c(1,51))])
  colnames(member_47df)[2]<-"G"
  member_48df=as.data.frame(ortho_para[,(c(1,52))])
  colnames(member_48df)[2]<-"G"
  member_49df=as.data.frame(ortho_para[,(c(1,53))])
  colnames(member_49df)[2]<-"G"
  member_50df=as.data.frame(ortho_para[,(c(1,54))])
  colnames(member_50df)[2]<-"G"
  member_51df=as.data.frame(ortho_para[,(c(1,55))])
  colnames(member_51df)[2]<-"G"
  member_52df=as.data.frame(ortho_para[,(c(1,56))])
  colnames(member_52df)[2]<-"G"
  member_53df=as.data.frame(ortho_para[,(c(1,57))])
  colnames(member_53df)[2]<-"G"
  member_54df=as.data.frame(ortho_para[,(c(1,58))])
  colnames(member_54df)[2]<-"G"
  member_55df=as.data.frame(ortho_para[,(c(1,59))])
  colnames(member_55df)[2]<-"G"
  member_56df=as.data.frame(ortho_para[,(c(1,60))])
  colnames(member_56df)[2]<-"G"
  member_57df=as.data.frame(ortho_para[,(c(1,61))])
  colnames(member_57df)[2]<-"G"
  member_58df=as.data.frame(ortho_para[,(c(1,62))])
  colnames(member_58df)[2]<-"G"
  member_59df=as.data.frame(ortho_para[(c(1,63))])
  colnames(member_59df)[2]<-"G"
  ortho_para_comb<- rbind(member_1df,member_2df,member_3df,member_4df,member_5df,member_6df,member_7df,member_8df,member_9df,member_10df,member_11df,member_12df,member_13df,member_14df,member_15df,member_16df,member_17df,member_18df,member_19df,member_20df,member_21df,member_22df,member_23df,member_24df,member_25df,member_26df,member_27df,member_28df,member_29df,member_30df,member_31df,member_32df,member_33df,member_34df,member_35df,member_36df,member_37df,member_38df,member_39df,member_40df,member_41df,member_42df,member_43df,member_44df,member_45df,member_46df,member_47df,member_48df,member_49df,member_50df,member_51df,member_52df,member_53df,member_54df,member_55df,member_56df,member_57df,member_58df,member_59df)
  colnames(ortho_para_comb)<-c("MCL_Group","G")
  #remove blanks
  ortho_para_members= as.data.frame(ortho_para_comb[(ortho_para_comb$G !=""& ortho_para_comb$G !="NA"),])
  
  #need to split speces before removing species identifiers
  ortho_para_members<- as.data.frame(ortho_para_members[order(ortho_para_members$G),])
  #manchurian genes are rows 1:4827  White are 4828:9392
  manch_ortho_para_members<- as.data.frame(ortho_para_members[c(1:4827),])
  white_ortho_para_members<- as.data.frame(ortho_para_members[c(4828:9392),])
  
  #trim species and isoform information from gene identifiers
  #trim species and isoform information from gene identifiers in both columns
  #W|c112197_g1_i1:1-351(+) to c115197_g1_i1
  #first the front "W|"  | is a metacharacter so must be preceded by \\
  manch_ortho_para_members<- as.data.frame(apply(manch_ortho_para_members,2,function(y) gsub(".\\|c", "c", y)))
  white_ortho_para_members<- as.data.frame(apply(white_ortho_para_members,2,function(y) gsub(".\\|c", "c", y)))
  #then the back ")" is also a metacharacter andany  character is denoted as "." with any number of items denoted as "*" so any number of any characters is ".*" 
  manch_ortho_para_members<- as.data.frame(apply(manch_ortho_para_members,2,function(y) gsub(":.*\\)", "", y)))
  white_ortho_para_members<- as.data.frame(apply(white_ortho_para_members,2,function(y) gsub(":.*\\)", "", y)))

  #Add flag if gene ID appears in ortho-para list. example a$flag[a$id %in% temp] <- 1
  keep_mgID$ortho_para<-0
  keep_mgID$ortho_para[keep_mgID$Manch_ID %in% manch_ortho_para_members$G]<-1
  keep_wgID$ortho_para<-0
  keep_wgID$ortho_para[keep_wgID$white_ID %in% white_ortho_para_members$G]<-1
 
  #add MCL_group information for ortho_para members.
  keep_mgID <- merge(keep_mgID, manch_ortho_para_members, by.x = "Manch_ID", by.y = "G", all.x=TRUE)
  #combine MCL_Group columns into one
  #example TEST$UNIT[is.na(TEST$UNIT)] <- TEST$STATUS[is.na(TEST$UNIT)]
  #however, the columns cannot be factors, so I am changing them to characters.  May need to change them back to factors later
  keep_mgID$MCL_Group.x<- as.character(keep_mgID$MCL_Group.x)
  keep_mgID$MCL_Group.y<- as.character(keep_mgID$MCL_Group.y)
  keep_mgID$MCL_Group.x[is.na(keep_mgID$MCL_Group.x)]<- keep_mgID$MCL_Group.y[is.na(keep_mgID$MCL_Group.x)]
  
  #add MCL_group information for ortho_para members.
  keep_wgID <- merge(keep_wgID, white_ortho_para_members, by.x = "white_ID", by.y = "G", all.x=TRUE)
  #combine MCL_Group columns into one
  #example TEST$UNIT[is.na(TEST$UNIT)] <- TEST$STATUS[is.na(TEST$UNIT)]
  #however, the columns cannot be factors, so I am changing them to characters.  May need to change them back to factors later
  keep_wgID$MCL_Group.x<- as.character(keep_wgID$MCL_Group.x)
  keep_wgID$MCL_Group.y<- as.character(keep_wgID$MCL_Group.y)
  keep_wgID$MCL_Group.x[is.na(keep_wgID$MCL_Group.x)]<- keep_wgID$MCL_Group.y[is.na(keep_wgID$MCL_Group.x)]
  

  #create subset of orthoMCL groups with just paralog groups
  #white
  paralogs_W= groupsMW[(groupsMW$Manch_count == 0),]
  #combine all columns with gene identifiers into a single column
  member_1df=as.data.frame(paralogs_W[,(c(1,5))])
  colnames(member_1df)[2]<-"G"
  member_2df=as.data.frame(paralogs_W[,(c(1,6))])
  colnames(member_2df)[2]<-"G"
  member_3df=as.data.frame(paralogs_W[,(c(1,7))])
  colnames(member_3df)[2]<-"G"
  member_4df=as.data.frame(paralogs_W[,(c(1,8))])
  colnames(member_4df)[2]<-"G"
  member_5df=as.data.frame(paralogs_W[,(c(1,9))])
  colnames(member_5df)[2]<-"G"
  member_6df=as.data.frame(paralogs_W[,(c(1,10))])
  colnames(member_6df)[2]<-"G"
  member_7df=as.data.frame(paralogs_W[,(c(1,11))])
  colnames(member_7df)[2]<-"G"
  paralogs_W_comb<- rbind(member_1df,member_2df,member_3df,member_4df,member_5df,member_6df,member_7df)

  #Manchurian 
  paralogs_M= groupsMW[(groupsMW$White_count == 0),]
  #combine all columns with gene identifiers into a single column
  member_1df=as.data.frame(paralogs_M[,(c(1,5))])
  colnames(member_1df)[2]<-"G"
  member_2df=as.data.frame(paralogs_M[,(c(1,6))])
  colnames(member_2df)[2]<-"G"
  member_3df=as.data.frame(paralogs_M[,(c(1,7))])
  colnames(member_3df)[2]<-"G"
  member_4df=as.data.frame(paralogs_M[,(c(1,8))])
  colnames(member_4df)[2]<-"G"
  member_5df=as.data.frame(paralogs_M[,(c(1,9))])
  colnames(member_5df)[2]<-"G"
  member_6df=as.data.frame(paralogs_M[,(c(1,10))])
  colnames(member_6df)[2]<-"G"
  paralogs_M_comb<- rbind(member_1df,member_2df,member_3df,member_4df,member_5df,member_6df)
  
  #trim species and isoform information from gene identifiers
  #trim species and isoform information from gene identifiers in both columns
  #W|c112197_g1_i1:1-351(+) to c115197_g1_i1
  #first the front "W|"  | is a metacharacter so must be preceded by \\
  paralogs_M_comb<- as.data.frame(apply(paralogs_M_comb,2,function(y) gsub(".\\|c", "c", y)))
  paralogs_W_comb<- as.data.frame(apply(paralogs_W_comb,2,function(y) gsub(".\\|c", "c", y)))
  #then the back ")" is also a metacharacter andany  character is denoted as "." with any number of items denoted as "*" so any number of any characters is ".*" 
  paralogs_M_comb<- as.data.frame(apply(paralogs_M_comb,2,function(y) gsub(":.*\\)", "", y)))
  paralogs_W_comb<- as.data.frame(apply(paralogs_W_comb,2,function(y) gsub(":.*\\)", "", y)))

  #Add flag if gene ID appears in para list. example a$flag[a$id %in% temp] <- 1
  keep_mgID$para<-0
  keep_mgID$para[keep_mgID$Manch_ID %in% paralogs_M_comb$G]<-1
  keep_wgID$para<-0
  keep_wgID$para[keep_wgID$white_ID %in% paralogs_W_comb$G]<-1
  
  #add MCL_group information for para members.
  keep_mgID <- merge(keep_mgID, paralogs_M_comb, by.x = "Manch_ID", by.y = "G", all.x=TRUE)
  #combine MCL_Group columns into one
  #example TEST$UNIT[is.na(TEST$UNIT)] <- TEST$STATUS[is.na(TEST$UNIT)]
  #however, the columns cannot be factors, so I am changing them to characters.  May need to change them back to factors later
  keep_mgID$MCL_group<- as.character(keep_mgID$MCL_group)
  keep_mgID$MCL_Group.x[is.na(keep_mgID$MCL_Group.x)]<- keep_mgID$MCL_group[is.na(keep_mgID$MCL_Group.x)]
  
  #add MCL_group information for ortho_para members.
  keep_wgID <- merge(keep_wgID, paralogs_W_comb, by.x = "white_ID", by.y = "G", all.x=TRUE)
  #combine MCL_Group columns into one
  #example TEST$UNIT[is.na(TEST$UNIT)] <- TEST$STATUS[is.na(TEST$UNIT)]
  #however, the columns cannot be factors, so I am changing them to characters.  May need to change them back to factors later
  keep_wgID$MCL_group<- as.character(keep_wgID$MCL_group)
  keep_wgID$MCL_Group.x[is.na(keep_wgID$MCL_Group.x)]<- keep_wgID$MCL_group[is.na(keep_wgID$MCL_Group.x)]
  
  #remove extra MCL_group columns and sort columns
  keep_mgID <-keep_mgID[,c(4,1,5,3,6,8,2)]
  colnames(keep_mgID)[1]<-"MCL_Group"
  keep_wgID <-keep_wgID[,c(4,1,5,3,6,8,2)]
  colnames(keep_wgID)[1]<-"MCL_Group"
  
  ####step 3: extract expression information from edgeR outputs if available####
  manch_QLF_annot <- read.delim("~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/EBSeq/manch_QLF_annot.txt", na.strings=c("NA","-","---NA---"))
  white_QLF_annot <- read.delim("~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/EBSeq/white_QLF_annot.txt", na.strings=c("NA","-","---NA---"))
  
  #identifying columns to keep from edgeR outputs
  keep_col<-c(1,2,6,7,11,12,16,17:31,37)
  
  manch_QLF_annot<-manch_QLF_annot[,keep_col]
  #rename column headers of one species before importing the other
  #example colnames(df) <- paste("prefix", colnames(df), sep = "_") would give "prefix_<existing_name>"
  colnames(manch_QLF_annot) <- paste("M", colnames(manch_QLF_annot), sep = "")
  colnames(manch_QLF_annot)[6:20]<-c("MlogFC_WE","MFDR_WE","MCH1","MCL1","MCL2","MCH2","MWL1","MWH1","MWL2","MWH2","MWH3","MWL3","MEL1","MEL2","MEH3")
  
  #combine manch expression and groups MW dataset by manch gene ID columns in each dataframe, keep all entries for keep_mgID
  query_m_expression_QLF <- merge(keep_mgID, manch_QLF_annot, by.x = "Manch_ID", by.y = "MSeqName", all.x=TRUE)
  #reorder columns
  query_m_expression_QLF <-query_m_expression_QLF[,c(2,1,4:7,9:29,3)]
  
  #repeat for white
  white_QLF_annot<-white_QLF_annot[,keep_col]
  #rename column headers of one species before importing the other
  #example colnames(df) <- paste("prefix", colnames(df), sep = "_") would give "prefix_<existing_name>"
  colnames(white_QLF_annot) <- paste("W", colnames(white_QLF_annot), sep = "")
  colnames(white_QLF_annot)[6:20]<-c("WlogFC_WE","WFDR_WE","WCL1","WCH1","WCH2","WCL2","WWH1","WWL2","WWH2","WWH3","WWL3","WEL1","WEH1","WEL3","WEH2")
  
  #combine manch expression and groups MW dataset by white gene ID columns in each dataframe, keep all entries for keep_wgID
  query_w_expression_QLF <- merge(keep_wgID, white_QLF_annot, by.x = "white_ID", by.y = "WSeqName", all.x=TRUE)
  #reorder columns
  query_w_expression_QLF <-query_w_expression_QLF[,c(2,1,4:7,9:29,3)]
  
  #####step 3: combine lines for Manch and white 1-to-1 ortholog members and output data
  query_mw_expression_QLF <- merge(query_m_expression_QLF, query_w_expression_QLF, by.x = "white_ID", by.y = "white_ID", all.x=TRUE, all.y=TRUE)
   
  #collapse MCL_group columns, ortho, ortho-para,para,singleton,columns
  #combine MCL_Group columns into one
  #example TEST$UNIT[is.na(TEST$UNIT)] <- TEST$STATUS[is.na(TEST$UNIT)]
  #however, the columns cannot be factors, so I am changing them to characters.  May need to change them back to factors later
  query_mw_expression_QLF$MCL_Group.x<- as.character(query_mw_expression_QLF$MCL_Group.x)
  query_mw_expression_QLF$MCL_Group.y<- as.character(query_mw_expression_QLF$MCL_Group.y)
  query_mw_expression_QLF$MCL_Group.x[is.na(query_mw_expression_QLF$MCL_Group.x)]<- query_mw_expression_QLF$MCL_Group.y[is.na(query_mw_expression_QLF$MCL_Group.x)]
  #combine ortholog flag columns into one
  query_mw_expression_QLF$ortholog.x<- as.character(query_mw_expression_QLF$ortholog.x)
  query_mw_expression_QLF$ortholog.y<- as.character(query_mw_expression_QLF$ortholog.y)
  query_mw_expression_QLF$ortholog.x[is.na(query_mw_expression_QLF$ortholog.x)]<- query_mw_expression_QLF$ortholog.y[is.na(query_mw_expression_QLF$ortholog.x)]
  #combine ortho_para flag columns into one
  query_mw_expression_QLF$ortho_para.x<- as.character(query_mw_expression_QLF$ortho_para.x)
  query_mw_expression_QLF$ortho_para.y<- as.character(query_mw_expression_QLF$ortho_para.y)
  query_mw_expression_QLF$ortho_para.x[is.na(query_mw_expression_QLF$ortho_para.x)]<- query_mw_expression_QLF$ortho_para.y[is.na(query_mw_expression_QLF$ortho_para.x)]
  #combine para flag columns into one
  query_mw_expression_QLF$para.x<- as.character(query_mw_expression_QLF$para.x)
  query_mw_expression_QLF$para.y<- as.character(query_mw_expression_QLF$para.y)
  query_mw_expression_QLF$para.x[is.na(query_mw_expression_QLF$para.x)]<- query_mw_expression_QLF$para.y[is.na(query_mw_expression_QLF$para.x)]
  #combine singleton flag columns into one
  query_mw_expression_QLF$singleton.x<- as.character(query_mw_expression_QLF$singleton.x)
  query_mw_expression_QLF$singleton.y<- as.character(query_mw_expression_QLF$singleton.y)
  query_mw_expression_QLF$singleton.x[is.na(query_mw_expression_QLF$singleton.x)]<- query_mw_expression_QLF$singleton.y[is.na(query_mw_expression_QLF$singleton.x)]
  
  #reorder columns
  query_mw_expression_QLF <-query_mw_expression_QLF[,c(2:28,1,34:54)]
  #remove suffix
  colnames(query_mw_expression_QLF)[1:6]<-c("MCL_Group","Manch_ID","ortholog","ortho_para","para","singleton")
  
  write.csv(query_mw_expression_QLF, file = file.path(output_folder,"query_mw_expression_QLF.csv"), row.names = FALSE)
  
  ##step 4: collapse table to genes with at least some expression information
  query_mw_expression_QLF_collapse= as.data.frame(query_mw_expression_QLF[(!is.na(query_mw_expression_QLF$MlogFC_CW)|(!is.na(query_mw_expression_QLF$WlogFC_CW))),])
  
  write.csv(query_mw_expression_QLF_collapse, file = file.path(output_folder,"query_mw_expression_QLF_collapse.csv"), row.names = FALSE)
  
print("~~open in Excel to double-check appropriatness of query,sort by ortho(ascending):ortho_para:para:singleton:MCL_group, add ortho_para prefixes to MCL_group names, add seqDesc suffix to MCL_group names, trim unnecessary columns, transpose data, change 'MCL_group' to 'Sample', and add 'TRT' column. Save in output_folder as 'query_mw_expression_transpose.csv' Then use boxplot script to create boxplots~~")
}