###################################################################################
###################### Used for WOE and variable selection ########################
###################################################################################

#Read in libraries
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(gridExtra)){install.packages("gridExtra")}
if(!require(captioner)){install.packages("captioner")}
if(!require(GGally)){install.packages("GGally")}
if(!require(pander)){install.packages("pander")}
if(!require(Rprofet)){install.packages("Rprofet")}
if(!require(ROCit)){install.packages("ROCit")}
if(!require(kableExtra)){install.packages("kableExtra")}
if(!require(readxl)){install.packages("readxl")}
if(!require(Information)){install.packages("Information")}
if(!require(ClustOfVar)){install.packages("ClustOfVar")}
if(!require(corrplot)){install.packages("corrplot")}
if(!require(cluster)){install.packages("cluster")}
if(!require(factoextra)){install.packages("factoextra")}

library(tidyverse)
library(gridExtra)
library(captioner)
library(pander)
library(kableExtra)
library(Rprofet)
library(ROCit)
library(GGally)
library(readxl)
library(Information)
library(grid)
library(ClustOfVar)
library(corrplot)
library(cluster)
library(factoextra)
library(NbClust)
library(plyr)


###################################################################################
###################### Useful Functions  ##########################################
###################################################################################

#this is to transform data into WOE values based on bins

new_data_to_woe <- function(WOE_lists, new_data){
  j=0
  woe_data = new_data
  pb <-  winProgressBar(title = "converting data...", min = 0,
                        max = length(WOE_lists$Tables), width = 300)
  for(i in WOE_lists$Tables){
    j = j + 1
    setWinProgressBar(pb, j, title=paste0("converting is ", round(j/length(WOE_lists$Tables)*100, 0),
                                          "% done"))
    if(colnames(i)[1] %in% colnames(woe_data)){
      if(grepl("^\\[", as.character(i[1,1]))){
        for(t in 1:nrow(i)){
          
          unlisted_split_bin <- unlist(strsplit(as.character(i[t,1]), ",", fixed=T))
          no_par <- gsub("]","", unlisted_split_bin, fixed=T)
          no_brak <- gsub("[","", no_par, fixed=T)
          bin_numeric = as.numeric(no_brak)
          woe_data[woe_data[,colnames(i)[1]] >= bin_numeric[1] & woe_data[,colnames(i)[1]] <= bin_numeric[2], colnames(i)[1]]  = i[t,4]
        }
      }
      else if(i[1,1] == "NA"){
        for(t in 1:nrow(i)){
          if(t == 1){
            woe_data[is.na(woe_data[,colnames(i)[1]]),colnames(i)[1]] = i[t,4]
          }
          else{
            unlisted_split_bin <- unlist(strsplit(as.character(i[t,1]), ",", fixed=T))
            no_par <- gsub("]","", unlisted_split_bin, fixed=T)
            no_brak <- gsub("[","", no_par, fixed=T)
            bin_numeric = as.numeric(no_brak)
            woe_data[woe_data[,colnames(i)[1]] >= bin_numeric[1] & woe_data[,colnames(i)[1]] <= bin_numeric[2], colnames(i)[1]]  = i[t,4]
          }
        }
      }
      else{
        for(t in 1:nrow(i)){
          woe_data[woe_data[,colnames(i)[1]]==i[t,1],colnames(i)[1]] = as.character(i[t,4])
        }
        woe_data[,colnames(i)[1]] = as.numeric(unlist(woe_data[,colnames(i)[1]]))
      }  
    }
    else{
      message("No Vars Found")
    }
  }
  close(pb)
  return(woe_data)
  
}

#perform kmeans clustering on WOE data

WOEClust_kmeans<-function(WOE_df,id,target,num_clusts){
  
  dat = WOE_df
  
  ## Adjusting ID, target inputs to be column numbers
  if(is.character(id) == TRUE){id = match(id,colnames(dat))}
  if(is.character(target) == TRUE){target = match(target,colnames(dat))}
  
  ## Removes ID, target columns. Performs variable clustering and generates data frame containing variable and their cluster number.
  dat2<-dat[,-c(1,2)]
  Groups<-ClustOfVar::kmeansvar(X.quanti =dat2[,1:ncol(dat2)], X.quali = NULL, num_clusts,
                                iter.max = 150, nstart = 1, matsim = FALSE)
  Memberships<-data.frame(Groups$cluster)
  Memberships<-cbind(rownames(Memberships),Memberships)
  rownames(Memberships) <- NULL
  names(Memberships)<-c("Variable","Group")
  
  ## Storing the column numbers
  varcol = colnames(dat2)
  varcol = match(varcol,colnames(dat))
  
  ## Function from WOEFun3 that calculates IV for each variables
  IVFun1<-function(dat,idcol,targetcol,varcol){
    dat2<-dat[,c(idcol,targetcol,varcol)]
    dat2$Bins<-as.character(dat2[,3])
    
    NumBad<-aggregate(dat2[,2]~Bins,data=dat2,FUN=sum)
    NumGood<-aggregate((ifelse(dat2[,2]==1,0,1))~Bins,data=dat2,FUN=sum)
    
    
    IVF1_i<-(NumBad[,2]/sum(NumBad[,2]))-(NumGood[,2]/sum(NumGood[,2]))
    IVF2_i<-log((NumBad[,2]/sum(NumBad[,2]))/(NumGood[,2]/sum(NumGood[,2])))
    IVF3_i<-log(((NumBad[,2]+.5)/(sum(NumBad[,2])+.5))/((NumGood[,2]+.5)/(sum(NumGood[,2])+.5)))
    IVF23_i<-ifelse(IVF2_i==-Inf|IVF2_i==Inf,IVF3_i,IVF2_i)
    IVF_i<-IVF1_i*IVF23_i
    
    IVF<-sum(IVF_i)
    IVF2<-c(colnames(dat2[3]),IVF)
    
    return(c(IVF2))
  }
  
  ## Wrapper function for IVFun1
  IVFun2<-function(varcol){IVFun1(dat=dat,idcol=id,targetcol=target,varcol)}
  
  ## Calling the IVFun functions to calculate the IV for each variable
  IV_and_Var<-lapply(varcol, IVFun2)
  IV_and_Var2<-data.frame(Variable=sapply(IV_and_Var,"[",1),
                          IV=as.numeric(sapply(IV_and_Var,"[",2)))
  
  ## Managing and organizing the data frame to order it by cluster number and then IV
  Memberships2 = merge(Memberships,IV_and_Var2, By = "Variable")
  
  KMEANSGroups = Memberships2[order(Memberships2$Group,-Memberships2$IV),]
  
  KMEANSGroups$Variable = as.character(KMEANSGroups$Variable)
  
  return(KMEANSGroups)
}

#perform heirarical clustering with a default of wards method

WOEClust_hclust<-function(WOE_df,id,target,num_clusts,method='ward.D', IV){
  
  dat = WOE_df
  
  ## Adjusting ID, target inputs to be column numbers
  if(is.character(id) == TRUE){id = match(id,colnames(dat))}
  if(is.character(target) == TRUE){target = match(target,colnames(dat))}
  
  woe_trans = t(dat[,-c(id,target)])
  dist.probes = dist(woe_trans)
  
  probes.complete=hclust(dist.probes, method = method)
  
  groups = cutree(probes.complete, k=num_clusts)
  
  Memberships<-data.frame(groups)
  Memberships<-cbind(rownames(Memberships),Memberships)
  rownames(Memberships) <- NULL
  names(Memberships)<-c("Variable","Group")
  Memberships$Variable = as.character(Memberships$Variable)
  
  infoVal = IV
  
  Memberships = merge(x = Memberships, y = infoVal, by = "Variable", all.x = TRUE)
  
  Memberships = Memberships[order(Memberships$Group,-Memberships$IV),]
  
  Memberships$Variable = as.character(Memberships$Variable)
  
  
  return(Memberships)
}

###################################################################################
###################### Beginning of Analysis  #####################################
###################################################################################


#read in our data
data = read_excel("final_model_data.xlsx", sheet="Row 1 Observations")

#turn 1/0 code into yes and no
data$`Closure Reason Precense` <- ifelse(data$`Closure Reason Precense` == 1, "Yes", "No")
data$`Delinquent and Overlimit` <- ifelse(data$`Delinquent and Overlimit` == 1, "Yes", "No")

#remove the id and undeeded variables to be binned
model_data = data[,-c(1,28,29)]

#Use the information package function to create binned, WOE and IV values

infoTables <- create_infotables(data = model_data,
                                y = "Bad",
                                bins = 10,
                                parallel = T, ncore = 15)

#get the IV value from the object
IV_table = infoTables$Summary

#removing the variables with IV < 0.1
cut_data = model_data[,!colnames(model_data) %in% IV_table$Variable[21:25]]

#testing out various binning sizes on the cut data

infoTables_20 <- create_infotables(data = cut_data,
                                   y = "Bad",
                                   bins = 20,
                                   parallel = T, ncore = 8)

infoTables_10 <- create_infotables(data = cut_data,
                                   y = "Bad",
                                   bins = 10,
                                   parallel = T, ncore = 8)

infoTables_5 <- create_infotables(data = cut_data,
                                  y = "Bad",
                                  bins = 5,
                                  parallel = T, ncore = 8)

#make an example plot of the WOE transformed data
plot_infotables(infoTables_10, infoTables_10$Summary$Variable[1])+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

#create a list of all of the plots of WOE for all variables
bin_10 = vector(mode = "list", length = 20)
bin_20 = vector(mode = "list", length = 20)
bin_5 = vector(mode = "list", length = 20)

for (i in 1:20){
  bin_10[[i]] = plot_infotables(infoTables_10, infoTables_10$Summary$Variable[i])+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  bin_20[[i]] = plot_infotables(infoTables_20, infoTables_20$Summary$Variable[i])+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  bin_5[[i]] = plot_infotables(infoTables_5, infoTables_5$Summary$Variable[i])+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

#plot all of the 10 binned data into one plot
grid.arrange(bin_10[[1]], bin_10[[2]], bin_10[[3]], bin_10[[4]],bin_10[[5]], bin_10[[6]], bin_10[[7]], bin_10[[8]],bin_10[[9]], bin_10[[10]], bin_10[[11]], bin_10[[12]],bin_10[[13]], bin_10[[14]], bin_10[[15]], bin_10[[16]],bin_10[[17]], bin_10[[18]], bin_10[[19]], bin_10[[20]])

#use the above defined function to transform our data into WOE data
woe_data = new_data_to_woe(infoTables_10, cut_data)
woe_data_ided = cbind(data$DebtDimId, woe_data)
colnames(woe_data_ided)[1] = "ID"

#make a correlation plot
corrplot(cor(woe_data_ided[,-c(1,2)]), order="hclust", hclust.method = "ward")

#perform heirarical clustering on the woe data
hclust_woe = WOEClust_hclust(woe_data_ided, id="ID", target="Bad", 10,IV=infoTables_10$Summary, method="ward.D")

#do some summaries of the summary output variables

sum(data$`Profit Per Account`) #3511931

sum(abs(data$`Profit Per Account`)) #5646207

mean_1 = mean(data$`Net Purchases Total`[data$Bad==1])
mean_0 = mean(data$`Net Purchases Total`[data$Bad==0])
mean_df = data.frame(Bad = c(1,0), mean = c(mean_1, mean_0))

ggplot(data=data,aes(x=data$`Net Purchases Total`, color=as.factor(Bad)))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, fill="white")+
  geom_density(alpha=0.6)+scale_color_brewer(palette="Dark2") + 
  geom_vline(data=mean_df, aes(xintercept=mean, color=as.factor(Bad)),
             linetype="dashed")+
  labs(title = "Net Purchases by Group\n", x = "Net Purchases over 10 months", y = "density", color = "Bad")

#create an index table of what the WOE values mean and what bin proportion of the population each on is in

Index_WOE = data.frame(matrix(nrow = 152, ncol = 6))
t = 0
for (i in infoTables_10$Tables){
  for (j in 1:nrow(i)){
    t=t+1
    Index_WOE[t,1] <- names(i[1])
    Index_WOE[t,2:6] <- i[j,]
  }
}

colnames(Index_WOE) <- c("Variable","Bin", names(infoTables_10$Tables[[1]])[2:5])
