#Function to quantify the association between fluorescence channels using mutual information (entropy method)
#In function it is encoded as NRF2 and p53 stains, but can be any pair of immunofluorescence stains
#data is the dataframe with CellProfiler quantification of NRF2 and p53 staining to be analyzed, the 
#first column has ImageNumber in it (this is how my IF data was always formatted, can just add a first column
#of ones if data is not formatted like this), second and third columns are stain values 
#nbins is the number of bins you want to break distribution up into (if a low and a high state- then 2, 
#if low, medium and high state- then 3, etc)
#lowstatesplit is the proportion of data in low state provided as a fraction/decimal
#Output is mutual information in bits
#12-31-19: function now returns MI so it can be saved to a variable

bivariatedist_MI <- function(data,nbins,lowstatesplit){
  
  library(tidyverse)
  library(purrr)
  library(mclust)
  library(infotheo)
  
  #MI calculation
  MIcalc<-function(data,nbins,lowstatesplit){
    
    #Create function to calculate cdf of variable
    P_NRF2<-ecdf(data[,2]) 
    P_p53<-ecdf(data[,3])
    
    #Transform random variable by its own CDF --> gives variable with a uniform distribution
    NRF2_uniform<-P_NRF2(data[,2])
    p53_uniform<-P_p53(data[,3])
    
    #Classify each point as being in low or high NRF2 and p53 states
    NRF2_df<-data.frame(NRF2_uniform)
    p53_df<-data.frame(p53_uniform)
    
    NRF2_df$State[NRF2_df$NRF2_uniform<((1/nbins)*2*lowstatesplit)]<-"low"
    NRF2_df$State[NRF2_df$NRF2_uniform>=((1/nbins)*2*lowstatesplit)]<-"high"
    
    p53_df$State[p53_df$p53_uniform<((1/nbins)*2*lowstatesplit)]<-"low"
    p53_df$State[p53_df$p53_uniform>=((1/nbins)*2*lowstatesplit)]<-"high"
  
  
    #Calculate NRF2 and p53 cutoffs in terms of original units
    NRF2cutoff<-quantile(data[,2],(1/nbins)*2*lowstatesplit)
    p53cutoff<-quantile(data[,3],(1/nbins)*2*lowstatesplit)
    
    #For clinical cases with only 2 images per case, use hypergeometric test to see
    #if the cells from one image are statistically enriched in either the low or high state 
    #of EITHER stains, if so exclude case
    if((length(unique(data$ImageNumber))==2) & (exists("shuffledata")==FALSE)){
      NRF2_df$ImageNumber<-data$ImageNumber
      p53_df$ImageNumber<-data$ImageNumber
      popsize<-length(data[,1])
      
      #NRF2 staining 
      breakdown<-NRF2_df %>% group_by(State,ImageNumber) %>% tally()
      state_totals<-NRF2_df %>% group_by(State) %>% tally()
      cellsperimage_totals<-NRF2_df %>% group_by(ImageNumber) %>% tally()
      
      #Testing if Image #1 is enriched in low state for NRF2 staining
      pop_success<-(state_totals %>% filter(State=='low'))$n
      sample_size<-(cellsperimage_totals %>% filter(ImageNumber==min(NRF2_df$ImageNumber)))$n
      sample_success<-(breakdown %>% filter(State=='low' & ImageNumber==min(NRF2_df$ImageNumber)))$n
      
      #Hypothesis test, lower.tail= FALSE to test for overrepresentation, we subtract sample successes by 1, when P[X ≥ x] is needed.
      #(P(X >x-1) = P(X ≥x))
      pval_Image1_low_NRF2<-phyper(sample_success-1,pop_success,popsize-pop_success,sample_size,lower.tail=FALSE) 
      
      #Testing if Image #1 is enriched in high state for NRF2 staining
      pop_success<-(state_totals %>% filter(State=='high'))$n
      sample_size<-(cellsperimage_totals %>% filter(ImageNumber==min(NRF2_df$ImageNumber)))$n
      sample_success<-(breakdown %>% filter(State=='high' & ImageNumber==min(NRF2_df$ImageNumber)))$n
      
      #Hypothesis test, lower.tail= FALSE to test for overrepresentation, we subtract sample successes by 1, when P[X ≥ x] is needed.
      #(P(X >x-1) = P(X ≥x))
      pval_Image1_high_NRF2<-phyper(sample_success-1,pop_success,popsize-pop_success,sample_size,lower.tail=FALSE)
      
      #Don't have to check for Image2 b/c it will produce the same p-values (but flipped for high and low state)
      
      #p53 staining
      breakdown<-p53_df %>% group_by(State,ImageNumber) %>% tally()
      state_totals<-p53_df %>% group_by(State) %>% tally()
      cellsperimage_totals<-p53_df %>% group_by(ImageNumber) %>% tally()
      
      #Testing if Image #1 is enriched in low state for p53 staining
      pop_success<-(state_totals %>% filter(State=='low'))$n
      sample_size<-(cellsperimage_totals %>% filter(ImageNumber==min(NRF2_df$ImageNumber)))$n
      sample_success<-(breakdown %>% filter(State=='low' & ImageNumber==min(NRF2_df$ImageNumber)))$n
      
      #Hypothesis test, lower.tail= FALSE to test for overrepresentation, we subtract sample successes by 1, when P[X ≥ x] is needed.
      #(P(X >x-1) = P(X ≥x))
      pval_Image1_low_p53<-phyper(sample_success-1,pop_success,popsize-pop_success,sample_size,lower.tail=FALSE) 
      
      #Testing if Image #1 is enriched in high state for p53 staining
      pop_success<-(state_totals %>% filter(State=='high'))$n
      sample_size<-(cellsperimage_totals %>% filter(ImageNumber==min(NRF2_df$ImageNumber)))$n
      sample_success<-(breakdown %>% filter(State=='high' & ImageNumber==min(NRF2_df$ImageNumber)))$n
      
      #Hypothesis test, lower.tail= FALSE to test for overrepresentation, we subtract sample successes by 1, when P[X ≥ x] is needed.
      #(P(X >x-1) = P(X ≥x))
      pval_Image1_high_p53<-phyper(sample_success-1,pop_success,popsize-pop_success,sample_size,lower.tail=FALSE)
      
      #Exclude sample if Image1 (same p vals for Image2 but flipped so only have to consider 1) is enriched 
      # for one state in NRF2 AND p53 staining
      if((pval_Image1_high_p53<0.05 | pval_Image1_low_p53<0.05) & (pval_Image1_high_NRF2<0.05 | pval_Image1_low_NRF2<0.05)){
        stop("This case has image enrichment in either the high or low state for NRF2 and p53 staining, therefore cannot be analyzed")
      }
      else{
        message("Although this case only has 2 images, the images do not have enrichment in either the high or low state for NRF2 and p53 staining and can therefore be analyzed")
      }
      
    }
    
    #Calculate mutual information
    MI<-mutinformation(NRF2_df$State,p53_df$State,method="emp")
    MI_bits<-log2(exp(1))*MI #Convert output from mutinformation (in nats) to bits (more conventional MI units)
    
    info=list(MI_bits,NRF2cutoff,p53cutoff)
    
    return(info)
    
  }
  
  #Process data and calculate MI 
  data=data.frame(data) #convert 3 column data to dataframe
  data<-data[ , purrr::map_lgl(data, is.numeric)] #select only numeric columns from dataframe
  colnames(data)<-c("ImageNumber","NRF2","p53") #name of columns
  MI<-MIcalc(data,nbins,lowstatesplit) #call MI calculation (function above)
  
  #Randomly shuffle first vector of staining (bootstrap) to generate MI of null distribution
  r<-sample(1:length(data[,1]), length(data[,1]), replace=FALSE)
  shuffledata<-transform(data,NRF2=data[r,2])
  null_MI<-MIcalc(shuffledata,nbins,lowstatesplit)
  
  #Bootstrap actual dataset and null distribution 100 times with replacement to generate 90% CIs
  data_bootstrap<-matrix(nrow=100,ncol=2)
  null_bootstrap<-matrix(nrow=100,ncol=2)
  for(i in 1:100){
    ran<-sample(1:length(data[,2]), length(data[,2]), replace=TRUE)
    data_ran<-data[ran,]
    null_ran<-shuffledata[ran,]
    data_bootstrap[i,1]<-unlist(MIcalc(data_ran,nbins,lowstatesplit)[1])
    data_bootstrap[i,2]<-data_bootstrap[i,1]-MI[[1]]
    null_bootstrap[i,1]<-unlist(MIcalc(null_ran,nbins,lowstatesplit)[1])
    null_bootstrap[i,2]<-null_bootstrap[i,1]-null_MI[[1]]
  }
  
  data_CI_lower<-unname(quantile(data_bootstrap[,1],0.05))
  data_CI_upper<-unname(quantile(data_bootstrap[,1],0.95))
  
  null_CI_lower<-unname(quantile(null_bootstrap[,1],0.05))
  null_CI_upper<-unname(quantile(null_bootstrap[,1],0.95))
  
  cat("Mutual information for data +-90%CIs,NRF2 cutoff, p53 cutoff:\n")
  MIprint<-c(MI[1],data_CI_lower,data_CI_upper,MI[2],MI[3])
  print(cat(unlist(MIprint),sep='\n'))
  
  cat("\n")
  
  cat("Mutual information for null +-90%CIs:\n")
  nullMI<-c(null_MI[1],null_CI_lower,null_CI_upper)
  print(cat(unlist(nullMI),sep='\n'))
  
  return(MI[1])
  
}




