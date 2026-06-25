library(readxl)
library(dplyr)
library(roahd)


########################## AGE EXPERIMENT ####################################
age_participants = setdiff(1:21,c(4,13,15))
n_age=length(age_participants)
elec = c("C3","C4","O1","O2")

################### pairwise analysis age ##########################
cond_age=c( "ambiguous", "younger", "older")

pdf("Pairwise_AGE_AlphaBeta_4_electrodes.pdf",width=8,height=8)
 
lt = 1250 #length of post-stimulus interval
m <- 1
for (j in 1:2){ #conditions   
  for (h in (j+1):3){ #conditions
         
         mean_alpha = matrix(0, nrow = lt, ncol= 4)
         mean_beta = matrix(0, nrow = lt, ncol= 4)
         
    for (l in 1:4){ #electrodes 
         
           alpha = matrix(0, nrow = n_age, ncol= lt)  
           beta = matrix(0, nrow = n_age, ncol= lt)  
       
       for (k in 1:n_age){  #B stands for baseline, #P stands for post-stimulus
         i= age_participants[k]
         
         file = paste0("CSV_files/time_Age_B_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
         tB = read.csv(file,sep=",",dec=".",header = F)    # time baseline
         file = paste0("CSV_files/time_Age_P_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
         tP = as.numeric( read.csv(file,sep=",",dec=".",header = F)  ) # time post-stimulus
         #First condition j alpha 
         file = paste0("CSV_files/alpha_band_Age_B_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
         abB = read.csv(file,sep=",",dec=".",header = F)    # alpha-band baseline
         abB = mean(apply(abB, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
         abB = rep(abB, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
         file = paste0("CSV_files/alpha_band_Age_P_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
         abP = read.csv(file,sep=",",dec=".",header = F)    # alpha-band post-stimulus
         abP = apply(abP, 2, mean)  #Function of time: mean over all frequencies  
         #Second condition h  alpha
         file = paste0("CSV_files/alpha_band_Age_B_",cond_age[[h]],"_P",i,"_",elec[l],".csv")
         abB2 = read.csv(file,sep=",",dec=".",header = F)    # alpha-band baseline
         abB2 = mean(apply(abB2, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
         abB2 = rep(abB2, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
         file = paste0("CSV_files/alpha_band_Age_P_",cond_age[[h]],"_P",i,"_",elec[l],".csv")
         abP2 = read.csv(file,sep=",",dec=".",header = F)    # alpha-band post-stimulus
         abP2 = apply(abP2, 2, mean)  #Function of time: mean over all frequencies  
         
         #First condition j beta 
         file = paste0("CSV_files/beta_band_Age_B_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
         bbB = read.csv(file,sep=",",dec=".",header = F)    # beta-band baseline
         bbB = mean(apply(bbB, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
         bbB = rep(bbB, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
         file = paste0("CSV_files/beta_band_Age_P_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
         bbP = read.csv(file,sep=",",dec=".",header = F)    # beta-band post-stimulus
         bbP = apply(bbP, 2, mean)  #Function of time: mean over all frequencies  
         
         #Second condition h  beta
         file = paste0("CSV_files/beta_band_Age_B_",cond_age[[h]],"_P",i,"_",elec[l],".csv")
         bbB2 = read.csv(file,sep=",",dec=".",header = F)    # beta-band baseline
         bbB2 = mean(apply(bbB2, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
         bbB2 = rep(bbB2, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
         file = paste0("CSV_files/beta_band_Age_P_",cond_age[[h]],"_P",i,"_",elec[l],".csv")
         bbP2 = read.csv(file,sep=",",dec=".",header = F)    # beta-band post-stimulus
         bbP2 = apply(bbP2, 2, mean)  #Function of time: mean over all frequencies  
         
         
         alpha[k, ] = abP/abB - abP2/abB2
         beta[k, ] = bbP/bbB - bbP2/bbB2
         
       }
       
       mean_alpha[,l] <-apply(alpha,2,mean) 
       mean_beta[,l] <-apply(beta,2,mean) 
       
    }

         
         df_alpha = data.frame( t=rep(tP,4), electrode = rep(c("C3","C4","O1","O2"),each=lt),
                                position = rep(c("Central","Occipital"),each=2*lt),
                                Mean = as.numeric(mean_alpha))
         df_alpha$Electrode = paste(df_alpha$position,"-",df_alpha$electrode)
       
         
         p1 <- ggplot(df_alpha, aes(x=t, y= Mean, group=electrode, color=Electrode,linetype=Electrode)) + 
           geom_line(linewidth=1) +
           geom_hline(yintercept = 0,col="black")+
           # ggtitle(paste0("Mean difference of CWT Post-stimulus/Baseline Ratios\nAge experiment. Alpha band \n",
           #                toupper(cond_age[[j]])," vs ", toupper(cond_age[[h]])))+ 
           ggtitle("ALPHA band")+ 
           xlab("time (s)") + ylab("Mean CWT amplitude ratios difference") +
           #   coord_cartesian( ylim=c(-0.13,0.2)) +
           scale_color_manual(values=c("magenta3", "magenta3", "slateblue", "slateblue")) +
           scale_linetype_manual(values=c(1,4,1,4)) +
           theme_minimal()
         
         df_beta = data.frame( t=rep(tP,4), electrode = rep(c("C3","C4","O1","O2"),each=lt),
                                position = rep(c("Central","Occipital"),each=2*lt),
                                Mean = as.numeric(mean_beta))
         df_beta$Electrode = paste(df_beta$position,"-",df_beta$electrode)
         
         p2<- ggplot(df_beta, aes(x=t, y= Mean, group=electrode, color=Electrode,linetype=Electrode)) + 
           geom_line(linewidth=1) +
           geom_hline(yintercept = 0,col="black")+
           # ggtitle(paste0("Mean difference of CWT Post-stimulus/Baseline Ratios\nAge experiment. Beta band \n",
           #                toupper(cond_age[[j]])," vs ", toupper(cond_age[[h]])))+ 
           ggtitle("BETA band")+ 
           xlab("time (s)") + ylab("Mean CWT amplitude ratios difference") +
           #   coord_cartesian( ylim=c(-0.13,0.2)) +
           scale_color_manual(values=c("magenta3", "magenta3", "slateblue", "slateblue")) +
           scale_linetype_manual(values=c(1,4,1,4)) +
           theme_minimal()
         
         # gridExtra::grid.arrange(p1, p2)
         

         tit = grid::textGrob(paste("AGE task:",toupper(cond_age[[j]])," vs ", toupper(cond_age[[h]]) ), gp=grid::gpar(fontsize=16,font=1))
         gridExtra::grid.arrange(p1, p2,top=tit, nrow=2)
     }
}

dev.off()
 
 ########################## EMOTION EXPERIMENT ####################################
 emo_participants = setdiff(1:21,c(5,10,13,15))
 n_emo=length(emo_participants)
 elec = c("C3","C4","O1","O2")
 
 
######################### pairwise analysis emotion ################################
cond_emo=c("ambiguous","positive",  "negative")
  
 pdf("Pairwise_EMOTION_AlphaBeta_4_electrodes.pdf",width=8,height=8)
 
 lt = 1250 #length of post-stimulus interval
 m <- 1
 for (j in 1:2){ #conditions   
   for (h in (j+1):3){ #conditions
     
     mean_alpha = matrix(0, nrow = lt, ncol= 4)
     mean_beta = matrix(0, nrow = lt, ncol= 4)
     
     for (l in 1:4){ #electrodes 
       
          alpha = matrix(0, nrow = n_age, ncol= lt)  #1250 is time_P length
          beta = matrix(0, nrow = n_age, ncol= lt)  #1250 is time_P length
       
       
       for (k in 1:n_emo){  #B stands for baseline, #P stands for post-stimulus
         i= emo_participants[k]
         
         file = paste0("CSV_files/time_Emo_B_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
         tB = read.csv(file,sep=",",dec=".",header = F)    # time baseline
         file = paste0("CSV_files/time_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
         tP = as.numeric( read.csv(file,sep=",",dec=".",header = F)  ) # time post-stimulus
         #First condition j alpha 
         file = paste0("CSV_files/alpha_band_Emo_B_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
         abB = read.csv(file,sep=",",dec=".",header = F)    # alpha-band baseline
         abB = mean(apply(abB, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
         abB = rep(abB, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
         file = paste0("CSV_files/alpha_band_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
         abP = read.csv(file,sep=",",dec=".",header = F)    # alpha-band post-stimulus
         abP = apply(abP, 2, mean)  #Function of time: mean over all frequencies  
         #Second condition h  alpha
         file = paste0("CSV_files/alpha_band_Emo_B_",cond_emo[[h]],"_P",i,"_",elec[l],".csv")
         abB2 = read.csv(file,sep=",",dec=".",header = F)    # alpha-band baseline
         abB2 = mean(apply(abB2, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
         abB2 = rep(abB2, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
         file = paste0("CSV_files/alpha_band_Emo_P_",cond_emo[[h]],"_P",i,"_",elec[l],".csv")
         abP2 = read.csv(file,sep=",",dec=".",header = F)    # alpha-band post-stimulus
         abP2 = apply(abP2, 2, mean)  #Function of time: mean over all frequencies  
         
         #First condition j beta 
         file = paste0("CSV_files/beta_band_Emo_B_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
         bbB = read.csv(file,sep=",",dec=".",header = F)    # beta-band baseline
         bbB = mean(apply(bbB, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
         bbB = rep(bbB, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
         file = paste0("CSV_files/beta_band_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
         bbP = read.csv(file,sep=",",dec=".",header = F)    # beta-band post-stimulus
         bbP = apply(bbP, 2, mean)  #Function of time: mean over all frequencies  
         
         #Second condition h  beta
         file = paste0("CSV_files/beta_band_Emo_B_",cond_emo[[h]],"_P",i,"_",elec[l],".csv")
         bbB2 = read.csv(file,sep=",",dec=".",header = F)    # beta-band baseline
         bbB2 = mean(apply(bbB2, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
         bbB2 = rep(bbB2, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
         file = paste0("CSV_files/beta_band_Emo_P_",cond_emo[[h]],"_P",i,"_",elec[l],".csv")
         bbP2 = read.csv(file,sep=",",dec=".",header = F)    # beta-band post-stimulus
         bbP2 = apply(bbP2, 2, mean)  #Function of time: mean over all frequencies  
         
         
         alpha[k, ] = abP/abB - abP2/abB2
         beta[k, ] = bbP/bbB - bbP2/bbB2
         
       }
       
          mean_alpha[,l] <-apply(alpha,2,mean) 
          mean_beta[,l] <-apply(beta,2,mean) 
          
     }
     
     
     df_alpha = data.frame( t=rep(tP,4), electrode = rep(c("C3","C4","O1","O2"),each=lt),
                            position = rep(c("Central","Occipital"),each=2*lt),
                            Mean = as.numeric(mean_alpha))
     df_alpha$Electrode = paste(df_alpha$position,"-",df_alpha$electrode)
     
     
     p1 <- ggplot(df_alpha, aes(x=t, y= Mean, group=electrode, color=Electrode,linetype=Electrode)) + 
       geom_line(linewidth=1) +
       geom_hline(yintercept = 0,col="black")+
       # ggtitle(paste0("Mean difference of CWT Post-stimulus/Baseline Ratios\nEmotion experiment. Alpha band \n",
       #                toupper(cond_emo[[j]])," vs ", toupper(cond_emo[[h]])))+ 
       ggtitle("ALPHA band")+ 
       xlab("time (s)") + ylab("Mean CWT amplitude ratios difference") +
       #   coord_cartesian( ylim=c(-0.13,0.2)) +
       scale_color_manual(values=c("magenta3", "magenta3", "slateblue", "slateblue")) +
       scale_linetype_manual(values=c(1,4,1,4)) +
       theme_minimal()
     
     df_beta = data.frame( t=rep(tP,4), electrode = rep(c("C3","C4","O1","O2"),each=lt),
                           position = rep(c("Central","Occipital"),each=2*lt),
                           Mean = as.numeric(mean_beta))
     
     df_beta$Electrode = paste(df_beta$position,"-",df_beta$electrode)
    
     
     p2<- ggplot(df_beta, aes(x=t, y= Mean, group=electrode, color=Electrode,linetype=Electrode)) + 
       geom_line(linewidth=1) +
       geom_hline(yintercept = 0,col="black")+
       # ggtitle(paste0("Mean difference of CWT Post-stimulus/Baseline Ratios\nEmotion experiment. Beta band \n",
       #                toupper(cond_emo[[j]])," vs ", toupper(cond_emoe[[h]])))+ 
       ggtitle("BETA band")+ 
       xlab("time (s)") + ylab("Mean CWT amplitude ratios difference") +
       #   coord_cartesian( ylim=c(-0.13,0.2)) +
       scale_color_manual(values=c("magenta3", "magenta3", "slateblue", "slateblue")) +
       scale_linetype_manual(values=c(1,4,1,4)) +
       theme_minimal()
     
     # gridExtra::grid.arrange(p1, p2)
     
     tit = grid::textGrob(paste("EMOTION task:",toupper(cond_emo[[j]])," vs ", toupper(cond_emo[[h]]) ), gp=grid::gpar(fontsize=16,font=1))
     gridExtra::grid.arrange(p1, p2,top=tit, nrow=2)
     
   }
 }
 
 dev.off()
