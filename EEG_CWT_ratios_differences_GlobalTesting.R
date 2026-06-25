library(readxl)
library(dplyr)
library(roahd)
library(ggplot2)

source("visualization_real_data.R")
source("aux_fun_registration.R")

# Fast and Fair Simultaneous Confidence Bands for Functional Parameters
# https://www.dliebl.com/ffscb/
# devtools::install_github("lidom/ffscb")

library("ffscb")

#### EXAMPLE
# # Generate a sample
# p          <- 200
# N          <- 80
# grid       <- make_grid(p, rangevals=c(0,1))
# mu         <- meanf_poly(grid,c(0,.25))
# names(mu)  <- grid
# cov.m      <- make_cov_m(cov.f = covf_nonst_matern, grid=grid, cov.f.params=c(2, 1/4, 1/4))
# sample     <- make_sample(mu,cov.m,N)
# 
# # Compute the estimate, hat.mu, and its covariance, hat.cov.mu
# hat.mu     <- rowMeans(sample)
# hat.cov    <- crossprod(t(sample - hat.mu)) / N
# hat.cov.mu <- hat.cov / N
# # Compute the tau-parameter
# # I.e., the 'roughness parameter function needed for the FFSCB-bands
# hat.tau    <- tau_fun(sample)
# 
# # Make and plot confidence bands
# b <- confidence_band(x          = hat.mu,
#                      cov.x      = hat.cov.mu,
#                      tau        = hat.tau,
#                      df         = N-1,
#                      type       = c("FFSCB.t", "Bs", "BEc", "naive.t"),
#                      conf.level = 0.95,
#                      n_int      = 4)
# plot(b)

#### FUNCIÓN NUESTRA PARA LLAMAR A confidence_band
conf.bands <- function(x, t, n_int=n_int, tit ="",xl="", yl="", H0=T, col.line=1, col.band=gray(.75),ylims=NULL,mu0=0, sig.level = 0.05, alternative="ts", H1text=""){ 
  # x is an array with individuals in rows and time points in columns
  # t is the vector of recorded time points
  # n_int is the nb of subintervals to consider for the conf. bands
  # tit, xl, yl are graphical parameters
  # H0: logic, wether to plot the horizontal line H0: mu(t)=0
  
  ### mis opciones
  # mu0 : value of mu under H0
  # alternative: "ts", "l", "g"
  
  sample <- t(x)
  rownames(sample) <- as.character(t)
  
  N=nrow(x)
  # Compute the estimate, hat.mu, and its covariance, hat.cov.mu
  hat.mu     <- rowMeans(sample)
  hat.cov    <- crossprod(t(sample - hat.mu)) / N
  hat.cov.mu <- hat.cov / N
  # Compute the tau-parameter 
  # I.e., the 'roughness parameter function needed for the FFSCB-bands
  hat.tau    <- tau_fun(sample)
  
  # Make and plot confidence bands
  if ( alternative!= "ts"){cl=1 - 2*sig.level}else{cl=1 - sig.level}# si sig.level = 0.05 => si one-sided cl=0.9, pq el cálculo de la banda está hecha para two-sided, 
                                                                    # y de cada lado se dejaría 0.1/2=0.05, por lo que la banda unilateral contiene el resto, i.e., 0.95
  b <- confidence_band(x          = hat.mu, 
                       cov.x      = hat.cov.mu, 
                       tau        = hat.tau, 
                       df         = N-1,
                       # type       = c("FFSCB.t", "Bs", "BEc", "naive.t"),  #todos los métodos
                       type       = "FFSCB.t",  #sólo el que nos interesa
                       conf.level = cl,  #confidence level, see line 67
                       n_int      = n_int)
  # plot(b, main=tit, xlab=xl, ylab=yl)   #Plot básico con varios tipos de ICs (simultaneos, naive y demás)
  
  #Plot más sofisticado, sólo con el FFSCB al 95% global y por subintervalos
  
  FF8_t_band        <- b[,-1] 
  FF8_crit_value_u  <- (FF8_t_band[,1] - hat.mu)/sqrt(diag(hat.cov.mu))
  
  grid=t
  
  ## Plots
  # width     <- 7
  # height    <- 3.125
  cex       <- 1.1
  cexs      <- 0.95
  ##
  # layout(mat = matrix(c(1,2,3,3), nrow=2, ncol=2),
  #        heights = c(1, 1),     # Heights of the two rows
  #        widths =  c(1, 1.25))  # Widths of the two columns
  # 
  # par(family = "serif", ps=13, cex.axis=1.05, font.main = 1)
  # par(mar=c(3.1, 2.1, 2, 1.1))
  # plot(y=hat.tau,x=grid*100, type="l", main = "", xlab = "", ylab="")
  # mtext(text = expression(paste("Roughness estimate ",hat(tau))), 3, line = 0.4, adj = 0, cex=cex)
  # ##
  # par(mar=c(3.1, 2.1, 2, 1.1))
  # matplot(y=cbind(FF8_crit_value_u), x=grid*100, col="black", lty=c(2,1), type = "l", 
  #         main="", xlab = "", ylab="", lwd=c(1.5,1), ylim=c(min(FF8_crit_value_u),4.1))
  # abline(v=c( 0:8 * 100/8), col="gray", lwd=1.3)
  # mtext(text = expression(paste("Fair adaptive critical value function ", hat(u)[alpha/2]^"*")), 3, line = 0.4, adj = 0, cex=cex)
  # mtext(text = "% of Stance Phase", 1, line = 2.25, cex=cexs)
  ## 
  par(mar=c(3.1, 4.1, 3.1, 2.1))
  if (is.null(ylims)){ylims=c(min(FF8_t_band),max(FF8_t_band))}
  matplot( y = b, x = grid, type="n", ylim=ylims, 
           ylab=yl, xlab =xl, main=tit)
  if (alternative =="l"){
    polygon(x=c(grid,rev(grid)), y=c(rep(ylims[1],length(t)),rev(FF8_t_band[,1])), col = col.band , border = col.band)
    
  }else{
    if (alternative=="g"){
      polygon(x=c(grid,rev(grid)), y=c(FF8_t_band[,2],rep(ylims[2],length(t))), col = col.band , border = col.band)
      
    }else{
      polygon(x=c(grid,rev(grid)), y=c(FF8_t_band[,2],rev(FF8_t_band[,1])), col = col.band , border = col.band)
      
    }
  }
  abline(  h = mu0, lwd=0.7)
  abline(v=seq(grid[1],grid[length(grid)],,n_int+1) , col="gray")
  lines(   y = hat.mu,  x = grid, col=col.line, lty=1)
  if (H0==T){
    c = parse(text=H1text)[[1]]
    axis(4, at = mu0, labels = bquote(H[0]:~mu==.(mu0)~.(c)))
  }
  legend(x=1500, y=5, legend = c(expression(paste("Estimated mean")), 
                                 expression(FF[t]^8)), y.intersp=1.05, bg="white", box.col = "white",
         lty=c(1,1), lwd = c(1.5,10), col=c(col.line, col.band), cex =cexs, seg.len=2)
  #mtext(text = xl, 1, line = 1.75, cex=cexs, font=2)
  #mtext(text = yl, 2, line = 1.85, cex=cexs, font=2)
  # mtext(text = paste(tit,"\n 95% Simultaneous Confidence Band (SCB)\n99.375% SCB over each subinterval"), 3, line = 0.4, cex=cex)
  #mtext(text = tit, 3, line = 0.4, font=2, cex=cex, adj=0)
  #text(x = -1250, y = 3.75, labels="99.375% SCB", srt=90)
  box()
  return(FF8_t_band)
}
#### SAME WITH GGPLOT GRAPHICS
conf.bands_ggplot <- function(x, t, n_int=n_int, tit ="",xl="", yl="", H0=T, col.line=1, col.band=gray(.75),ylims=NULL,mu0=0, sig.level = 0.05, alternative="ts", H1text=""){ 
  # x is an array with individuals in rows and time points in columns
  # t is the vector of recorded time points
  # n_int is the nb of subintervals to consider for the conf. bands
  # tit, xl, yl are graphical parameters
  # H0: logic, wether to plot the horizontal line H0: mu(t)=0
  
  ### mis opciones
  # mu0 : value of mu under H0
  # alternative: "ts", "l", "g"
  
  sample <- t(x)
  rownames(sample) <- as.character(t)
  
  N=nrow(x)
  # Compute the estimate, hat.mu, and its covariance, hat.cov.mu
  hat.mu     <- rowMeans(sample)
  hat.cov    <- crossprod(t(sample - hat.mu)) / N
  hat.cov.mu <- hat.cov / N
  # Compute the tau-parameter 
  # I.e., the 'roughness parameter function needed for the FFSCB-bands
  hat.tau    <- tau_fun(sample)
  
  # Make and plot confidence bands
  if ( alternative!= "ts"){cl=1 - 2*sig.level}else{cl=1 - sig.level}# si sig.level = 0.05 => si one-sided cl=0.9, pq el cálculo de la banda está hecha para two-sided, 
  # y de cada lado se dejaría 0.1/2=0.05, por lo que la banda unilateral contiene el resto, i.e., 0.95
  b <- confidence_band(x          = hat.mu, 
                       cov.x      = hat.cov.mu, 
                       tau        = hat.tau, 
                       df         = N-1,
                       # type       = c("FFSCB.t", "Bs", "BEc", "naive.t"),  #todos los métodos
                       type       = "FFSCB.t",  #sólo el que nos interesa
                       conf.level = cl,  #confidence level, see line 67
                       n_int      = n_int)
  # plot(b, main=tit, xlab=xl, ylab=yl)   #Plot básico con varios tipos de ICs (simultaneos, naive y demás)
  
  #Plot más sofisticado, sólo con el FFSCB al 95% global y por subintervalos
  
  FF8_t_band        <- b[,-1] 
  FF8_crit_value_u  <- (FF8_t_band[,1] - hat.mu)/sqrt(diag(hat.cov.mu))
  
  grid=t
  
  # Build a data frame
  df <- data.frame(
    x = grid,
    mu_hat = hat.mu,
    lower = FF8_t_band[,1],
    upper = FF8_t_band[,2]
  )
  if (is.null(ylims)){ylims=c(min(FF8_t_band),max(FF8_t_band))}
  # Base plot
  p <- ggplot(df, aes(x = x)) +
    
    # Confidence band (polygon equivalent)
    {
      if (alternative == "l") {
        geom_ribbon(aes(ymin = ylims[1], ymax = lower),
                    fill = col.band, color = col.band)
      } else if (alternative == "g") {
        geom_ribbon(aes(ymin = upper, ymax = ylims[2]),
                    fill = col.band, color = col.band)
      } else {
        geom_ribbon(aes(ymin = lower, ymax = upper),
                    fill = col.band, color = col.band)
      }
    } +
    
    # Estimated mean line
    geom_line(aes(y = mu_hat), color = col.line, linetype = 1) +
    
    # Horizontal line (mu0)
    # geom_hline(yintercept = mu0, size = 0.3) +
    geom_segment(aes(x = min(x), xend = max(x), y = mu0, yend = mu0), size=0.5) +
    
    
    # Vertical grid lines
    geom_vline(xintercept = seq(grid[1], grid[length(grid)], length.out = n_int + 1),
               color = "gray") +
    
    # Labels
    labs(
      x = xl,
      y = yl,
      title = tit
    ) +  theme(plot.margin = unit(c(2,2,2,2), "cm")) +
    
    theme_minimal()
  
  # Add annotation for H0 if needed
  if (H0 == TRUE) {
    f <-paste0('H[0] * ":"~  mu * "=" * ', mu0,'   *  ~~~', H1text)
    p <- p +
      annotate("text", x = max(grid), y = mu0,
               label = f,
               parse=T,
               vjust = 1, angle = 90)
  }
  
  # Legend (manual)
  # p <- p +
  #   scale_color_manual(
  #     values = c("mean" = col.line, "band" = col.band),
  #     labels = c("Estimated mean", "FF[t]^8")
  #   )
  
  # Print
 #print(p)
  return(list(FF8_t_band=FF8_t_band, p=p))
}


############################## AGE EXPERIMENT ################################
age_participants = setdiff(1:21,c(4,13,15))
n_age=length(age_participants)
cond_age=c( "ambiguous", "younger","older")
elec = c("C3","C4","O1","O2")
col = gray(1:20/20)[c(11,14,16,19)]


alpha_plots = vector(mode="list",length=3)
beta_plots = vector(mode="list",length=3)
alpha_plots_h1_g = vector(mode="list",length=3) #for H1: mu>1
beta_plots_h1_g = vector(mode="list",length=3)



for (j in 1:3){ #conditions 
  
  alpha_plots[[j]] = vector(mode="list",length=4)
  beta_plots[[j]]  = vector(mode="list",length=4)
  alpha_plots_h1_g[[j]] = vector(mode="list",length=4) #for H1: mu>1
  beta_plots_h1_g[[j]] = vector(mode="list",length=4)
  
  for (l in 1:4){ #electrodes 
    
    alpha = matrix(0, nrow = n_age, ncol= 1250)  #1250 is time_P length
    beta = matrix(0, nrow = n_age, ncol= 1250)  #1250 is time_P length
    
    
    for (k in 1:n_age){  #B stands for baseline, #P stands for post-stimulus
      i= age_participants[k]
      
      file = paste0("CSV_files/time_Age_B_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
      tB = read.csv(file,sep=",",dec=".",header = F)    # time baseline
      file = paste0("CSV_files/time_Age_P_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
      tP = as.numeric( read.csv(file,sep=",",dec=".",header = F)  ) # time post-stimulus
      
      file = paste0("CSV_files/alpha_band_Age_B_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
      abB = read.csv(file,sep=",",dec=".",header = F)    # alpha-band baseline
      abB = mean(apply(abB, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
      abB = rep(abB, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
      file = paste0("CSV_files/alpha_band_Age_P_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
      abP = read.csv(file,sep=",",dec=".",header = F)    # alpha-band post-stimulus
      abP = apply(abP, 2, mean)  #Function of time: mean over all frequencies  
      
      # par(mfrow=c(2,2))
      # plot(tP,abP, type="l", main = paste0("alpha_band_Age_P_",cond_age[[j]],"_P",i,"_",elec[l]))
      # lines(tP, abB,col="red")
      # plot(tP,abP/abB, type="l", main = paste0("Ratio_alpha_band_Age_P_",cond_age[[j]],"_P",i,"_",elec[l]))
      # lines(tP, rep(1,length(tP)),col="red")
      
      file = paste0("CSV_files/beta_band_Age_B_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
      bbB = read.csv(file,sep=",",dec=".",header = F)    # beta-band baseline
      bbB = mean(apply(bbB, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
      bbB = rep(bbB, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
      file = paste0("CSV_files/beta_band_Age_P_",cond_age[[j]],"_P",i,"_",elec[l],".csv")
      bbP = read.csv(file,sep=",",dec=".",header = F)    # beta-band post-stimulus
      bbP = apply(bbP, 2, mean)  #Function of time: mean over all frequencies  
      
      # plot(tP,bbP, type="l", main = paste0("beta_band_Age_P_",cond_age[[j]],"_P",i,"_",elec[l]))
      # lines(tP, bbB,col="red")
      # plot(tP,bbP/bbB, type="l", main = paste0("Ratio_beta_band_Age_P_",cond_age[[j]],"_P",i,"_",elec[l]))
      # lines(tP, rep(1,length(tP)),col="red")
      
      
      alpha[k, ] = abP/abB
      beta[k, ] = bbP/bbB
      
    }
    
    ####### ONE-SIDED TEST OR ONE-SIDED CI: H0: bbP/bbB = 1     H1: bbP/bbB < 1 (and bbP/bbB >1).  (same with alpha band)  #########
    
    
    if(l==1){
      tit = paste0("ALPHA band. Mean and 95% CB\n", elec[l])
      yl = "CWT amplitude ratio"
    }
    if(l==2){
      tit = paste0("\n", elec[l])
      yl=""
    }
    if(l==3){
      tit = paste0("\n", elec[l])
      yl=""
    }
    if(l==4){
      tit = paste0("\n", elec[l])
      yl=""
    }
    p <- conf.bands_ggplot(alpha, tP, n_int=1,
                           tit=tit,
                           xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                           col.band = col[1],#ylims = c(-5,5),
                           mu0=1, alternative="l", H1text = "H[1]:~mu<1")$p     #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
    
    alpha_plots[[j]][[l]] <- p
    
    p <- conf.bands_ggplot(alpha, tP, n_int=1,
                           tit=tit,
                           xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                           col.band = col[1],#ylims = c(-5,5),
                           mu0=1, alternative="g", H1text = "H[1]:~mu>1")$p     #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
    
    alpha_plots_h1_g[[j]][[l]] <- p
    
    

    if(l==1){
      tit = paste0("BETA band. Mean and 95% CB\n", elec[l])
      yl = "CWT amplitude ratio"
    }else{
      yl=""
    }
    
    p <- conf.bands_ggplot(beta, tP, n_int=1,
                           tit=tit,
                           xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                           col.band = col[2],,#ylims = c(-5,5),
                           mu0=1, alternative="l", H1text = "H[1]:~mu<1")$p    #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
    
    beta_plots[[j]][[l]]  <- p
    
    p <- conf.bands_ggplot(beta, tP, n_int=1,
                           tit=tit,
                           xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                           col.band = col[2],,#ylims = c(-5,5),
                           mu0=1, alternative="g", H1text = "H[1]:~mu>1")$p    #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
    
    beta_plots_h1_g[[j]][[l]]  <- p
    
  }
}


pdf("AGE_r<1_AND_r>1_global.pdf",width=13,height=8)

for (j in 1:3){ #conditions 
  
  pl = c(alpha_plots[[j]][1:4],beta_plots[[j]][1:4])
  
  tit = grid::textGrob(bquote("AGE task:"~.(toupper(cond_age[[j]])) ~~~~~H[1]:~mu(t)<1), gp=grid::gpar(fontsize=16,font=1))
  gridExtra::grid.arrange(grobs= pl,top=tit, nrow=2)

  pl = c(alpha_plots_h1_g[[j]][1:4],beta_plots_h1_g[[j]][1:4])
  
  tit = grid::textGrob(bquote("AGE task:"~.(toupper(cond_age[[j]])) ~~~~~H[1]:~mu(t)>1), gp=grid::gpar(fontsize=16,font=1))
  gridExtra::grid.arrange(grobs= pl,top=tit, nrow=2)
  
}
dev.off()

 
 ######## pairwise analysis age ONE_SIDED TESTS ALL COMPARISONS #########
 
 alphaPair_plots = vector(mode="list",length=3)
 betaPair_plots = vector(mode="list",length=3)
 

 for (j in 1:3){ #conditions   #TODOS CONTRA TODOS PQ EL TEST NO ES SIMÉTRICO
   alphaPair_plots[[j]] = vector(mode="list",length=3)
   betaPair_plots[[j]] = vector(mode="list",length=3)
     for (h in 1:3){ #conditions
       alphaPair_plots[[j]][[h]] = vector(mode="list",length=4)
       betaPair_plots[[j]][[h]] = vector(mode="list",length=4)
       
       if(j!=h){
         
         for (l in 1:4){ #electrodes 
           
       alpha = matrix(0, nrow = n_age, ncol= 1250)  #1250 is time_P length
       beta = matrix(0, nrow = n_age, ncol= 1250)  #1250 is time_P length
       
       
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
       
       ####### GLOBAL (1 single sub-interval) ONE-SIDED TEST OR ONE-SIDED CI: H0: bbP/bbB - bbP2/bbB2 = 0     H1: bbP/bbB - bbP2/bbB2 <0.  (same with alpha band)  #########

       if(l==1){
         tit = paste0("ALPHA band. Mean and 95% CB\n", elec[l])
         yl = "Difference of CWT amplitude ratios"
       }
       if(l==2){
         tit = paste0("\n", elec[l])
         yl=""
       }
       if(l==3){
         tit = paste0("\n", elec[l])
         yl=""
       }
       if(l==4){
         tit = paste0("\n", elec[l])
         yl=""
       }
       p <- conf.bands_ggplot(alpha, tP, n_int=1,
                  tit=tit,
                  xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                  col.band = col[1],#ylims = c(-5,5),
                  mu0=0, alternative="l", H1text = "H[1]:~mu<0")$p     #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
       
       alphaPair_plots[[j]][[h]][[l]] <- p

 
       if(l==1){
         tit = paste0("BETA band. Mean and 95% CB\n", elec[l])
         yl = "Difference of CWT amplitude ratios"
       }else{
         yl=""
       }
       
       p <- conf.bands_ggplot(beta, tP, n_int=1,
                  tit=tit,
                  xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                  col.band = col[2],,#ylims = c(-5,5),
                  mu0=0, alternative="l", H1text = "H[1]:~mu<0")$p    #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
       
       betaPair_plots[[j]][[h]][[l]]  <- p
       
       }
     }
   }
 }
 


pdf("Pairwise_AGE_Alpha_band_H1_diff<0_ALL_COMP_global_test.pdf",width=13,height=8)
#all pairs of conditions, two comparisons for each

for (j in 1:2){ #conditions   
  for (h in (j+1):3){ #conditions
    
    #cond j - cond h
    pl = c(alphaPair_plots[[j]][[h]][1:4],betaPair_plots[[j]][[h]][1:4])
    

    tit = grid::textGrob(bquote("AGE task:"~.(toupper(cond_age[[j]]))~vs ~.(toupper(cond_age[[h]])) ~~~~~H[1]:~mu[i](t)-mu[j](t)<0), gp=grid::gpar(fontsize=16,font=1))
    gridExtra::grid.arrange(grobs= pl,top=tit, nrow=2)

    #cond h - cond j
    pl = c(alphaPair_plots[[h]][[j]][1:4],betaPair_plots[[h]][[j]][1:4])
    
    tit = grid::textGrob(bquote("AGE task:"~.(toupper(cond_age[[h]]))~vs ~.(toupper(cond_age[[j]])) ~~~~~H[1]:~mu[i](t)-mu[j](t)<0), gp=grid::gpar(fontsize=16,font=1))
    gridExtra::grid.arrange(grobs= pl,top=tit, nrow=2)

    
  }
}
dev.off()
 

########################## EMOTION EXPERIMENT ################################
emo_participants = setdiff(1:21,c(5,10,13,15))
n_emo=length(emo_participants)

cond_emo = c("ambiguous", "positive","negative")
elec = c("C3","C4","O1","O2")

alpha_plots = vector(mode="list",length=3)
beta_plots = vector(mode="list",length=3)
alpha_plots_h1_g = vector(mode="list",length=3) #for H1: mu>1
beta_plots_h1_g = vector(mode="list",length=3)


for (j in 1:3){ #conditions 
  
  alpha_plots[[j]] = vector(mode="list",length=4)
  beta_plots[[j]]  = vector(mode="list",length=4)
  alpha_plots_h1_g[[j]] = vector(mode="list",length=4) #for H1: mu>1
  beta_plots_h1_g[[j]] = vector(mode="list",length=4)
  
  for (l in 1:4){ #electrodes 
    
    alpha = matrix(0, nrow = n_emo, ncol= 1250)  #1250 is time_P length
    beta = matrix(0, nrow = n_emo, ncol= 1250)  #1250 is time_P length
    
    
    for (k in 1:n_emo){  #B stands for baseline, #P stands for post-stimulus
      i= emo_participants[k]
      
      file = paste0("CSV_files/time_Emo_B_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
      tB = read.csv(file,sep=",",dec=".",header = F)    # time baseline
      file = paste0("CSV_files/time_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
      tP = as.numeric( read.csv(file,sep=",",dec=".",header = F)  ) # time post-stimulus
      
      file = paste0("CSV_files/alpha_band_Emo_B_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
      abB = read.csv(file,sep=",",dec=".",header = F)    # alpha-band baseline
      abB = mean(apply(abB, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
      abB = rep(abB, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
      file = paste0("CSV_files/alpha_band_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
      abP = read.csv(file,sep=",",dec=".",header = F)    # alpha-band post-stimulus
      abP = apply(abP, 2, mean)  #Function of time: mean over all frequencies  
      
      # par(mfrow=c(2,2))
      # plot(tP,abP, type="l", main = paste0("alpha_band_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l]))
      # lines(tP, abB,col="red")
      # plot(tP,abP/abB, type="l", main = paste0("Ratio_alpha_band_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l]))
      # lines(tP, rep(1,length(tP)),col="red")
      
      file = paste0("CSV_files/beta_band_Emo_B_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
      bbB = read.csv(file,sep=",",dec=".",header = F)    # beta-band baseline
      bbB = mean(apply(bbB, 2, mean) )  #Baseline value: mean over all freqs and timepoints!!!!
      bbB = rep(bbB, length(tP))  #Constant curve with baseline value, length equal to post-stimulus length
      file = paste0("CSV_files/beta_band_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l],".csv")
      bbP = read.csv(file,sep=",",dec=".",header = F)    # beta-band post-stimulus
      bbP = apply(bbP, 2, mean)  #Function of time: mean over all frequencies  
      
      # plot(tP,bbP, type="l", main = paste0("beta_band_Emo_P_",cond_emo[[j]],"_P",i,"_",elec[l]))
      # lines(tP, bbB,col="red")
      # plot(tP,bbP/bbB, type="l", main = paste0("Ratio_beta_band_Age_P_",cond_emo[[j]],"_P",i,"_",elec[l]))
      # lines(tP, rep(1,length(tP)),col="red")
      
      
      alpha[k, ] = abP/abB
      beta[k, ] = bbP/bbB
      
    }
    
    ####### ONE-SIDED TEST OR ONE-SIDED CI: H0: bbP/bbB = 1     H1: bbP/bbB < 1 (and bbP/bbB >1).  (same with alpha band)  #########
    
    
    if(l==1){
      tit = paste0("ALPHA band. Mean and 95% CB\n", elec[l])
      yl = "CWT amplitude ratio"
    }
    if(l==2){
      tit = paste0("\n", elec[l])
      yl=""
    }
    if(l==3){
      tit = paste0("\n", elec[l])
      yl=""
    }
    if(l==4){
      tit = paste0("\n", elec[l])
      yl=""
    }
    p <- conf.bands_ggplot(alpha, tP, n_int=1,
                           tit=tit,
                           xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                           col.band = col[1],#ylims = c(-5,5),
                           mu0=1, alternative="l", H1text = "H[1]:~mu<1")$p     #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
    
    alpha_plots[[j]][[l]] <- p
    
    p <- conf.bands_ggplot(alpha, tP, n_int=1,
                           tit=tit,
                           xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                           col.band = col[1],#ylims = c(-5,5),
                           mu0=1, alternative="g", H1text = "H[1]:~mu>1")$p     #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
    
    alpha_plots_h1_g[[j]][[l]] <- p
    
    
    
    if(l==1){
      tit = paste0("BETA band. Mean and 95% CB\n", elec[l])
      yl = "CWT amplitude ratio"
    }else{
      yl=""
    }
    
    p <- conf.bands_ggplot(beta, tP, n_int=1,
                           tit=tit,
                           xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                           col.band = col[2],,#ylims = c(-5,5),
                           mu0=1, alternative="l", H1text = "H[1]:~mu<1")$p    #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
    
    beta_plots[[j]][[l]]  <- p
    
    p <- conf.bands_ggplot(beta, tP, n_int=1,
                           tit=tit,
                           xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                           col.band = col[2],,#ylims = c(-5,5),
                           mu0=1, alternative="g", H1text = "H[1]:~mu>1")$p    #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
    
    beta_plots_h1_g[[j]][[l]]  <- p
    
    
  }
}


pdf("EMOTION_r<1_AND_r>1_global.pdf",width=13,height=8)

for (j in 1:3){ #conditions 
  
  pl = c(alpha_plots[[j]][1:4],beta_plots[[j]][1:4])
  
  tit = grid::textGrob(bquote("EMOTION task:"~.(toupper(cond_emo[[j]])) ~~~~~H[1]:~mu(t)<1), gp=grid::gpar(fontsize=16,font=1))
  gridExtra::grid.arrange(grobs= pl,top=tit, nrow=2)
  
  pl = c(alpha_plots_h1_g[[j]][1:4],beta_plots_h1_g[[j]][1:4])
  
  tit = grid::textGrob(bquote("EMOTION task:"~.(toupper(cond_emo[[j]])) ~~~~~H[1]:~mu(t)>1), gp=grid::gpar(fontsize=16,font=1))
  gridExtra::grid.arrange(grobs= pl,top=tit, nrow=2)
  
}
dev.off()



########### pairwise analysis emo ONE_SIDED TESTS ALL COMPARISONS##############

alphaPair_plots = vector(mode="list",length=3)
betaPair_plots = vector(mode="list",length=3)

for (j in 1:3){ #conditions   #TODOS CONTRA TODOS PQ EL TEST NO ES SIMÉTRICO
  alphaPair_plots[[j]] = vector(mode="list",length=3)
  betaPair_plots[[j]] = vector(mode="list",length=3)
  for (h in 1:3){ #conditions
    alphaPair_plots[[j]][[h]] = vector(mode="list",length=4)
    betaPair_plots[[j]][[h]] = vector(mode="list",length=4)
    
    if(j!=h){
      
         for (l in 1:4){ #electrodes 
           
         alpha = matrix(0, nrow = n_emo, ncol= 1250)  #1250 is time_P length
         beta = matrix(0, nrow = n_emo, ncol= 1250)  #1250 is time_P length
         
         
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
         
         ####### GLOBAL (1 single sub-interval) ONE-SIDED TEST OR ONE-SIDED CI: H0: bbP/bbB - bbP2/bbB2 = 0     H1: bbP/bbB - bbP2/bbB2 <0.  (same with alpha band)  #########
         # 
         # conf.bands(alpha, tP, n_int=1,
         #            tit=paste0("ALPHA BAND Difference of CWT Ratios: Post_stimulus/Baseline. Mean and 95% Conf. band  \nEMOTION experiment ",
         #                       toupper(cond_emo[[j]])," vs ", toupper(cond_emo[[h]]),", Elec. ",elec[l]),
         #            xl = "time (s)", yl = "Difference of CWT amplitude ratios", H0=T, col.line = 1, col.band = col[l],#ylims = c(-5,5),
         #            mu0=0, alternative="l", H1text = "H[1]:~mu<0")#alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
         # 
         # alphaPair_plots[[m]] <- recordPlot()
         # 
         # conf.bands(beta, tP, n_int=1,
         #            tit=paste0("BETA BAND Difference of CWT Ratios: Post_stimulus/Baseline. Mean and 95% Conf. band  \nEMOTION experiment ",
         #                       toupper(cond_emo[[j]])," vs ", toupper(cond_emo[[h]]),", Elec. ",elec[l]),
         #            xl = "time (s)", yl = "Difference of CWT amplitude ratios", H0=T, col.line = 1, col.band = col[l],,#ylims = c(-5,5),
         #            mu0=0, alternative="l", H1text = "H[1]:~mu<0")#alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
         # 
         # betaPair_plots[[m]] <- recordPlot()
         
         if(l==1){
           tit = paste0("ALPHA band. Mean and 95% CB\n", elec[l])
           yl = "Difference of CWT amplitude ratios"
         }
         if(l==2){
           tit = paste0("\n", elec[l])
           yl=""
         }
         if(l==3){
           tit = paste0("\n", elec[l])
           yl=""
         }
         if(l==4){
           tit = paste0("\n", elec[l])
           yl=""
         }
         p <- conf.bands_ggplot(alpha, tP, n_int=1,
                                tit=tit,
                                xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                                col.band = col[1],#ylims = c(-5,5),
                                mu0=0, alternative="l", H1text = "H[1]:~mu<0")$p     #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
         
         alphaPair_plots[[j]][[h]][[l]] <- p
         
         # conf.bands(beta, tP, n_int=1,
         #            tit=paste0("BETA BAND Difference of CWT Ratios: Post_stimulus/Baseline. Mean and 95% Conf. band  \nAGE experiment ",
         #                       toupper(cond_age[[j]])," vs ", toupper(cond_age[[h]]),", Elec. ",elec[l]),
         #            xl = "time (s)", yl = "Difference of CWT amplitude ratios", H0=T, col.line = 1, col.band = col[l],,#ylims = c(-5,5),
         #            mu0=0, alternative="l", H1text = "H[1]:~mu<0")#alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
         
         if(l==1){
           tit = paste0("BETA band. Mean and 95% CB\n", elec[l])
           yl = "Difference of CWT amplitude ratios"
         }else{
           yl=""
         }
         
         p <- conf.bands_ggplot(beta, tP, n_int=1,
                                tit=tit,
                                xl = "time (s)", yl = yl, H0=F, col.line = 1, 
                                col.band = col[2],,#ylims = c(-5,5),
                                mu0=0, alternative="l", H1text = "H[1]:~mu<0")$p    #alternative="ts", H1text = "H[1]:~mu!=0")  #alternative="g", H1text = "H[1]:~mu>1")
         
         betaPair_plots[[j]][[h]][[l]]  <- p
         

       }
     }
   }
 }
 
 
 
pdf("Pairwise_EMO_Alpha_band_H1_diff<0_ALL_COMP_global_test.pdf",width=13,height=8)
#all pairs of conditions, two comparisons for each

 for (j in 1:2){ #conditions   
   for (h in (j+1):3){ #conditions
     
     #cond j - cond h
     pl = c(alphaPair_plots[[j]][[h]][1:4],betaPair_plots[[j]][[h]][1:4])
     

     tit = grid::textGrob(bquote("EMOTION task:"~.(toupper(cond_emo[[j]]))~vs ~.(toupper(cond_emo[[h]])) ~~~~~H[1]:~mu[i](t)-mu[j](t)<0), gp=grid::gpar(fontsize=16,font=1))
     gridExtra::grid.arrange(grobs= pl,top=tit, nrow=2)   

     #cond h - cond j
     pl = c(alphaPair_plots[[h]][[j]][1:4],betaPair_plots[[h]][[j]][1:4])
     

     tit = grid::textGrob(bquote("EMOTION task:"~.(toupper(cond_emo[[h]]))~vs ~.(toupper(cond_emo[[j]])) ~~~~~H[1]:~mu[i](t)-mu[j](t)<0), gp=grid::gpar(fontsize=16,font=1))
     gridExtra::grid.arrange(grobs= pl,top=tit, nrow=2)     

   
 }
 }
dev.off()

 