################################
################################# 
# Dataset: ENIGMA_Subsample_Thickness.xlsx
# Analysis covariance model by Lara M. Wierenga
# contact details l.m.wierenga@fsw.leidenuniv.nl
# created june 7 2018 
# Updated feb 14th 2020
################################
################################ 

# This script is developed to run on the ENIGMA dataset, according to project proposal 
# Entitled: Sex differences in variability in brain structure across the lifespan: Results from the ENIGMA Lifespan working group 

# First and senior authors and email addresses: 
#   First: Lara M. Wierenga (l.m.wierenga@fsw.leidenuniv.nl), Leiden University
# Seniors: Eveline A. Crone (ECrone@FSW.leidenuniv.nl), Leiden University and Christian K. Tamnes (c.k.tamnes@psykologi.uio.no), University of Oslo
# 
# Co-author names and e-mail addresses:
#   Sophia Frangou (sophia.frangou@mssm.edu), Icahn School of Medicine at Mount Sinai
# Gaelle Doucet (gaelle.doucet@mssm.edu), Icahn School of Medicine at Mount Sinai
# Danai Dima (danai.dima@kcl.ac.uk), King's College London
# Paul Thompson (pthomp@usc.edu), USC
# Christian K. Tamnes (c.k.tamnes@psykologi.uio.no), Oslo University

# In addition, all members of the ENIGMA Lifespan working group who contribute data for this project and edit the manuscript, will be coauthors.

# The folowing research questions will be tested: 
# - 1) Do males show greater variability in brain structure than females?
# - 2) Is greater male variability in brain structure observed at both lower and upper extremities? 
# - 3) Does the sex difference in variability in brain structure differ across the lifespan? 
# - 4) Do females show greater diversity than males in regional volumes?



# ============= GLOBAL SETTINGS ===============
## dependencies: this is part of my standard library, most of these are used in this code
packages <- c("lattice", "lme4","reshape2", "ggplot2","reshape2","plyr","abind", "Hmisc", "RColorBrewer","mgcv","gtools","scales","grid","methods","igraph","devtools","data.table","Rmpfr","circlize","qvalue", "pgirmess","boot", "R.matlab","MatchIt", "lawstat", "car","randomForest","GGally","plotly","quantregForest","data.table","reprtree","lsr","rogme")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)

## Aestetics
colour_sex <-c("#fcbf16","#1596a9")
font_size <- 8
text_col <- "black"
B=100#00 # nr of permuations, for testing set to 100, for running set 10000
set <- 7:14 # set the columns with variables

setwd("~/Desktop/ENIGMA/")
ifelse(dir.exists("Tables/"),FALSE,dir.create("tables/"))
ifelse(dir.exists("Figures/"),FALSE,dir.create("Figures/"))
ifelse(dir.exists("Figures/Variance_age/"),FALSE,dir.create("Figures/Variance_age/"))
# ifelse(dir.exists("Figures/Plot_interactions/"),FALSE,dir.create("Figures/Plot_interactions/"))
ifelse(dir.exists("Figures/Plot_correlations/"),FALSE,dir.create("Figures/Plot_correlations/"))
ifelse(dir.exists("Figures/Variance_rousellet/"),FALSE,dir.create("Figures/Variance_rousellet/"))
ifelse(dir.exists("Figures/Age_by_sex/"),FALSE,dir.create("Figures/Age_by_sex/"))


# ============= FUNCTIONS ===============
cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  
  cd  <- md/csd                        ## cohen's d
}

var_ratio_sign_M_F <- function(xxx,vname,B){
  # variance ratio significance M > F
  perm=NULL
  n=nrow(xxx)
  sex=xxx$sex
  
  for(b in 1:B)
  {
    sex_permute=sex[sample(n,n)]
    aa=aggregate(xx[[vname]],by=list(gr=sex_permute),function(x) var(x, na.rm = TRUE))
    ratio=log(aa[aa$gr=="M",-1]/aa[aa$gr=="F",-1])
    perm=rbind(perm,ratio)
  }
  
  aa=aggregate(xx[[vname]],by=list(gr=xx$sex),function(x) var(x, na.rm = TRUE))
  var_ratio <- round(log(aa[aa$gr=="M",-1]/aa[aa$gr=="F",-1]),3)
  
  ## p.value
  p_value_var=mean(perm>=as.numeric(var_ratio),na.rm = TRUE)
  return(p_value_var)
}

var_ratio_sign_F_M <- function(xxx,vname,B){
  # variance ratio significance F > M
  perm=NULL
  n=nrow(xxx)
  sex=xxx$sex
  
  for(b in 1:B)
  {
    sex_permute=sex[sample(n,n)]
    aa=aggregate(xx[[vname]],by=list(gr=sex_permute),function(x) var(x, na.rm = TRUE))
    ratio=log(aa[aa$gr=="F",-1]/aa[aa$gr=="M",-1])
    perm=rbind(perm,ratio)
  }
  
  aa=aggregate(xx[[vname]],by=list(gr=xx$sex),function(x) var(x, na.rm = TRUE))
  var_ratio <- round(log(aa[aa$gr=="F",-1]/aa[aa$gr=="M",-1]),3)
  
  ## p.value
  p_value_var=mean(perm>=as.numeric(var_ratio),na.rm = TRUE)
  return(p_value_var)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# ============= Load dataset  ===============
xx <- data.frame(read.table('~/ENIGMA_Subsample_Thickness.csv',sep=",",header=TRUE))


for (i in c(1,3:5)){
  xx[,i] <- as.factor(xx[,i])
}
xx$sex <- ifelse(xx$sex==1,"M","F") # change binary variable of males and famles


# ---------------------------------------------------
# # 1) Do males show greater variability in brain structure than females?
# ---------------------------------------------------

## A. Create data frame with residuals using random forest analysis (controlling for sex (mean diff), age, FS_Version and Scanner_Strengh) 

#select brain variables
vnames=names(xx)[c(set)] 

# estimate sex difference in mean and variance
table_ttest <- NULL
for(vname in vnames){
  
  #use linear model to covary for cohort
  df=data.frame(y=xx[[vname]],age=xx$age,Cohort=xx$Cohort,FS_Version=xx$FS_Version,Scanner_Strength=xx$Scanner_Strengh)
  lm_cohort <- lm(y ~ poly(age,3,raw=T) + Cohort, data=df)
  X1 <- as.data.frame(predict(lm_cohort, newdata=df,type="terms"))
  df$y <- df$y-X1$Cohort

  df <- subset(df, select=-c(Cohort))
  # remove effects of other covariates
  ff=randomForest(y~.,data=df,ntrees=500,proximity=TRUE,importance=TRUE, do.trace=100)
  df$residuals=ff$y -predict(ff) 
  df=cbind(df,sex=xx$sex)
  
  #plot males and females by age
  
  p2 <- ggplot(data=df,aes(x=age,y=residuals,col=sex)) +  
    geom_point() + 
    ggtitle(vname) +
    scale_colour_manual(values=(colour_sex)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.ticks=element_blank(),
          plot.title = element_text (size = font_size+4,color="black",hjust=.5)#family="Arial"), 
    )
  
  print(p2)
  
  filename<- paste("Figures/Age_by_sex/Age_by_sex_",vname,".pdf",sep="")
  dev.copy(pdf,filename,width=30/2.54, height=22/2.54)
  dev.off()
  
  
  #estiamte mean sex differences
  M_female <- t.test(residuals ~ sex, data=df)[5]$estimate[1]
  M_male <- t.test(residuals ~ sex, data=df)[5]$estimate[2]
  p <- round(t.test(residuals ~ sex, data=df)$p.value,3)
  D <- round(cohensD(residuals ~ sex, data=df),3)
  
  # remove effects of covariates and mean sex difference
  
  df2=data.frame(y=xx[[vname]],sex=xx$sex,age=xx$age,Cohort=xx$Cohort,FS_Version=xx$FS_Version,Scanner_Strength=xx$Scanner_Strengh) 
  lm_cohort <- lm(y ~ poly(age,3,raw=T) + Cohort, data=df2)
  X1 <- as.data.frame(predict(lm_cohort, newdata=df2,type="terms"))
  df2$y <- df2$y-X1$Cohort
  
  df2 <- subset(df2, select=-c(Cohort))
  

  ff=randomForest(y~.,data=df2,ntrees=500,proximity=TRUE,importance=TRUE, do.trace=100)
  df2$residuals=ff$y -predict(ff) 

  # estiatme variance ratio on mean corrected values
  aa=aggregate(df2$residuals,by=list(gr=df2$sex),function(x) var(x, na.rm = TRUE))
  var_ratio <- round(log(aa[aa$gr=="M",-1]/aa[aa$gr=="F",-1]),3)
  p_value_var_M_F <- var_ratio_sign_M_F(df2,vname,B)  
  p_value_var_F_M <- var_ratio_sign_F_M(df2,vname,B)  
  
  table_ttest <- rbind(table_ttest,cbind(vname,M_female,M_male,p,D,var_ratio,p_value_var_M_F,p_value_var_F_M))
  
}


write.table(table_ttest,file=paste("Tables/table_variance_diff_",vnames[length(vnames)],".csv",sep=""),sep=",",row.names=F)


# plot variance ratio
table_ttest_plot <- table_ttest
rownames(table_ttest_plot) <- table_ttest_plot[,1]
table_ttest_plot <- as.data.frame(table_ttest_plot)
table_ttest_plot$var_ratio <- as.numeric(levels(table_ttest_plot$var_ratio))[table_ttest_plot$var_ratio]
table_ttest_plot$p_value_var_M_F <- as.numeric(levels(table_ttest_plot$p_value_var_M_F))[table_ttest_plot$p_value_var_M_F]

table_ttest_plot$col_label <- ifelse(table_ttest_plot$var_ratio>0,0,1)

for ( i in c(7)){
  table_ttest_plot[,i] <- sapply(table_ttest_plot[,i], function(x) ifelse(x < .01,"**",ifelse(x<.05 && x >.01,"*","")))
}

p1 <- ggplot(table_ttest_plot,aes(x=reorder(vname,var_ratio),y=table_ttest_plot$var_ratio)) 
p1 <- p1 + geom_bar(stat="identity",aes(fill=as.factor(col_label))) + 
  scale_fill_manual(values=rev(colour_sex)) +
  coord_flip() + 
  xlab("")+
  ylab("variance ratio")+
  geom_text(aes(reorder(vname,var_ratio), (var_ratio) + .02 , label = p_value_var_M_F),size=font_size-5,color=text_col, data = table_ttest_plot)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.text.x = element_text (size = font_size,color=text_col),#,family="Arial"), 
        axis.text.y = element_text (size = font_size,color=text_col),#,family="Arial"), 
        axis.title.y = element_text (size = font_size,color=text_col),#,family="Arial"), 
        axis.title.x = element_text (size = font_size,color=text_col),#,family="Arial"), 
        axis.ticks=element_blank(),
        legend.position="none")
print(p1)

filename<- paste("Figures/Variance_diff_",vnames[length(vnames)],".pdf",sep="")
dev.copy(pdf,filename,width=30/2.54, height=22/2.54)
dev.off()



# ---------------------------------------------------
# - 2) Is greater male variability in brain structure observed at both lower and upper extremities? 
# ---------------------------------------------------

# plot only significant effects
vnames_sign <- as.character(table_ttest_plot[table_ttest_plot$p_value_var_M_F!="n.s.",][[1]])

#plot significants in plot code by rouselette et al. 

for(vname in vnames_sign){

  #use linear model to covary for cohort
  resid_DF_all <- data.frame(y=xx[[vname]],sex=xx$sex,age=xx$age,Cohort=xx$Cohort,FS_Version=xx$FS_Version,Scanner_Strength=xx$Scanner_Strengh)
  lm_cohort <- lm(y ~ poly(age,3,raw=T) + Cohort, data=resid_DF_all)
  X1 <- as.data.frame(predict(lm_cohort, newdata=resid_DF_all,type="terms"))
  resid_DF_all$y <- resid_DF_all$y-X1$Cohort
  
  resid_DF_all <- subset(resid_DF_all, select=-c(Cohort))

    # remove mean sex difference and covariates
    ff=randomForest(y~.,data=resid_DF_all,ntrees=500,proximity=TRUE,importance=TRUE, do.trace=100)
    residuals=ff$y -predict(ff) 
    resid_DF_all$y <- residuals
    resid_DF_all <- resid_DF_all[c(2,1)]
    
    scaleFUN <- function(x) sprintf("%.0f", x) # te use round labels on x axis
    
    title_name <-vname
    
    # compute shift function
    sf <- shifthd(data=resid_DF_all, formula = y ~ sex, nboot = 200)
    sf <- round(sf,4)
    sf$difference <- -sf$difference
    sf$ci_lower <- -sf$ci_lower
    sf$ci_upper <- -sf$ci_upper
    
    sf_rev <- sf*-1
    
    # plot shift function
    psf <- plot_sf(sf, plot_theme = 2,symb_fill=(colour_sex))
    
    # change axis labels
    psf <- psf +
      labs(x = "Quantiles of scores",
           y = "Differences M - F")
    
    psf <- add_sf_lab(psf, sf, y_lab_nudge = .1,link_col = c(colour_sex)) + 
      scale_y_continuous(labels = scaleFUN) +
      ggtitle("")+
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.text.x = element_text (size = font_size,color=text_col),#,family="Arial"),
            axis.text.y = element_text (size = font_size,color=text_col),#,family="Arial"), 
            axis.title.y = element_text (size = font_size,color=text_col),#,family="Arial"), 
            axis.title.x = element_text (size = font_size,color=text_col),#,family="Arial"), 
            axis.ticks=element_blank(),
            plot.title=element_text (size =16,color=text_col,hjust=.5),
            # legend.text = element_text (size = 12,color="#888888",family="Helvetica"),
            legend.position="none")
    # print(psf)
    
    # scatterplots
    p <- plot_scat2(resid_DF_all,
                    xlabel = vname,
                    ylabel = "Adjusted scores",
                    alpha = .9,
                    shape = 21,
                    size = 1
                    #                   colour = "grey10",
                    #                   fill = "grey90"
    ) 
    
    p <- plot_dec_links(p, sf,
                        dec_size = 1,
                        md_size = 1.5,
                        add_rect = FALSE,
                        rect_alpha = 0.2,
                        rect_col = "grey50",
                        link_col = c(colour_sex),
                        add_lab = TRUE) # superimposed deciles + rectangle
    
    p <- p + coord_flip() +  scale_colour_manual(values=(colour_sex)) +  scale_fill_manual(values=(colour_sex)) +
      ggtitle(title_name) + 
      scale_y_continuous(labels = scaleFUN) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.text.x = element_text (size = font_size,color=text_col),#,family="Arial"), 
            axis.text.y = element_blank(),#element_text (size = font_size,color="#888888"),#,family="Arial"), 
            axis.title.y = element_blank(),#element_text (size = font_size,color="#888888"),#,family="Arial"), 
            axis.title.x = element_text (size = font_size,color=text_col),#,family="Arial"), 
            axis.ticks=element_blank(),
            plot.title=element_text (size =16,color=text_col,hjust=.4),
            legend.position="none")# flip axes
    
    # print(p)    
    print(cowplot::plot_grid(p, psf, labels=c("A", "B"),ncol = 1, nrow = 2, rel_heights = c(1, 1, 1), label_size = 10, hjust = -0.5, scale=.95))
    

  filename<- paste("Figures/Variance_rousellet/Variance_rousellet_",vname,".pdf",sep="")
  dev.copy(pdf,filename,width=30/2.54, height=18/2.54)
  dev.off()
  
}






# ---------------------------------------------------
# - 3) Does the sex difference in variability in brain structure differs across the lifespan? 
# ---------------------------------------------------


Final_table <- matrix(ncol=17,nrow=length(vnames))
colnames(Final_table) <- c("anatomical_structure","best_age_model","Intercept","Intercept_SE","Intercept_sign","bage","bage_SE","bage_sign","bage_sign2","bsex","bsex_SE","bsex_sign","bsex_sign2","bagesex","bagesex_SE","bagesex_sign","bagesex_sign2")
k=0

for (j in c(set)){ 
  # select data
  k <- k+1
  data1 <- subset(xx,age<30)
  temp <- NULL
  data1$age=data1$age
  data1$x=data1$sex
  data1$y=data1[,j] #raw datapoints
  data1=data1[!is.na(data1$y) & !is.na(data1$x),]
  data1$x=as.factor(data1$x)
  
  #use linear model to covary for cohort
  lm_cohort <- lm(y ~ poly(age,3,raw=T) + Cohort, data=data1)
  X1 <- as.data.frame(predict(lm_cohort, newdata=data1,type="terms"))
  data1$y <- data1$y-X1$Cohort
  
  data1 <- subset(data1, select=-c(Cohort))
 
  
  frmla = y ~ x + age  + FS_Version + Scanner_Strengh
  fit.rf = randomForest(frmla, data=data1, ntree=500, importance=TRUE)
  data1$residuals_randomforest <- abs(data1$y -predict(fit.rf))
  
  # calculate lm on residuals 
  model <- lm(residuals_randomforest ~ poly(age,1) * x, data=data1)

  plot_graph <- function(model) {
    data1.pred <- data1[,c("residuals_randomforest","age","x")]
    X1 <- as.data.frame(predict(model, newdata=data1.pred,type="terms"))
    data1.pred$pred <- predict(model, newdata=data1.pred) #
    
    # get confidence intervals
    preds <- predict(model, newdata=data1.pred, type="terms", se.fit = TRUE)
    fit<-preds$fit
    data1.pred$lwr <- 1.96* rowSums(preds$se.fit[,c(1:2)]) 
    data1.pred$upr <- 1.96* rowSums(preds$se.fit[,c(1:2)])
    
    p1 <- ggplot(data=data1.pred,aes(x=age,y=pred,fill=x)) +  
      scale_colour_manual(values=(colour_sex)) +
      scale_fill_manual(values=(colour_sex)) +
      geom_point(aes(y=residuals_randomforest,col=x)) + 
      geom_line(data=data1.pred,aes(y=pred,x=age,col=x),lwd=1.1) +
      geom_ribbon(data=data1.pred,aes(ymin=pred-lwr,ymax=pred+upr),alpha=0.2)+
      ggtitle(names(data1)[j]) +
      ylab("Residuals")+
      xlab("Age") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.ticks=element_blank(),
            plot.title = element_text (size = font_size+4,color="black",hjust=.5)#family="Arial"), 
      )
    
    print(p1)
    
    filename<- paste("Figures/Variance_age/Linear_Age_effect_noage_corr_",names(data1)[j],"allSites.pdf",sep="")
    dev.copy(pdf,filename)
    dev.off()
    
  }
  
  
 
# ---------------------------------------------------
# - 4) Do females show greater diversity than males in regional volumes?
# ---------------------------------------------------

plot_corr_matrix <- function(xxx,group_name,vnames,analysis_name,p_valueSign){
  
  xxx <- data.frame(xxx[vnames],sex=xxx$sex,age=xxx$age,Cohort=xxx$Cohort,FS_Version=xxx$FS_Version,Scanner_Strength=xxx$Scanner_Strengh)
  xxx <- xxx[complete.cases(xxx),]
  
  # Remove mean sex effects and covariates
  resid_DF=NULL
  for(vname in vnames){
    
    #use linear model to covary for cohort
    df=data.frame(y=xxx[[vname]],sex=xxx$sex,age=xxx$age,Cohort=xxx$Cohort,FS_Version=xxx$FS_Version,Scanner_Strength=xxx$Scanner_Strength)
    lm_cohort <- lm(y ~ poly(age,3,raw=T) + Cohort, data=df)
    X1 <- as.data.frame(predict(lm_cohort, newdata=df,type="terms"))
    df$y <- df$y-X1$Cohort
    
    df <- subset(df, select=-c(Cohort))
    # remove effects of other covariates
    ff=randomForest(y~.,data=df,ntrees=500,proximity=TRUE,importance=TRUE, do.trace=100)
    residuals=ff$y -predict(ff) 
    resid_DF=cbind(resid_DF,residuals)
    }
  
  resid_DF=data.frame(resid_DF)
  names(resid_DF)=vnames
  
  xxx <- cbind(xxx[,c(group_name)],resid_DF)
  
  
  names(xxx)[1] <- "group"
  xxx <- xxx[complete.cases(xxx),]
  
  N_subjects <- nrow(xxx)
  N_F <- nrow(xxx[xxx$group=="F",])
  N_M <- nrow(xxx[xxx$group=="M",])
  
  # Standardized ROIs 
  resid_Ms=scale(xxx[xxx$group=="M",c(2:ncol(xxx))])
  resid_Fs=scale(xxx[xxx$group=="F",c(2:ncol(xxx))])
  
  n_ROIs <- dim(resid_Ms)[2]
  
  # bind scaled residuals males and females
  zi=rbind(resid_Ms,resid_Fs)
  group=c(rep("M",nrow(resid_Ms)),rep("F",nrow(resid_Fs)))
  
  # estiamte correlation matrix 
  cm=(cor(zi[group=="M",]))
  cf=(cor(zi[group=="F",]))
  
  dd=cm-cf # differ correlation matrix (M - F)
  write.table(dd,file=paste0("tables/corr_males_minus_females_",vnames[length(vnames)],".csv"),sep=",",row.names=T)
  
  # dd=cf-cm # M - F
  vv=as.numeric(dd)
  
  ### Permutations to calculate p-value for group difference in correlation
  ha=NULL
  for(b in 1:B)
  {
    GG=sample(group,size=length(group))
    cmb=cor(zi[GG=="M",])
    cfb=cor(zi[GG=="F",])
    hh=as.numeric(cmb-cfb)
    ha=c(ha,hh)
  }
  
  mm=matrix(ha,ncol=length(vv),byrow=TRUE)
  mm=rbind(mm,vv)
  
  # calculate p-values
  pvalue_correlation_stat_all <- NULL
  for (i in 1:(n_ROIs*n_ROIs)){
    temp=mean(mm[,i]>=vv[i])
    pvalue_correlation_stat_all <- cbind(pvalue_correlation_stat_all,temp)
  }
  
  matrix_p_values <- matrix(pvalue_correlation_stat_all,nrow=n_ROIs,ncol=n_ROIs)
  row.names(matrix_p_values) <- vnames
  colnames(matrix_p_values) <- vnames
  write.table(matrix_p_values,file=paste0("tables/corr_males_minus_females_p_",vnames[length(vnames)],".csv"),sep=",",row.names=T)
  
  # plot correlations
  matrix_p_plot <- matrix_p_values
  matrix_p_plot[upper.tri(matrix_p_plot)] <- 1
  n_m_f <- length(matrix_p_plot[matrix_p_plot<.05]==TRUE)
  
  for(i in 1:ncol(matrix_p_plot)){
    for (j in 1:nrow(matrix_p_plot)){
      if(matrix_p_plot[i,j]<p_valueSign){
        print(paste(round(matrix_p_plot[i,j],3),names(xxx[i+1]),names(xxx[j+1]),sep="_"))
        p_sign <- ggplot(xxx,aes(x=xxx[,i+1],y=xxx[,j+1],colour=xxx$group))  
        p_sign <- p_sign + geom_point() + 
          scale_colour_manual(values=colour_sex)+
          xlab(names(xxx[i+1]))+
          ylab(names(xxx[j+1]))+
          ggtitle(group_name) +
          geom_smooth(method="lm") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(),
                axis.ticks=element_blank(),
                plot.title = element_text (size = font_size+4,color="black",hjust=.5)#family="Arial"), 
          )
        
        print(p_sign)
        
        filename<- paste("Figures/Plot_correlations/",group_name,"_",names(xxx[i+1]),"_",names(xxx[j+1]),".pdf",sep="")
        dev.copy(pdf,filename)#,width=8.5/2.54, height=5/2.54)
        dev.off()
        
      }
    }
  }
  
  #######
  #### P_value Fs > Ms (reversed of above)
  ######
  cm=cor(zi[group=="M",])
  cf=cor(zi[group=="F",])
  dd_f=cf-cm
  vv_f=as.numeric(dd_f)
  write.table(dd_f,file=paste0("tables/corr_females_min_males_",vnames[length(vnames)],".csv"),sep=",",row.names=T)
  
  
  ### Permutations 
  ha_f=NULL
  for(b in 1:B){
    GG_f=sample(group,size=length(group))
    cmb_f=cor(zi[GG_f=="M",])
    cfb_f=cor(zi[GG_f=="F",])
    hh_f=as.numeric(cfb_f-cmb_f)
    ha_f=c(ha_f,hh_f)
    }
  
  mm_f=matrix(ha_f,ncol=length(vv_f),byrow=TRUE)
  mm_f=rbind(mm_f,vv_f)
  
  # calculate p-values
  pvalue_correlation_stat_all_f <- NULL
  for (i in 1:(n_ROIs*n_ROIs)){
    temp=mean(mm_f[,i]>=vv_f[i])
    pvalue_correlation_stat_all_f <- cbind(pvalue_correlation_stat_all_f,temp)
  }
  
  matrix_p_values_f <- matrix(pvalue_correlation_stat_all_f,nrow=n_ROIs,ncol=n_ROIs)
  row.names(matrix_p_values_f) <- vnames
  colnames(matrix_p_values_f) <- vnames
  write.table(matrix_p_values_f,file=paste0("tables/corr_females_min_males_p_",vnames[length(vnames)],".csv"),sep=",",row.names=T)
  
  
  # plot correlations
  matrix_p_plot_f <- matrix_p_values_f
  matrix_p_plot_f[upper.tri(matrix_p_plot_f)] <- 1
  n_f_m <- length(matrix_p_plot_f[matrix_p_plot_f<.05]==TRUE)
  
  for(i in 1:ncol(matrix_p_plot_f)){
    for (j in 1:nrow(matrix_p_plot_f)){
      if(matrix_p_plot_f[i,j]<p_valueSign){
        print(paste(round(matrix_p_plot_f[i,j],3),names(xxx[i+1]),names(xxx[j+1]),sep="_"))
        p_sign <- ggplot(xxx,aes(x=xxx[,i+1],y=xxx[,j+1],colour=xxx$group))  
        p_sign <- p_sign + geom_point() + 
          scale_colour_manual(values=colour_sex)+
          xlab(names(xxx[i+1]))+
          ylab(names(xxx[j+1]))+
          geom_smooth(method="lm") +
          theme_bw() +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(),
                axis.ticks=element_blank(),
                plot.title = element_text (size = font_size+4,color="black",hjust=.5)#family="Arial"), 
          )
        
        print(p_sign)
        
        filename<- paste("Figures/Plot_correlations/",group_name,"_",names(xxx[i+1]),"_",names(xxx[j+1]),".pdf",sep="")
        dev.copy(pdf,filename)#,width=8.5/2.54, height=5/2.54)
        dev.off()
        
      }
    }
  }
  
  
  # save number of males>females
  temp_m_f <- data.frame(t(c(n_m_f,n_f_m)))
  names(temp_m_f) <- c("N_males_>_females","N_females_>_males")
  write.table(temp_m_f,file=paste("Tables/sign_p_m_f_",vnames[length(vnames)],".csv",sep=""),sep=",",row.names=F)
  
  
  # pvalues associated with each permutation
  pval_perm=nrow(mm)-apply(mm,2,rank)
  pval_perm=(pval_perm+.5)/(nrow(mm)+1)

  # pvalues associated with each permutation
  pval_perm_f=nrow(mm_f)-apply(mm_f,2,rank)
  pval_perm_f=(pval_perm_f+.5)/(nrow(mm_f)+1)
  
  
  ############
  #### Display levelplot
  ############
  order_var <- 1:n_ROIs #rev(order_var)
  
  corr_difference_order <- dd
  names2 <- vnames 
  
  corr_difference_order <- corr_difference_order[order_var,order_var]
  corr_difference_order <- corr_difference_order[,rev(1:n_ROIs)]
  
  # plot difference matrix + significance
  # select colors for difference plot male females
  matrix_p_values_inv <- matrix_p_values_f
  diag(matrix_p_values_inv) <- 1
  p_table_order_inv_sign <- ifelse(matrix_p_values_inv<=p_valueSign,1,0)
  p_table_order_sign <- ifelse(matrix_p_values<=p_valueSign,1,0)
  p_table_order_new <- p_table_order_sign+(p_table_order_inv_sign)
  p_table_order_new <- p_table_order_new[order_var,(order_var)]
  p_table_order_new[lower.tri(p_table_order_new)] <- 1
  p_table_order_new <- p_table_order_new[,rev(1:n_ROIs)]
  corr_difference_order <- corr_difference_order*(p_table_order_new)
  
  cols = (colour_sex) # Fs
  
  rgb.palette <- colorRampPalette(c(cols[1],"white",cols[2]), space = "rgb") #blue to purple 
  
  # make sure 0 is in middle
  min <- min(corr_difference_order,na.rm=T)
  max <- max(corr_difference_order,na.rm=T)
  ifelse(max>-min, min <--max, max <--min)
  steps.all <- as.numeric((max-min)/99)
  col.seq.all <- seq(min, max , steps.all)
  
  p_table_order_new <- round(corr_difference_order,2)
  
  myPanel <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  p_table_order_new[cbind(x,y)],col="white",cex=.25)
  }
  
  title <- paste0(group_name,"\n",analysis_name,"\nN(Total,F,M) = ",N_subjects,",",N_F,",",N_M)
  p <- levelplot(corr_difference_order,
                 at=col.seq.all,
                 scales=list(tck=0, x=list(rot=90,alternating=2),cex=.3,col=text_col),
                 col.regions=(rgb.palette(120)),
                 main=paste(title, sep=""),
                 xlab="",
                 ylab=(list("",col="#888888")),
                 panel=myPanel,
                 par.settings=list(axis.line=list(col="transparent")),
                 shrink = c(0.9, 0.9)
  )
  
  
  
  print(p)
  
  filename<- paste("Figures/Plot_correlations/Matrix_",group_name,"_",analysis_name,"_",vnames[length(vnames)],"_p_",p_valueSign,".pdf",sep="")
  
  dev.copy(pdf,filename,width=27.4/2.54, height=27.4/2.54)
  dev.off()
  
}


#test two sided
plot_corr_matrix(xx,"sex",names(xx)[c(set)] ,"Anatomical correlation males vs females",.025)





