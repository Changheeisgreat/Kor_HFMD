### set directory

files<-list.files('09MAY/')
no_files=length(files)/2
data_directory=paste0(getwd(),'/09MAY/')

### automatical regression by years
years<-list(2014,2015,2016,2017,2018,2019)
for (year in years){
  x<-sprintf('dataset_%s_epidemic_period.csv',year)
  x<-paste0(data_directory,x)
  df<-read.csv(x)
  model1 <- lm(log(Rt) ~ Cu_hfmd+PHSMs, data =df)
  a<-summary(model1)
  
  model2 <- lm(log(Rt) ~ Cu_hfmd+Rainyseason, data =df)
  b<-summary(model2)
  
  model3 <- lm(log(Rt) ~ Cu_hfmd+PHSMs+Rainyseason, data =df)
  c<-summary(model3)
  
  models<-list(a,b,c)
  
  file<-sprintf('reg_result_%s.txt',year)
  for (model in models){
    cat(year,"\n",toString(model[1]),"\n",file=file,append=TRUE)
    for (no in c(4,8,9,10,11)){
      write.table(model[no],file=file,append=TRUE)
      cat("\n",file=file,append=TRUE)
    }
  }
  
  est<-as.data.frame(c[4])[3,1]
  est_ci<-as.data.frame(confint(model3))[3,1]
  est_CI<-as.data.frame(confint(model3))[3,2]
  cat(
    "ESTIMATE \n",
    est,
    "\n 95CI under \n",
    est_ci,
    "\n 95CI up \n",
    est_CI,
    "\n",
    file=file,append=TRUE
  )
  y<-sprintf('dataset_%s.csv',year)
  y<-paste0(data_directory,y)
  df2<-read.csv(y)
  "S0"<-0.999
  "E0"<-0.00015
  "I0"<-0.00015
  "R0"<-1-(S0+E0+I0)
  
  "sigma"=0.575
  "gamma"=0.6
  
  "S_0R_0"=2.24775
  beta0=(S_0R_0/S0)*gamma
  beta=beta0
  
  reg_coef_mat_phsm=c(0.0,est,est_ci,est_CI)
  np=length(reg_coef_mat_phsm)
  
  phsm=df2$PHSMs
  n=length(phsm)
  
  rain_est<-as.data.frame(c[4])[4,1]
  rain_est_ci<-as.data.frame(confint(model3))[4,1]
  rain_est_CI<-as.data.frame(confint(model3))[4,2]
  
  reg_coef_mat_rain=c(0.0, rain_est, rain_est_ci, rain_est_CI)
  rain=df2$Rainyseason
  
  #Setting for model output 
  S = c()
  E = c()
  I = c()
  R = c()
  t = c()
  N = c()
  Inc  = c()
  Rt = c()
  betaT=c()
  #SEIR model and save output
  for(j in 1:np){
    dt=1
    S[1]<-S0
    E[1]=E0
    I[1]=I0
    R[1]=R0
    t[1]=0
    N[1]=S[1]+E[1]+I[1]+R[1] 
    betaT[1]=beta
    Rt[1]=beta/gamma
    for(i in 2:n){
      #betaT[i]=beta*(exp(reg_coef_mat_hdsc[j] *hday[i]))
      
      #betaT[i]= beta + reg_coef_mat_Hum[j]*Hum[i] - reg_coef_mat_hdsc[j]*hday[i]
      betaT[i]= beta*(exp(reg_coef_mat_rain[1]*(rain[i])))*(exp(reg_coef_mat_phsm[j]*phsm[i]))
    }
    for(i in 2:n-1){
      dS=-betaT[i]*I[i]*S[i]
      dE=betaT[i]*I[i]*S[i]-sigma*E[i]
      dI=sigma*E[i]-gamma*I[i]
      dR=gamma*I[i]
      S[i+1]=S[i]+dt*dS
      E[i+1]=E[i]+dt*dE
      I[i+1]=I[i]+dt*dI
      R[i+1]=R[i]+dt*dR
      t[i+1]=t[i]+dt
    }
    for(i in 2:n){
      Inc[i]=betaT[i]*I[i]*S[i]
      Rt[i]=betaT[i]*S[i]/gamma
    }
    
    matrix=data.frame(time=t,s=S,e=E,i=I,r=R,inc=Inc,rt=Rt)
    form=sprintf('/SEIR_Rt_%s_%s.csv',year,j)
    form=paste0(getwd(),form)
    write.csv(matrix,file=form)
  }
  file1=sprintf("SEIR_Rt_%s_1.csv",year)
  file2=sprintf("SEIR_Rt_%s_2.csv",year)
  file3=sprintf("SEIR_Rt_%s_3.csv",year)
  file4=sprintf("SEIR_Rt_%s_4.csv",year)
  SEIR_without_phsm_rain1 <-read.csv(file1)
  SEIR_with_phsm_rain2 <-read.csv(file2)
  SEIR_with_phsm_rain3 <-read.csv(file3)
  SEIR_with_phsm_rain4 <-read.csv(file4)
  
  final=sprintf("Result_%s.csv",year)
  result=data.frame(df2$days,
                    df2$PHSMs,
                    df2$Rainyseason,
                    SEIR_without_phsm_rain1$i,
                    SEIR_with_phsm_rain2$i,
                    SEIR_with_phsm_rain3$i,
                    SEIR_with_phsm_rain4$i
                    )
  write.csv(result,file=final)
}



  