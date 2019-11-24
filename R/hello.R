#'Linear Regression Fitting
#'
#' lr is used for fitting linear regression models with numeric or categorical covariates.
#' It can be used to estimate regression coefficients with least square method and make
#' inference based on relevant t-test and F-test (sequential F-test can be realized by F_test).
#'
#'@import stats
#'
#'@param formula an object of class \code{formula}: a symbolic form of the model to be fitted
#'which should contain both respond and covariate variables. The details of model
#'specification form are given under "Details".
#'@param data an optional data frame,list or environment containing the variables in the model.
#'If not found in data or by default, the variables would be taken from environment where ls is
#'called.
#'@param coding optional,the method to be used for fitting categorical covaraites in the model.
#'There are two options: coding = 'reference' (by default) or 'means' (see 'Details').
#'@param intercept logical. If TRUE (by default), the corresponding fitting model will contain
#'intercept term.
#'@param reference an optional categorical covariate group which would be considered as reference
#'group when \code{coding} = 'reference' (by default). If not specified, model fitting will take
#'first group occured among the corresponding covariate as reference.
#'
#'@keyword lr linearregression
#'
#'@details Models for lr are specified symbolically like 'y ~ x1 + x2'. A typical model has the form
#''response ~ covariates' where response is the numeric response vector and covariates are a
#'series of numeric or categorical terms which specifies a linear predictor for response. A
#'covariate specification of the form 'x1 + x2' indicates covariates will contain all the
#'observations in 'x1' and 'x2'.
#'
#'In addition, regarding to the variable in a data frame contianed in model formula, they can be
#'either expressed as mtcars$mpg or mpg with \code{data} = mtcars where 'mtcars' is the name of a
#'data frame, and 'mpg' is the variable name in 'mtcars'.
#'
#'\code{coding} method is used to cope with model containing categorical covariates. Two common methods:
#''cell reference coding' and 'cell means coding' are supported by setting \code{coding} = 'reference'
#'(by default) or 'means'. Former one takes one group of the corresponding categorical covariate as a
#'reference group and reserve intercept in the model, while latter one just eliminates intercept in the
#'model.
#'
#'@examples
#'## A example for numeric respond and covariate variables
#'y = c(23,24,26,37,38,25,36,40)
#'x1 = c(1,2,3,4,5,6,7,8)
#'x2 = c(23,32,34,20,24,56,34,24)
#'result = lr(y~x1+x2)
#'
#'## A example for numeric respond and categorical variables
#'x3 = c("M","F","F","U","M","F","F","U")
#'result = lr(y~x1+x2+x3)
#'
#'## A example for variables in data frame
#'result = lr(mpg~cyl+disp+hp, data=mtcars)
#'result = lr(mtcars$mpg~mtcars$cyl+mtcars$disp+mtcars$hp)
#'
#'@usage
#'lr(formula,data,coding='reference',intercept=TRUE,reference=1)
#'
#'@aliases factorize linear_regression lr.fit
#'
#'@export lr
#'@export lr.fit
#'@export factorize
#'
#'@return linear regression fitting coefficients and reference results

lr=function(formula,data,coding='reference',intercept=TRUE, reference=1){
  if(length(formula)!=3) stop("The model form is incorrect.")
  y=try(get(as.character(formula[[2]])),silent=TRUE)
  if(any(strsplit(as.character(formula[[2]]), "")[[1]]=="$")){
    y=try(get(as.character(formula[[2]])[[3]],get(as.character(formula[[2]])[[2]])),silent=TRUE)
  }
  if(class(y)=="try-error"){
    if(length(data)>0){
      a=as.character(formula[[2]])
      y=try(get(a,data),silent=TRUE)
    }
    if(class(y)=="try-error") stop("Cannot find data to fit model.")
  }
  x_names=labels(terms(formula))
  x_num=length(x_names)
  ###add intercept
  if(intercept){
    x=matrix(1,length(y),1)
    covariate_name=c("intercept")
  }else{
    x=vector()
  }
  ###form x_matrix
    for (i in 1:x_num){
      x1=try(get(x_names[i]),silent=TRUE)
      if(any(strsplit(x_names[i], "")[[1]]=="$")){
        position=which(strsplit(x_names[i], "")[[1]]=="$")
        variable=substring(x_names[i],position+1,nchar(x_names[i]))
        environment=substring(x_names[i],1,position-1)
        x1=try(get(variable,get(environment)),silent=TRUE)
      }
      if(class(x1)=="try-error"){
        if(length(data)>0){
          x1=try(get(x_names[i],data),silent=TRUE)
        }
        if(class(x1)=="try-error") stop("Cannot find data to fit model.")
      }
      if(typeof(x1)!="double"){
        names=unique(x1)
        x1=factorize(x1)
        if(coding=='reference'){
          if(reference){
            x1=x1[,-1]
          }else if(any(colnames(x1)==reference)){
            x1=x1[,-which(colnames(x1)==reference)]
          }
        }
        covariate_name=c(covariate_name,colnames(x1))
      }else{
        covariate_name=c(covariate_name,x_names[i])
      }
      x=cbind(x,x1)
    }
    colnames(x)=covariate_name
  if(coding=='means'){
    x=x[,-'intercept']
  }
  result=lr.fit(x,y)
  return(result)
}



lr.fit=function(x,y){
  if(length(dim(x))==0){
    n_row=length(x)
    n_col=1
  }else{
    n_row=nrow(x)
    n_col=ncol(x)
  }
  ### 1 estimate
  ###estimating the coefficients using Least Squares Method
  ####1.1 determine x is full-column-rank
  x_matrix=t(x)%*%x
  x_matrix_inverse=try(solve(t(x)%*%x),silent = TRUE)
  if(class(x_matrix_inverse)=="try-error") stop("Multiple covariates 'x's can't be linear dependent.")
  ###1.2 estimated coefficients,i.e. beta_hat
  beta_estimates=x_matrix_inverse%*%t(x)%*%y
  ###1.3 fitted-values,i.e. Y_hat
  hat_matrix=x%*%x_matrix_inverse%*%t(x)
  fitted_values=hat_matrix%*%y
  ###1.4 residual,i.e. Y-Y_hat
  residual=y-fitted_values
  ###1.5 SSY, SSR and SSE
  SSE=sum(residual^2)
  SSR=sum((fitted_values-mean(y))^2)
  SSY=sum((y-mean(y))^2)
  ###1.6 rank
  df_SSE=n_row-n_col
  df_SSR=n_col-1
  df_SSY=n_row-1
  ###1.7 R-squared and adjusted R-squared
  R_square=SSR/SSY
  adj_R_square=1-SSE/df_SSE/(SSY/df_SSY)
  #####2 inference
  ###2.0 overall F-test
  F_value=SSR/df_SSR/(SSE/df_SSE)
  F_p_value=pf(F_value, df_SSR, df_SSE, lower.tail=F)
  ###2.1 estimated-variance of estimated-coefficients,i.e. Var_hat(beta_hat)
  beta_variance=SSE/df_SSE*x_matrix_inverse
  beta_SE=sqrt(diag(beta_variance))
  T_value=beta_estimates/beta_SE
  T_p_value=2*pt(-abs(T_value),df_SSE)
  beta_CI_lower=beta_estimates-qt(0.05/2,df_SSE)*beta_SE
  beta_CI_upper=beta_estimates+qt(0.05/2,df_SSE)*beta_SE
  ###3 formatting results
  beta_result=cbind(beta_estimates,beta_SE,T_value,T_p_value,beta_CI_lower,beta_CI_upper)
  colnames(beta_result)=c("Estimates","SE", "T value", "p value", "95%CI"," ")
  Sum_of_squares=matrix(c(SSE,df_SSE,SSR,df_SSR,SSY,df_SSY),3,2)
  rownames(Sum_of_squares)=c("SSE","SSR","SSY")
  colnames(Sum_of_squares)=c("SS","df")
  F_test_and_R_square=matrix(c(F_value,F_p_value,R_square,adj_R_square),1,4)
  colnames(F_test_and_R_square)=c("F value","p value","R^2","adjusted R^2")
  colnames(fitted_values)="y_hat"
  colnames(residual)="y-y_hat"
  result=list(beta_result,F_test_and_R_square,fitted_values,residual,
              Sum_of_squares,x_matrix_inverse,beta_variance)
  names(result)=c("Coefficients","F test and R square","Fitted values","Residuals",
                  "Sum of Squares","X inverse matrix","Coefficients Variance")
  return(result)
}

factorize=function(x){
  n=length(x)
  q=unique(x)
  group_num=length(q)
  x_matrix=matrix(0,n,group_num)
  for(i in 1:group_num){
    x_matrix[,i]=x==q[i]
  }
  colnames(x_matrix)=q
  return(x_matrix)
}

