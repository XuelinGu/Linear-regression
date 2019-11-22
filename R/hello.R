#'Hello, world!
#'
#' This is an example function named 'hello'
#' which prints 'Hello, world!'.
#'
#' Some useful keyboard shortcuts for package authoring:
#'
#'   Install Package:           'Ctrl + Shift + B'
#'   Check Package:             'Ctrl + Shift + E'
#'   Test Package:              'Ctrl + Shift + T'
#'
#'@param x,y,formula input value
#'@keyword linear_regression
#'@details value all components of linear_regression results are listed below:
#'x matrix
#'formula y ~ x1+x2
#'@import stats
#'@examples
#'x1 <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#'y <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#'x2 <- c(2, 10, 20, 2, 10, 20, 2, 10, 20, 2)
#'x <- c(x1, x2)
#'formula=y~x1+x2
#'lm.D9 <- linear_regression(formula)
#'lm.D90 <- estimate(x1,y)
#'lm.D99 = factorize(x)
#'square(x)

#'@usage
#'linear_regression(formula=y~x1+x2)
#'estimate(x,y)
#'factorize(x)
#'square(x)
#'@aliases hello factorize linear_regression estimate square
#'@export factorize linear_regression estimate hello square
#'@return the matrix_x

hello <- function() {
  print("Hello, world!")
}
square=function(x){
  return(x^2)
}

linear_regression=function(formula=y~x1+x2){
  if(length(formula)!=3) stop("The model form is incorrect.")
  dependent=try(get(as.character(formula[[2]])),silent=TRUE)
  if(class(dependent)=="try-error") stop("Could not find data to fit model.")
  covariate_name=labels(terms(formula))
  covariate_num=length(covariate_name)
  x=matrix(1,length(dependent),1)
  for (i in 1:covariate_num){
    x1=get(covariate_name[i])
    x=cbind(x,x1)
  }
  return(x)
}

estimate=function(x,y){
  n_row=nrow(x)
  n_col=ncol(x)
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
  rank_SSE=n_row-n_col
  rank_SSR=n_col-1
  rank_SSY=n_row-1
  ###1.7 R-squared and adjusted R-squared
  R_square=SSR/SSY
  adj_R_square=1-SSE/rank_SSE/(SSY/rank_SSY)
  #####2 inference
  ###2.0 overall F-test
  F_value=SSR/rank_SSR/(SSE/rank_SSE)
  F_p_value=pf(F_value, rank_SSR, rank_SSE, lower.tail=F)
  ###2.1 estimated-variance of estimated-coefficients,i.e. Var_hat(beta_hat)
  beta_variance=SSE/rank_SSE*x_matrix_inverse
  beta_SE=sqrt(diag(beta_variance))
  T_value=beta_estimates/beta_SE
  T_p_value=2*pt(-abs(T_value),rank_SSE)
  beta_CI_lower=beta_estimates-qt(0.05/2,rank_SSE)*beta_SE
  beta_CI_lower=beta_estimates+qt(0.05/2,rank_SSE)*beta_SE
  ###3 formatting results
  return(list(beta_estimates,beta_SE,T_value,T_p_value,beta_CI_lower,beta_CI_lower,fitted_values,residual,
              SSE,SSR,SSY,rank_SSE,rank_SSR,rank_SSY,F_p_value,R_square,adj_R_square, x_matrix_inverse))
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
