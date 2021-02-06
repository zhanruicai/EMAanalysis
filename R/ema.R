

#' Get the data summarize for the EMA variable
#'
#' This function returns the summary of the EMA variable, which include the mean, variance, autocorrelation,
#' \eqn{\lambda_o} and \eqn{\lambda_s} for each subject. The function will return a n*6 matrix, where n is the
#' number of subjects. In each row, we record the subject's ID, and his/her five summary variables listed above.
#'
#' @param ID a vector contain all the IDs for each measurement in x.
#' @param x a vector contain all the measurements for all the subjects.
#'
#' @return Matrix of dimension n*6, where n is the number of subject. The six columns are: ID,
#' mean, variance, autocorrelation, \eqn{\lambda_o} and \eqn{\lambda_s}, respectively.
#'
#' @examples
#' ID = rep(seq(1:10), each = 8)
#' x = rnorm(80)
#' summarize_EMA_data(ID, x)
#'
#' @references
#' \insertRef{buu2020association}{EMAanalysis}
#'
#' @export
summarize_EMA_data = function(ID, x)
{
  m1 = tapply(x, INDEX = ID, function(aa){if(length(na.omit(aa))<1) return(0) else return(mean(aa, na.rm = T))})
  v1 = tapply(x, INDEX = ID, function(aa){if(length(na.omit(aa))<=1) return(NA) else return(var(aa, na.rm = T))})
  lambda_o = tapply(x, INDEX = ID, function(aa){if(length(na.omit(aa))<=1) return(NA) else return(get_lambda_o(aa))})

  lambda_s = tapply(x, INDEX = ID, function(aa){if(length(na.omit(aa))<=1) return(NA) else return(get_lambda_s(aa))})
  atcor = tapply(x, INDEX = ID, function(aa){if(length(na.omit(aa))<=1) return(NA) else return(autocorrelation(aa))})
  data.frame(ID = names(m1), mean=m1, variance=v1, autocorrelation = atcor, lambda_o, lambda_s)
}


#' Calculate the first order auto-correlation
#'
#' This function calculate the first order auto-correlation, i.e., the correlation between consecutive observations, for the EMA data in one subject.
#' A higher value indicates that measured variable (for example, emotion) is sustained over time and shows less homeostatic recovery.
#'
#' @export
#' @param x numeric vector of measurements for one subject
#'
#' @examples
#' x <- 1:10
#' autocorrelation(x)
#'
#' @references
#' \insertRef{buu2020association}{EMAanalysis}
autocorrelation = function(x)
{
  x = na.omit(x)
  if(length(x)<=1) return(0)
  xbar = mean(x)
  if(sum((x-xbar)^2)==0) return(0)
  a1 = x[-1]
  a2 = x[-length(x)]
  return(sum((a1-xbar)*(a2-xbar))/sum((x-xbar)^2))
}

#' Calculate \eqn{\lambda_o}
#'
#' This function calculate \eqn{\lambda_o}, an instability measure of the EMA data in one subject.The data vector is
#' first standardized. Then we calculate \eqn{\lambda_o} as
#' \deqn{\lambda_o = n_o/t_o}
#' where \eqn{n_o} is the number of transitions from outside to inside the standard range, and \eqn{t_o}
#' is the number of assessments with a score outside the standard range. \eqn{\lambda_o} indicates the conditional probability of transitioning to inside the standard range, given
#' being outside.
#'
#' @export
#' @param x numeric vector of measurements for one subject
#'
#' @examples
#' x <- 1:10
#' get_lambda_o(x)
#'
#' @references
#' \insertRef{buu2020association}{EMAanalysis}
get_lambda_o = function(x)
{
  x = na.omit(x)
  if(length(x)<=1) return(0)
  x = x-mean(x)
  a1 = abs(x)>sd(x)
  if(sum(a1)==0) return(0)
  l1 = length(x)
  if(a1[l1]==0){
    n0 = sum(a1[2:l1]-a1[1:(l1-1)]<0)
  }else{
    n0 = sum(a1[2:l1]-a1[1:(l1-1)]<0)+1
  }
  n0/sum(a1)
}


#' Calculate \eqn{\lambda_s}
#'
#' This function calculate \eqn{\lambda_s}, an instability measure of the EMA data in one subject.The data vector is
#' first standardized. Then we calculate \eqn{\lambda_s} as
#' \deqn{\lambda_s = n_s/t_s}
#' where \eqn{n_s} is the number of transitions from inside to outside the standard range, and \eqn{t_s}
#' is the number of assessments with a score inside the standard range. \eqn{\lambda_s} measures the conditional probability of transitioning from inside to outside.
#'
#' @export
#' @param x numeric vector of measurements for one subject
#'
#' @examples
#' x <- 1:10
#' get_lambda_s(x)
#'
#' @references
#' \insertRef{buu2020association}{EMAanalysis}
get_lambda_s = function(x)
{
  x = na.omit(x)
  if(length(x)<=1) return(0)
  x = x-mean(x)
  a1 = abs(x)<=sd(x)
  if(sum(a1)==0) return(0)
  l1 = length(x)
  if(a1[l1]==0){
    ns = sum(a1[2:l1]-a1[1:(l1-1)]<0)
  }else{
    ns = sum(a1[2:l1]-a1[1:(l1-1)]<0)+1
  }
  ns/sum(a1)
}


