n_i=4;m=100;
coef_pi.true=c(0.5,-1,0.5,-0.25,-0.1);pi_model='center'

n=n_i*m
set.seed(12)
Z1=rnorm(n,0,1)
Z2=rnorm(n,0,1)
Z3=rnorm(n,0,1)
Z4=rnorm(n,0,1)	
X1=exp(Z1/2)	
X2=Z2/(1+exp(Z1))+10	
X3=(Z1*Z3/25+0.6)^3
X4=(Z1+Z4+20)^2
inter.info.raw=cbind(rep(1:n,each=n_i),rep(1:n_i,m)+n_i*rep(0:(m-1),each=n_i^2))
inter.info=data.frame(inter.info.raw[inter.info.raw[,1]!=inter.info.raw[,2],])
colnames(inter.info)=c('id','inter_id')
last_child_each_house_id=cumsum(rep(n_i,times=m))
first_child_each_house_id=c(1,last_child_each_house_id+1)[1:m]

##----- Real propensity and outcome model -------##
logit.f_propen=coef_pi.true[1]+coef_pi.true[2]*Z1+coef_pi.true[3]*Z2+coef_pi.true[4]*Z3+
  coef_pi.true[5]*Z4		
#if (sim_number<200){hist(p_A1.true,probability=T)}
b=rnorm(m,0,1)
b.v=rep(b,each=n_i)
f_propen_with_b.v=exp(logit.f_propen+b.v)/(1+exp(logit.f_propen+b.v))
A=rbinom(n,1,prob=f_propen_with_b.v)
A_inter_sum=apply(as.matrix(1:n),1,cal_sum_A_inter.fun,inter.info,A)
# if(inter_model==1){
#   A.m=matrix(A,nrow=n_i,byrow=F);A_friend.m=A.m[n_i:1,]
#   A_friend=as.vector(A_friend.m) # group 1 first, then group 2 etc..
#   A_inter=A_friend
# }else if (inter_model==2){
#   #			A_inter=as.numeric(A_inter_sum==(n_i-1))
#   #			A_inter=A_inter_sum
  A_inter=A_inter_sum/n_i
# }
  
# coef_y.true <- coef_y_inter.true <- c(2,-1.5,-2.7,3,-1,0.5,6,1,8) # order is 1,Z1,Z2,Z3,Z4,A,A_friend,A*Z1,A_friend*Z2
coef_y.true <- coef_y_inter.true <- c(2,-1,-2.7,3,-1,0.5,6,1,8) # order is 1,Z1,Z2,Z3,Z4,A,A_friend,A*Z1,A_friend*Z2

mu.true=coef_y.true[1]+coef_y.true[2]*Z1+coef_y.true[3]*Z2+
  coef_y.true[4]*Z3+coef_y.true[5]*Z4+
  coef_y.true[6]*A+coef_y.true[7]*A_inter+
  coef_y.true[8]*A*Z1+coef_y.true[9]*A_inter*Z2
epsilon=rnorm(n,0,1)
y=mu.true+epsilon
mydata=data.frame(y,A,A_inter,Z1,Z2,Z3,Z4,X1,X2,X3,X4,f_propen_with_b.v)
mydata$group <- rep(1:m, each = n_i)

library(inferference)
compare_this <- interference(formula =  y | A ~ Z1 + Z2 + Z3 + Z4 + (1|group) | group,
             allocations = c(.1, .7),
             data = mydata,
#              model_method = 'oracle',
#              model_options = list(fixed.effects = coef_pi.true,
#                                   random.effects = 1),
             causal_estimation_options = list(variance_estimation ='naive'),
             method = 'simple')
compare_this$estimates[3:4, ]
