###########################
# Doubly robust estimator #
###########################
#rm(list=ls())
library(lme4)
setwd('/Users/Lan/Desktop/LAN/study/statistics/biostatistics/research/my paper/hajek')

#########################################
# Function to sum a section of a vector #
#########################################

sum_part_fun=function(from_to,v) {
	sum(v[from_to[1]:from_to[2]])
}

########################################
# Function to repeat n_i for n_i times #
########################################

rep_n_i_fun=function(n_i){
	return(rep(1/n_i,times=n_i))
}

########################################
# Function to cauculate the propensity #
#         of interference set for      #
#             each individual          #
########################################
cal_f_inter.fun=function(target.id,inter.info,A,f.v){
	f_inter.long=
		((f.v^A)*((1-f.v)^(1-A)))[inter.info$inter_id]
	ind=(inter.info$id==target.id)
	f_inter=ifelse(sum(ind)==0,1,prod(f_inter.long[ind]))
	return(f_inter)
}

########################################
# Function to cauculate the total trt  #
#         of interference set for      #
#             each individual          #
########################################

cal_sum_A_inter.fun=function(target.id,inter.info,A){
	sum_A_inter=sum(A[inter.info[inter.info$id==target.id,2]])
	return(sum_A_inter)
}

##########################################
# Functions to find the interference set #
##########################################

find_inter_fun=function(target_id,new_id,new_housecode){
	target_housecode=new_housecode[new_id==target_id]
	kids_id_same_house=new_id[new_housecode==target_housecode]
	temp=cbind(target_id,kids_id_same_house)
	index=(temp[,1]!=temp[,2])
	if(is.null(index)){
		break
	}else{
		return(temp[index==T,])	
	}
}

############################################
# Functions to integrate the random effect #
#        b out of the propensity score     #
############################################

f_a_inter_integrand.f=function(b,A,vec,target.id,ind,sd_b){
	dnorm(b,sd=sd_b)*prod((exp(A*(vec+b))/(1+exp(vec+b)))[c(ind,target.id)] )
}

partial_l_v_gamma.f=function(b,A,vec,target.id,ind,sd_b,L){
	dnorm(b,sd=sd_b)*prod((exp(A*(vec+b))/(1+exp(vec+b)))[c(ind,target.id)] )*
		(sum(  (L*exp(vec+b)/(1+exp(vec+b)))[c(ind,target.id)] ))
}

partial_l_v_sigmab2.f=function(b,A,vec,target.id,ind,sd_b){
	dnorm(b,sd=sd_b)*prod((exp(A*(vec+b))/(1+exp(vec+b)))[c(ind,target.id)] )*
		(b^2)
}

to_integ_fun_for_d_lv_dgamma=function(dummy,target.id,inter.info,A,vec,ind,sd_b,covariates.m,f_a_inter){
	L=covariates.m[,dummy]
	d_lv_d_gamma=sum((A*L)[c(ind,target.id)])-integrate(Vectorize(partial_l_v_gamma.f,vectorize.args='b'),lower =-4,upper =4,
			vec =vec,A=A,target.id=target.id,
			ind=ind,sd_b=sd_b,L=covariates.m[,dummy])$value/f_a_inter
	return(d_lv_d_gamma)
}

to_integ_fun=function(target.id,inter.info,A,logit.f_propen,sd_b,num_covariates=num_covariates,
				covariates.m=NULL){
	ind=inter.info$inter_id[inter.info$id==target.id]
	f_a_inter=integrate(Vectorize(f_a_inter_integrand.f,vectorize.args='b'),lower =-4,upper =4,
			vec =logit.f_propen,A=A,target.id=target.id,ind=ind,sd_b=sd_b)$value
	d_lv_d_sigmab2=integrate(Vectorize(partial_l_v_sigmab2.f,vectorize.args='b'),lower =-4,upper =4,
		vec =logit.f_propen,A=A,target.id=target.id,ind=ind,sd_b=sd_b)$value/(2*(sd_b)^4*f_a_inter)-
		1/(2*(sd_b)^2)
	d_lv_d_gamma.v=apply(as.matrix(1:(num_covariates+1)),1,to_integ_fun_for_d_lv_dgamma,
			vec =logit.f_propen,A=A,target.id=target.id,ind=ind,
			sd_b=sd_b,covariates.m=covariates.m,f_a_inter=f_a_inter)
					# the '1' in 'num_covariates+1' is for the intercept
	return(c(f_a_inter,d_lv_d_sigmab2,d_lv_d_gamma.v))
}

######################################################
# Function to fit a regression for the outcome model #
#  Depending on the covariates, this model would be  #
#       correctly specified or mis specified         #
######################################################

fit_lm_fun=function(Iflm.tru,n_i,m,n,covariates.m,alpha,y,A,A_inter,inter_size,dia,mydata,cal_CE){

	if(dia==T){
		temp.J.list=apply(as.matrix(n_i),1,rep_n_i_fun)
		J.vec=matrix(unlist(temp.J.list),ncol=1)
		y.lm= lm(formula = y ~A+A_inter+agecat1_2+agecat2_3+ newmom_edu + newfloor + 
      		season+sanitation +watersourcecat+ breastfed + sexcat,
			data=mydata)
		hat_beta=matrix(summary(y.lm)$coefficient[,1],ncol=1)
		Z.m=cbind(1,A,A_inter,covariates.m)
		Z_A1.m=cbind(1,1,alpha^(1/J.vec-1),covariates.m)
		Z_A0.m=cbind(1,0,alpha^(1/J.vec-1),covariates.m)
		Z_A.m=cbind(1,alpha,alpha^(1/J.vec-1),covariates.m)
#		Z_A1.m=cbind(1,1,alpha*inter_size*J.vec,covariates.m)
#		Z_A0.m=cbind(1,0,alpha*inter_size*J.vec,covariates.m)
#		Z_A.m=cbind(1,alpha,alpha*inter_size*J.vec,covariates.m)
		Z_obsA1_minus_ij.m=cbind(1,1,A_inter,covariates.m)
		Z_obsA0_minus_ij.m=cbind(1,0,A_inter,covariates.m)
	}else{
		J.vec=matrix(rep(1/n_i,times=n),ncol=1)
		L1=covariates.m[,1];L2=covariates.m[,2];L3=covariates.m[,3];L4=covariates.m[,4]
		y.lm=lm(y~L1+L2+L3+L4+A+A_inter+A*L1+A_inter*L2)
		hat_beta=matrix(summary(y.lm)$coefficient[,1],ncol=1)
		Z_A.m =cbind(1,L1,L2,L3,L4,alpha,alpha*inter_size*J.vec,alpha*L1,(alpha*inter_size*J.vec)*L2)
		Z_A1.m=cbind(1,L1,L2,L3,L4,1,alpha*inter_size*J.vec,1*L1,(alpha*inter_size*J.vec)*L2)
		Z_A0.m=cbind(1,L1,L2,L3,L4,0,alpha*inter_size*J.vec,0*L1,(alpha*inter_size*J.vec)*L2)
		Z.m=cbind(1,L1,L2,L3,L4,A,A_inter,A*L1,A_inter*L2)
		Z_obsA1_minus_ij.m=cbind(1,L1,L2,L3,L4,1,A_inter,1*L1,A_inter*L2)
		Z_obsA0_minus_ij.m=cbind(1,L1,L2,L3,L4,0,A_inter,0*L1,A_inter*L2)
	}
	temp=t(Z.m)*(summary(y.lm)$residuals)
	hat_M=summary(y.lm)$cov.unscaled%*%(temp%*%t(temp))%*%summary(y.lm)$cov.unscaled
		# hat var(X'(Y-X\beta))

	hat_yij_reg =Z_A.m%*%hat_beta
	hat_y1ij_reg=Z_A1.m%*%hat_beta
	hat_y0ij_reg=Z_A0.m%*%hat_beta
	hat_y1_obsA_.ij=Z_obsA1_minus_ij.m%*%hat_beta
	hat_y0_obsA_.ij=Z_obsA0_minus_ij.m%*%hat_beta
	hat_y_obsA.ij=Z.m%*%hat_beta
	hat_var_y_reg_part1=cal_hat_var_reg_part1_fun(Z_withalpha.m=Z_A.m,J.vec=J.vec,hat_M=hat_M,m=m)
	hat_var_y1_reg_part1=cal_hat_var_reg_part1_fun(Z_withalpha.m=Z_A1.m,J.vec=J.vec,hat_M=hat_M,m=m)
	hat_var_y0_reg_part1=cal_hat_var_reg_part1_fun(Z_withalpha.m=Z_A0.m,J.vec=J.vec,hat_M=hat_M,m=m)
#	hat_var_y_reg_part1=t(J.vec)%*%Z_A.m%*%hat_M%*%t(Z_A.m)%*%J.vec/(m^2)
#	hat_var_y1_reg_part1=t(J.vec)%*%Z_A1.m%*%hat_M%*%t(Z_A1.m)%*%J.vec/(m^2)
#	hat_var_y0_reg_part1=t(J.vec)%*%Z_A0.m%*%hat_M%*%t(Z_A0.m)%*%J.vec/(m^2)
#	hat_var_y_reg_part1=t(J.vec)%*%Z.m%*%hat_M%*%t(Z.m)%*%J.vec/(m^2)
#	hat_var_y1_reg_part1=t(J.vec)%*%Z_obsA1_minus_ij.m%*%hat_M%*%t(Z_obsA1_minus_ij.m)%*%J.vec/(m^2)
#	hat_var_y0_reg_part1=t(J.vec)%*%Z_obsA0_minus_ij.m%*%hat_M%*%t(Z_obsA0_minus_ij.m)%*%J.vec/(m^2)
	if(cal_CE==T){
		Z_withalpha.m.list=list(Z_A.m=Z_A.m,Z_A1.m=Z_A1.m,Z_A0.m=Z_A0.m)
	}else{
		Z_withalpha.m.list=NULL
	}
	return(list(hat_yij_reg=hat_yij_reg,hat_y1ij_reg=hat_y1ij_reg,hat_y0ij_reg=hat_y0ij_reg,
		hat_y_obsA.ij=hat_y_obsA.ij,hat_y1_obsA_.ij=hat_y1_obsA_.ij,hat_y0_obsA_.ij=hat_y0_obsA_.ij,
		hat_var_y_reg_part1=hat_var_y_reg_part1,hat_var_y1_reg_part1=hat_var_y1_reg_part1,
		hat_var_y0_reg_part1=hat_var_y0_reg_part1,hat_M=hat_M,J.vec=J.vec,Z_withalpha.m.list=Z_withalpha.m.list))

}

######################################
# Function to calculate the part1 of #  
# hat var for regression estimators  #
######################################

cal_hat_var_reg_part1_fun=function(Z_withalpha.m,J.vec,hat_M,m){
	hat_var_reg_part1=t(J.vec)%*%Z_withalpha.m%*%hat_M%*%t(Z_withalpha.m)%*%J.vec/(m^2)
	return(hat_var_reg_part1)
}


##############################################
# Function to fit the propensity score model #
##############################################

fit_propen_fun=function(Ifpropen.tru,covariates.m,A,A_inter,inter_model,inter_size,n_i,m,n,dia,mydata){
	if(dia==T){
		p.lm= glmer(formula = A ~agecat1_2+agecat2_3+ newmom_edu + newfloor + 
      		season+sanitation +watersourcecat+ breastfed + sexcat+(1|housecode),
			family=binomial,data=mydata)
	}else{
		L1=covariates.m[,1];L2=covariates.m[,2];L3=covariates.m[,3];L4=covariates.m[,4]
	      household_id=rep(1:m,each=n_i)
	      p.lm=glmer(A~L1+L2+L3+L4+(1|household_id),family=binomial)
	}
	hat_gamma=matrix(fixef(p.lm),ncol=1)
	hat_logit_p=cbind(1,covariates.m)%*%matrix(hat_gamma,ncol=1)
#	temp=summary(p.lm)@REmat
#	est_hat_ranef_sd=as.numeric(temp[,dim(summary(p.lm)@REmat)[2]])
	est_hat_ranef_sd=sqrt(as.numeric(VarCorr(p.lm)))
	return(list(hat_logit_p,est_hat_ranef_sd))
}

###############################################
#   Function to calculate the sandwich part   # 
# in the variance estimator for IPW estimator #
###############################################

cal_ipwvar_w_sandwich_fun=function(hat_psi_v_ipw,d_lv.m,hat_inv_V.m,m,num_covariates){

	hat_E_psi_v_ipw_d_lv=rowSums(matrix(rep(hat_psi_v_ipw,
		times=num_covariates+2),nrow=num_covariates+2,byrow=T)*d_lv.m)/m
	hat_HV_inv_H_ipw=sum((hat_E_psi_v_ipw_d_lv%*%hat_inv_V.m)*
		hat_E_psi_v_ipw_d_lv)
	hat_var_y_ipw=sum(hat_psi_v_ipw^2)/(m^2)-hat_HV_inv_H_ipw/m
	return(hat_var_y_ipw)
}

####################################################
# Function to calculate hat Y Y1 Y0, one at a time #
####################################################
cal_Y_Y1_Y0_est_fun=function(alpha,Y_Y1_Y0,A,A_inter,y,n_i,m,n,first_child_each_house_id,
			last_child_each_house_id,lm_info.mu_tru,lm_info.mu_mis,
			propen_info.pi_tru,propen_info.pi_mis,mean.temp.propen.pi.tru.m,
			mean.temp.propen.pi.mis.m,pi_inter,num_covariates,dia){

	pi_A=alpha^A*(1-alpha)^(1-A)
	hat_f_a_inter.tru=mean.temp.propen.pi.tru.m[1,]

	if(dia==T){
		d_lv.tru.m=mean.temp.propen.pi.tru.m[-1,last_child_each_house_id]
	}else{
		d_lv.tru.m=mean.temp.propen.pi.tru.m[-1,4*(0:(m-1))+1]
	}
	#-----------------------------------------------#
	#      estimators only using correct model      #
	#-----------------------------------------------#
	hat_E_d_lv_d_lv.tru.m=d_lv.tru.m%*%t(d_lv.tru.m)/m
	hat_E_d_lv_d_lv_inv.tru.m=solve(hat_E_d_lv_d_lv.tru.m)

	if(Y_Y1_Y0=='Y'){
		A_pi_over_f.tru=pi_A*pi_inter/hat_f_a_inter.tru;
		hat_yij_reg.mu_tru=lm_info.mu_tru[[1]];hat_y_obsA.ij.tru=lm_info.mu_tru[[4]];
			hat_var_y_reg_part1.mu_tru=lm_info.mu_tru[[7]]
	}else if (Y_Y1_Y0=='Y1'){
		A_pi_over_f.tru=pi_inter*A/hat_f_a_inter.tru;
		hat_yij_reg.mu_tru=lm_info.mu_tru[[2]];hat_y_obsA.ij.tru=lm_info.mu_tru[[5]]
			hat_var_y_reg_part1.mu_tru=lm_info.mu_tru[[8]]
	}else if (Y_Y1_Y0=='Y0'){
		A_pi_over_f.tru=pi_inter*(1-A)/hat_f_a_inter.tru;
		hat_yij_reg.mu_tru=lm_info.mu_tru[[3]];hat_y_obsA.ij.tru=lm_info.mu_tru[[6]]
			hat_var_y_reg_part1.mu_tru=lm_info.mu_tru[[9]]
	}
	##---------- Calculate estimators ---------##

	hat_yij_ipw.pi_tru=y*A_pi_over_f.tru
	hat_yij_DR.pi_tru_mu_tru=hat_yij_reg.mu_tru+hat_yij_ipw.pi_tru-
						hat_y_obsA.ij.tru*A_pi_over_f.tru
	first_last_child_each_house_id.m=cbind(first_child_each_house_id,last_child_each_house_id)
	hat_yi_ipw.pi_tru=apply(first_last_child_each_house_id.m,1,sum_part_fun,hat_yij_ipw.pi_tru)/n_i
	hat_yi_reg.mu_tru=apply(first_last_child_each_house_id.m,1,sum_part_fun,hat_yij_reg.mu_tru)/n_i
	hat_yi_DR.pi_tru_mu_tru=apply(first_last_child_each_house_id.m,1,sum_part_fun,hat_yij_DR.pi_tru_mu_tru)/n_i

	hat_y_ipw.pi_tru=mean(hat_yi_ipw.pi_tru)
	hat_y_reg.mu_tru=mean(hat_yi_reg.mu_tru)
	hat_y_DR.pi_tru_mu_tru=mean(hat_yi_DR.pi_tru_mu_tru)
	hat_psi_v_ipw.pi_tru=hat_yi_ipw.pi_tru-hat_y_ipw.pi_tru
	hat_var_y_ipw.pi_tru=cal_ipwvar_w_sandwich_fun(hat_psi_v_ipw=hat_psi_v_ipw.pi_tru,d_lv.m=d_lv.tru.m,
						hat_inv_V.m=hat_E_d_lv_d_lv_inv.tru.m,m=m,num_covariates=num_covariates)
	hat_var_y_reg.mu_tru=hat_var_y_reg_part1.mu_tru+var(hat_yi_reg.mu_tru)/m
	hat_var_y_DR.pi_tru_mu_tru=var(hat_yi_DR.pi_tru_mu_tru)/m

	if(dia==T){
		estimates=c(hat_y_ipw.pi_tru,hat_y_reg.mu_tru,hat_y_DR.pi_tru_mu_tru,		
			hat_var_y_ipw.pi_tru,hat_var_y_reg.mu_tru,hat_var_y_DR.pi_tru_mu_tru)
		DR.i.data=data.frame(pi_tru_mu_tru=hat_yi_DR.pi_tru_mu_tru)
		return(list(estimates=estimates,psi.tru=hat_psi_v_ipw.pi_tru,DR.i.data=DR.i.data,
			hat_yi_reg.mu_tru=hat_yi_reg.mu_tru))
	}else{
		#----------------------------------------------------------#
		#      estimators using one or more misspecified model     #
		#----------------------------------------------------------#
		hat_f_a_inter.mis=mean.temp.propen.pi.mis.m[1,]
		d_lv.mis.m=mean.temp.propen.pi.mis.m[-1,4*(0:(m-1))+1]
		hat_E_d_lv_d_lv.mis.m=d_lv.mis.m%*%t(d_lv.mis.m)/m
		hat_E_d_lv_d_lv_inv.mis.m=solve(hat_E_d_lv_d_lv.mis.m)
		if(Y_Y1_Y0=='Y'){
			A_pi_over_f.mis=pi_A*pi_inter/hat_f_a_inter.mis
			hat_yij_reg.mu_mis=lm_info.mu_mis[[1]];hat_y_obsA.ij.mis=lm_info.mu_mis[[4]];
				hat_var_y_reg_part1.mu_mis=lm_info.mu_mis[[7]]
		}else if (Y_Y1_Y0=='Y1'){
			A_pi_over_f.mis=pi_inter*A/hat_f_a_inter.mis			
			hat_yij_reg.mu_mis=lm_info.mu_mis[[2]];hat_y_obsA.ij.mis=lm_info.mu_mis[[5]]
				hat_var_y_reg_part1.mu_mis=lm_info.mu_mis[[8]]
		}else if (Y_Y1_Y0=='Y0'){
			A_pi_over_f.mis=pi_inter*(1-A)/hat_f_a_inter.mis
			hat_yij_reg.mu_mis=lm_info.mu_mis[[3]];hat_y_obsA.ij.mis=lm_info.mu_mis[[6]]
				hat_var_y_reg_part1.mu_mis=lm_info.mu_mis[[9]]
		}
		##---------- Calculate estimators ---------##
		hat_yij_ipw.pi_mis=y*A_pi_over_f.mis
		hat_yij_DR.pi_mis_mu_tru=hat_yij_reg.mu_tru+hat_yij_ipw.pi_mis-
							hat_y_obsA.ij.tru*A_pi_over_f.mis
		hat_yij_DR.pi_tru_mu_mis=hat_yij_reg.mu_mis+hat_yij_ipw.pi_tru-
							hat_y_obsA.ij.mis*A_pi_over_f.tru
		hat_yij_DR.pi_mis_mu_mis=hat_yij_reg.mu_mis+hat_yij_ipw.pi_mis-
							hat_y_obsA.ij.mis*A_pi_over_f.mis
		hat_yi_ipw.pi_mis=apply(first_last_child_each_house_id.m,1,sum_part_fun,hat_yij_ipw.pi_mis)/n_i
		hat_yi_reg.mu_mis=apply(first_last_child_each_house_id.m,1,sum_part_fun,hat_yij_reg.mu_mis)/n_i
		hat_yi_DR.pi_mis_mu_tru=
			apply(first_last_child_each_house_id.m,1,sum_part_fun,hat_yij_DR.pi_mis_mu_tru)/n_i
		hat_yi_DR.pi_tru_mu_mis=
			apply(first_last_child_each_house_id.m,1,sum_part_fun,hat_yij_DR.pi_tru_mu_mis)/n_i
		hat_yi_DR.pi_mis_mu_mis=
			apply(first_last_child_each_house_id.m,1,sum_part_fun,hat_yij_DR.pi_mis_mu_mis)/n_i

		hat_y_ipw.pi_mis=mean(hat_yi_ipw.pi_mis)
		hat_y_reg.mu_mis=mean(hat_yi_reg.mu_mis)
		hat_y_DR.pi_mis_mu_tru=mean(hat_yi_DR.pi_mis_mu_tru)
		hat_y_DR.pi_tru_mu_mis=mean(hat_yi_DR.pi_tru_mu_mis)
		hat_y_DR.pi_mis_mu_mis=mean(hat_yi_DR.pi_mis_mu_mis)
		hat_psi_v_ipw.pi_mis=hat_yi_ipw.pi_mis-hat_y_ipw.pi_mis

		hat_var_y_reg.mu_mis=hat_var_y_reg_part1.mu_mis+var(hat_yi_reg.mu_mis)/m
		hat_var_y_ipw.pi_mis=cal_ipwvar_w_sandwich_fun(hat_psi_v_ipw=hat_psi_v_ipw.pi_mis,d_lv.m=d_lv.mis.m,
					hat_inv_V.m=hat_E_d_lv_d_lv_inv.mis.m,m=m,num_covariates=num_covariates)
		hat_var_y_DR.pi_mis_mu_tru=var(hat_yi_DR.pi_mis_mu_tru)/m
		hat_var_y_DR.pi_tru_mu_mis=var(hat_yi_DR.pi_tru_mu_mis)/m
		hat_var_y_DR.pi_mis_mu_mis=var(hat_yi_DR.pi_mis_mu_mis)/m

		estimates=c(hat_y_ipw.pi_tru,hat_y_reg.mu_tru,	
			hat_y_DR.pi_tru_mu_mis,hat_y_DR.pi_mis_mu_tru,hat_y_DR.pi_tru_mu_tru,
			hat_y_ipw.pi_mis,hat_y_reg.mu_mis,hat_y_DR.pi_mis_mu_mis,
			hat_var_y_ipw.pi_tru,hat_var_y_reg.mu_tru,
			hat_var_y_DR.pi_tru_mu_mis,hat_var_y_DR.pi_mis_mu_tru,hat_var_y_DR.pi_tru_mu_tru,
			hat_var_y_ipw.pi_mis,hat_var_y_reg.mu_mis,hat_var_y_DR.pi_mis_mu_mis)
		DR.i.data=data.frame(pi_mis_mi_tru=hat_yi_DR.pi_mis_mu_tru,pi_tru_mu_mis=hat_yi_DR.pi_tru_mu_mis,
			pi_tru_mu_tru=hat_yi_DR.pi_tru_mu_tru,pi_mis_mu_mis=hat_yi_DR.pi_mis_mu_mis)
		return(list(estimates=estimates,psi.tru=hat_psi_v_ipw.pi_tru,
				psi.mis=hat_psi_v_ipw.pi_mis,DR.i.data=DR.i.data,
				hat_yi_reg.mu_tru=hat_yi_reg.mu_tru,hat_yi_reg.mu_mis=hat_yi_reg.mu_mis))
	}
}

#########################################################
# Function to calculate hat Y,Y1,Y0 for different alpha #
#########################################################
 
sub_cal_Y_Y1_Y0_est_fun=function(alpha,A,A_inter,y,n_i,m,n,
		first_child_each_house_id,last_child_each_house_id,lm_info.mu_tru,
		lm_info.mu_mis,propen_info.pi_tru,propen_info.pi_mis,mean.temp.propen.pi.tru.m,
		mean.temp.propen.pi.mis.m,Z_covariates.m,X_covariates.m,inter.info,num_covariates,dia){

	pi_inter=apply(as.matrix(1:n),1,cal_f_inter.fun,inter.info,A,rep(alpha,n))

	result_y=cal_Y_Y1_Y0_est_fun(alpha,Y_Y1_Y0='Y',A=A,A_inter=A_inter,y=y,n_i=n_i,m=m,n=n,
		first_child_each_house_id=first_child_each_house_id,last_child_each_house_id=last_child_each_house_id,
		lm_info.mu_tru=lm_info.mu_tru,lm_info.mu_mis=lm_info.mu_mis,
		propen_info.pi_tru=propen_info.pi_tru,propen_info.pi_mis=propen_info.pi_mis,
		mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
		mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,
		pi_inter=pi_inter,num_covariates=num_covariates,dia=dia)	
	result_y1=cal_Y_Y1_Y0_est_fun(alpha,Y_Y1_Y0='Y1',A=A,A_inter=A_inter,y=y,n_i=n_i,m=m,n=n,
		first_child_each_house_id=first_child_each_house_id,last_child_each_house_id=last_child_each_house_id,
		lm_info.mu_tru=lm_info.mu_tru,lm_info.mu_mis=lm_info.mu_mis,
		propen_info.pi_tru=propen_info.pi_tru,propen_info.pi_mis=propen_info.pi_mis,
		mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
		mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,
		pi_inter=pi_inter,num_covariates=num_covariates,dia=dia)	
	result_y0=cal_Y_Y1_Y0_est_fun(alpha,Y_Y1_Y0='Y0',A=A,A_inter=A_inter,y=y,n_i=n_i,m=m,n=n,
		first_child_each_house_id=first_child_each_house_id,last_child_each_house_id=last_child_each_house_id,
		lm_info.mu_tru=lm_info.mu_tru,lm_info.mu_mis=lm_info.mu_mis,
		propen_info.pi_tru=propen_info.pi_tru,propen_info.pi_mis=propen_info.pi_mis,
		mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
		mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,
		pi_inter=pi_inter,num_covariates=num_covariates,dia=dia)	

	if(dia==T){
		temp.list=list(y.info=c(result_y$estimates,result_y1$estimates,result_y0$estimates),
			psi.tru=cbind(result_y$psi.tru,result_y1$psi.tru,result_y0$psi.tru),
			DR.i.data=list(yi_DR=result_y$DR.i.data,y1i_DR=result_y1$DR.i.data,y0i_DR=result_y0$DR.i.data),
			hat_y_i_reg.tru=list(hat_yi_reg.tru=result_y$hat_yi_reg.mu_tru,hat_y1i_reg.tru=result_y1$hat_yi_reg.mu_tru,
				hat_y0i_reg.tru=result_y0$hat_yi_reg.mu_tru),hat_M.tru=lm_info.mu_tru$hat_M,
			Z_withalpha.tru.m.list=lm_info.mu_tru$Z_withalpha.m.list,J.vec=lm_info.mu_tru$J.vec)
	}else{
		temp.list=list(y.info=c(result_y$estimates,result_y1$estimates,result_y0$estimates),
			psi.tru=cbind(result_y$psi.tru,result_y1$psi.tru,result_y0$psi.tru),
			psi.mis=cbind(result_y$psi.mis,result_y1$psi.mis,result_y0$psi.mis),
			DR.i.data=list(yi_DR=result_y$DR.i.data,y1i_DR=result_y1$DR.i.data,y0i_DR=result_y0$DR.i.data),
			hat_y_i_reg.tru=list(hat_yi_reg.tru=result_y$hat_yi_reg.mu_tru,
				hat_y1i_reg.tru=result_y1$hat_yi_reg.mu_tru,hat_y0i_reg.tru=result_y0$hat_yi_reg.mu_tru),
			hat_y_i_reg.mis=list(hat_yi_reg.mis=result_y$hat_yi_reg.mu_mis,
				hat_y1i_reg.mis=result_y1$hat_yi_reg.mu_mis,hat_y0i_reg.mis=result_y0$hat_yi_reg.mu_mis),
			hat_M.tru=lm_info.mu_tru$hat_M,Z_withalpha.tru.m.list=lm_info.mu_tru$Z_withalpha.m.list,
			hat_M.mis=lm_info.mu_mis$hat_M,Z_withalpha.mis.m.list=lm_info.mu_mis$Z_withalpha.m.list,
			J.vec=lm_info.mu_tru$J.vec)

	}
	return(temp.list)
}

##################################################
# Function to calculate causal effect estimators #
#          and their variance estimators         #
##################################################

cal_CE_fun=function(CE,hat_y_y1_y0_alpha_1.info,hat_y_y1_y0_alpha_0.info,
		mean.temp.propen.pi.tru.m,mean.temp.propen.pi.mis.m,m,
		last_child_each_house_id,num_covariates,dia){

	J.vec=hat_y_y1_y0_alpha_1.info$J.vec  # J.vec doesn't depend on alpha
	if(dia==T){
		est_num=3
		d_lv.tru.m=mean.temp.propen.pi.tru.m[-1,last_child_each_house_id]
	}else{
		est_num=8
		Z_withalpha1.mis.m.list=hat_y_y1_y0_alpha_1.info$Z_withalpha.mis.m.list
		Z_withalpha0.mis.m.list=hat_y_y1_y0_alpha_0.info$Z_withalpha.mis.m.list
		hat_M.mis=hat_y_y1_y0_alpha_1.info$hat_M.mis   # hat_M just depends on observed covariates
									     # it doesn't depend on alpha
		if(CE=='DE'){
			hat_CE_i.mis=hat_y_y1_y0_alpha_1.info$hat_y_i_reg.mis$hat_y1i_reg.mis-
					 hat_y_y1_y0_alpha_1.info$hat_y_i_reg.mis$hat_y0i_reg.mis
			psi_CE.mis=	hat_y_y1_y0_alpha_1.info$psi.mis[,2]-hat_y_y1_y0_alpha_1.info$psi.mis[,3]
			CE_Z_withalpha.mis.m=Z_withalpha1.mis.m.list$Z_A1.m-Z_withalpha1.mis.m.list$Z_A0.m
		}else if (CE=='IE'){
			hat_CE_i.mis=hat_y_y1_y0_alpha_1.info$hat_y_i_reg.mis$hat_y0i_reg.mis-
					 hat_y_y1_y0_alpha_0.info$hat_y_i_reg.mis$hat_y0i_reg.mis
			psi_CE.mis=	hat_y_y1_y0_alpha_1.info$psi.mis[,3]-hat_y_y1_y0_alpha_0.info$psi.mis[,3]
			CE_Z_withalpha.mis.m=Z_withalpha1.mis.m.list$Z_A0.m-Z_withalpha0.mis.m.list$Z_A0.m
		}else if (CE=='TE'){
			hat_CE_i.mis=hat_y_y1_y0_alpha_1.info$hat_y_i_reg.mis$hat_y1i_reg.mis-
					 hat_y_y1_y0_alpha_0.info$hat_y_i_reg.mis$hat_y0i_reg.mis
			psi_CE.mis=	hat_y_y1_y0_alpha_1.info$psi.mis[,2]-hat_y_y1_y0_alpha_0.info$psi.mis[,3]
			CE_Z_withalpha.mis.m=Z_withalpha1.mis.m.list$Z_A1.m-Z_withalpha0.mis.m.list$Z_A0.m
		}else if (CE=='OE'){
			hat_CE_i.mis=hat_y_y1_y0_alpha_1.info$hat_y_i_reg.mis$hat_yi_reg.mis-
					 hat_y_y1_y0_alpha_0.info$hat_y_i_reg.mis$hat_yi_reg.mis
			psi_CE.mis=	hat_y_y1_y0_alpha_1.info$psi.mis[,1]-hat_y_y1_y0_alpha_0.info$psi.mis[,1]
			CE_Z_withalpha.mis.m=Z_withalpha1.mis.m.list$Z_A.m-Z_withalpha0.mis.m.list$Z_A.m
		}
		d_lv.tru.m=mean.temp.propen.pi.tru.m[-1,4*(0:(m-1))+1]
		d_lv.mis.m=mean.temp.propen.pi.mis.m[-1,4*(0:(m-1))+1]
		hat_E_d_lv_d_lv.mis.m=d_lv.mis.m%*%t(d_lv.mis.m)/m
		hat_E_d_lv_d_lv_inv.mis.m=solve(hat_E_d_lv_d_lv.mis.m)

		hat_var_CE_ipw.pi_mis=cal_ipwvar_w_sandwich_fun(hat_psi_v_ipw=psi_CE.mis,d_lv.m=d_lv.mis.m,
							hat_inv_V.m=hat_E_d_lv_d_lv_inv.mis.m,m=m,num_covariates=num_covariates)	
		hat_var_CE_reg_part1.mis=cal_hat_var_reg_part1_fun(
				Z_withalpha.m=CE_Z_withalpha.mis.m,J.vec=J.vec,hat_M=hat_M.mis,m=m)	
		hat_var_CE_reg.mu_mis=hat_var_CE_reg_part1.mis+var(hat_CE_i.mis)/m
	}
	hat_M.tru=hat_y_y1_y0_alpha_1.info$hat_M.tru
	Z_withalpha1.tru.m.list=hat_y_y1_y0_alpha_1.info$Z_withalpha.tru.m.list
	Z_withalpha0.tru.m.list=hat_y_y1_y0_alpha_0.info$Z_withalpha.tru.m.list
	if(CE=='DE'){
		hat_CE_i.tru=hat_y_y1_y0_alpha_1.info$hat_y_i_reg.tru$hat_y1i_reg.tru-
				 hat_y_y1_y0_alpha_1.info$hat_y_i_reg.tru$hat_y0i_reg.tru
		hat_CE=hat_y_y1_y0_alpha_1.info$y.info[(2*est_num+1):(3*est_num)]-
			 hat_y_y1_y0_alpha_1.info$y.info[(4*est_num+1):(5*est_num)]
		psi_CE.tru=	hat_y_y1_y0_alpha_1.info$psi.tru[,2]-hat_y_y1_y0_alpha_1.info$psi.tru[,3]
		hat_DR_CE_i=hat_y_y1_y0_alpha_1.info$DR.i.data$y1i_DR-hat_y_y1_y0_alpha_1.info$DR.i.data$y0i_DR
		CE_Z_withalpha.tru.m=Z_withalpha1.tru.m.list$Z_A1.m-Z_withalpha1.tru.m.list$Z_A0.m
	}else if (CE=='IE'){
		hat_CE_i.tru=hat_y_y1_y0_alpha_1.info$hat_y_i_reg.tru$hat_y0i_reg.tru-
				 hat_y_y1_y0_alpha_0.info$hat_y_i_reg.tru$hat_y0i_reg.tru
		hat_CE=hat_y_y1_y0_alpha_1.info$y.info[(4*est_num+1):(5*est_num)]-
			 hat_y_y1_y0_alpha_0.info$y.info[(4*est_num+1):(5*est_num)]
		psi_CE.tru=	hat_y_y1_y0_alpha_1.info$psi.tru[,3]-hat_y_y1_y0_alpha_0.info$psi.tru[,3]
		hat_DR_CE_i=hat_y_y1_y0_alpha_1.info$DR.i.data$y0i_DR-hat_y_y1_y0_alpha_0.info$DR.i.data$y0i_DR
		CE_Z_withalpha.tru.m=Z_withalpha1.tru.m.list$Z_A0.m-Z_withalpha0.tru.m.list$Z_A0.m
	}else if (CE=='TE'){
		hat_CE_i.tru=hat_y_y1_y0_alpha_1.info$hat_y_i_reg.tru$hat_y1i_reg.tru-
				 hat_y_y1_y0_alpha_0.info$hat_y_i_reg.tru$hat_y0i_reg.tru
		hat_CE=hat_y_y1_y0_alpha_1.info$y.info[(2*est_num+1):(3*est_num)]-
			 hat_y_y1_y0_alpha_0.info$y.info[(4*est_num+1):(5*est_num)]
		psi_CE.tru=	hat_y_y1_y0_alpha_1.info$psi.tru[,2]-hat_y_y1_y0_alpha_0.info$psi.tru[,3]
		hat_DR_CE_i=hat_y_y1_y0_alpha_1.info$DR.i.data$y1i_DR-hat_y_y1_y0_alpha_0.info$DR.i.data$y0i_DR
		CE_Z_withalpha.tru.m=Z_withalpha1.tru.m.list$Z_A1.m-Z_withalpha0.tru.m.list$Z_A0.m
	}else if (CE=='OE'){
		hat_CE_i.tru=hat_y_y1_y0_alpha_1.info$hat_y_i_reg.tru$hat_yi_reg.tru-
				 hat_y_y1_y0_alpha_0.info$hat_y_i_reg.tru$hat_yi_reg.tru
		hat_CE=hat_y_y1_y0_alpha_1.info$y.info[1:est_num]-hat_y_y1_y0_alpha_0.info$y.info[1:est_num]
		psi_CE.tru=	hat_y_y1_y0_alpha_1.info$psi.tru[,1]-hat_y_y1_y0_alpha_0.info$psi.tru[,1]
		hat_DR_CE_i=hat_y_y1_y0_alpha_1.info$DR.i.data$yi_DR-hat_y_y1_y0_alpha_0.info$DR.i.data$yi_DR
		CE_Z_withalpha.tru.m=Z_withalpha1.tru.m.list$Z_A.m-Z_withalpha0.tru.m.list$Z_A.m
	}
	hat_E_d_lv_d_lv.tru.m=d_lv.tru.m%*%t(d_lv.tru.m)/m
	hat_E_d_lv_d_lv_inv.tru.m=solve(hat_E_d_lv_d_lv.tru.m)


	hat_var_CE_ipw.pi_tru=cal_ipwvar_w_sandwich_fun(hat_psi_v_ipw=psi_CE.tru,d_lv.m=d_lv.tru.m,
						hat_inv_V.m=hat_E_d_lv_d_lv_inv.tru.m,m=m,num_covariates=num_covariates)
	hat_var_CE_reg_part1.tru=cal_hat_var_reg_part1_fun(
			Z_withalpha.m=CE_Z_withalpha.tru.m,J.vec=J.vec,hat_M=hat_M.tru,m=m)	
	hat_var_CE_reg.mu_tru=hat_var_CE_reg_part1.tru+var(hat_CE_i.tru)/m

	if(dia==T){
		hat_var_CE_DR=var(hat_DR_CE_i)/m
		temp.list=list(hat_CE=hat_CE,hat_var_CE_ipw.pi_tru=hat_var_CE_ipw.pi_tru,
				hat_var_CE_reg.mu_tru=hat_var_CE_reg.mu_tru,hat_var_CE_DR=hat_var_CE_DR)
	}else{
		hat_var_CE_DR=apply(hat_DR_CE_i,2,var)/m
		temp.list=list(hat_CE=hat_CE,hat_var_CE_ipw.pi_tru=hat_var_CE_ipw.pi_tru,hat_var_CE_ipw.pi_mis=hat_var_CE_ipw.pi_mis,
				hat_var_CE_reg.mu_tru=hat_var_CE_reg.mu_tru,hat_var_CE_reg.mu_mis=hat_var_CE_reg.mu_mis,
				hat_var_CE_DR=hat_var_CE_DR)
	}
	return(temp.list)
}
               
####################################
# Function to calculate estimators #
####################################

cal_est_fun=function(n_i,m,n,first_child_each_house_id,last_child_each_house_id,alpha_0,alpha_1,
		inter.info,inter_model,inter_size,mydata,num_covariates,dia,cal_CE){

	y=mydata$y;A=mydata$A;A_inter=mydata$A_inter;

	if(dia==T){
		var.name=c('agecat1_2','agecat2_3','newmom_edu','newfloor','season',
			'sanitation','watersourcecat','breastfed','sexcat')
		Z_covariates.m=matrix(unlist(mydata[,var.name]),ncol=ncol(mydata[,var.name]))
		colnames(Z_covariates.m)=var.name
		lm_info_alpha_0.mu_mis=NULL;lm_info_alpha_1.mu_mis=NULL
		propen_info.pi_mis=NULL;mean.temp.propen.pi.mis.m=NULL
	}else{
		var.name=c('Z1','Z2','Z3','Z4')
		Z_covariates.m=matrix(unlist(mydata[,var.name]),ncol=ncol(mydata[,var.name]))
		colnames(Z_covariates.m)=var.name
		var.name=c('X1','X2','X3','X4')
		X_covariates.m=matrix(unlist(mydata[,var.name]),ncol=ncol(mydata[,var.name]))
		colnames(X_covariates.m)=var.name

		##----- Outcome ----##
		##------ Wrong model -------##
		lm_info_alpha_0.mu_mis=fit_lm_fun(Iflm.tru=F,n_i=n_i,m=m,n=n,X_covariates.m,alpha=alpha_0,
				y=y,A=A,A_inter=A_inter,inter_size=inter_size,dia=dia,mydata=mydata,cal_CE=cal_CE)
		lm_info_alpha_1.mu_mis=fit_lm_fun(Iflm.tru=F,n_i=n_i,m=m,n=n,X_covariates.m,alpha=alpha_1,
				y=y,A=A,A_inter=A_inter,inter_size=inter_size,dia=dia,mydata=mydata,cal_CE=cal_CE)

		##------ Propensity ------##
		##------ Wrong model -------##
		propen_info.pi_mis=fit_propen_fun(Ifpropen.tru=F,X_covariates.m,A=A,A_inter=A_inter,
					inter_model=inter_model,inter_size=inter_size,n_i=n_i,m=m,n=n,dia=dia,mydata=mydata)
		hat_logit_p.mis=propen_info.pi_mis[[1]]
		est_hat_ranef_sd.mis=propen_info.pi_mis[[2]]
		mean.temp.propen.pi.mis.m=apply(as.matrix(1:n),1,to_integ_fun,inter.info,A,
			hat_logit_p.mis,sd_b=est_hat_ranef_sd.mis,covariates.m=cbind(1,X_covariates.m),
			num_covariates=num_covariates)
	}

	##----- Outcome ----##
	##------ Right model -------##
	lm_info_alpha_0.mu_tru=fit_lm_fun(Iflm.tru=T,n_i=n_i,m=m,n=n,Z_covariates.m,alpha=alpha_0,
				y=y,A=A,A_inter=A_inter,inter_size=inter_size,dia=dia,mydata=mydata,cal_CE=cal_CE)
	lm_info_alpha_1.mu_tru=fit_lm_fun(Iflm.tru=T,n_i=n_i,m=m,n=n,Z_covariates.m,alpha=alpha_1,
				y=y,A=A,A_inter=A_inter,inter_size=inter_size,dia=dia,mydata=mydata,cal_CE=cal_CE)

	##------ Propensity ------##
	##------ Right model -------##
	propen_info.pi_tru=fit_propen_fun(Ifpropen.tru=T,Z_covariates.m,A=A,A_inter=A_inter,
				inter_model=inter_model,inter_size=inter_size,n_i=n_i,m=m,n=n,dia=dia,mydata=mydata)
	hat_logit_p.tru=propen_info.pi_tru[[1]]
	est_hat_ranef_sd.tru=propen_info.pi_tru[[2]]
	mean.temp.propen.pi.tru.m=apply(as.matrix(1:n),1,to_integ_fun,inter.info,A,
		hat_logit_p.tru,sd_b=est_hat_ranef_sd.tru,covariates.m=cbind(1,Z_covariates.m),
		num_covariates=num_covariates)


	##---------- Calculate estimators ---------##
	#--- Ys ---#
	hat_y_y1_y0_alpha_0.info=sub_cal_Y_Y1_Y0_est_fun(alpha=alpha_0,A=A,A_inter=A_inter,y=y,n_i=n_i,m=m,n=n,
		first_child_each_house_id=first_child_each_house_id,last_child_each_house_id=last_child_each_house_id,
		lm_info.mu_tru=lm_info_alpha_0.mu_tru,lm_info.mu_mis=lm_info_alpha_0.mu_mis,propen_info.pi_tru=propen_info.pi_tru,
		propen_info.pi_mis=propen_info.pi_mis,mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
		mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,Z_covariates.m=Z_covariates.m,
		X_covariates.m=X_covariates.m,inter.info=inter.info,num_covariates=num_covariates,dia=dia)

	hat_y_y1_y0_alpha_1.info=sub_cal_Y_Y1_Y0_est_fun(alpha=alpha_1,A=A,A_inter=A_inter,y=y,n_i=n_i,m=m,n=n,
		first_child_each_house_id=first_child_each_house_id,last_child_each_house_id=last_child_each_house_id,
		lm_info.mu_tru=lm_info_alpha_1.mu_tru,lm_info.mu_mis=lm_info_alpha_1.mu_mis,propen_info.pi_tru=propen_info.pi_tru,
		propen_info.pi_mis=propen_info.pi_mis,mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
		mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,Z_covariates.m=Z_covariates.m,
		X_covariates.m=X_covariates.m,inter.info=inter.info,num_covariates=num_covariates,dia=dia)

	if(cal_CE==T){
		#--- CEs ---#
		DE_info=cal_CE_fun(CE='DE',hat_y_y1_y0_alpha_1.info=hat_y_y1_y0_alpha_1.info,
			hat_y_y1_y0_alpha_0.info=hat_y_y1_y0_alpha_0.info,
			mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
			mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,
			m=m,last_child_each_house_id=last_child_each_house_id,
			num_covariates=num_covariates,dia=dia)
		IE_info=cal_CE_fun(CE='IE',hat_y_y1_y0_alpha_1.info=hat_y_y1_y0_alpha_1.info,
			hat_y_y1_y0_alpha_0.info=hat_y_y1_y0_alpha_0.info,
			mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
			mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,
			m=m,last_child_each_house_id=last_child_each_house_id,
			num_covariates=num_covariates,dia=dia)
		TE_info=cal_CE_fun(CE='TE',hat_y_y1_y0_alpha_1.info=hat_y_y1_y0_alpha_1.info,
			hat_y_y1_y0_alpha_0.info=hat_y_y1_y0_alpha_0.info,
			mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
			mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,
			m=m,last_child_each_house_id=last_child_each_house_id,
			num_covariates=num_covariates,dia=dia)
		OE_info=cal_CE_fun(CE='OE',hat_y_y1_y0_alpha_1.info=hat_y_y1_y0_alpha_1.info,
			hat_y_y1_y0_alpha_0.info=hat_y_y1_y0_alpha_0.info,
			mean.temp.propen.pi.tru.m=mean.temp.propen.pi.tru.m,
			mean.temp.propen.pi.mis.m=mean.temp.propen.pi.mis.m,
			m=m,last_child_each_house_id=last_child_each_house_id,
			num_covariates=num_covariates,dia=dia)
		return(c(unlist(DE_info),unlist(IE_info),unlist(TE_info),unlist(OE_info)))
	}else{
		return(hat_y_y1_y0_alpha_1.info[[1]])
	}
}

#############################
# Function to simulate data #
#############################

sim_data_fun=function(dummy,n_i,m,alpha_0,alpha_1,coef_pi.true,
			coef_y.true,inter_model,inter_size,dia,cal_CE){
	if(dia==T){
#		source('Hajek_diarrhea_rcode_for_inclusion_local_pretoo.r')
		entry_data=entry_data_post
		y=(as.vector(table(id_post,new_dia_post)[,2]>0))+0
		A=as.numeric(entry_data$dose3)
		n=length(y)
		housecode_original=entry_data$housecode
#		m=length(unique(housecode_original))
		original_house_kid_num=table(entry_data$housecode)
		original_house_kid_num.dataframe=data.frame(housecode=names(original_house_kid_num),
			kids_num=as.numeric(original_house_kid_num))
		groupsize.v=original_house_kid_num.dataframe[,'kids_num']
		last_child_each_house_id=cumsum(groupsize.v)
		first_child_each_house_id=c(1,last_child_each_house_id+1)[1:m]
		entry_data$raw_id=as.numeric(as.factor((entry_data$id)))
		entry_data$raw_housecode=as.numeric(as.factor((entry_data$housecode)))
		temp_inter_info=apply(as.matrix(1:n),1,find_inter_fun,entry_data$raw_id,
				entry_data$raw_housecode)
		temp_inter=apply(as.matrix(1:length(temp_inter_info)),1,
			list_matrix_2_vector_byrow_fun,temp_inter_info)
		inter.info=data.frame(matrix(unlist(temp_inter),ncol=2,byrow=T))
		colnames(inter.info)=c('id','inter_id')
		A_inter_sum=apply(as.matrix(1:n),1,cal_sum_A_inter.fun,inter.info,A)
#		A_inter=A_inter_sum
		temp.J.list=apply(as.matrix(n_i),1,rep_n_i_fun)
		J.vec=matrix(unlist(temp.J.list),ncol=1)
#		A_inter=A_inter_sum*J.vec
		A_inter=as.numeric((A_inter_sum+1)*J.vec==1)

		mydata=data.frame(y,A,A_inter,
					entry_data[,c('raw_id','raw_housecode','agecat','agecat1_2',
					'agecat2_3','newmom_edu','newfloor','season','sanitation',
					'watersourcecat','breastfed','sexcat')])
		names(mydata)[names(mydata)=="raw_id"] <- "id"
		names(mydata)[names(mydata)=="raw_housecode"] <- "housecode"
		num_covariates=9
	}else{
		##------------------ covariates ----------------------##
		##-- Z is the real covariates (unobserved directly) --##
		##----------- X are the observed covariates ----------##
		n=n_i*m
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
		if(inter_model==1){
			A.m=matrix(A,nrow=n_i,byrow=F);A_friend.m=A.m[n_i:1,]
			A_friend=as.vector(A_friend.m) # group 1 first, then group 2 etc..
			A_inter=A_friend
		}else if (inter_model==2){
#			A_inter=as.numeric(A_inter_sum==(n_i-1))
#			A_inter=A_inter_sum
			A_inter=A_inter_sum/n_i
		}
		mu.true=coef_y.true[1]+coef_y.true[2]*Z1+coef_y.true[3]*Z2+
				coef_y.true[4]*Z3+coef_y.true[5]*Z4+
				coef_y.true[6]*A+coef_y.true[7]*A_inter+
				coef_y.true[8]*A*Z1+coef_y.true[9]*A_inter*Z2
		epsilon=rnorm(n,0,1)
		y=mu.true+epsilon
		mydata=data.frame(y,A,A_inter,Z1,Z2,Z3,Z4,X1,X2,X3,X4)
		num_covariates=4
		}
		result=cal_est_fun(n_i=n_i,m=m,n=n,first_child_each_house_id=first_child_each_house_id,
			last_child_each_house_id=last_child_each_house_id,alpha_0=alpha_0,alpha_1=alpha_1,
			inter.info=inter.info,inter_model=inter_model,inter_size=inter_size,
			mydata=mydata,num_covariates=num_covariates,dia=dia,cal_CE=cal_CE)
		return(result)
}

#################
# Main function #
#################

main_fun=function(n_i,m,alpha_0,alpha_1,coef_pi.true,coef_y.true,
			rnm_result,inter_model=1,inter_size,sim_number,dia,cal_CE){

	result=apply(as.matrix(1:sim_number),1,sim_data_fun,n_i,m,alpha_0=alpha_0,
			alpha_1=alpha_1,coef_pi.true=coef_pi.true,coef_y.true=coef_y.true,
			inter_model=inter_model,inter_size=inter_size,dia=dia,cal_CE=cal_CE)
	if(cal_CE==F){		rownames(result)=rnm_result}
	if(dia==T){
		return(result)
	}else{
		cat('\n','n_i:',n_i,'\n','m:',m,'\n')
		cat('alpha_1:',alpha_1,'\n')
	##------ Calculate the truth -------##
		bar_y_alpha_1=2+0.5*alpha_1+6*(alpha_1*inter_size/n_i)*inter
		bar_y1_alpha_1=2+0.5*1+6*(alpha_1*inter_size/n_i)*inter
		bar_y0_alpha_1=2+0.5*0+6*(alpha_1*inter_size/n_i)*inter
		mean_emp=apply(result,1,mean)
		bias_emp=abs(rep(c(bar_y_alpha_1,bar_y1_alpha_1,bar_y0_alpha_1),each=8)-mean_emp[c(1:8,17:24,33:40)])
		var_emp=apply(result,1,var)
		#print(mean_emp);cat('\n')
		print(bias_emp);cat('\n','Var:\n')
		print(var_emp[c(1:8,17:24,33:40)])
		var_hat_emp=cbind(mean_emp[c(1:8,17:24,33:40)+8],var_emp[c(1:8,17:24,33:40)])
		colnames(var_hat_emp)=c('est_var','emp_var')
		print(var_hat_emp)
		save(alpha_1,coef_y_inter.true,coef_pi.true,inter,inter_size,result,file=
		sprintf('DR_inter%d_ni%d_m%d_alpha_%fpi_%ssim%dintersize%d.RData',inter,n_i,m,alpha_1,pi_model,sim_number,inter_size))
		cat('##########################\n')
	#	return(bias_emp)
		return(result)
	}
}



#sim_number=10000
sim_number=2
#coef_pi.true=c(.5,-1,2,-1.25,-0.1);pi_model='skewed'
#coef_pi.true=c(1,-1,0.5,-0.25,-0.1);pi_model='center'
coef_pi.true=c(0,-1,0.5,-0.25,-0.1);pi_model='center'
coef_pi.true=c(0.5,-1,0.5,-0.25,-0.1);pi_model='center'
#n_alpha_0.m=cbind(2,500,c(1,3,5,7,9)/10)
#n_alpha_0.m=cbind(2,100,c(1,3,5,7,9)/10)
#n_alpha_0.m=cbind(c(50,100,200,500,1000,5000,1e4),.3)
#n_alpha_0.m=cbind(c(100),.3)


##----- With Interference -----##

#set.seed(12)
#inter=1
coef_y_inter.true=c(2,-1.5,-2.7,3,-1,0.5,6,1,8) # order is 1,Z1,Z2,Z3,Z4,A,A_friend,A*Z1,A_friend*Z2

#apply(as.matrix(n_alpha_0.m),1,main_fun,coef_pi.true=coef_pi.true,
#			coef_y.true=coef_y_inter.true,
#			rnm_result=rnm_result,inter_size=1)

set.seed(12)
inter=1
rnm_result=paste(rep(c('hat','hat_var'),each=8),rep(c('y','y1','y0'),each=16),c('ipw.pi_tru','reg.mu_tru',
		'DR.pi_tru_mu_mis','DR.pi_mis_mu_tru','DR.pi_tru_mu_tru',
		'ipw.pi_mis','reg.mu_mis','DR.pi_mis_mu_mis'),sep='_')
n_i=4;m=50;
alpha_0=0.1;alpha_1=0.7
#alpha_1=0.9

#n_alpha_0.m=cbind(n_i,m,alpha)

#debug(sim_data_fun)

##---------- Sims just on Y --------##
final_result=main_fun(n_i,m,alpha_0,alpha_1,coef_pi.true=coef_pi.true,
			coef_y.true=coef_y_inter.true,rnm_result=rnm_result,
			inter_model=2,inter_size=(n_i-1),sim_number,dia=F,cal_CE=F)

##---------- Sims on CE --------##
#final_result=main_fun(n_i,m,alpha_0,alpha_1,coef_pi.true=coef_pi.true,
#			coef_y.true=coef_y_inter.true,rnm_result=rnm_result,
#			inter_model=2,inter_size=(n_i-1),sim_number,dia=F,cal_CE=T)


##---- Diarrhea data ----##
#source('Hajek_diarrhea_rcode_for_inclusion_local_pretoo.r')
#n_i.v=as.numeric(house_kid_num_post);m=length(unique(entry_data_post$housecode))
#final_result=main_fun(n_i=n_i.v,m,alpha_0,alpha_1,coef_pi.true=coef_pi.true,
#			coef_y.true=coef_y_inter.true,rnm_result=rnm_result,
#			inter_model=2,inter_size=(n_i-1),sim_number=1,dia=T,cal_CE=T)

