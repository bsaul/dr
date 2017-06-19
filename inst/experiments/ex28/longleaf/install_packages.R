## Install packages necssary for running simulations on cluster ####

print(getwd())
packages <- c('geepack_1.2-1', 'rootSolve_1.7', 'geex_0.1.0', 'dr_0.1.0')

basedir <- '/nas/longleaf/home/saulb/'

installed <- installed.packages()
for(i in seq_along(packages)){
  pkg_name <- gsub('(_|\\.|-|[0-9])', '', packages[i])
  pkg_vers <- gsub('([a-z]|_)', '', packages[i] )
  
  print(paste('Attempting to install and load', pkg_name, pkg_vers))
  
  # check if correct version is installed
  if((pkg_name %in% rownames(installed) == FALSE) 
     # | (installed[rownames(installed) == pkg_name, 'Version'] != pkg_vers)
  ){
    
    install.packages(paste0(basedir ,"Rpackages/", packages[i], ".tar.gz"),
                     lib=paste0(basedir, "Rlibs/"), repos=NULL)
    lib_try <- library(pkg_name, lib.loc=paste0(basedir, "Rlibs/"), character.only = TRUE,
                       quietly = TRUE, logical.return = TRUE, warn.conflicts = FALSE,
                       verbose = FALSE)
    if(lib_try){
      print(paste(pkg_name, "installed"))
    } else {
      print(paste(pkg_name, "NOT installed"))
    }
  }
}

