library(RODBC)

RShowDoc("RODBC", package="RODBC")

ch <- odbcConnect('nautilus-vm.mathstat.dal.ca',uid = "srdbuser", pwd = "srd6us3r!")
