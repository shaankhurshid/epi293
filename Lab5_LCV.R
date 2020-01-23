
#### Load the R functions for LCV

source("RunLCV.R")

#### Read and analyze HDL - CAD
data_hdl = read.table("~/Documents/MPH/EPI293/hdl_cad_lcv_data.txt",header=T)

N1=data_hdl$N.x
N2=data_hdl$N.y
LCV_hdl = RunLCV(data_hdl$baseL2,data_hdl$Z.x,data_hdl$Z.y,ldsc.intercept=1,n.1=N1,n.2=N2)

print(LCV_hdl)

#### Read and analyze LDL - CAD

data_ldl = read.table("~/Documents/MPH/EPI293/ldl_cad_lcv_data.txt",header=T)

N1=data_ldl$N.x
N2=data_ldl$N.y
LCV_ldl = RunLCV(data_ldl$baseL2,data_ldl$Z.x,data_ldl$Z.y,ldsc.intercept=1,n.1=N1,n.2=N2)

print(LCV_ldl)


