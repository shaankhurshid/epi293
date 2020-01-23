# Dependencies
library("devtools")
library(MRPRESSO)

#Load MR input data for LDL
mrdata=read.table("~/Documents/MPH/EPI293/ldl_cad_data.txt",sep=" ",row.names=1, header=T)

#Run MR-PRESSO
mr_ldl=mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, 
             DISTORTIONtest = TRUE, data = mrdata, NbDistribution = 5000,  SignifThreshold = 0.05)

mr_ldl$"Main MR results"
mr_ldl$"MR-PRESSO results"

#Load MR input data for HDL
mrdata=read.table("~/Documents/MPH/EPI293/hdl_cad_data.txt",sep=" ",row.names=1, header=T)

#Run MR-PRESSO
mr_hdl=mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = mrdata, NbDistribution = 5000,  SignifThreshold = 0.05)

mr_hdl$"Main MR results"
mr_hdl$"MR-PRESSO results"


