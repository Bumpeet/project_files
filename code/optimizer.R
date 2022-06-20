# Loading libraries
library(nloptr) # For COBYLA, Nelder-Mead and Subplex
library(data.table)

# Reading the original data

original_data <- read.csv("original_data.csv")
original_data = original_data$mol_frac
original_data = c(original_data[1:8], original_data[10:17])

objective <- function(k)
{
  
  system(paste0("isom.PFR.exe -override=E[1]=",k[1],
                ",E[2]=",k[2],",E[3]=",k[3],",E[4]=",k[4],",E[5]=",k[5],
                ",E[6]=",k[6],",E[7]=",k[7],",E[8]=",k[8],",E[9]=",k[9],
                ",E[10]=",k[10],",E[11]=",k[11],",E[12]=",k[12],",E[13]=",k[13],
                ",E[14]=",k[14],",E[15]=",k[15],",E[16]=",k[16],",E[17]=",k[17],
                ",E[18]=",k[18],",E[19]=",k[19],",E[20]=",k[20],",E[21]=",k[21],
                ",E[22]=",k[22],",E[23]=",k[23],",E[24]=",k[24],",E[25]=",k[25],
                ",E[26]=",k[26],",E[27]=",k[27],",E[28]=",k[28],",E[29]=",k[29],
                ",E[30]=",k[30],",E[31]=",k[31],",E[32]=",k[32],",E[33]=",k[33],
                ",E[34]=",k[34],",E[35]=",k[35],",E[36]=",k[36],",E[37]=",k[37],
                ",E[38]=",k[38],",E[39]=",k[39],",E[40]=",k[40],",E[41]=",k[41],
                ",E[42]=",k[42],",E[43]=",k[43],",E[44]=",k[44],",E[45]=",k[45],
                ",E[46]=",k[46],",E[47]=",k[47],",E[48]=",k[48],",E[49]=",k[49],
                ",E[50]=",k[50],",E[51]=",k[51],",E[52]=",k[52],",E[53]=",k[53],
                ",E[54]=",k[54]))
  
  # Reading data
  dt <- fread("ISOM.PFR_res.csv", select = c("y_i[1]","y_i[2]","y_i[3]",
                                                           "y_i[4]","y_i[5]","y_i[6]",
                                                           "y_i[7]","y_i[8]","y_i[10]","y_i[11]"
                                                           ,"y_i[12]","y_i[13]","y_i[14]"
                                                           ,"y_i[15]","y_i[16]","y_i[17]"))[6]
  simulated_data = unlist(dt[1])
  
  error = sqrt(sum(((simulated_data-original_data)/original_data)^2))*(1/17)
  
  cat("Error value is: ",error,"\n")
  return(error)
}


initial_par = c(148.93, 154.28, 143.17, 151.41,
                150.98, 155.92, 152.96, 149.95,
                127.28, 139.07, 64.5, 77.06,
                146.14, 160.28, 98.28,  105.4,
                3.51,  4.79, 180.2,  400.43,
                187.05, 300.79,  51.08, 341.89,
                129.75, 88.64, 135.45, 129.29,
                154.54, 98.63, 150.29, 102.35,
                168.13, 90.7, 177.32, 222.9,
                59.8, 54.08, 330.28, 329.06,
                284.97, 166.32, 165.82, 112.05,
                265, 264, 263.5, 265,
                295.62, 295, 295.19, 294,
                294.22, 278.81)
# initial_par = opt$par
lower_bound = initial_par-0.15*initial_par
upper_bound = initial_par+0.15*initial_par

opt=cobyla(x0 = initial_par,
           fn = objective,
           lower = lower_bound,
           upper = upper_bound,
           nl.info=TRUE,
           control = list(xtol_rel= 1e-9, maxeval=2000))

