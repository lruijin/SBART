setwd("//jackdaw.biostat.wusm.wustl.edu/r.lu/simBART/results_sparse")
funcs = c("f1") #, "f2")
P = 5
sigmas = c("01", "1", "3", "5")
ns = c(100)#, 100)

sum_table = function(x, ndigits = 4){
    sum.table = 
        c(paste0(format(colMeans(x),digits = ndigits),"(", 
                 format(apply(x,2,sd),digits = ndigits),")"),
          paste0(format(apply(x,2,median),digits = ndigits), "(",
                 format(apply(x,2,quantile,0.025),digits = ndigits), ",",
                 format(apply(x,2,quantile,0.975),digits = ndigits), ")"))
    sum.table = as.data.frame(matrix(sum.table,nrow = 1))
    return(sum.table)
}

sum_table_binary = function(x, ndigits = 2){
    sum.table = 
        paste0(format(colMeans(x)*100, digits = ndigits), "%")
    sum.table = as.data.frame(matrix(sum.table,nrow = 1))
    return(sum.table)
}

infoTable <- table1 <- table2 <- table3 <- table4 <- table5 <- table6 <- NULL
infoTable_sum <- table1_sum <- table2_sum <- table3_sum <- table4_sum <- 
    table5_sum <- table6_sum <- NULL


for(func in funcs){
    for(sigma in sigmas){
        for(n in ns){
            for (s in 1:200){
                load(paste0(func,"_P",P,"_sigma",sigma,"_n",n,
                            "/result",s,".Rdata"))
                infoTable = rbind(infoTable, 
                                  c(Time,acc.rate))
                table1 = rbind(table1,
                               c(theta >= theta.CI.lower & 
                                     theta <= theta.CI.upper,
                                 eta >= eta.CI.lower & eta <= eta.CI.upper))
                table2 = rbind(table2,
                               c(theta.CI.upper-theta.CI.lower,
                                 eta.CI.upper - eta.CI.lower))
                table3 = rbind(table3,
                               c(eta.est, sse_eta))
                table4 = rbind(table4,
                               c(theta.est, sse_theta))
                table5 = rbind(table5,
                               c(rmse_simBart_test, mae_simBart_test,
                                 rmse_simBart_train, mae_simBart_train))
                table6 = rbind(table6, sel)
                if(s == 200){
                    sum_tmp = sum_table(infoTable)
                    names(sum_tmp) = c("Time.mean", "acc.rate.mean",
                                       "Time.median", "acc.rate.median")
                    infoTable_sum = 
                        rbind(infoTable_sum,
                              data.frame(
                                  fx = func,
                                  sigma = as.numeric(paste0("0.",sigma)),
                                  size = n,
                                  sum_tmp))
                    sum_tmp = sum_table_binary(table1)
                    names(sum_tmp) = c(paste0("theta", 1:(P-1)),
                                       paste0("eta", 1:P))
                    table1_sum = rbind(table1_sum,
                                       data.frame(
                                           fx = func,
                                           sigma = as.numeric(paste0("0.",sigma)),
                                           size = n,
                                           sum_tmp))
                    sum_tmp = sum_table(table2)
                    names(sum_tmp) = c(paste0("theta", 1:(P-1) ,".mean"),
                                       paste0("eta", 1:P, ".mean"),
                                       paste0("theta", 1:(P-1), ".median"),
                                       paste0("eta", 1:P, ".median"))
                    
                    table2_sum = rbind(table2_sum,
                                       data.frame(
                                           fx = func,
                                           sigma = as.numeric(paste0("0.",sigma)),
                                           size = n,
                                           sum_tmp))
                    
                    sum_tmp = sum_table(table3)
                    names(sum_tmp) = c(paste0("eta", 1:P, ".mean"),
                                       "SSE_eta.mean",
                                       paste0("eta", 1:P, ".median"),
                                       "SSE_eta.median")
                    
                    table3_sum = rbind(table3_sum,
                                       data.frame(
                                           fx = func,
                                           sigma = as.numeric(paste0("0.",sigma)),
                                           size = n,
                                           sum_tmp))
                    sum_tmp = sum_table(table4)
                    names(sum_tmp) = c(paste0("theta", 1:(P-1), ".mean"),
                                       "SSE_theta.mean",
                                       paste0("theta", 1:(P-1), ".median"),
                                       "SSE_theta.median")
                    
                    table4_sum = rbind(table4_sum,
                                       data.frame(
                                           fx = func,
                                           sigma = as.numeric(paste0("0.",sigma)),
                                           size = n,
                                           sum_tmp))
                    
                    sum_tmp = sum_table(table5)
                    names(sum_tmp) = c("RMSE_test.mean", "MAE_test.mean",
                                       "RMSE_train.mean", "MAE_train.mean",
                                       "RMSE_test.median", "MAE_test.median",
                                       "RMSE_train.median", "MAE_train.median")
                    
                    table5_sum = rbind(table5_sum,
                                       data.frame(
                                           fx = func,
                                           sigma = as.numeric(paste0("0.",sigma)),
                                           size = n,
                                           sum_tmp))
                   
                    sum_tmp = sum_table(table6)
                    names(sum_tmp) = c(paste0("theta", 1:(P-1), ".sel"),
                                       paste0("theta", 1:(P-1), ".sel"))
                    table6_sum = rbind(table6_sum,
                                       data.frame(
                                           fx = func,
                                           sigma = as.numeric(paste0("0.",sigma)),
                                           size = n,
                                           sum_tmp))
                    infoTable <- table1 <- table2 <- table3 <- table4 <- table5 <- table6 <-NULL
                }
            }
        }
    }
}


infoTable_sum
theta
table1_sum[,1:5]
table3_sum[,c(3,2,1,11)] %>% arrange(size, sigma, fx)
# theta
# table4_sum[,1:5] %>% inner_join(table1_sum[,1:5]) %>%
#     select(c(3,2,1,4,6,5,7)) %>%
#     arrange(size, sigma, fx)
table5_sum[,c(3,2,1,8,9)] %>% arrange(size, sigma, fx)
save(infoTable_sum, table1_sum, table2_sum, table3_sum,table5_sum,
     table6_sum,
     file = "f1_P3_n100_sparse_sum.Rdata")
