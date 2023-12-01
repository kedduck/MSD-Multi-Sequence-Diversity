
library(Biostrings)
library(corrplot)

########定义函数NVI
cal_NVI <- function(coll) {
  NVI_data <- apply(coll, 2, function(col){
    col_without_gap <- col[col != "-"]
    diff_aa_num <- length(unique(col_without_gap))
    max_aa_freq <- max(table(col_without_gap))
    length_col <- length(col_without_gap)
    NVI <- ifelse(length_col > 1, log(diff_aa_num / max_aa_freq) / log(length_col), 0)
    return(NVI)
  })
}


########定义带滑窗的函数NVI
cal_NVI_average <- function(data, length,win_size) {
  average_data <- sapply(1:length, function(i) {
    if (i <= win_size) {
      mean(data[1:(i + win_size)])
    } else if (length - i <= win_size) {
      mean(data[(i - win_size):length
                ])
    } else {
      mean(data[(i - win_size):(i + win_size)])
    }
  })
  return(average_data)
}


#######定义一致性比较
cal_iden <- function(seq){
  seq_num <- nrow(seq)
  seq_iden<- matrix(NA, seq_num, seq_num)
  
  for (i in 1:seq_num) {
    for (j in 1:i) {
      seq_data <- seq[c(i, j), , drop = FALSE]
      seq_no_gap <- seq_data[, !apply(seq_data == "-", 2, all)]
      seq_length <- ncol(seq_no_gap)
      iden_aa = 0
      for (k in 1:seq_length) {
        if (seq_no_gap[1,k] == seq_no_gap[2,k]){
          iden_aa <- iden_aa + 1
        }
      }
      seq_iden[j, i] <- iden_aa / seq_length
      seq_iden[i, j] <- iden_aa / seq_length
    }
  }
  colnames(seq_iden) <- rownames(seq)
  rownames(seq_iden) <- rownames(seq)
  return(seq_iden)
}


#######定义相似性比较
cal_sim <- function(seq) {
  aa_group1 <- c("F","Y","W")
  aa_group2 <- c("V","I","L")
  aa_group3 <- c("R","K","H")
  aa_group4 <- c("D","E")
  aa_group5 <- c("S","T")
  aa_group6 <- c("N","Q")
  group_list <- list(aa_group1,aa_group2,aa_group3,aa_group4,aa_group5,aa_group6)
  seq_num <- nrow(seq)
  seq_sim <- matrix(NA, seq_num, seq_num)
  
  for (i in 1:seq_num) {
    for (j in 1:i) {
      seq_data <- seq[c(i, j), , drop = FALSE]
      seq_no_gap <- seq_data[, !apply(seq_data == "-", 2, all)]
      seq_length <- ncol(seq_no_gap)
      sim_aa = 0
      for (k in 1:seq_length) {
        if (seq_no_gap[1,k] == seq_no_gap[2,k]){
          sim_aa <- sim_aa + 1
        } else{
          for(l in group_list){
            if (seq_no_gap[1,k] %in% l && seq_no_gap[2,k] %in% l ){
              sim_aa <-  sim_aa + 1
              break
            }
          }
        }
      }

      seq_sim[j, i] <- sim_aa / seq_length
      seq_sim[i, j] <- sim_aa / seq_length
    }
  }
  colnames(seq_sim) <- rownames(seq)
  rownames(seq_sim) <- rownames(seq)
  return(seq_sim)
}

####读入数据
align_file <- readAAStringSet("RNase.fas")
df <- as.matrix(align_file)

####开始运行NVI
NVI_data <- cal_NVI(df)

####求平均值方差等
NVI_data_mean <- mean(NVI_data)
NVI_data_var <- var(NVI_data)
NVI_data_sd <- sd(NVI_data)
NVI_static <- paste("mean", NVI_data_mean, 
                    "Var", NVI_data_var, 
                    "sd", NVI_data_sd,
                    sep="\t")

####求滑窗后的值
NVI_average_data <- cal_NVI_average(NVI_data, ncol(df), 2)

####保存结果文件
write.csv(NVI_data, file = "NVI.csv", row.names = FALSE)
write.table(NVI_static, file = "NVI_mean.csv", row.names = FALSE)
write.csv(NVI_average_data, file = "NVI_average.csv", row.names = FALSE)

#####NVI画图
pdf(file="NVI_average.pdf")
plot(NVI_average_data, las=1, xlab="", ylab="NVI",
     type="l", lwd=1.5, col="#007EA7", ylim=c(-1,1))
dev.off()

####开始运行iden_sim
iden_matrix <- cal_iden(df)

sim_matrix <- cal_sim(df)

#####iden或者sim 单独画图
# cor.plot <- corrplot(corr = iden_matrix, method="color",
#                      is.corr = FALSE, col.lim = c(0,1),
#                      tl.srt = -30,
#                      type="lower",tl.pos = 'lt',outline=FALSE, 
#                      number.cex = 0.8, number.font = 1,
#                      addgrid.col = "white",cl.pos='r',addCoef.col = "white")

#####iden或者sim 合并画图
mix_matrix <- iden_matrix 
mix_matrix[upper.tri(mix_matrix)] <- sim_matrix[upper.tri(sim_matrix)]

####保存结果文件
write.csv(mix_matrix, file = "iden_sim.csv")
#####画图
pdf(file="siden_sim.pdf")
cor.plot <- corrplot(corr = mix_matrix, method="color",
                     is.corr = FALSE, col.lim = c(0,1),
                     tl.srt = 30, tl.pos = 'lt', tl.col = "black",
                     type="full",col = COL1('YlOrRd'),
                     number.cex = 0.8, number.font = 1,
                     addgrid.col = "white",cl.pos='r',addCoef.col = "black")
dev.off()

