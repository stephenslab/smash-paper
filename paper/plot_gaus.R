library(dplyr)
library(reshape2)
library(ggplot2)


##barplots in main body of paper
mean.table = melt(table.1.v1)
names(mean.table) = c("method", "testfunction", "mse")

pdf("paper/smash_gaus_1_v1.pdf",width = 10, height = 5)
ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
                             & method != "SMASH, truth, Symm8"
                             & method != "TI-thresh, truth, Symm8"
                             & method != "TI-thresh, SMASH, Symm8"
                             & method != "SMASH, joint, Haar"), aes(x = testfunction, y = mse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Test Function") + 
  ylab("Relative MSE") + 
  ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 1") +
  geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
dev.off()


mean.table = melt(table.3.v1)
names(mean.table) = c("method", "testfunction", "mse")

pdf("paper/smash_gaus_3_v1.pdf",width = 10, height = 5)
ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
                             & method != "SMASH, truth, Symm8"
                             & method != "TI-thresh, truth, Symm8"
                             & method != "TI-thresh, SMASH, Symm8"
                             & method != "SMASH, joint, Haar"), aes(x = testfunction, y = mse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Test Function") + 
  ylab("Relative MSE") + 
  ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 3") +
  geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
dev.off()



# qplot(x = testfunction, y = mse, fill = method, xlab = "Test Function", ylab = "Relative MSE", main = "MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 1",
#                        data = mean.table %>% filter(method != "SURE, homo, Symm8"
#                                                     & method != "SMASH, truth, Symm8"
#                                                     & method != "TI-thresh, truth, Symm8"
#                                                     & method != "TI-thresh, SMASH, Symm8"
#                                                     & method != "SMASH, joint, Haar"
#                        )
#                        , geom = "bar", stat = "identity",
#                        position = "dodge")



mean.table = melt(table.1.v5)
names(mean.table) = c("method", "testfunction", "mse")


pdf("paper/smash_gaus_1_v5.pdf",width = 10, height = 5)
ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
                             & method != "SMASH, truth, Symm8"
                             & method != "TI-thresh, truth, Symm8"
                             & method != "TI-thresh, SMASH, Symm8"
                             & method != "SMASH, joint, Haar"
                             & method != "Ebayes, homo, Symm8"
                             & method != "PostMean, homo, Symm8"), aes(x = testfunction, y = mse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Test Function") + 
  ylab("Relative MSE") + 
  ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for Clipped Blocks variance function, SNR = 1") +
  theme(plot.title = element_text(size = rel(0.9))) +
  geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
dev.off()



mean.table = melt(table.3.v5)
names(mean.table) = c("method", "testfunction", "mse")


pdf("paper/smash_gaus_3_v5.pdf",width = 10, height = 5)
ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
                             & method != "SMASH, truth, Symm8"
                             & method != "TI-thresh, truth, Symm8"
                             & method != "TI-thresh, SMASH, Symm8"
                             & method != "SMASH, joint, Haar"
                             & method != "Ebayes, homo, Symm8"
                             & method != "PostMean, homo, Symm8"), aes(x = testfunction, y = mse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Test Function") + 
  ylab("Relative MSE") + 
  ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for Clipped Blocks variance function, SNR = 3") +
  theme(plot.title = element_text(size = rel(0.9))) +
  geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
dev.off()



##error bars?



##boxplot?


##lineplot
ggplot(mean.table, aes(x = testfunction, y = mse, color = method, group = method)) + 
       geom_line(linetype = 2) + 
       geom_point(shape = 19) + 
       xlab("Test Function") + 
       ylab("Relative MSE")


