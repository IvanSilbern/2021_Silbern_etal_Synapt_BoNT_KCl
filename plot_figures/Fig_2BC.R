library(data.table)

CaEGTA <- fread("Figures\\Fig_2BC\\Fig_2B_source.txt")
BTX    <- fread("Figures\\Fig_2BC\\Fig_2C_source.txt")

pvalue.threshold <- 0.01
fc.threshold <- 1.2

pdf("Figures\\Fig_2BC\\Fig_2B.pdf")

plot(y = -log10(CaEGTA$q.val),
     x = CaEGTA$log2FC,
     # pch = 21,
     # bg = "lightgrey",
     #xlim = c(-2, 2),
     type = "n",
     #main = "-log10 Q-value vs log2 Ratio",
     xlab = "log2-Ratio (Ca / EGTA)",
     ylab = "-log10 q-value", cex = 2)
points(y = -log10(CaEGTA$q.val[!CaEGTA$Candidate]),
       x =  CaEGTA$log2FC[!CaEGTA$Candidate],
       pch = 21,
       bg  = "lightgrey")
points(y = -log10(CaEGTA$q.val[CaEGTA$Candidate]),
       x =  CaEGTA$log2FC[CaEGTA$Candidate],
       bg = "orange",
       pch = 21)

abline(v = log2(1/fc.threshold), col = "blue", lty = 3, lwd = 3)
abline(v = log2(fc.threshold),   col = "blue", lty = 3, lwd = 3)
abline(h = -log10(pvalue.threshold), col = "blue", lty = 3, lwd = 3)

dev.off()

pdf("Figures\\Fig_2BC\\Fig_2C.pdf")

plot(y = -log10(CaEGTA$q.val),
     x = CaEGTA$log2FC,
     # pch = 21,
     # bg = "lightgrey",
     #xlim = c(-2, 2),
     type = "n",
     #main = "-log10 Q-value vs log2 Ratio",
     xlab = "log2-Ratio (Mock / BTX)",
     ylab = "-log10 q-value")
points(y = -log10(BTX$q.val[!BTX$Candidate]),
       x =  BTX$log2FC[!BTX$Candidate],
       pch = 21,
       bg  = "lightgrey")
points(y = -log10(BTX$q.val[BTX$Candidate]),
       x =  BTX$log2FC[BTX$Candidate],
       bg = "orange",
       pch = 21)

abline(v = log2(1/fc.threshold), col = "blue", lty = 3, lwd = 3)
abline(v = log2(fc.threshold),   col = "blue", lty = 3, lwd = 3)
abline(h = -log10(pvalue.threshold), col = "blue", lty = 3, lwd = 3)

dev.off()