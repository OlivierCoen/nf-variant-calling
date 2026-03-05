source("tools.R")

# Import tubular files produced by BCFTools We have many because we called freebayes on separate
# chromosome regions to spread the work.
files = list.files("tab/", full.names = T)
snps = lapply(files, fread, sep = "\t", header = T, na.strings = ".")
snps = rbindlist(snps)
n = colnames(snps)
setnames(snps, splitToColumns(n,"]", 2)) # trims the left part of BCFTools column names

# Import a tabular file indicating chromosome lengths (= the fasta index)
chroms = fread("/path/to/fa.fai", header = F, select = 1:2, col.names = c("chrom", "len"))

# Extract columns for read counts
RO = snps[,grep(":RO", n), with = F] # read count for the reference allele
AO = snps[,grep(":AO", n), with = F] # read count for the alt allele
snps = snps[,1:5, with = F] # so we can only retain other columns in the original table
setnames(AO, splitToColumns(colnames(AO), ":", 1))
setnames(RO, splitToColumns(colnames(RO), ":", 1))
AO = AO[,colnames(RO), with = F] # makes sure AO columns are in the same order as RO's
AO[is.na(AO)] = 0L
RO[is.na(RO)] = 0L

# total read depth for both alleles, per pool
dp = RO + AO

# we impose a minimum depth of 30X for each sample
dp30 = rowSums(dp < 30L) == 0L

# total depth must not be too high (< 90th percentile)
quantiles = sapply(dp, quantile, probs = 0.9)
tooHigh = Map(function(d, q) d > q, dp, quantiles )
dpOK = Reduce(`+`, tooHigh) == 0L

# we retain SNPs for which both alleles are not too rare
# = between 10% and 90% in at least one pool
freq = AO / dp
notRare = rowSums(freq >= 0.1 & freq <= 0.9) > 0L
notRare[is.na(notRare)] = F

goodQual = snps$QUAL > 10 # filter on SNP quality (based on visual inspection of its distribution)

# we evaluate the effect of filters on SNP density, to check that these filters
# do not penalize regions in particular. We compute the number of SNPs per 1 Mb window
snps[,genomePOS := scaffToGenomeCoord(CHROM, POS, chroms$len, chroms$chrom)]
snps[,window := genomePOS %/% 1e6 * 1e6]

nSNPs = snps[,.N, by = window]

nSNPs2 = snps[dp30 & dpOK & notRare,.N, by = window]
setorder(nSNPs, window)
nSNPs = merge(nSNPs, nSNPs2, by = "window")
with(nSNPs, plot(N.x, N.y))

retained = dp30 & dpOK & notRare
snps = snps[retained, 1:5, with = F]
RO = RO[retained,]
AO = AO[retained,]
dp = dp[retained,]

rm(dpOK, dp30, tooHigh, notRare, freq, goodQual, retained)


mantelro2<-function (ro, ao) {
  # returns the p-value of the Cochran-Mantel-Haenszel test of H0: read counts are independent of the trait (here, sex)
  # ro and ao = 4-column matrices: pop1_trait1, pop1_trait2, pop2_trait1, pop2_trait2
  # this function only works for two populations. Must be modified to account for more.
  s.x = ao + ro
  s.y = cbind(ro[,1]+ro[,2], ao[,1]+ao[,2], ro[,3]+ro[,4], ao[,3]+ao[,4])
  n = cbind(s.y[,1] + s.y[,2],  s.y[,3] +  s.y[,4])

  DELTA <- rowSums(ro[,c(1, 3)] - s.x[,c(1, 3) ] * s.y[, c(1, 3)]/n)
  YATES <- ifelse(abs(DELTA) >= 0.5, 0.5, 0)

  STATISTIC <- ((abs(DELTA) - YATES)^2/rowSums(cbind(s.x[,1]*s.x[,2]*s.y[,1]*s.y[,2], s.x[,3]*s.x[,4]*s.y[,3]*s.y[,4])/(n^2 * (n - 1))))

  pchisq(STATISTIC, 1, lower.tail = FALSE)
}

# we compute the p-values, first for uninfected populations
pVal = mantelro2(RO[,cbind(BVR_M, BVR_F, BUX_M, BUX_F)], AO[,cbind(BVR_M, BVR_F, BUX_M, BUX_F)])

# then for populations infected by the virus (non differentiation expected between sexes)
pValSTLC = mantelro2(RO[,cbind(STM_M, STM_A+STM_T, LC_M, LC_A+LC_T)], AO[,cbind(STM_M, STM_A+STM_T, LC_M, LC_A+LC_T)])

snps[,p := p.adjust(pVal, method = "fdr")]
snps[,pSTLC := p.adjust(pValSTLC, method = "fdr")]

# we prepare sliding-windows manhattan plots
melted = melt(snps, id.vars = c("CHROM","POS"), measure.vars = c("p","pSTLC"), variable.name = "type", value.name = "p")

melted[,window := as.integer(POS %/% 2e4 * 2e4 + 1e4)] # mid-position of contigous 20-kb windows

# we compute the 5th percentile of the p-value, for each window
perWindow = melted[,.(p = quantile(p, 0.05, na.rm = T)), by = .(type, window, CHROM)]

perWindow[,genomePOS := scaffToGenomeCoord(CHROM, window, chroms$len, chroms$chrom)]
perWindow[,chrom := stri_extract_first(CHROM, regex = "[0-9]+")]
perWindow[!grepl("SUPER",CHROM), chrom := "rest"]
setorder(perWindow, genomePOS)

# we compute the chromosome mid position in genome coordinates,
# to show on the x axis.
midChr = perWindow[,mean(range(genomePOS)), by = chrom]

library(ggplot2)

ggplot() +
  geom_point(perWindow, mapping = aes(y = -log10(p), x = genomePOS/1e6, col = chrom), size = 0.3) +
  xlab("Position (Mb)") +
  scale_color_manual(values = rep(c("grey","black"), length.out = length(unique(perWindow$chrom)))) +
  scale_x_continuous(label = midChr$chrom, breaks= midChr$V1/1e6) +
  scale_y_continuous(expand = c(0, 0.005)) +     # remove space between plot area and x axis

  theme_bw() +
  theme(
    axis.text.x=element_text(size=6),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


snps[, Wref := AO$BVR_M > RO$BVR_M]
snps[, maleFreq := RO[,BVR_M + BUX_M] / dp[,BVR_M + BUX_M]]

f = snps[,CHROM == "SUPER_19" & pZW > 1e-3 & POS > 32e6]
with(snps[f], plot(POS, maleFreq, type = "b"))
sel = data.table(snps[f, .(CHROM, POS, REF, ALT, pZW)], maleFreq = maleF[f])
writeT(snps[f, .(CHROM, POS, REF, ALT, maleFreq)], "snps.txt")

# plotting allele frequencies in males and females for SNPs with significant differentiation between sexes ------------

invert = sample(c(T, F), nrow(AO), T) # randomly picks the reference or alt allele, to avoid mapping biases (or reference biases, as the reference appears to be mostly composed of the Z allele)

poolList = list(`Unbiased Females` = c("BVR_F", "BUX_F"),
                `Unbiased Males` = c("BVR_M", "BUX_M"),
                `Biased Females` = grep("_A|_T", names(AO), value = T),
                `Biased Males` = c("LC_M","STM_M"))

types = names(poolList)

for(type in types) {
  pools = poolList[[type]]
  counts = rowSums(RO[, pools, with = F])
  counts[invert] = rowSums(AO[invert, pools, with = F])
  snps[,(type) := counts / rowSums(dp[, pools, with = F])]
}

snps[,logp := log10(p)]
chrom = "SUPER_19" # the super scaffold of sex chromsomes
quant = snps[CHROM == chrom, quantile(-logp, 0.99, na.rm = T)]

melted = melt(snps[CHROM == chrom & -logp >= quant], id.vars = "logp", measure.vars = types, variable.name = "type", value.name = "Allele frequency")
melted[,c("type", "sex") := data.table(splitToColumns(type, " "))]
melted[,type := factor(type, levels = c("Unbiased", "Biased"))] # to order the types properly in the plot


library(scales)
library(gridExtra)
ggplot(melted) +
  scale_y_continuous(labels = label_number(drop0trailing = TRUE)) +
  geom_bin_2d(aes(x = -logp, y = `Allele frequency`), bins = 50) +
  scale_fill_gradientn(colors = c("grey95", "dodgerblue4")) +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  facet_grid(type ~ sex)
