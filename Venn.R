# Venn and intersect plotting Dominikas proteomic data

# load required packages
library(UpSetR)
library(limma)
library(venneuler)

# load raw data
data <- read.csv("Dominika_proteomics.csv")

# extract biorep data
colnames(data)
bioreps <- data[,c(1,213:236)]
colnames(bioreps)

# simplify by combining bioreps of all anther lengths
bioreps$B1 <- rowSums(bioreps[,c(2:7)])
bioreps$B2 <- rowSums(bioreps[,c(8:13)])
bioreps$B3 <- rowSums(bioreps[,c(14:19)])
bioreps$B5 <- rowSums(bioreps[,c(20:25)])
colnames(bioreps)

# simplify by combining anther length reps
bioreps$six <- rowSums(bioreps[,c(2,8,14,20)])
bioreps$seven <- rowSums(bioreps[,c(3,9,15,21)])
bioreps$eight <- rowSums(bioreps[,c(4,10,16,22)])
bioreps$nine <- rowSums(bioreps[,c(5,11,17,23)])
bioreps$ten <- rowSums(bioreps[,c(6,12,18,24)])
bioreps$eleven <- rowSums(bioreps[,c(7,13,19,25)])
colnames(bioreps)

simple_reps <- bioreps[,c(1,26:31)]
colnames(simple_reps)

# convert values >0 to 1 to make a binary count matrix
for (x in 2:ncol(simple_reps)) {
  simple_reps[,x] <- ifelse(simple_reps[,x]>0,1,0)
}

# plot intersects in intersect plot
setnames <- colnames(simple_reps[,2:7])
setnames
upset(simple_reps,
      sets=setnames,
      keep.order = TRUE,
      sets.x.label = c("0.6mm", "0.7mm", "0.8mm", "0.9mm", "1.0mm", "1.1mm"),
      order.by = "freq",
      nintersects = 20)

simple_reps$total <- rowSums(simple_reps[,2:7])

onlypointsix <- subset(simple_reps, six==1 &
                         total==1)

# Plot a very basic venn diagram with limma
vennDiagram(simple_reps[,2:5], names=setnames, cex=1, circle.col = c(1:5))

# Plot with venneuler
venneuler(simple_reps[,2:5])

