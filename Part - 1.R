# Question 1 - A) Downloading the date in the file "gene_expression.tsv"

download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
              destfile = "gene_expression.tsv")
# B) read into R-script, making gene accession number as the row number

a <- read.table("gene_expression.tsv")
head(a)
a <- read.table("gene_expression.tsv", header = TRUE, row.names = 1)
head(a)
str(a)


#Question 2 - Making mean the other column and showing the values of first six genes

a$mean <- rowMeans(a)
head(a)
str(a)


#Question 3 - List of 10 genes with highest mean expression

order(-a$mean)
ord <- a[order(-a$mean),]
head(ord,10)


#Question 4 - genes with mean less than 10

subset(a, mean<10)


#Question 5 - histogram of mean values 

hist(a$mean)
hist(a$mean, main = "Histogram of mean values", xlab = "Mean")


#Question 6 -  Downloading "growth_data.csv" file into R-script

download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv",
              destfile = "growth_data.csv")

 #read table into R-script and specify columns names of data
b <- read.csv("growth_data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE )
head(b)
str(b)
colnames(b)


#Question 7 - Calculating mean and standard deviation of tree circumference at the start and end at both Sites

# subset of northeast (NE) Site in trees growth
subset(b, Site == "northeast")
NE <- subset(b, Site == "northeast")
head(NE)
str(NE)

#Calculation of mean and standard deviation of tree circumference at start of NE
mean(NE$Circumf_2004_cm)
sd(NE$Circumf_2004_cm)

#Calculation of mean and standard deviaion of tree circumference at end of NE
mean(NE$Circumf_2019_cm)
sd(NE$Circumf_2019_cm)

#subset of southwest(SW) Site in trees growth
subset(b, Site =="southwest")
SW <- subset(b, Site == "southwest")
head(SW)
str(SW)

#Calculation of mean and standard deviation of tree circumference at start of SW
mean(SW$Circumf_2004_cm)
sd(SW$Circumf_2004_cm)

#Calculation of mean and standard deviation of tree circumference at end of SW
mean(SW$Circumf_2019_cm)
sd(SW$Circumf_2019_cm)


#Question 8 - Boxplot of tree circumference at start and end at both sites

#at start and end of northeast(NE) Site
boxplot(NE$Circumf_2004_cm, NE$Circumf_2019_cm, main = "Boxplot of Circumference of tree at northeast Site",
      ylab = "Circumference of tree", names = c("2004_NE", "2019_NE"))
#at start and end of southwest(SW) Site
boxplot(SW$Circumf_2004_cm, SW$Circumf_2019_cm, main = "Boxplot of Circumference of tree at southwest Site",
        ylab = "Circumference of tree", names = c("2004_SW", "2019_SW"))


#Question 9 - growth of mean at each site over past 10 years

#growth of mean at NE Site
NE$Circumf_2019_cm - NE$Circumf_2009_cm
x <- NE$Circumf_2019_cm - NE$Circumf_2009_cm
mean(x)
head(NE)
str(NE)

#growth of mean at SW Site
SW$Circumf_2019_cm - SW$Circumf_2009_cm
y <- SW$Circumf_2019_cm - SW$Circumf_2009_cm
mean(y)
head(SW)
str(SW)


#Question 10 - Performing t.test and wilcox.test to estimate p-value at both Sites

#t.test
?t.test
t.test(x , y)
#wilcox.test
?wilcox.test
wilcox.test(x, y)

