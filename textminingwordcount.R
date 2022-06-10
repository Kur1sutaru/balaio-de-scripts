# Text mining and word cloud fundamentals in R
# Retrieved from http://www.sthda.com/english/wiki/text-mining-and-word-cloud-fundamentals-in-r-5-simple-steps-you-should-know

setwd("C:/Users/crist/Downloads/Pam revisao preditores")

# Paste the text in a plain text file (e.g : ml.txt. Save the file
combined_txt_towordcloud <- read.delim("C:/Users/crist/Downloads/Pam revisao preditores/test.txt")
View(combined_txt_towordcloud)

#Install and load the required packages
#Type the R code below, to install and load the required packages

# Install
# for text mining
install.packages("tm")
# for text stemming
install.packages("SnowballC")
# word-cloud generator
install.packages("wordcloud")
# color palettes
install.packages("RColorBrewer")

# Load
library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")

# load the text
# The text is loaded using Corpus() function from text mining (tm) package. 
# Corpus is a list of a document (in our case, we only have one document).
docs <- Corpus(VectorSource(combined_txt_towordcloud))

# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)
# Text stemming
docs <- tm_map(docs, stemDocument)

dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 20)

# Step 5 : Generate the Word cloud

set.seed(666)
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))










