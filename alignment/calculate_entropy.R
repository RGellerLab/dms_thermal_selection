shannon.entropy <- function(p)
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log(p.norm)*p.norm)
}

library(Biostrings)

# get entropy
# read alignment
res=readAAMultipleAlignment("suplementary data 1.fa")
res=as.matrix(res)
ent=apply(res,2,function(x) shannon.entropy(table(x)))

