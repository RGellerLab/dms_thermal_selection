library(Biostrings)

codons=names(GENETIC_CODE)
cod.df=data.frame(matrix(data=NA,ncol=64,nrow=64))
rownames(cod.df)=codons
colnames(cod.df)=codons
rm(codons)
# Generate df of all codons by all codons
# indicate the number of differences
for (i in 1:64){
  # i=1
  ref=rownames(cod.df)[i]
  for (c in 1:64){
    mut=colnames(cod.df)[c]
    r.let=unlist(strsplit(ref,""))
    m.let=unlist(strsplit(mut,""))
    cod.df[ref,mut]=sum(sapply(1:3, function(x) (r.let[x]!=m.let[x])+0))
    rm(m.let,c,r.let)
  }
  rm(i,ref,mut)
}


########################
## remove stops since not relevant for amino acid mutations
stps=c("TAA","TGA","TAG")
cod.df=cod.df[!row.names(cod.df)%in%stps, !colnames(cod.df)%in%stps ]
rm(stps)


### table of all res as mutations with AA translation
sum.res=data.frame(codon=row.names(cod.df),
                   AA=as.character(translate(DNAStringSet(row.names(cod.df)),no.init.codon = T)),
                   stringsAsFactors = F)
# codons reached by 1,2,3
sum.res$codons_by_single=NA
sum.res$codons_by_double=NA
sum.res$codons_by_triple=NA
## AA reached by 1,2,3
sum.res$AA_reached_by_single_mutation=NA
sum.res$AA_reached_by_double_mutation=NA
sum.res$AA_reached_by_triple_mutation=NA
## sums of aa reached
sum.res$number_of_AA_reached_by_single=NA
sum.res$number_of_AA_reached_by_single_double=NA
sum.res$number_of_AA_reached_by_triple_only=NA


## loop over
for (r in sum.res$codon){
  # r=sum.res$codon[1]
  
  # translate reference codon
  aa=as.character(translate(DNAString(r),no.init.codon = T))
  # Get codons that are 1, 2, or 3 amino acids away from reference codond
  cods1=colnames(cod.df)[cod.df[r,]==1]
  cods2=colnames(cod.df)[cod.df[r,]==2]
  cods3=colnames(cod.df)[cod.df[r,]==3]
  
  #Translate the codons that are 1,2,3 away into AA
  aa1=as.character(translate(DNAStringSet(cods1),no.init.codon = T))
  aa2=as.character(translate(DNAStringSet(cods2),no.init.codon = T))
  aa3=as.character(translate(DNAStringSet(cods3),no.init.codon = T))
  
  
  ## remove syn muts and keep unique changes
  aa1=unique(aa1[!aa1%in%aa])
  aa2=unique(aa2[!aa2%in%aa])
  aa3=unique(aa3[!aa3%in%aa])
  
  # count number of aa reached
  # by 1 and 2 and unique
  aa12=unique(c(aa1,aa2))
  # by 3 and unique
  aa3only=aa3[!aa3%in%c(aa12)]
  
  ## introduce data into table for each mutation type
  sum.res$codons_by_single[sum.res$codon==r]=paste0(cods1,collapse = ";")
  sum.res$AA_reached_by_single_mutation[sum.res$codon==r]=paste0(unique(aa1),collapse = ";")
  sum.res$codons_by_double[sum.res$codon==r]=paste0(cods2,collapse = ";")
  sum.res$AA_reached_by_double_mutation[sum.res$codon==r]=paste0(unique(aa2),collapse = ";")
  sum.res$codons_by_triple[sum.res$codon==r]=paste0(cods3,collapse = ";")
  sum.res$AA_reached_by_triple_mutation[sum.res$codon==r]=paste0(unique(aa3),collapse = ";")
  
  ## sum up the number reached by each type of mutation
  sum.res$number_of_AA_reached_by_single[sum.res$codon==r]=length(aa1)
  sum.res$number_of_AA_reached_by_single_double[sum.res$codon==r]=length(aa12)
  sum.res$number_of_AA_reached_by_triple_only[sum.res$codon==r]=length(aa3only)
  rm(aa,aa1,aa2,aa3,aa3only,cods1,cods2,cods3,r)
}
rm(cod.df)
write.csv(sum.res,"aa_reached_by_different_number_of_mutations.csv",row.names = F)
