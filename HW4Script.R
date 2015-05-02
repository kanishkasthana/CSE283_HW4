input="SLAMMER"
#Each amino acid is stored
atomdata=read.csv("atomicfrequencies.txt", header=FALSE)
#Extracting Relevant Data
atomdata=atomdata[,c(1,3,12,13,14,15,16)]
#Labeling Columns
names(atomdata)=c("Symbol","Name","Carbons","Hydrogens","Nitrogens","Oxygens","Sulphurs")
#Converting NA values in Sulphur to Zeros
atomdata$Sulphurs[is.na(atomdata$Sulphurs)]=0
atomdata$Nitrogens=as.numeric(atomdata$Nitrogens)
atomdata$Oxygens=as.numeric(atomdata$Oxygens)

totalatoms=c(0,0,0,0,0)
#Vector of Amino acid elements
aminoAcidVector=unlist(strsplit(input,""));

#All data gives a table of number of atoms in each Amino Acid for input
alldata=sapply(aminoAcidVector,function(aminoAcid){
  data=unlist(data.frame(atomdata[atomdata$Symbol==aminoAcid,c(3,4,5,6,7)]))
  return(data)
});

print(alldata)
#Sum of Atomic Composition of all Amino Acids
sumofallAminoAcids=rowSums(alldata)

#Calculating number of peptide bonds
numberofPeptideBonds=nchar(input)-1
#Number of Hydrogens to Subtract
hydrogenstoSubtract=2*numberofPeptideBonds
#Number of Oxygens to Subtract
oxygenstoSubtract=numberofPeptideBonds
#Correcting Values
sumofallAminoAcids["Hydrogens"]=sumofallAminoAcids["Hydrogens"]-hydrogenstoSubtract
sumofallAminoAcids["Oxygens"]=sumofallAminoAcids["Oxygens"]-oxygenstoSubtract
#Printing values
print(sumofallAminoAcids)

#Now that we know the relative abundance of the isotope of carbon is 0.01 we can do a biomial expansion
#to get the different probabilities of getting different combinations of number of c-12 and c-13 in amino acids

carbonDistribution=dbinom(0:sumofallAminoAcids["Carbons"],sumofallAminoAcids["Carbons"],0.01)
#Plotting Carbon Distribution
plot(0:sumofallAminoAcids["Carbons"],carbonDistribution,type="h",main="Carbon Bionomial Probability Distribution", ylab="Probability",xlab="Number of C-13s")

oxygenDistribution=dbinom(0:sumofallAminoAcids["Oxygens"],sumofallAminoAcids["Oxygens"],0.1)
#Plotting Oxygen Distribution
plot(0:sumofallAminoAcids["Oxygens"],oxygenDistribution,type="h",main="Oxygen Bionomial Probability Distribution", ylab="Probability", xlab="Number of O-18s")



