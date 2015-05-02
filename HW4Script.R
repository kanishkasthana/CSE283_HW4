#This Script was written by Kanishk Asthana kasthana@eng.ucsd.edu.
#The only thing you need to do is change the input peptide and everything else will be calculated accordingly
input="SLAMMER"
#Source of data: http://www.bmrb.wisc.edu/ref_info/aadata.dat
atomdata=read.csv("atomicfrequencies.txt", header=FALSE)
#Extracting Relevant Data
atomdata=atomdata[,c(1,3,12,13,14,15,16)]
#Labeling Columns
names(atomdata)=c("Symbol","Name","Carbons","Hydrogens","Nitrogens","Oxygens","Sulphurs")
#Converting NA values in Sulphur to Zeros
atomdata$Sulphurs[is.na(atomdata$Sulphurs)]=0
#Converting Data types to Numeric
atomdata$Nitrogens=as.numeric(atomdata$Nitrogens)
atomdata$Oxygens=as.numeric(atomdata$Oxygens)

#Generating Vector of Amino acid elements
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

#Part A of Problem 1


#Now that we know the relative abundance of the isotope of carbon is 0.01 we can do a biomial expansion
#to get the different probabilities of getting different combinations of number of c-12 and c-13 in amino acids

carbonDistribution=dbinom(0:sumofallAminoAcids["Carbons"],sumofallAminoAcids["Carbons"],0.01)
pdf("CarbonBinomiaDistribution.pdf")
#Plotting Carbon Distribution
plot(0:sumofallAminoAcids["Carbons"],carbonDistribution,type="h",main="Carbon Binomial Probability Distribution", ylab="Probability",xlab="Number of C-13s")
dev.off()
pdf("OxygenBinomialDistrbution.pdf")
oxygenDistribution=dbinom(0:sumofallAminoAcids["Oxygens"],sumofallAminoAcids["Oxygens"],0.1)
#Plotting Oxygen Distribution
plot(0:sumofallAminoAcids["Oxygens"],oxygenDistribution,type="h",main="Oxygen Binomial Probability Distribution", ylab="Probability", xlab="Number of O-18s")
dev.off()
#Taking the outer product of the two distributions to get the joint distribution for all possible combinations of oxygen and carbon atoms
combinedMatrix=outer(carbonDistribution,oxygenDistribution,FUN="*")

#Total Mass of Carbon atoms in the input peptide can also vary, the variation can be found out as follows:

carbonMasses=seq(12*sumofallAminoAcids["Carbons"],13*sumofallAminoAcids["Carbons"],by=1)
oxygenMasses=seq(16*sumofallAminoAcids["Oxygens"],18*sumofallAminoAcids["Oxygens"],by=2)

#Taking outer sum of the total possible carbon and oxygen Masses
combinedMasses=outer(carbonMasses,oxygenMasses,FUN="+")

#Adding masses for Hydrogen, Nitrogen and Sulphur

massToAdd=1*sumofallAminoAcids["Hydrogens"]+14*sumofallAminoAcids["Nitrogens"]+32*sumofallAminoAcids["Sulphurs"]

#New adjusted Masses:
combinedMasses=combinedMasses+massToAdd

uniqueMasses=unique(as.numeric(combinedMasses))
print(uniqueMasses)
#The objective now is calculate the probability of occurance of these unique masses from the combinedProbability Distribution

probabilitiesForUniqueMasses=sapply(uniqueMasses,function(mass){
  #Getting Logical indexes from the Mass table for using with the probability Distribution table
  logicalIndexesForDistribution=(combinedMasses==mass);
  massProbabilities=combinedMatrix[logicalIndexesForDistribution]
  
  #Returing total probabilty for each mass
  return(sum(massProbabilities))
});

print("Ordered Probabilties for Masses:")
print(probabilitiesForUniqueMasses)
pdf("Problem1AFigure.pdf")
plot(uniqueMasses,probabilitiesForUniqueMasses,type='h',xlab="Mass",ylab="Probability of Occurance",main="Mass Spectrum for Part A")
sink("Problem1AIsotopeProfile.txt")
print("Isotope Profile for 5 lowest masses P0,P1,...P4");
print(probabilitiesForUniqueMasses[1:5])
unlink("Problem1AIsotopeProfile.txt")
dev.off()

#Part B of Problem 1

#For the case of Deuterium having 0.6 percent of hydrogens replaced is similar to having a relative abundance of 0.6 for deuterium
hydrogenDistribution=dbinom(0:sumofallAminoAcids["Hydrogens"],sumofallAminoAcids["Hydrogens"],0.6)

#Plotting Hydrogen Distribution
pdf("HydrogenBinomialDistribution.pdf")
plot(0:sumofallAminoAcids["Hydrogens"],hydrogenDistribution,type="h",main="Hydrogen Binomial Probability Distribution", ylab="Probability",xlab="Number of H-2s")
dev.off()

#Taking the outerproduct of hydrogenDistribution with combinedMatrix to get the joint distribution for hydrogen, oxygen and carbon
completeMatrix=outer(combinedMatrix,hydrogenDistribution,FUN="*")

#Similarly the distribution of masses can be found out like in the previous case

#Taking outer sum of the total possible carbon and oxygen Masses
combinedMasses=outer(carbonMasses,oxygenMasses,FUN="+")

#Computing Hydrogen Masses
hydrogenMasses=seq(1*sumofallAminoAcids["Hydrogens"],2*sumofallAminoAcids["Hydrogens"],by=1)
#Taking the outer sum of above matrix with hydrogen masses
completeMasses=outer(combinedMasses,hydrogenMasses,FUN="+")
#Adding masses for Nitrogen and Sulphur
massToAdd=14*sumofallAminoAcids["Nitrogens"]+32*sumofallAminoAcids["Sulphurs"]
#New adjusted Masses:
completeMasses=completeMasses+massToAdd

#Unique Masses:
uniqueCompleteMasses=unique(as.numeric(completeMasses))

#Calculating summed probabilities
probabilitiesForUniqueCompleteMasses=sapply(uniqueCompleteMasses, function(mass){
  #Getting Logical indexes from the Mass table for using with the probability Distribution table
  logicalIndexesForDistribution=(completeMasses==mass);
  massProbabilities=completeMatrix[logicalIndexesForDistribution]
  #Returing total probabilty for each mass
  return(sum(massProbabilities))
});

#print("Ordered Probabilties for Masses after including Hydrogen:")
#print(probabilitiesForUniqueCompleteMasses)
pdf("Problem1BFigure.pdf")
plot(uniqueCompleteMasses,probabilitiesForUniqueCompleteMasses,type='h',xlab="Mass",ylab="Probability of Occurance",main="Mass Spectrum for Part B")
sink("Problem1BIsotopeProfile.txt")
print("Isotope Profile for 5 lowest masses P0,P1,...P4");
print(probabilitiesForUniqueCompleteMasses[1:5])
sink()
unlink("Problem1AIsotopeProfile.txt")
dev.off()
