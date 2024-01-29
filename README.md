# Seasonal-polyandry-in-conch
R code and data for Hooks and Burgess, Variation in polyandry, reproductive output, and within-brood genetic diversity in a marine snail population across seasons and years. Marine Ecology Progress Series

This README file was generated on 2024-Jan-29 by Scott Burgess

GENERAL INFORMATION

1. Title of Dataset: Variation in polyandry, reproductive output, and within-brood genetic diversity in a marine snail population across seasons and years.

2. Author Information
	A. Principal Investigator Contact Information
		Name: Scott Burgess
		Institution: Florida State University
		Address: 319 Stadium Drive, Tallahassee, FL, USA 32306
		Email: sburgess@bio.fsu.edu

	B. Lead author Contact Information
		Name: Alexandra Hooks
		Institution: Florida State University
		Address: 319 Stadium Drive, Tallahassee, FL, USA 32306
		Email: hooksap@gmail.com


3. Date of data collection (single date, range, approximate date): 2018-2020 

4. Geographic location of data collection: Turkey Point, Florida, USA

5. Information about funding sources that supported the collection of the data: American Museum of Natural History, Florida State University Department of Biological Science, Florida State University Coastal and Marine Laboratory, Association of University Women, National Science Foundation (Bio Oce 18-29867).



DATA & FILE OVERVIEW

1. File List: 
Figure 2.R
Figure 3.R
Figure 4.R
Figure 5.R

BestCluster_2018.csv
BestCluster_2020.csv
Data for Figure 3 and 4.csv
Data A for Figure 5.csv
Data B for Figure 5.csv



2. Relationship between files: 
Figure 2.R uses BestCluster_2018.csv and BestCluster_2020.csv
Figure 3.R uses Data for Figure 3 and 4.csv
Figure 4.R uses Data for Figure 3 and 4.csv
Figure 5.R uses Data A for Figure 5.csv and Data B for Figure 5.csv




3. Metadata

BestCluster_2018.csv
ClusterIndex: A unique number for each ancestry cluster produced by the program Colony. 
Probability: The probability an individual (OffspringID) belongs to the assigned ClusterIndex
OffspringID: Unique identifier for each hatchling in 2018
FatherID: The inferred father identification estimated from the program Colony.
MotherID: Unique identifier for the known mother of each hatchling

BestCluster_2020.csv
ClusterIndex: A unique number for each ancestry cluster produced by the program Colony. 
Probability: The probability an individual (OffspringID) belongs to the assigned ClusterIndex
OffspringID: Unique identifier for each hatchling in 2020
FatherID: The inferred father identification estimated from the program Colony.
MotherID: Unique identifier for the known mother of each hatchling
 

Data for Figure 3 and 4.csv
mother.id: Unique identifier for each mother
year: Year that mothers and offspring were sampled.
time.season: Time in the reproductive season when egg strings were laid and offspring were sampled (E=Early=April; L=Late=July)
number.sire: The number of unique reconstructed paternal genotypes estimated from the program Colony.
shannon: Shannon diversity of the paternal genotypes within a mother.
simpson: Simpsons diversity of the paternal genotypes within a mother.
effective.num.sires.shannon:
effective.num.sires.simpson:
N: The number of hatchlings that were genotyped per mother
Hexp: Expected heterozygosity of hatchlings within a brood, estimated from microsatellite markers
Ar: Allelic richness of hatchlings within a brood, estimated from microsatellite markers
capsule.num: The number of egg capsules laid per mother
offspring.num: The average number of hatchlings per egg capsule 
capsule.size: The average capsule size in each mother's egg string
total.output: The estimated total number of hatchlings in the brood


Data A for Figure 5.csv
mother.id: Unique identifier for each mother
year: Year that mothers and offspring were sampled.
time.season: Time in the reproductive season when egg strings were laid and offspring were sampled (E=Early=April; L=Late=July)
capsule.num: The number of egg capsules laid per mother 
offspring.num: The average number of hatchlings per egg capsule 
capsule.num: The number of egg capsules laid per mother
total.output: The estimated total number of hatchlings in the brood

Data B for Figure 5.csv
date.hatched: The date that hatchlings emerged from the egg capsule
date.laid: The date that the egg string was laid
mother.id: Unique identifier for each mother
capsule.id: Unique identifier for each capsule, within each mother
offspring.num: The average number of hatchlings per egg capsule 
capsule.num: The number of egg capsules laid per mother
capsule.size: The average capsule size in each mother's egg string
time.season: Time in the reproductive season when egg strings were laid and offspring were sampled (E=Early=April; L=Late=July)
year: Year that mothers and offspring were sampled.
survive: Indicated if hatchling were alive (1=alive, 0=dead) on the hatching day.
dev.time: The number of days between date that the egg string was laid and the date that hatchlings emerged from the egg capsule

