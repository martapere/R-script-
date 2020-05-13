rm(list = ls())  # Clear environment
dev.off()  # Clear plots
cat("\014")  # Clear console



library("tidyverse")
library("pacman")
pacman::p_load(pacman, rio)

#Exercise 1
library("tidyverse")
read_tsv("small_file.txt") -> small.file

library("pacman")
pacman::p_load(pacman, rio)

rio_txt <- import("~/Desktop/small_file.txt")
head(rio_txt)
small.file <- rio_txt
class(small.file$Category)
class(small.file$Length)
small.file$Category <- as.factor(small.file$Category)
class(small.file$Category)

rio_csv <- import("~/Desktop/Child_Variants.csv")
child <- rio_csv
head(child)
ggplot(child, aes(x = CHR, y = COVERAGE)) + geom_point()

#Exercise 2
rio_csv <- import("~/Desktop/genomes.csv")
genomes <- rio_csv
genomes %>%
  filter(Chromosomes > 40)%>%
  select(Organism, Chromosomes)

genomes %>%
  filter(Plasmids > 0) %>%
  filter(Chromosomes > 1) %>%
  select(Organism, Chromosomes, Plasmids)

genomes %>%
  arrange(desc(Size)) %>%
  slice(1:10) %>%
  select(Organism, Size)

genomes %>%
  select(-Groups)

genomes %>%
  select(starts_with("O"))

genomes %>%
  distinct(Groups) %>%
  nrow()

genomes %>%
  arrange(Size)%>%
  filter(Chromosomes > 0)%>%
  distinct(Chromosomes,.keep_all = TRUE) -> genomes

ggplot(genomes, aes(x=Chromosomes, y=Size)) +geom_point()

ggplot(genomes, aes(x=Chromosomes, y=log(Size))) +geom_point()

#Exercise 3
#Exercise 3
rio_csv<- import ("~/Desktop/cancer_stats.csv")
cancer <- rio_csv

cancer %>%  
  filter(Class=="Digestive System", `Female Cases`> `Male Cases`)

cancer %>%
  filter(is.na('Male Cases'))
cancer %>%
  filter(is.na('Female Cases'))

cancer %>%
  mutate(SurvivalRate =  `Male Deaths`/`Male Cases`)%>%
  arrange(desc(SurvivalRate))%>%
  select(Class, Site, SurvivalRate)

cancer %>%
  filter(str_detect(Site,"acute"))

cancer %>%
  filter(Class == "soft tissue") %>%
  select(Class, Site )    

cancer %>%
  filter(str_length(Site) <=4 )

cancer %>%
  filter(endsWith(Site,"y"))

#Exercise 4
rio_csv<- import ("~/Desktop/tidy_data1.csv")
tidy_data1 <- rio_csv

rio_csv<- import ("~/Desktop/tidy_data2.csv")
tidy_data2 <- rio_csv

rio_csv<- import ("~/Desktop/tidy_data3.csv")
tidy_data3 <- rio_csv

#tidy1
tidy_data1n<- tidy_data1 %>%
  gather(key, value)

tidy_data1n %>%
  filter(!is.na(value))

#tidy2
tidy_data2 %>%
  pivot_longer(cols=A:E, names_to = "sample", values_to = "value")-> tidy_data2 

tidy_data2 %>%
  filter(!is.na(value))

#tidy3
tidy_data3n <-tidy_data3 %>%
  pivot_longer(cols=WT_1:KO_3, names_to = "sample", values_to = "value")

tidy_data3n %>%
  select(Probe_ID:Symbol)%>%
  distinct(Probe_ID, .keep_all = TRUE)

#back to genomes data
genomes %>%
  separate(Groups, c("Domain", "Kingdom", "Class"), ";") %>%
  filter(! str_detect(Organism, "'")) %>%
  filter(! (Kingdom=="Other")) %>%
  filter(! (Class=="Other"))

#Exercise 5
cancer %>%
  mutate(Deaths=`Male Deaths` + `Female Deaths`, Cases = `Male Cases` + `Female Cases`)

child.1<-child %>%
  mutate(Type = if_else(str_length(REF)==1 & str_length(ALT)==1, "SNP", "INDEL")) 

small.file %>%
  group_by(Category)%>%
  summarise(meanLength= mean(Length), sdLength= sd(Length))

child.2 <-child.1 %>% 
  filter(dbSNP==".") %>%
  filter(Type=="SNP") %>%
  group_by(GENE) %>%
  summarise(COVERAGE=mean(COVERAGE),count=n()) %>%
  filter(count>2) %>%
  arrange(desc(COVERAGE))


tidy_data2 %>%
  select(Chr, sample)%>%
  group_by(sample, Chr) %>%
  summarise(valuemean=mean(value))

cancer %>%
  mutate(difMF=abs(`Male Cases` - `Female Cases`)) %>%
  select(Class,difMF) %>%
  arrange(difMF) 

cancer %>%
  filter(!(is.na(`Male Cases`) | is.na(`Female Cases`)))%>%
  mutate(MaleRate=`Male Deaths`/ `Male Cases` , FemaleRate= `Female Deaths`/ `Female Cases` , difFM=abs(MaleRate-FemaleRate)) %>%
  select(Class,Site,difFM) %>%
  arrange(desc(difFM))

cancer %>%
  mutate(Cases = `Male Cases` + `Female Cases`, Deaths = `Male Deaths` + `Female Deaths`, survival= Cases - Deaths/Cases)%>%
  group_by(Class) %>%
  arrange(desc(survival))%>%
  select(Class, Site, survival)


child %>%
  group_by(CHR, GENE) %>%
  summarise(count=n()) %>%
  arrange(desc(count))

#Exercise 9
rio_csv <-import ("~/Desktop/dna_methylation.csv")
dna_methylation <- rio_csv

dna_methylation2<- dna_methylation %>%
  pivot_wider(names_from = State, values_from = Count) %>%
  mutate(metpercent=(100*Meth/(Meth+Unmeth)))

dna_methylation2 %>%
  group_by(Gene, Group)%>%
  summarise(meanMet= mean(metpercent))

rio_txt<-import ("~/Desktop/methylation_annotation.txt")
rio_txt -> methylation_annotation

dna_methylation %>%
  rename(Gene_name=Gene) %>%
  left_join(methylation_annotation)

child %>% 
  select(ALT,REF) %>%
  mutate(allvariants=str_c(ALT,">",REF)) %>%
  group_by(allvariants) %>%
  summarise(count=n()) %>%
  arrange(desc(count))

small.file1 <- small.file%>%
  group_by(Category)%>%
  summarise(minLength = min(Length))%>%
  right_join(small.file)%>%
  mutate(normalizedLength= Length - minLength)

ggplot(small.file1, aes (x=Category, y=normalizedLength)) + geom_point()

#Exercise 7

lowestqual <- function(genename) {
  child %>%
    filter(GENE == genename) %>%
    arrange(QUAL) %>%
    select(QUAL, GENE)%>%
    slice(1)
}

lowestqual("AGRN")


lowestqual2 <- function(tbl) {
  tbl %>%
    arrange(QUAL) %>%
    slice(1)
}

child %>% lowestqual2()


child %>%
  group_by(GENE) %>%
  lowestqual2() %>%
  arrange(QUAL)


