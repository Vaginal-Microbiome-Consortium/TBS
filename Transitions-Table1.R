#CST Transition Modeling
#Lacto, Other, L_iners

cstDF <- merge(vagitypeDF,
myClinData[myClinData$BodySite=="Vagina" &
myClinData$subjectType=="Mom",
c("true_ga",  "PID",  "VisitNum",
   "trimester", "SampleID")], by="SampleID")
cstDF <- cstDF %>% filter(true_ga > 0)
cstDF <- cstDF %>%
mutate(msmType = ifelse(cstDF$vagitypes %in% c(
"Atopobium_vaginae",
"Gardnerella_vaginalis", "Lactobacillus_crispatus_cluster", "Lactobacillus_gasseri_cluster", "Lactobacillus_iners", "Lactobacillus_jensenii", "Lachnospiraceae_BVAB1",
"No Type", "Prevotella_cluster2", "Sneathia_amnii", "TM7_OTU-H1"
),
as.character(vagitypes), "Other"))

cstDF$msmType <- factor(cstDF$msmType)
# levels are in alphabetical order
levels(cstDF$msmType)  <-  c(
"Other",
"Other",
"Other",
"Lacto",
"Lacto",
"Liners",
"Lacto",
"Other",
"Other",
"Other",
"Other", "Other"
)
cstDF$msmType  <-  factor(cstDF$msmType,  levels=c("Lacto",  "Liners",  "Other"))
cbind(1:length(levels(cstDF$msmType)),  levels(cstDF$msmType))

##	[,1] [,2]
## [1,] "1"	"Lacto"
## [2,] "2"	"Liners"
## [3,] "3"	"Other"




#Fit Markov chain model. Transition probabilities for a 90-day horizon.

library(msm)
myqmatrix  <-  matrix(rep(1,length(levels(cstDF$msmType))^2),
ncol=length(levels(cstDF$msmType)))
##  Encode  states  as  integers. cstDF$msmTypeInt <- as.integer(cstDF$msmType) cstDF$msmTrue_ga <- as.integer(cstDF$true_ga)
cstDF <- cstDF %>%
arrange(ParticipantID, true_ga) %>% mutate(msmSubject=as.integer(as.factor(as.character(ParticipantID))))
##msmType <- cstDF$msmTypeInt
msmCST <- msm(msmTypeInt ~ msmTrue_ga, subject=msmSubject, data=cstDF,
                         qmatrix=myqmatrix,  gen.inits=TRUE) cstpmatrix  <-  pmatrix.msm(msmCST,  t=90,  ci="normal") rownames(cstpmatrix$estimates) <- levels(cstDF$msmType) colnames(cstpmatrix$estimates) <- levels(cstDF$msmType) cstpmatrix


##	Lacto	Liners
## Lacto	0.75355   (0.62822,0.8312)   0.09197 (0.05509,0.1860)
## Liners 0.06631 (0.03323,0.1964) 0.70920  (0.55153,0.8068)
## Other	0.13035 (0.07669,0.2371) 0.34115 (0.23055,0.4401)
##	Other
## Lacto	0.15448	(0.09499,0.2452)
## Liners	0.22449	(0.14166,0.3307)
## Other	0.52850	(0.41762,0.6378)

#Lacto, Other, L_iners. Include race as a covariate. The differences are not significant.
AAsubjects <- as.character(
unique(mySurveyData$SubjectID[
!is.na(mySurveyData$african_american)  & mySurveyData$african_american == "Yes" & 
mySurveyData$caucasian == "No" & mySurveyData$american_indian_or_alaska_native == "No" & mySurveyData$asian == "No" & 
mySurveyData$native_hawaiian == "No"]))
Csubjects <- as.character(
unique(mySurveyData$SubjectID[
!is.na(mySurveyData$caucasian)  & 
mySurveyData$african_american == "No" & 
mySurveyData$caucasian == "Yes" & mySurveyData$american_indian_or_alaska_native == "No" & mySurveyData$asian == "No" & 
mySurveyData$native_hawaiian == "No"]))
Csubjects <- Csubjects[!is.na(Csubjects)]
cstRaceDF <- cstDF[cstDF$ParticipantID %in% c(AAsubjects, Csubjects),] 
cstRaceDF$ethnicity <- NA
cstRaceDF$ethnicity[cstRaceDF$ParticipantID %in% AAsubjects] <- "AA" cstRaceDF$ethnicity[cstRaceDF$ParticipantID %in% Csubjects] <- "C" cstRaceDF$ethnicity <- factor(cstRaceDF$ethnicity)
msmCST <- msm(msmTypeInt ~ msmTrue_ga,
subject=msmSubject, data=cstRaceDF, qmatrix=myqmatrix, gen.inits=TRUE, covariates =~ ethnicity)
AApmatrix  <-  pmatrix.msm(msmCST,covariates=list(ethnicity="AA"),t=90,  ci="normal") 
Cpmatrix  <-  pmatrix.msm(msmCST,covariates=list(ethnicity="C"),t=90,  ci="normal") rownames(AApmatrix$estimates) <- levels(cstDF$msmType) 
colnames(AApmatrix$estimates) <- levels(cstDF$msmType)
AApmatrix

##	Lacto	Liners
##	Lacto	0.760357	( 5.495e-01,1.1603) 0.120908 (2.229e-109,0.2695)
##	Liners	0.008111	( 2.298e-03,1.1603) 0.757802 (2.229e-109,0.8648)
##	Other	0.044507	( 1.397e-02,0.9583) 0.336228 (5.928e-109,0.4615)
##		Other	
##	Lacto	0.118735	( 4.742e-02,0.4367)
##	Liners	0.234086	( 5.565e-02,0.4567)
##	Other	0.619266	( 4.507e-01,1.1603)

rownames(Cpmatrix$estimates) <- levels(cstDF$msmType) 
colnames(Cpmatrix$estimates) <- levels(cstDF$msmType) 
Cpmatrix

##		Lacto	Liners
##	Lacto	0.82561	(0.19957,0.9103) 0.02509 (0.00739,0.6716)
##	Liners	0.09302	(0.01850,0.4837) 0.83070 (0.34341,0.9463)
##	Other	0.27175	(0.08964,0.5366) 0.13773 (0.04039,0.4015)
##		Other	
##	Lacto	0.14930	(0.06452,0.2754)
##	Liners	0.07628	(0.01582,0.3768)
##	Other	0.59052	(0.28805,0.7788)


#Lcrisp/Ljens/Lgass, Liners, Gvag, BVAB1/Sneathia/Pcluster2/TM7, Other

cstDF <- merge(vagitypeDF,
myClinData[myClinData$BodySite=="Vagina" &
myClinData$subjectType=="Mom",
c("true_ga",  "PID",  "VisitNum",
   "trimester", "SampleID")], by="SampleID")
cstDF <- cstDF %>% filter(true_ga > 0)
cstDF <- cstDF %>%
mutate(msmType = ifelse(cstDF$vagitypes %in% c(
"Atopobium_vaginae",
"Gardnerella_vaginalis", "Lachnospiraceae_BVAB1", "Lactobacillus_crispatus_cluster", "Lactobacillus_gasseri_cluster", "Lactobacillus_iners", "Lactobacillus_jensenii",
"No Type", "Prevotella_cluster2", "Sneathia_amnii", "TM7_OTU-H1"
),
as.character(vagitypes), "Other"))
cstDF$msmType <- factor(cstDF$msmType)
# levels are in alphabetical order
levels(cstDF$msmType)  <-  c(
"BVAB1",
"Gvag",
"BVAB1", "Lcr_je_ga", "Lcr_je_ga", "Liners", "Lcr_je_ga", "Other",
"Other",
"BVAB1",
"BVAB1", "BVAB1"
)
cstDF$msmType  <-  factor(cstDF$msmType,  levels=c("Lcr_je_ga",  "Liners",  "BVAB1",  "Gvag",  "Other"))
cbind(1:length(levels(cstDF$msmType)),  levels(cstDF$msmType))
##	[,1] [,2]
## [1,] "1"	"Lcr_je_ga"
## [2,] "2"	"Liners"
## [3,] "3"	"BVAB1"
## [4,] "4"	"Gvag"
## [5,] "5"	"Other"

#Fit Markov chain model. Transition probabilities for a 90-day horizon.

library(msm)
myqmatrix  <-  matrix(rep(1,length(levels(cstDF$msmType))^2),
ncol=length(levels(cstDF$msmType)))
##  Encode  states  as  integers. cstDF$msmTypeInt <- as.integer(cstDF$msmType) cstDF$msmTrue_ga <- as.integer(cstDF$true_ga)
cstDF <- cstDF %>%
arrange(ParticipantID, true_ga) %>% mutate(msmSubject=as.integer(as.factor(as.character(ParticipantID))))
msmCST <- msm(msmTypeInt ~ msmTrue_ga, subject=msmSubject, data=cstDF,
                         qmatrix=myqmatrix,  gen.inits=TRUE) cstpmatrix  <-  pmatrix.msm(msmCST,  t=90,  ci="normal") rownames(cstpmatrix$estimates) <- levels(cstDF$msmType) colnames(cstpmatrix$estimates) <- levels(cstDF$msmType) 
cstpmatrix

##	Lcr_je_ga	Liners
##	Lcr_je_ga	0.75883	(9.175e-96,Inf)	0.09295	(4.144e-02,Inf)
##	Liners	0.07279	(6.920e-95,Inf)	0.70997	(4.556e-01,Inf)
##	BVAB1	0.06928	(9.175e-96,Inf)	0.34373	(5.704e-02,Inf)
##	Gvag	0.03632	(3.744e-94,Inf)	0.34710	(1.732e-01,Inf)
##	Other	0.26039	(3.014e-94,Inf)	0.32296	(1.918e-01,Inf)
##		BVAB1		Gvag	
##	Lcr_je_ga	0.03546	(3.024e-62,Inf)	0.02573	(5.577e-08,Inf)
##	Liners	0.10847	(3.429e-58,Inf)	0.02585	(1.304e-06,Inf)
##	BVAB1	0.44245	(3.024e-62,Inf)	0.03101	(6.406e-08,Inf)
##	Gvag	0.33228	(7.871e-57,Inf)	0.20930	(2.950e-06,Inf)
##	Other	0.18313	(1.975e-57,Inf)	0.08672	(2.950e-06,Inf)
##		Other			
##	Lcr_je_ga	0.08703	(2.091e-02,Inf)
##	Liners	0.08292	(4.888e-02,Inf)
##	BVAB1	0.11353	(1.232e-02,Inf)
##	Gvag	0.07500	(2.255e-02,Inf)
##	Other	0.14681	(7.592e-02,Inf)


#Without confidence intervals:


cstpmatrix <- pmatrix.msm(msmCST, t=90) rownames(cstpmatrix) <- levels(cstDF$msmType) colnames(cstpmatrix) <- levels(cstDF$msmType) cstpmatrix


##	Lcr_je_ga	Liners	BVAB1	Gvag	Other
## Lcr_je_ga	0.75883120	0.09295199	0.0354602	0.02572613	0.08703047
## Liners	0.07278871	0.70997215	0.1084702	0.02584629	0.08292270
## BVAB1	0.06927681	0.34372788	0.4424534	0.03101058	0.11353136
## Gvag	0.03632308	0.34710271	0.3322757	0.20929913	0.07499934
## Other	0.26038982	0.32295834	0.1831266	0.08671559	0.14680964





