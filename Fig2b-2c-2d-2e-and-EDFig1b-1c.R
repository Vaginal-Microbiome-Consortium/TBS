#VaHMP: Pregnant vs. Non-Pregnant Matched

#Get 300x2 matched screen IDs from file

#-----------------------------------------------------------------------------------------------------------------
#pregnant	african_american		caucasian	hispanic_or_latino	Total
#No		156			61		83			300
#Yes		156			61		83			300
#-----------------------------------------------------------------------------------------------------------------

#Calculate vagitypes

vagitypes <- apply(mydata[,1:numuclusttaxa],1,which.max) 
maxprop <- mydata[matrix(c(1:nrow(mydata),vagitypes),ncol=2)] 
vagitypes <- colnames(mydata)[vagitypes]
vagitypes[maxprop < 0.30] <- "No Type" 
mydata$vagitypes <- vagitypes

#Barplot of pregnant subjects

mybarplot(mydata[mydata$pregnant == "Yes",1:numuclusttaxa])

#Barplot of nonpregnant subjects

mybarplot(mydata[mydata$pregnant == "No",1:numuclusttaxa])

#Barplot of pregnant AA subjects.

mybarplot(mydata[mydata$pregnant== "Yes" &
       mydata$ethnicity == "african_american",1:numuclusttaxa])

#Barplot of nonpregnant AA subjects.

mybarplot(mydata[mydata$pregnant == "No" &
       mydata$ethnicity == "african_american",1:numuclusttaxa])

#Barplot of pregnant EA subjects

mybarplot(mydata[mydata$pregnant == "Yes" &
       mydata$ethnicity == "caucasian",1:numuclusttaxa])

#Barplot of nonpregnant EA subjects

mybarplot(mydata[mydata$pregnant == "No" &
       mydata$ethnicity == "caucasian",1:numuclusttaxa])

#Barplot of pregnant H subjects

mybarplot(mydata[mydata$pregnant == "Yes" &
       mydata$ethnicity == "hispanic_or_latino",1:numuclusttaxa])

#Barplot of nonpregnant H subjects

mybarplot(mydata[mydata$pregnant== "No" &
       mydata$ethnicity == "hispanic_or_latino",1:numuclusttaxa])




#Table and chi-square test for ethnicity and income. There is a relationship between income and ethnicity.

myTableData <- mydata %>%
dplyr::select(ethnicity, income) %>%
mutate(under20K = (!is.na(income) & income %in% c("less than 15k", "15k to 20k"))) %>%
dplyr::select(ethnicity, under20K) (mytable <- table(myTableData))

#-------------------------------------------------
##under20K
## ethnicity		FALSE TRUE
##african_american	86	226
##caucasian		84	38
##hispanic_or_latino	48	118

chisq.test(mytable)

##	Pearson's Chi-squared test
## data:	mytable
## X-squared = 70.093, df = 2, p-value = 6.019e-16

#Diversity by pregnancy status.

library(dplyr) mydata <- mydata %>%
      mutate(alphaDiv = 1/rowSums(.[,1:numuclusttaxa]^2))

mydata %>%
      select(pregnant, alphaDiv) %>% group_by(pregnant) %>%
      summarise(meanAlpha = mean(alphaDiv), sdAlpha=sd(alphaDiv), nAlpha = n()) %>% 
      mutate(seAlpha = sdAlpha/sqrt(nAlpha),
      widthCI = qt(1-0.025, nAlpha -1)*seAlpha,
      lowerCI = meanAlpha - qt(1-0.025, nAlpha -1)*seAlpha, upperCI = meanAlpha + qt(1-0.025, nAlpha -1)*seAlpha) %>%
kable()

#------------------------------------------------------------------------------------------------------------------------------------------
#pregnant	meanAlpha	sdAlpha	nAlpha	seAlpha	widthCI	lowerCI	upperCI
#No	2.342535	1.8937719	300	0.109337	0.2151675	2.127367	2.557702
#Yes	1.658340	0.9229061	300	0.053284	0.1048592	1.553481	1.763199
#------------------------------------------------------------------------------------------------------------------------------------------

#Diversity by ethnicity and pregnancy status

mydata %>%
      select(pregnant, ethnicity, alphaDiv) %>% 
      group_by(pregnant, ethnicity) %>%
summarise(meanAlpha = mean(alphaDiv), sdAlpha=sd(alphaDiv), nAlpha = n()) %>% mutate(seAlpha = sdAlpha/sqrt(nAlpha),
      widthCI = qt(1-0.025, nAlpha -1)*seAlpha,
       lowerCI = meanAlpha - qt(1-0.025, nAlpha -1)*seAlpha, upperCI = meanAlpha + qt(1-0.025, nAlpha -1)*seAlpha) %>%
kable()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pregnant	ethnicity	meanAlpha	sdAlpha	nAlpha	seAlpha	widthCI	lowerCI	upperCI
No	african_american	2.604279	1.8200963	156	0.1457243	0.2878620	2.316417	2.892141
No	caucasian	2.084248	2.2879328	61	0.2929398	0.5859668	1.498281	2.670215
No	hispanic_or_latino	2.040407	1.6470130	83	0.1807832	0.3596353	1.680772	2.400043
Yes	african_american	1.794544	1.0250597	156	0.0820705	0.1621209	1.632423	1.956665
Yes	caucasian	1.555093	0.9117687	61	0.1167400	0.2335148	1.321578	1.788608
Yes	hispanic_or_latino	1.478221	0.6597555	83	0.0724176	0.1440616	1.334160	1.622283
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Two-sample t-test for difference in diversity between preg/non-preg.

t.test(alphaDiv ~ pregnant, data=mydata)

##
##	Welch Two Sample t-test
##
## data:	alphaDiv by pregnant
## t = 5.6252, df = 433.44, p-value = 3.327e-08
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##	0.4451379 0.9232522
## sample estimates:
##	mean in group No mean in group Yes
##	2.342535	1.658340

#Two-sample t-tests for differences in diversity between preg/non-preg by ethnicity. FDR correction applied across three groups

mydata %>%
      select(pregnant, ethnicity, alphaDiv) %>% 
      group_by(ethnicity, pregnant)	%>% 
      summarise(value=list(alphaDiv)) %>% spread(pregnant, value) %>% 
      group_by(ethnicity) %>%
mutate(p_value	=t.test(unlist(No), unlist(Yes))$p.value, test_stat=t.test(unlist(No), unlist(Yes))$statistic) %>%
      ungroup() %>%
mutate(q_value=p.adjust(p_value, method="fdr")) %>% dplyr::select(ethnicity,p_value, test_stat, q_value) %>% 
kable()

#----------------------------------------------------------------------------------------
#ethnicity	p_value	test_stat	q_value
#african_american	0.0000023	4.841586	0.0000069
#caucasian	0.0973152	1.678024	0.0973152
#hispanic_or_latino	0.0047050	2.886734	0.0070575
#-----------------------------------------------------------------------------------------

#Two-sample t-tests for differences in diversity between preg/non-preg by ethnicity for subjects with household incomes of $20k or less. FDR correction applied across three groups

mydata %>%
       dplyr::filter(income %in% c("15k to 20k", "less than 15k")) %>% 
       dplyr:: select(pregnant, ethnicity, alphaDiv) %>% 
       group_by(ethnicity, pregnant)	%>% 
       summarise(value=list(alphaDiv)) %>%
       spread(pregnant, value) %>% group_by(ethnicity) %>%
mutate(p_value	=t.test(unlist(No), unlist(Yes))$p.value, test_stat=t.test(unlist(No), unlist(Yes))$statistic) %>%
       ungroup() %>%
mutate(q_value=p.adjust(p_value, method="fdr")) %>% dplyr::select(ethnicity,p_value, test_stat, q_value) %>% 
kable()

#-------------------------------------------------------------------------------------------
#ethnicity	p_value	test_stat	q_value
#african_american	0.0001220	3.926096	0.0003659
#caucasian	0.1373284	1.552434	0.1373284
#hispanic_or_latino	0.0409553	2.081511	0.0614329
#--------------------------------------------------------------------------------------------

## # A tibble: 2 x 2
##	ethnicity	n
##	<chr>	<int>
## 1 african_american	19
## 2 caucasian	63

library(ggplot2) alphaEthBoxData <- mydata %>%
		select(pregnant, ethnicity, alphaDiv) %>% 
		melt()

## Using pregnant, ethnicity as id variables

alphaEthBoxDataCopy <- alphaEthBoxData alphaEthBoxDataCopy$ethnicity <- "All"
alphaEthBoxData <- rbind(alphaEthBoxData, alphaEthBoxDataCopy) alphaEthBoxData$ethnicity <- factor(alphaEthBoxData$ethnicity,
levels=c("african_american",
"caucasian", "hispanic_or_latino", "All"),
labels=c("AA", "EA", "SA", "All"))

alphaEthBoxPlot <- ggplot(alphaEthBoxData,
aes(x=ethnicity, y=value), group=pregnant) +
geom_boxplot(aes(fill=pregnant)) + ylab("Alpha Diversity") +
xlab("Race/Ethnicity") print(alphaEthBoxPlot

#Pregnant vs Non-Pregnant
#Vagitype Distribution
#Check if vagitypes are different between pregnant and non-pregnant. FDR q-value reported for each vagitype.

checkTypes <- c("crispatus", "iners", "jensenii", "gasseri", "vaginalis", "BVAB1") 
mychisq <- function(type, mydata, mycol) {
      mytable <- table(mydata[,c(mycol)], grepl(type, mydata$vagitypes)) 
      mychisq <- chisq.test(mytable)
      mychisq$p.value
}
myfisher <- function(type, mydata, mycol) {
      mytable <- table(mydata[,c(mycol)], grepl(type, mydata$vagitypes)) 
      myfisher <- fisher.test(mytable)
      myfisher$p.value
}

mypvalues <- sapply(checkTypes, myfisher, mydata, "pregnant") myqvalues <- p.adjust(mypvalues, method="fdr") kable(myqvalues, col.names=NA)

#------------------------------------------------
#crispatus	0.4986673
#iners		0.0046451
#jensenii	0.3333552
#gasseri	0.2185131
#vaginalis	0.0142673
#BVAB1	0.3519569
#------------------------------------------------

#Lacto vs. non-Lacto

lactoNonLacto <- rep("NonLacto", nrow(mydata)) 
lactoNonLacto[grep("Lactobacillus", mydata$vagitypes)] <- "Lacto" 
mytable <- table(mydata$pregnant, lactoNonLacto)
kable(mytable)

#----------------------------
#Lacto	NonLacto
#No	153	147
#Yes	204	96
#----------------------------

fisher.test(mytable)


##
##	Fisher's Exact Test for Count Data
##
## data:	mytable
## p-value = 3.057e-05
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##	0.3465539 0.6917677
## sample estimates:
## odds ratio
##	0.4903866

#Pregnant versus non pregnant AA.

myttest <- t.test(alphaDiv ~ pregnant, data=mydata[mydata$ethnicity =="african_american",]) myttest

##
##	Welch Two Sample t-test
##
## data:	alphaDiv by pregnant
## t = 4.8416, df = 244.34, p-value = 2.286e-06
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##	0.4803073 1.1391620
## sample estimates:
##	mean in group No mean in group Yes
##	2.604279	1.794544

#Vagitype Distribution

mypvalues <- sapply(checkTypes,
	myfisher,
	mydata[mydata$ethnicity == "african_american",] , "pregnant")
myqvalues <- p.adjust(mypvalues, method="fdr")
kable(myqvalues, col.names=NA)

#--------------------------------
#crispatus	0.2257246
#iners	0.0832681
#jensenii	0.2257246
#gasseri	1.0000000
#vaginalis	0.0761553
#BVAB1	0.2813142
#----------------------------------


#Lacto vs. non-Lacto

mydataAA <- mydata[mydata$ethnicity == "african_american",] 
lactoNonLacto <- rep("NonLacto", nrow(mydataAA)) lactoNonLacto[grep("Lactobacillus", mydataAA$vagitypes)] <- "Lacto" 
mytable <- table(mydataAA$pregnant, lactoNonLacto)
kable(mytable)

#--------------------------
#Lacto	NonLacto
#No	66	90
#Yes	91	65
#--------------------------


fisher.test(mytable)


##
##	Fisher's Exact Test for Count Data
##
## data:	mytable
## p-value = 0.006485
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##	0.3257702 0.8416434
## sample estimates:
## odds ratio
##	0.5249109

#Pregnant versus non pregnant EA.

myttest <- t.test(alphaDiv ~ pregnant, data=mydata[mydata$ethnic=="caucasian",]) myttest


##
##	Welch Two Sample t-test
##
## data:	alphaDiv by pregnant
## t = 1.678, df = 78.589, p-value = 0.09732
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##	-0.09857298	1.15688294
## sample estimates:
##	mean in group No mean in group Yes
##	2.084248	1.555093

#Vagitype Distribution

mypvalues <- sapply(checkTypes,
       myfisher,
       mydata[mydata$ethnicity == "caucasian",], "pregnant")
myqvalues <- p.adjust(mypvalues, method="fdr")
kable(myqvalues, col.names=NA)

#----------------------------------
#crispatus	0.1890621
#iners		0.1890621
#jensenii	1.0000000
#gasseri		1.0000000
#vaginalis	1.0000000
#BVAB1		0.4876033
#------------------------------------


#Lacto vs. non-Lacto

mydataEA <- mydata[mydata$ethnicity == "caucasian",] 
lactoNonLacto <- rep("NonLacto", nrow(mydataEA))
lactoNonLacto[grep("Lactobacillus", mydataEA$vagitypes)] <- "Lacto" 
mytable <- table(mydataEA$pregnant, lactoNonLacto)
kable(mytable)

#---------------------------
#Lacto	NonLacto
#No	45	16
#Yes	45	16
#--------------------------

fisher.test(mytable)

##
##	Fisher's Exact Test for Count Data
##
## data:	mytable
## p-value = 1
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##	0.4119491 2.4274847
## sample estimates:
## odds ratio
##	1

#Pregnant versus non pregnant H.

myttest <- t.test(alphaDiv ~ pregnant,
       data=mydata[mydata$ethnic=="hispanic or latino",])
myttest


##
##	Welch Two Sample t-test
##
## data:	alphaDiv by pregnant
## t = 2.8867, df = 107.66, p-value = 0.004705
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##	0.1761474 0.9482248
## sample estimates:
##	mean in group No mean in group Yes
##	2.040407	1.478221

#---------------------------------------------------
#Lacto	NonLacto
#No	42	41
#Yes	68	15
#-----------------------------------------------------

fisher.test(mytable)

##
##	Fisher's Exact Test for Count Data
##
## data:	mytable
## p-value = 3.245e-05
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##	0.1037513 0.4808143
## sample estimates:
## odds ratio
##	0.2281236



