#VaHMP: Pregnant vs. Non-Pregnant All
#Table of pregnancy status and race/ethnicity sample sizes.
#-----------------------------------------------------------------------------------------------------------------------------
#pregnant_derived	african	american	caucasian	hispanic or latino	Total
#No			1184			682		103			1969
#Yes			304			111		198			613
#-----------------------------------------------------------------------------------------------------------------------------

#Calculate vagitypes.

vagitypes <- apply(mydata[,1:numuclusttaxa],1,which.max) maxprop  <-  mydata[matrix(c(1:nrow(mydata),vagitypes),ncol=2)] vagitypes <- colnames(mydata)[vagitypes]
vagitypes[maxprop < 0.30] <- "No Type" mydata$vagitypes <- vagitypes

#Barplot of pregnant subjects.
mybarplot(mydata[mydata$pregnant_derived  ==  "Yes",1:numuclusttaxa])

#Barplot of nonpregnant subjects.
mybarplot(mydata[mydata$pregnant_derived  ==  "No",1:numuclusttaxa])

#Barplot of pregnant AA subjects.
mybarplot(mydata[mydata$pregnant_derived == "Yes" &
       mydata$ethnic  ==  "african  american",1:numuclusttaxa])

#Barplot of nonpregnant AA subjects.
mybarplot(mydata[mydata$pregnant_derived == "No" &
       mydata$ethnic  ==  "african  american",1:numuclusttaxa])

#Barplot of pregnant EA subjects.
mybarplot(mydata[mydata$pregnant_derived == "Yes" &
       mydata$ethnic  ==  "caucasian",1:numuclusttaxa])

#Barplot of nonpregnant EA subjects.
mybarplot(mydata[mydata$pregnant_derived == "No" &
       mydata$ethnic  ==  "caucasian",1:numuclusttaxa])

#Barplot of pregnant H subjects.
mybarplot(mydata[mydata$pregnant_derived == "Yes" &
       mydata$ethnic  ==  "hispanic  or  latino",1:numuclusttaxa])

#Barplot of nonpregnant H subjects.
mybarplot(mydata[mydata$pregnant_derived == "No" &
       mydata$ethnic  ==  "hispanic  or  latino",1:numuclusttaxa])

#Diversity by pregnancy status.
library(dplyr)   mydata <- mydata %>%
       mutate(alphaDiv  =  1/rowSums(.[,1:numuclusttaxa]^2))

mydata %>%
       select(pregnant_derived, alphaDiv) %>% 
       group_by(pregnant_derived) %>%
summarise(meanAlpha  =  mean(alphaDiv),  sdAlpha=sd(alphaDiv),  nAlpha  =  n())  %>% 
mutate(seAlpha = sdAlpha/sqrt(nAlpha),
       widthCI  =  qt(1-0.025,  nAlpha  -1)*seAlpha,
       lowerCI  =  meanAlpha  -  qt(1-0.025,  nAlpha  -1)*seAlpha, 
       upperCI  =  meanAlpha  +  qt(1-0.025,  nAlpha  -1)*seAlpha)  %>%
       kable()
       
       
#--------------------------------------------------------------------------------------------------------------------------------------------
#pregnant_derived	meanAlpha	sdAlpha	nAlpha	seAlpha	widthCI	lowerCI	upperCI
#No	2.214214	1.825688	1969	0.0411437	0.0806899	2.133524	2.294904
#Yes	1.723586	1.077374	613	0.0435147	0.0854563	1.638130	1.809043
#---------------------------------------------------------------------------------------------------------------------------------------------

#Diversity by ethnicity and pregnancy status

mydata %>%
       select(pregnant_derived, ethnic, alphaDiv) %>% 
       group_by(pregnant_derived, ethnic) %>%
       summarise(meanAlpha  =  mean(alphaDiv),  sdAlpha=sd(alphaDiv),  nAlpha  =  n())  %>% 
       mutate(seAlpha = sdAlpha/sqrt(nAlpha),
       widthCI  =  qt(1-0.025,  nAlpha  -1)*seAlpha,
       lowerCI  =  meanAlpha  -  qt(1-0.025,  nAlpha  -1)*seAlpha, 
       upperCI  =  meanAlpha  +  qt(1-0.025,  nAlpha  -1)*seAlpha)  %>%
       kable()


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#pregnant_derived	ethnic	meanAlpha	sdAlpha	nAlpha	seAlpha	widthCI	lowerCI	upperCI
#No	african american	2.394169	1.8448900	1184	0.0536160	0.1051931	2.288976	2.499362
#No	caucasian	1.968356	1.8464544	682	0.0707044	0.1388248	1.829531	2.107181
#No	hispanic or latino	1.773524	1.0467624	103	0.1031406	0.2045788	1.568945	1.978103
#Yes	african american	1.898126	1.2238472	304	0.0701925	0.1381264	1.759999	2.036252
#Yes	caucasian	1.578837	1.0101964	111	0.0958836	0.1900188	1.388819	1.768856
#Yes	hispanic or latino	1.536754	0.7990555	198	0.0567864	0.1119873	1.424767	1.648741
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Two-sample t-test for difference in diversity between preg/non-preg.

t.test(alphaDiv ~ pregnant_derived, data=mydata)
##
##	Welch Two Sample t-test
##
## data:	alphaDiv by pregnant_derived
## t = 8.1927, df = 1758.4, p-value = 4.86e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##	0.3731726 0.6080832
## sample estimates:
##	mean in group No mean in group Yes
##	2.214214	1.723586

#Two-sample t-tests for differences in diversity between preg/non-preg by ethnicity. FDR correction applied across three groups.

mydata %>%
      dplyr::select(pregnant_derived, ethnic, alphaDiv) %>% 
      group_by(ethnic,   pregnant_derived)	%>% 
      summarise(value=list(alphaDiv))  %>% 
      spread(pregnant_derived, value) %>%
      group_by(ethnic) %>%
      mutate(p_value	=t.test(unlist(No), unlist(Yes))$p.value, 
      test_stat=t.test(unlist(No), unlist(Yes))$statistic) %>%
      ungroup() %>%
      mutate(q_value=p.adjust(p_value, method="fdr")) %>% 
      dplyr::select(ethnic,p_value, test_stat, q_value) %>% 
       kable()
       
       #-------------------------------------------------------------------------------------------
       #ethnic			p_value		test_stat	q_value
       #african american	0.0000000	5.615988	0.0000001
       #caucasian		0.0012282	3.269598	0.0018424
       #hispanic or latino	0.0459545	2.010961	0.0459545
       #-------------------------------------------------------------------------------------------
       
       
#Two-sample t-tests for differences in diversity between preg/non-preg by ethnicity for subjects with household incomes of $20k or less. FDR correction applied across three groups.

mydata %>%
      dplyr::filter(income  %in%  c("15k  to  20k",  "less  than  15k"))  %>% 
      dplyr:: select(pregnant_derived, ethnic, alphaDiv) %>% 
      group_by(ethnic, pregnant_derived)	%>% 
      summarise(value=list(alphaDiv))  %>%
      spread(pregnant_derived, value) %>% group_by(ethnic) %>%
      mutate(p_value	=t.test(unlist(No), unlist(Yes))$p.value, 
      test_stat=t.test(unlist(No), unlist(Yes))$statistic) %>%
      ungroup() %>%
      mutate(q_value=p.adjust(p_value, method="fdr")) %>% 
      dplyr::select(ethnic,p_value, test_stat, q_value) %>% 
      kable()
      
      
      
#---------------------------------------------------------------------------------------------------------
# ethnic	p_value	test_stat	q_value
# african american	0.0000100	4.468179	0.0000299
# caucasian	0.0000446	4.182512	0.0000669
# hispanic or latino	0.0444624	2.044865	0.0444624
#---------------------------------------------------------------------------------------------------------

#Two-sample t-tests for differences in diversity between preg/non-preg by ethnicity for subjects with household incomes of $40k or more. FDR correction applied across three groups.

mydata %>%
      dplyr::filter(!is.na(income)) %>%
      dplyr::filter(!(income  %in%  c("20k  to  40k",  "15k  to  20k",  "less  than  15k")))  %>%
      dplyr:: select(pregnant_derived, ethnic, alphaDiv) %>% 
group_by(ethnic,  pregnant_derived)	%>% 
summarise(value=list(alphaDiv))  %>% 
spread(pregnant_derived, value) %>%
      group_by(ethnic) %>%
      mutate(p_value	=t.test(unlist(No), unlist(Yes))$p.value, 
      test_stat=t.test(unlist(No), unlist(Yes))$statistic) %>%
      ungroup() %>%
      mutate(q_value=p.adjust(p_value, method="fdr")) %>% 
      dplyr::select(ethnic,p_value, test_stat, q_value) %>% 
      kable()

#---------------------------------------------------------------------------------------------------------
#ethnic	p_value	test_stat	q_value
#african american	0.0072612	2.7450302	0.0108917
#caucasian	0.0008737	3.4243938	0.0026211
#hispanic or latino	0.9688592	-0.0418078	0.9688592
#---------------------------------------------------------------------------------------------------------

mydata %>%
      dplyr::filter(!is.na(income)) %>%
      dplyr::filter(!(income  %in%  c("20k  to  40k",  "15k  to  20k",  "less  than  15k")))  %>% 
      group_by(ethnic) %>%
      summarise(n=  n())



## # A tibble: 3 x 2
##	ethnic	n
##	<chr>	<int>
##  1 african american	187
## 2 caucasian	467
##  3  hispanic or latino	11

library(ggplot2) 
alphaEthBoxData <- mydata %>%
	select(pregnant_derived, ethnic, alphaDiv) %>% 
melt()

## Using pregnant_derived, ethnic as id variables

alphaEthBoxDataCopy <- alphaEthBoxData alphaEthBoxDataCopy$ethnic <- "All"
alphaEthBoxData <- rbind(alphaEthBoxData, alphaEthBoxDataCopy) 
alphaEthBoxData$ethnic <- factor(alphaEthBoxData$ethnic,
		levels=c("african  american", "caucasian", "hispanic or latino", "All"),
		labels=c("AA",  "EA",  "SA",  "All"))
alphaEthBoxPlot <- ggplot(alphaEthBoxData,
		aes(x=ethnic, y=value), 
		group=pregnant_derived)  +
		geom_boxplot(aes(fill=pregnant_derived))  +
		ylab("Alpha   Diversity")   +
		xlab("Race/Ethnicity") 
print(alphaEthBoxPlot)

#Pregnant vs Non-Pregnant
#Vagitype Distribution
#Check if vagitypes are different between pregnant and non-pregnant. FDR q-value reported for each vagitype

checkTypes  <-  c("crispatus",  "iners",  "jensenii",  "gasseri",  "vaginalis",  "BVAB1") mychisq <- function(type, mydata, mycol) {
       mytable <- table(mydata[,c(mycol)], grepl(type, mydata$vagitypes)) mychisq <- chisq.test(mytable)
       mychisq$p.value
}
myfisher <- function(type, mydata, mycol) {
       mytable <- table(mydata[,c(mycol)], grepl(type, mydata$vagitypes)) myfisher <- fisher.test(mytable)
       myfisher$p.value
}
mypvalues <- sapply(checkTypes, myfisher, mydata, "pregnant_derived") myqvalues <- p.adjust(mypvalues, method="fdr")
kable(myqvalues,  col.names=NA)

#---------------------------------------
# crispatus	0.7728340
# iners 		0.0000000
#jensenii		0.3551782
#gasseri		0.7228959
#vaginalis	0.0901729
#BVAB1		0.3551782
#--------------------------------------

# Pregnant versus non pregnant AA

myttest <- t.test(alphaDiv ~ pregnant_derived, data=mydata[mydata$ethnic=="african american",]) myttest

##
##	Welch Two Sample t-test
##
## data:	alphaDiv by pregnant_derived
## t = 5.616, df = 698.79, p-value = 2.822e-08
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##	0.3226253 0.6694615
## sample estimates:
##	mean in group No mean in group Yes
##	2.394169	1.898126

# Vagitype Distribution

mypvalues <- sapply(checkTypes,
       myfisher,
       mydata[mydata$ethnic == "african american",], "pregnant_derived")
myqvalues <- p.adjust(mypvalues, method="fdr")
kable(myqvalues,  col.names=NA)

#--------------------------------------------
crispatus	0.8526456
iners	0.0000127
jensenii	0.8526456
gasseri	0.8526456
vaginalis	0.0071042
BVAB1	0.2812383
#-------------------------------------------


# Lacto vs. non-Lacto
mydataAA <- mydata[mydata$ethnic == "african	american",]
lactoNonLacto  <-  rep("NonLacto",  nrow(mydataAA)) 
lactoNonLacto[grep("Lactobacillus",  mydataAA$vagitypes)]  <-  "Lacto" 
mytable <- table(mydataAA$pregnant_derived, lactoNonLacto) 
kable(mytable)

#--------------------------------------------------
#	Lacto	NonLacto
#No	606	578
#Yes	200	104

fisher.test(mytable)

##
##Fisher's Exact Test for Count Data
##
## data:	mytable
## p-value = 5.83e-06
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##	0.4147587 0.7143061
## sample estimates:
##  odds ratio
##	0.5454155

# Pregnant versus non pregnant EA

myttest  <-  t.test(alphaDiv  ~  pregnant_derived,  data=mydata[mydata$ethnic=="caucasian",]) myttest

##
##Welch Two Sample t-test
##
## data:	alphaDiv by pregnant_derived
## t = 3.2696, df = 250.2, p-value = 0.001228
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##	0.1548863 0.6241507
## sample estimates:
##	mean in group No mean in group Yes
##	1.968356	1.578837

#Vagitype Distribution

mypvalues <- sapply(checkTypes,
       myfisher,
       mydata[mydata$ethnic == "caucasian",], "pregnant_derived")
myqvalues <- p.adjust(mypvalues, method="fdr")
kable(myqvalues,  col.names=NA)

#-----------------------------------------------------------
#crispatus	0.8918688
#iners		0.6403198
#jensenii		0.8918688
#gasseri		1.0000000
#vaginalis	0.7541106
#BVAB1		0.0678166
#-----------------------------------------------------------

# Lacto vs. non-Lacto

mydataEA <- mydata[mydata$ethnic == "caucasian",] 
lactoNonLacto  <-  rep("NonLacto",  nrow(mydataEA))
lactoNonLacto[grep("Lactobacillus",  mydataEA$vagitypes)]  <-  "Lacto" 
mytable <- table(mydataEA$pregnant_derived, lactoNonLacto) 
kable(mytable)

#--------------------------------------------
#Lacto	NonLacto
#No		479	203
#Yes		80	31
#-------------------------------------------

fisher.test(mytable)

##
##Fisher's Exact Test for Count Data
##
## data:	mytable
## p-value = 0.7374
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##	0.5649423 1.4514107
## sample estimates:
##  odds ratio
##	0.9144404

# Pregnant versus non pregnant H.

myttest <- t.test(alphaDiv ~ pregnant_derived,
data=mydata[mydata$ethnic=="hispanic or latino",])
myttest

##
##Welch Two Sample t-test
##
## data:	alphaDiv by pregnant_derived
## t = 2.011, df = 165.34, p-value = 0.04595
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##	0.004302883 0.469237591
## sample estimates:
##	mean in group No mean in group Yes
##	1.773524	1.536754

#Vagitype Distribution
mypvalues <- sapply(checkTypes,
       myfisher,
       mydata[mydata$ethnic == "hispanic or latino",], "pregnant_derived")
myqvalues <- p.adjust(mypvalues, method="fdr")
kable(myqvalues,  col.names=NA)


#-----------------------------------------------------
#crispatus	0.0537774
#iners		0.2147489
#jensenii		0.0537774
#gasseri		0.2147489
#vaginalis	0.4145109
#BVAB1		0.2147489
#-----------------------------------------------------

#Lacto vs. non-Lacto

mydataSA <- mydata[mydata$ethnic == "hispanic or latino",] 
lactoNonLacto  <-  rep("NonLacto",  nrow(mydataSA)) 
lactoNonLacto[grep("Lactobacillus",  mydataSA$vagitypes)]  <-  "Lacto" 
mytable <- table(mydataSA$pregnant_derived, lactoNonLacto) 
kable(mytable)

#---------------------------------------------------
#Lacto	NonLacto
#No	63	40
#Yes	155	43
#--------------------------------------------------

fisher.test(mytable)

##
##	Fisher's Exact Test for Count Data
##
## data:	mytable
## p-value = 0.002632
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##	0.2514720 0.7621022
## sample estimates:
##  odds ratio
##	0.4382236



