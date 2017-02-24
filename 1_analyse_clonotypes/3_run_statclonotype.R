
#
# IMGT/StatClonotype uses a generic statistical procedure [1] for identifying
# significant changes in IG and TR differences of proportions of IMGT clonotypes
# (AA) diversity and expression [5].
#
library(IMGTStatClonotype)

data(MID1)
data(MID2)

MID1<-clonRem(MID1)
MID2<-clonRem(MID2)

# ...the IMGTstandardized approach allows a clear distinction between 
# the clonotype diversity (numbers of IMGTclonotypes (AA) perV,Dor J gene), 
# and the clonotype expression (numbers of sequences assigned, unambiguously, to a given IMGT clonotype (AA) perV,Dor J gene) [16]...

# Numbers of IMGT clonotypes (AA) in the two compared sets from the
# IMGT/HighV-QUEST output for clonotype diversity
Ndiv<-clonNumDiv(MID1,MID2)
# Numbers of IMGT clonotypes (AA) in the two compared sets from the
# IMGT/HighV-QUEST output for clonotype expression
Nexp<-clonNumExp(MID1,MID2)

# Significance of the difference in proportions with 95% confidence in-
# terval (CI) for IMGT clonotype (AA) diversity between two sets from
# IMGT/HighV-QUEST output
div<-sigrepDiv(Ndiv,MID1,MID2)
# Significance of the difference in proportions with 95% confidence in-
# terval (CI) for IMGT clonotype (AA) expression between two sets from
# IMGT/HighV-QUEST output
exp<-sigrepExp(Nexp,MID1,MID2)

diffpropGph(div)$Vgenes
diffpropGph(exp)$Vgenes

