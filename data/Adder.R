


#Adder


M <- matrix(rnorm(16),4,4)
row.names(M) <- c("High","Mid1", "Mid2","Low")
colnames(M) <- row.names(M)

M


Parts <- c("Mid1", "Mid2","Low")
OutName <- "Middle"

Adder <- function(M,Parts,OutName) {
  
  ColSum <- apply(  M[,colnames(M) %in% Parts ], 1, sum  ) #Adds the relevant columns
  Mnew <- cbind(M,ColSum)
  colnames(Mnew) <- c( colnames(Mnew)[-ncol(Mnew)], OutName) # Replace name with the one provided
  
  RowSum <- apply(  Mnew[ row.names(Mnew) %in% Parts, ], 2, sum  ) #Adds the relevant rows
  Mnew <- rbind(Mnew,RowSum)
  row.names(Mnew) <- c( row.names(Mnew)[-nrow(Mnew)], OutName) # Replace name with the one provided 

  return( Mnew[!(row.names(Mnew) %in% Parts), !(colnames(Mnew) %in% Parts)] ) #return the matrix, minus the promted rows
}

Adder(M,Parts,OutName)

library("XLConnect")


setwd("C:/Users/ps583/Documents/GitHub/TRICE/data")

df1 <- readNamedRegionFromFile(file="ModelData_v2.xlsx",name="nontrade")
df2 <- readNamedRegionFromFile(file="ModelData_v2.xlsx",name="trademat")
df3 <- readNamedRegionFromFile(file="TradeMat.xlsx",name="trademat")

Groups <- unique(df1$Group)

DF <- aggregate(df1$GDP2011,by=list(factor(df1$Group)),FUN=function(x) sum(x,na.rm=TRUE))
colnames(DF) <- c("Group","GDP.2011")

DF <- cbind( DF, "Pop.2011"=aggregate(df1$POP2011,by=list(factor(df1$Group)),FUN=function(x) sum(x,na.rm=TRUE))[,2] )
DF <- cbind( DF, "Co2.2011"=aggregate(df1$CO2_2011,by=list(factor(df1$Group)),FUN=function(x) sum(x,na.rm=TRUE))[,2] )
DF 



#More clearly define the trade matrix
temp1 <- df2[,-1]
dim(temp1)
row.names(temp1) <- colnames(temp1)
NewTradeMat <- data.matrix(temp1)
NewTradeMat[is.na(NewTradeMat)] <- 0

#MatList <- vector(mode="list",15)

for(i in Groups) {
  print(i) #for each group
  Parts <- df1[df1$Group==i,"UnStatRName"] #get the matrix labels we'll need
  print(Parts)
  
  if( length(Parts)>1 ) NewTradeMat <- Adder(NewTradeMat,Parts,i)
}

NewTradeMat

RoW <- colnames(NewTradeMat) [! (colnames(NewTradeMat) %in% Groups) ]
FinalTradeMat <- Adder(NewTradeMat,RoW,"Rest.Of.World")

#RoW Info


#Custom Sort
Target <- c("Brazil","Japan", "EU", "Sub.Saharan.Africa", "Canada", "United.States", "Latin.America", "Rest.Of.World", "Southeast.Asia", "Mideast.OIL", "Russian.Federation", "India", "South.Africa", "China", "Ex.USSR")
FinalTradeMat <- FinalTradeMat[match(Target,row.names(FinalTradeMat)),]
FinalTradeMat <- FinalTradeMat[, match(row.names(FinalTradeMat),colnames(FinalTradeMat))]
DF <- DF[match(row.names(FinalTradeMat),as.character(DF$Group)),]

write.csv(FinalTradeMat,file="ReducedTradeMat.csv")
write.csv(DF,file="ReducedNonTradeData.csv")
write.csv(RoW,file="RestOfWorld.csv")
