the.data<-as.matrix(read.table("ENB18data.txt"))
View(the.data)
my.data<-the.data[sample(1:768,300),c(1:5,7)]
View(my.data)

plot(my.data[,1],my.data[,6],
     main=expression("Scatter Plots for X1 vs. Y2"),
     xlab="Relative compactness(X1)",
     ylab="Cooling Load(Y2)",
     col="red")


plot(my.data[,2],my.data[,6],
     main=expression("Scatter Plots for X2 vs. Y2"),
     xlab="Surface area(X2)",
     ylab="Cooling Load(Y2)",
     col="blue")

plot(my.data[,3],my.data[,6],
     main=expression("Scatter Plots for X3 vs. Y2"),
     xlab="Wall area(X3)",
     ylab="Cooling Load(Y2)",
     col="green")

plot(my.data[,4],my.data[,6],
     main=expression("Scatter Plots for X4 vs. Y2"),
     xlab="Roof area(X4)",
     ylab="Cooling Load(Y2)",
     col="yellow")

plot(my.data[,5],my.data[,6],
     main=expression("Scatter Plots for X5 vs. Y2"),
     xlab="Overall height(X5)",
     ylab="Cooling Load(Y2)",
     col="orange")

hist(my.data[,1],
     main=expression("Histogram for Relative compactness(X1)"),
     xlab="Relative compactness(X1)",
     ylab="Frequency of Relative compactness(X1)",
     col="red")

hist(my.data[,2],
     main=expression("Histogram for Surface area(X2)"),
     xlab="Surface area(X2)",
     ylab="Frequency of Surface area(X2)",
     col="blue")

hist(my.data[,3],
     main=expression("Histogram for Wall area(X3)"),
     xlab="Wall area(X3)",
     ylab="Frequency of Wall area(X3)",
     col="green" )

hist(my.data[,4],
     main=expression("Histogram for Roof area(X4)"),
     xlab="Roof area(X4)",
     ylab="Frequency of Roof area(X4)",
     col="yellow" )

hist(my.data[,5],
     main=expression("Histogram for Overall height(X5)"),
     xlab="Overall height(X5)",
     ylab="Frequency of Overall height(X5)",
     col="orange" )

hist(my.data[,6],
     main=expression("Histogram for Cooling Load(Y2)"),
     xlab="Cooling Load(Y2)",
     ylab="Frequency of Cooling Load(Y2)",
     col="brown")

#Copy the original data
original<-my.data
View(original)

#Linear Feature Scaling
my.data[,6]<-(my.data[,6]-min(my.data[,6]))/(max(my.data[,6])-min(my.data[,6]))
View(my.data[,6])
hist(my.data[,6])     

#Negation and Transformation
my.data[,5]<-(7-my.data[,5])/(max(my.data[,5])-min(my.data[,5]))
View(my.data[,5])
hist(my.data[,5])

#Log Transformation and Scaling
my.data[,4]<-log(my.data[,4])
my.data[,4]<-1-(my.data[,4]-min(my.data[,4]))/(max(my.data[,4])-min(my.data[,4]))
View(my.data[,4])

#Standardization
my.data[,3]<-(my.data[,3]-mean(my.data[,3]))/sd(my.data[,3])
View(my.data[,3])

#Polynomial Transformation and Scaling
my.data[,2]<-(my.data[,2])^(1/2)
my.data[,2]<-1-(my.data[,2]-min(my.data[,2]))/(max(my.data[,2])-min(my.data[,2]))
View(my.data[,2]) 

your.data<-cbind(my.data[,2],my.data[,3],my.data[,4],my.data[,5],my.data[,6])
View(your.data)

write.table(your.data,"Shantanu-transformed.txt")

source("AggWaFit7181.R")
View(my.data)

#WAM(p=1)
fit.QAM(my.data,"out_AM.txt","stat_AM.txt",g=AM,g.inv=invAM)
#WPM(p=2)
fit.QAM(my.data,"out_QM.txt","stat_QM.txt",g=QM,g.inv=invQM)
#OWA
fit.OWA(my.data,"out_OWA.txt","stat_OWA.txt")
#Choquet Integral
to.fit2<-my.data[,c(1,2,4,5,6)]
fit.choquet(to.fit2,"out_choq.txt","stat_choq.txt")
#WPM(p=0.5)
fit.QAM(to.fit2,"out_PM05.txt","stat_PM05.txt",g=PM05,g.inv=invPM05)
View(my.data)

#Prediction using Weighted Airthmetic mean
QAM(c(0.82,0.79,0.072,0.66,0),c(0,0.792,0.1630,0,0.044))
QAM
#Prediction using Weighted Power mean(p=2)
QAM(c(0.82,0.79,0.072,0.66,0),c(0,0.0873,0,0.8941,0.0184),QM,invQM)
QAM
#Prediction using Weighted Power mean(p=0.5)
QAM(c(0.82,0.79,0.072,0.66,0),c(0.2838,0,0,0.5939,0.1221),PM05,invPM05)
QAM
#Prediction using Ordered Weight Average
OWA(c(0.82,0.79,0.072,0.66,0),c(0.070,0.3594,0.5464,0,0.02376))
OWA
#Prediction using Choquet integral
choquet(c(0.66,0.82),c(0,0.189,0.774,1))
choquet
