#Estructura del bosque
setwd("C:/Users/David Esteban/Google Drive/Materias cursadas/Ecología Forestal II/Indices")

#Parcelas
p16<-read.csv2("parcela16.csv")
p17<-read.csv2("parcela17.csv")
p18<-read.csv2("parcela18.csv")

#Alturas de copa
p16$HC<-p16$H-p16$Hcom
p17$HC<-p17$H-p17$Hcom
p18$HC<-p18$H-p18$Hcom

#Estructura diametrica


with(p16, plot(D, H, xlab="Diametro (cm)", ylab="Altura (m)", pch=0))
with(p17, plot(D, H, xlab="Diametro (cm)", ylab="Altura (m)", pch=1))
with(p18, plot(D, H, xlab="Diametro (cm)", ylab="Altura (m)", pch=2))

## Función estimar n de Liocur-Meyer III N=k*exp(-rD^n)
estimar.n<-function(N, D,datos,x,y,z){
        n=seq(x,y,z)
        CME<-function(mod){
                cme= CME=anova(mod)$Mean[length(anova(mod)$Mean)]
        }
        cme=vector()
        for(i in 1:length(n)){
                mod<-lm(log(N)~I(D^n[i]),datos)
                cme[i]=CME(mod)
        } 
        plot(n,cme,type="l", xlab="m",ylab="CME del modelo lineal")
        
        data=data.frame(cme,n)
        print(data[which.min(data$cme),])
        
}

#Distribución diamétrica de cada parcela

h16<-with(p16, hist(D, xlab="Diametro",ylab="",breaks = 2, main="", las=2))
tab16<-data.frame(h16$counts,h16$mids)
t16<-tab16[-c(6:8),]
colnames(t16)<-c("N","xi")

h17<-with(p17, hist(D, xlab="Diametro",ylab="Frecuencia", main="", las=2))
tab17<-data.frame(h17$counts,h17$mids)
t17<-tab17[-c(3,5),]
colnames(t17)<-c("N","xi")

h18<-with(p18, hist(D, xlab="Diametro",ylab="Frecuencia", main="", las=2))
t18<-data.frame(h18$counts,h18$mids)
colnames(t18)<-c("N","xi")

            ## MODELAMIENTO DE LAS CLASES DIAMETRICAS

#Parcela 16
mod16<-lm(log(N)~I(xi^2),t16)
mod16<-lm(log(N)~I(xi^-2.318), t16) #mejor modelo
mod16<-lm(log(N)~xi, t16)

summary(mod16)

#Grafica del modelo y las distribuciones
h16<-with(p16, hist(D, xlab="Diametro",ylab="", main="", las=2, xlim=c(10,100)))
par(new=TRUE)
curve(exp(0.8785+(0.2221/2))*exp((127.2*x^-2.318)),from=10,to=100,xlim=c(10,100), axes=F, col="Red", ylab="",xlab="")


#Parcela 17
with(t17, estimar.n(N, xi,t17,10.0,10.3,0.001))

mod17<-lm(log(N)~I(xi^2),t17)
mod17<-lm(log(N)~I(xi^10.098), t17) #Mejor modelo
mod17<-lm(log(N)~xi, t17)

summary(mod17)
h17<-with(p17, hist(D, xlab="Diametro",ylab="Frecuencia", main="", las=2))
par(new=T)
curve(exp(2.209e+00+(0.10967/2))*exp((-1.053e-11*x^10.098)),from=10,to=13, axes=F, col="Red", ylab="",xlab="")

#Parcela 18
with(t18, estimar.n(N, xi,t18,-280,-100,1))

mod18<-lm(log(N)~I(xi^2),t18)
mod18<-lm(log(N)~I(xi^-100), t18)
mod18<-lm(log(N)~xi, t18)

summary(mod18)

#Teniendo en cuenta que ningun modelo se le ajusto, se opta por partir la base
#de datos en dos 

t18a<-t18[c(1:4),]
t18b<-t18[c(5:7),]
t18c<-t18[c(4:7),]

with(t18a, estimar.n(N, xi,t18a,-10,1,1))

mod18a<-lm(log(N)~I(xi^2),t18a)
mod18a<-lm(log(N)~I(xi^-8), t18a)  #Mejor ajuste R2
mod18a<-lm(log(N)~xi, t18a)

summary(mod18a)

with(t18b, estimar.n(N, xi,t18b,11,12,0.01))

mod18b<-lm(log(N)~I(xi^2),t18b)
mod18b<-lm(log(N)~I(xi^11.47), t18b) #Mejor modelo
mod18b<-lm(log(N)~xi, t18b)

summary(mod18b)

mod18c<-lm(log(N)~log(xi)+xi, t18c)

summary(mod18c)

#Grafica de la paracela 18

h18<-with(p18, hist(D, xlab="Diametro",ylab="Frecuencia", main="", las=2, xlim=c(10,45),ylim=c(0,40)))
par(new=T)
curve(exp(6.149e-01+(0.0125/2))*exp(1.842e+09*x^-8),from=10,to=30,xlim=c(10,45),ylim=c(0,40), axes=F, col="Red", ylab="",xlab="")
par(new=T)
curve(exp(-218.0351+(0.0241/2))*exp(-2.5293*x)*x^86.9937, from=20, to=45, xlim=c(10,45), ylim=c(0,40), axes=F, ylab="",xlab="",col="blue")


#Clases altimetricas 
a16<-with(p16, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2))
tab16<-data.frame(a16$counts,a16$mids)
colnames(tab16)<-c("N","h")
with(tab16, sum(N))
tab16$Por<-100*tab16$N/56
tab16

a17<-with(p17, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2))
tab17<-with(a17, data.frame(counts,mids))
tab17<-tab17[-c(5,6),]
colnames(tab17)<-c("N","h")
with(tab17, sum(N))
tab17$Por<-100*tab17$N/22
tab17

a18<-with(p18, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2))
tab18<-with(a18, data.frame(counts,mids))
colnames(tab18)<-c("N","h")
with(tab18, sum(N))
tab18$Por<-100*tab18$N/74
tab18
 
t16$por<-100*t16$N/sum(t16$N)
t16           
t17$por<-100*t17$N/sum(t17$N)
t17
t18$por<-100*t18$N/sum(t18$N)
t18


                  ## MODELAMIENTO DE LAS CLASES ALIMETRICAS 

#Parcela 16
with(tab16, estimar.n(N, h,tab16,0.9,1.0,0.001))
m16<-lm(log(N)~I(h^2),tab16)
m16<-lm(log(N)~I(h^0.979), tab16) #mejor modelo
m16<-lm(log(N)~h, tab16)

summary(m16)

#Grafica del modelo y las distribuciones
a16<-with(p16, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2, xlim=c(5,30)))
par(new=TRUE)
curve(exp(3.71461+(0.10211/2))*exp((-0.09046*x^0.979)),from=5,to=30,xlim=c(5,30), axes=F, col="Red", ylab="",xlab="")

#Parcela 17
with(tab17, estimar.n(N, h,tab17,0.6,0.7,0.001))
m17<-lm(log(N)~I(h^2),tab17)
m17<-lm(log(N)~I(h^0.666), tab17) #mejor modelo
m17<-lm(log(N)~h, tab17)

summary(m17)

#Grafica del modelo y las distribuciones
a17<-with(p17, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2,  xlim=c(4,11)))
par(new=TRUE)
curve(exp(4.9124+(0.19223/2))*exp((-1.0417*x^0.666)),from=4,to=11,xlim=c(4,11), axes=F, col="Red", ylab="",xlab="")

#Parcela 18
tab18
a18<-with(p18, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2, xlim=c(2,24)))

with(tab18, estimar.n(N, h,tab18,-10,100,1))
m18<-lm(log(N)~I(h^2),tab18)
m18<-lm(log(N)~I(h^0), tab18) #mejor modelo
m18<-lm(log(N)~h, tab18)

summary(mod18)
#Diametricas
{
        par(mfrow=c(1,3))
        h16<-with(p16, hist(D, xlab="Diametro",ylab="", main="", las=2, xlim=c(10,100)))
        par(new=TRUE)
        curve(exp(0.8785+(0.2221/2))*exp((127.2*x^-2.318)),from=10,to=100,
              xlim=c(10,100), axes=F, col="Red", ylab="",xlab="")
        h17<-with(p17, hist(D, xlab="Diametro",ylab="Frecuencia", main="", las=2))
        par(new=T)
        curve(exp(2.209e+00+(0.10967/2))*exp((-1.053e-11*x^10.098)),from=10,
              to=13, axes=F, col="Red", ylab="",xlab="")
        h18<-with(p18, hist(D, xlab="Diametro",ylab="Frecuencia", main="", 
                            las=2, xlim=c(10,45)))
        par(new=T)
        curve(exp(6.149e-01+(0.0125/2))*exp(1.842e+09*x^-8),from=10,to=30,
              xlim=c(10,45),ylim=c(0,40), axes=F, col="Red", ylab="",xlab="")
        par(new=T)
        curve(exp(-218.0351+(0.0241/2))*exp(-2.5293*x)*x^86.9937, from=20, to=45, 
              xlim=c(10,45), ylim=c(0,40), axes=F, ylab="",xlab="",col="blue")
}

#Alturas
{
        par(mfrow=c(1,3))
        a16<-with(p16, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2, xlim=c(5,30)))
        a17<-with(p17, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2,  xlim=c(4,11)))
        a18<-with(p18, hist(H, xlab="Altura",ylab="Frecuencia", main="", las=2, xlim=c(2,24)))
}

