n=4

pin=runif(n,min=0,max=1)
pout=runif(n,min=0,max=1)

link=matrix(nrow = 4,ncol=4)
link[1,]=c(0,0,1,0)
link[2,]=c(0,0,1,0)
link[3,]=c(0,0,0,1)
link[4,]=c(0,1,1,0)

x=pin
y=pout

convergence(link,x,y)
weights=pageweights(x,y)
f=0.5
pcm=constructpcm(link,f)
library("e1071")
cm=corr_matrix(link,pcm,weights,f)
cl=clustering(cm,4,weights,0.5)
