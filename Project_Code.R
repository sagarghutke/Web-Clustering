convergence<-function(link,pagein,pageout)
{
  pin<-pagein
  pout<-pageout
  for(i in 1:length(pin))
  {
    pin_sum=0
    pout_sum=0;
    for(j in 1:length(pin))
    {
      if(link[j,i]==1) pin_sum =pin_sum+ pout[j]
      if(link[i,j]==1) pout_sum =pout_sum+ pin[j]
    }
    pin[i]=pin_sum
    pout[i]=pout_sum
    normalise(pin)
    normalise(pout)
  }
  eval.parent(substitute(pagein<-pin))
  eval.parent(substitute(pageout<-pout))
}

normalise<-function(vector)
{
  val=0;
  v<-vector
  for(i in v)
    val=val+i*i;
  val=sqrt(val)
  
  for(i in 1:length(v))
  {
    v[i]=v[i]/val;
  }
  eval.parent(substitute(vector<-v))
}

average<-function(x)
{
 sum(x)/length(x)
}

avgin=average(x)
avgout=average(y)

pageweights<-function(pin,pout)
{
  w=vector(length=length(pin))
  avgin=average(x)
  avgout=average(y)
  maxin=max(pin)
  minin=min(pin)
  maxout=max(pout)
  minout=min(pout)
  for(i in 1:length(pin))
  {
    w[i]=1+max((pin[i]-avgin)/(maxin-minin),(pout[i]-avgout)/(maxout-minout))
  }
  w
}

constructpcm<-function(link,f)
{
  l=sqrt(length(link))
  pcm=matrix(nrow = l,ncol = l)
  for(i in 1:l)
  {
    for(j in 1:l)
    {
      if(i==j) pcm[i,j]=1
      else if(link[i,j]==1) pcm[i,j]=f
      else pcm[i,j]=0
    }
  }
  pcm
}

update_pm<-function(f,pcm)
{
  factor=0
  k=0
  for(i in 1:nrow(pcm))
  {
    factor=f
      for(j in 1:ncol(pcm))
      {
        if(pcm[i,j]==f)
        {
          k=1
          temp=0
          while(k<=ncol(pcm))
          {
            if(pcm[j,k]!=0&&pcm[j,k]!=1)
            { temp=factor*pcm[j,k];
            if(temp>pcm[i,k])
              pcm[i,k]=temp
            }
            k=k+1;
          }
          factor=factor*f
        }
       
      }
  }
  pcm
}

calc_sl<-function(pcm,i,j,f)
{
  log10(pcm[i,j])/log10(f)
}

dist_mat<-function(link)
{
  for(i in 1:nrow(link))
  {
    for(j in 1:ncol(link))
    {
      if(i==j)link[i,j]=0
      else if(link[i,j]!=1)link[i,j]=Inf
    }
  }
  link
}

corr_weight<-function(w1,w2)
{
  max(w1,w2)
}

corr_matrix<-function(link,pcm,weights,f)
{
  cm=matrix(nrow=nrow(link),ncol = ncol(link))
  distmat=dist_mat(link)
  distances=allShortestPaths(distmat)
  pcm_up=update_pm(f,pcm)
  for( i in 1:nrow(link))
  {
    for(j in 1:ncol(link))
    {
      vector=extractPath(distances,i,j)
      prod=1
      for(k in 1:length(vector)-1)
      {
        if(i==j) break
        prod=prod*corr_weight(weights[vector[k]],weights[vector[k+1]])
      }
      if(i==j) cm[i,j]=1
      else
      {
      sl=calc_sl(pcm_up,i,j,f)
      cm[i,j]=prod*(f^sl)
      }
    }
  }
  cm
}


modrow<-function(cm,i)
{
  prod=1
  for(j in 1:ncol(cm))
  {
    prod=prod+(cm[i,j]^2)
  }
  
  sqrt(prod)
}

modcol<-function(cm,i)
{
  prod=1
  for(j in 1:nrow(cm))
  {
    prod=prod+(cm[j,i]^2)
  }
  
  sqrt(prod)
}

rownumerator<-function(cm,i,j)
{
  prod=1
  for(k in 1:ncol(cm))
  {
    prod=prod+cm[i,k]*cm[j,k]
  }
  prod
}

colnumerator<-function(cm,i,j)
{
  prod=1
  for(k in 1:nrow(cm))
  {
    prod=prod+cm[k,i]*cm[k,j]
  }
  prod
}

simin<-function(cm,i,j)
{
  denom=modcol(cm,i)*modcol(cm,j)
  num=colnumerator(cm,i,j)
  num/denom
}

simout<-function(cm,i,j)
{
  denom=modrow(cm,i)*modrow(cm,j)
  num=rownumerator(cm,i,j)
  num/denom
}

outweight<-function(cm,i,j)
{
  num=modrow(cm,i)+modrow(cm,j)
  denom=modrow(cm,i)+modrow(cm,j)+modcol(cm,i)+modcol(cm,j)
  num/denom
}

inweight<-function(cm,i,j)
{
  num=modcol(cm,i)+modcol(cm,j)
  denom=modrow(cm,i)+modrow(cm,j)+modcol(cm,i)+modcol(cm,j)
  num/denom
}

sim<-function(cm,i,j)
{
  ins=inweight(cm,i,j)*simin(cm,i,j)
  outs=outweight(cm,i,j)*simout(cm,i,j)
  ins+outs
}

clustering<-function(cm,l,weights,tresh)
{
  t=tresh
  rownum=1
  clusters=list()
  centroid_val=vector()
  centroid=vector()
  cluster_size=vector()
  clusters[[rownum]]=1
  centroid_val[rownum]=weights[1]
  centroid[rownum]=1
  cluster_size[rownum]=1
  max=0;
  cluster_number=0;
  for(i in 2:l)
  {
    for(j in 1:rownum)
    {
      similarity=sim(cm,i,centroid[j])
      similarity
      if(similarity>t&&similarity>max)
      {
        max=similarity
        cluster_number=j
      }
    }
    
    if(cluster_number!=0)
    {
      print(cluster_number)
      flush.console()
      cluster_size[cluster_number]=cluster_size[cluster_number]+1
      clusters[[cluster_number]][cluster_size[cluster_number]]=i
      centroid_val[cluster_number]=calc_centroid(cluster_number,cluster_size[cluster_number],cm,clusters)
      centroid[cluster_number]=calc_cluster(cluster_number,cluster_size[cluster_number],cm,clusters)
      }
    else
    {
      print("new",cluster_number)
      flush.console()
      rownum=rownum+1
      clusters[[length(clusters)+1]]=c(i)
      cluster_size[rownum]=1
      centroid_val[rownum]=weights[i]
      centroid[rownum]=i
    }
  }
  print(centroid)
  clusters
}

calc_centroid<-function(cnum,csize,cm,clusters)
{
  cerow=0
  c=0
  for(i in 1:csize)
  {
    c=clusters[[cnum]][i]
    c=unlist(c)
    flush.console()
   cerow=cerow+modrow(cm,c)
  }
  cerow=cerow/csize
  
  cecol=0
  c=0
  for(i in 1:csize)
  {
    c=clusters[[cnum]][i]
    c=unlist(c)
    cecol=cecol+modcol(cm,c)
  }
  
  cecol=cecol/csize
  
  sqrt(cerow^2+cecol^2)
}

calc_cluster<-function(cnum,csize,cm,clusters)
{
  c1=0
  c2=0
  avgsim=0
  max=0
  ce=0
  for(i in 1:csize)
  {
    avgsim=0
    c1=clusters[[cnum]][i]
    c1=unlist(c1)
    for(j in 1:csize)
    {
      if(i==j) next
      c2=clusters[[cnum]][j]
      c2=unlist(c2)
      avgsim=avgsim+sim(cm,c1,c2)
    }
    
    avgsim=avgsim/(csize-1)
    if(avgsim>max)
    {
      max=avgsim
      ce=i
    }
    
  }
  ce
}

input.array<-read.csv(file="~/edges.csv",head=TRUE, sep=",")
adjMatrix.array<-array(0:0, dim=c(5000,5000))
tmax=5000
for(i in 1:tmax){
  if(input.array[i,1]<=5000 && input.array[i,2]<=5000)
    adjMatrix.array[input.array[i,1],input.array[i,2]]=1;
}




