#1st col of data for group; 2nd col of data for value.
ggcdf=function(data){
  group=unique(data[,1])
  i=1
  n=length(group)
  p=ggplot()
  while(i<=n){
    eval(parse(text=paste(sep="","
      temp=data[data[,1]=='",group[i],"',2]
      temp=as.numeric(temp)
      temp=temp[order(temp)]
      len=length(temp)[1]
      iid<-c()
      for(j in 1:len){iid[j]=j/len}
      ",gsub("-","_",group[i]),"<-data.frame(temp,iid) 
      p=p+geom_line(data=",gsub("-","_",group[i]),",aes(x=",gsub("-","_",group[i]),"[,1],y=",
gsub("-","_",group[i]),"[,2],color='",gsub("-","_",group[i]),"'),size=1.5)
        
                          ")))
    print(i)
    i=i+1
  }
  return(p)
}

