pie_plot=function(temp.data){
  temp.data=data.frame(table(temp.data))
  colnames(temp.data)=c('type','nums')
  
  temp.data$nums=as.numeric(temp.data$nums)
  label_value <- paste('(', round(temp.data$nums/sum(temp.data$nums) * 100, 1), '% n:',temp.data$nums,')', sep = '')  
  label <- paste(temp.data$type, label_value, sep = '')
  temp.data$type=label
  temp.data=temp.data[order(temp.data$nums,decreasing = T),]
  order.temp=temp.data$type
  order.temp=order.temp[grep("Unknown",order.temp,invert = T)]
  temp.data$type=ordered(temp.data$type,levels=c(order.temp,temp.data$type[grep("Unknown",temp.data$type)]))
  p <- ggplot(data = temp.data,mapping = aes(x = 'Content', y = nums, fill = type)) + 
    geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme_void()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),legend.title = element_blank())  
  return(p)
}

library(ggforce)

pie_plot2=function(temp.data){
  temp.data=data.frame(table(temp.data))
  colnames(temp.data)=c('x','y')
  
  temp.data$y=as.numeric(temp.data$y)
  label_value <- paste('(', round(temp.data$y/sum(temp.data$y) * 100, 1), '% n:',temp.data$y,')', sep = '')  
  label <- paste(temp.data$x, label_value, sep = '')
  temp.data$x=label
  temp.data=temp.data[order(temp.data$y,decreasing = T),]
  order.temp=temp.data$x
  order.temp=order.temp[grep("Unknown",order.temp,invert = T)]
  temp.data$x=ordered(temp.data$x,levels=c(order.temp,temp.data$x[grep("Unknown",temp.data$x)]))
  p <- ggplot(data = temp.data) + 
    geom_arc_bar(aes(x0=0, y0=0, r0=0, r=1, color=NULL,
                     amount=y, fill=x, explode=.1), 
                 stat='pie') +  
    coord_equal() + labs(x = '', y = '', title = '') + theme_void()+
    theme(axis.text = element_blank(),axis.ticks = element_blank(),legend.title = element_blank())  
  return(p)
}
