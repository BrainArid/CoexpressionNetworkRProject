plot2Groups <- function(GroupA, GroupB, byRow=TRUE, main="Mean Plot", xlab="Group A", ylab="Group B", file=NA, histA=FALSE, histB=FALSE,breaksA=100,breaksB=100)
{
  #if(!is.na(file))
    png(filename=file);
  
  a <- hist(GroupB,breaks=breaksB)
  
  if(histA && histB)
  {
    layout(matrix(c(3,1,1,1,3,1,1,1,3,1,1,1,4,2,2,2), 4, 4, byrow = TRUE));
  }
  else if(histA)
  {
    layout(matrix(c(1,1,1,2), 4, 1, byrow = TRUE));
  }
  else if(histB)
  {
    a <- hist(GroupB,breaks=breaksB)
    layout(matrix(c(2,1,1,1), 1, 4, byrow = FALSE))
  }

  
  plot(x=GroupA, y=GroupB, main=main,xlab=xlab,ylab=ylab, type="n");
  points(x=GroupA, y=GroupB, pch=20, cex=0.5);
  
  
  if(histA)
  {
    hist(GroupA,breaks=breaksA);
  }
  if(histB)
  {
    ##hist(GroupB);
    
    barplot(height=a$counts,width=1, space=0, horiz=TRUE)
    width <- a$breaks[2] - a$breaks[1]
    axis(2, at=(pretty(a$breaks) - a$breaks[1])/width,labels=pretty(a$breaks))
  }
    if(histA && histB)
    {
      dev.off()
    }
    
  #if(!is.na(file))
    dev.off();
}