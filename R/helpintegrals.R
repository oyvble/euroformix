#Derived integral expressions for different hypotheses.

#####
#K=1#
#####
.calcInt1p = function(LikFun,myInt,mxLims,thetaOther,...) {
  return(LikFun(thetaOther)) #simply call the likelihood function
}

#####
#K=2#
#####

#Only unknowns
.calcInt2p0K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fx = function(x) {
    return(LikFun(c(x,thetaOther)))
  }
  return(2*myInt(int_fx,mxLims[1,1],mxLims[1,2]))
}

#1 or 2 conditionals
.calcInt2p1K = function(LikFun,myInt,mxLims,thetaOther,...) { 
  int_fx = function(x) {
    return(LikFun(c(x,thetaOther)))
  }  
  return(myInt(int_fx,mxLims[1,1],mxLims[1,2]))
}

#####
#K=3#
#####

#Only unknowns
.calcInt3p0K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fx = function(x,part) { 
    int_fzx = function(z) { 
      return(LikFun(c(x,z,thetaOther)))
    } #end z
    fromZ = mxLims[2,1]
    upperZ = ifelse(part==1,x,1-2*x)
    toZ = min(upperZ,mxLims[2,2])
    val = 0
    if(fromZ < toZ) val = myInt(int_fzx, fromZ,toZ) #integrate on z
    return(val) #integrate on z
  }
  cc = 0 #evaluation counter
  xM = 1/3 #mid point for x
  fromX = mxLims[1,1]
  toX = mxLims[1,2]
  val1 <- val2 <- 0
  if(fromX < xM) val1 = myInt(int_fx, fromX  , xM,part=1)
  if(xM < toX)    val2 = myInt(int_fx, xM ,    toX,part=2)
  return(6*(val1 + val2)) #final integrated value
}

#1 conditional
.calcInt3p1K = function(LikFun,myInt,mxLims,thetaOther,...) { 
  int_fx = function(x) { 
    int_fzx = function(z) { 
      return(LikFun(c(x,z,thetaOther)))
    } #end z
    fromZ = mxLims[2,1]
    upperZ = (1-x)/2
    toZ = min(upperZ,mxLims[2,2])
    val = 0
    if(fromZ < toZ) val = myInt(int_fzx, fromZ,toZ) #integrate on z
    return(val) #integrate on z
  }
  return(2*myInt(int_fx,mxLims[1,1],mxLims[1,2]))
}

#2 or 3 conditional
.calcInt3p2K = function(LikFun,myInt,mxLims,thetaOther,...) { 
  int_fx = function(x) { 
    int_fzx = function(z) { 
      return(LikFun(c(x,z,thetaOther)))
    } #end z
    fromZ = mxLims[2,1]
    upperZ = (1-x)
    toZ = min(upperZ,mxLims[2,2])
    val = 0
    if(fromZ < toZ) val = myInt(int_fzx, fromZ,toZ) #integrate on z
    return(val) #integrate on z
  }
  return(myInt(int_fx,mxLims[1,1],mxLims[1,2]))
}

#####
#K=4#
#####

#Only unknowns
.calcInt4p0K = function(LikFun,myInt,mxLims,thetaOther) {
  int_fy = function(y) {  
    int_fxy = function(x,part) { 
      int_fzxy = function(z) { 
        return(LikFun(c(x,z,y,thetaOther)))
      }
      
      #set limits for Z
      fromZ = max(y,mxLims[2,1])
      upperZ = ifelse(part==1,x,1-2*x-y)
      toZ = min(upperZ,mxLims[2,2])
      val = 0
      if(fromZ < toZ) val = myInt(int_fzxy, fromZ, toZ) #integrate on z
      return(val)
    } #end x
    xM = (1-y)/3 #mid point
    xU = 0.5 - y #last point
    fromX = max(y,mxLims[1,1])
    toX = min(xU,mxLims[1,2])
    val1 <- val2 <- 0
    if(xM > fromX) val1 = myInt(int_fxy, fromX  , xM, part=1)
    if(toX > xM)    val2 = myInt(int_fxy, xM ,    toX, part=2)
    return(val1 + val2)
  } 
  return(24*myInt(int_fy,mxLims[3,1],mxLims[3,2]))
}

#3 or 4 knowns
.calcInt4p3K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fx = function(x) {  
    int_fyx = function(y) { 
      int_fzyx = function(z) { 
        return(LikFun(c(x,y,z,thetaOther)))
      }
      fromZ = mxLims[3,1]
      toZ = min(1-x-y,mxLims[3,2])
      val = 0
      if(fromZ < toZ) val = myInt(int_fzyx, fromZ, toZ) #integrate on z
      return(val)
    } #end x
    fromY = mxLims[2,1]
    toY = min(1-x,mxLims[2,2])
    val1 = 0
    if(fromY < toY) val1 = myInt(int_fyx, fromY  , toY)
    return(val1)
  } 
  return(myInt(int_fx,mxLims[1,1],mxLims[1,2]))
}

#2 knowns
.calcInt4p2K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fx = function(x) {  
    int_fyx = function(y) { 
      int_fzyx = function(z) { 
        return(LikFun(c(x,y,z,thetaOther)))
      }
      
      #set limits for Z
      fromZ = mxLims[3,1]
      toZ = min((1-x-y)/2,mxLims[3,2])
      val = 0
      if(fromZ < toZ) val = myInt(int_fzyx, fromZ, toZ) #integrate on z
      return(val)
    } #end x
    
    fromY = mxLims[2,1]
    toY = min(1-x,mxLims[2,2])
    val1 = 0
    if(fromY < toY) val1 = myInt(int_fyx, fromY  , toY)
    return(val1)
  } 
  return(2*myInt(int_fx,mxLims[1,1],mxLims[1,2]))
}

#1 knowns
.calcInt4p1K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fB1 = function(B1) {  
    int_fxB1 = function(x, part) { 
      int_fzxB1 = function(z) { 
        return(LikFun(c(B1,x,z,thetaOther)))
      }
      fromZ = mxLims[3,1]
      upperZ = ifelse(part==1,x,1-B1-2*x)
      toZ = min(upperZ,mxLims[3,2]) #cannot exceed 1/3
      val = 0
      if(fromZ < toZ) val = myInt(int_fzxB1, fromZ, toZ) #integrate on z
      return(val)
    } #end x
    fromX = mxLims[2,1]
    toX = min((1-B1)/2,mxLims[2,2]) #depends on value of B1
    xM = (1-B1)/3 #Middle X depends on value of B1
    val1 <- val2 <- 0
    if(fromX < xM) val1 = myInt(int_fxB1, fromX  , xM, part=1)
    if(xM < toX) val2 = myInt(int_fxB1, xM  , toX, part=2)
    return(val1 + val2)
  } 
  return(6*myInt(int_fB1,mxLims[1,1],mxLims[1,2]))
}

#####
#K=5#
#####

#Only unknowns
.calcInt5p0K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fw = function(w) { 
    int_fyw = function(y) {  
      int_fxyw = function(x,part) { 
        int_fzxyw = function(z) { 
          return(LikFun(c(x,z,y,w,thetaOther)))
        }
        val = 0
        fromZ = max(y,mxLims[2,1])
        upperZ = ifelse(part==1,x,1-2*x-y-w)
        toZ = min(upperZ,mxLims[2,2])
        if(fromZ < toZ) val = myInt(int_fzxyw, fromZ, toZ) #integrate on z
        return(val)
      } #end x
      xM = (1-y-w)/3 #mid point
      xU = 0.5 - y - 0.5*w #last point
      fromX = max(y,mxLims[1,1])
      toX = min(xU,mxLims[1,2])
      
      val1 <- val2 <- 0
      if(fromX < xM)  val1 = myInt(int_fxyw,fromX, xM, part=1)
      if(xM < toX)    val2 = myInt(int_fxyw,   xM, toX, part=2)
      return(val1 + val2)
    } #end y
    yU = (1-w)/4
    fromY = max(w,mxLims[3,1])
    toY = min(yU,mxLims[3,2])
    return( myInt(int_fyw, fromY,toY))
  } #end w
  return(120*myInt(int_fw,mxLims[4,1],mxLims[4,2]))
}

#4 or 5 knowns
.calcInt5p4K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fx = function(x) {  
    int_fyx = function(y) { 
      int_fzyx = function(z) { 
        int_fwzyx = function(w) { 
          return(LikFun(c(x,y,z,w,thetaOther)))
        }
        fromW = mxLims[4,1]
        toW = min(1-x-y-z,mxLims[4,2])
        val3 = 0
        if(fromW < toW) val3 = myInt(int_fwzyx, fromW, toW) #integrate on z
        return(val3)
      }
      fromZ = mxLims[3,1]
      toZ = min(1-x-y,mxLims[3,2])
      val2 = 0
      if(fromZ < toZ) val2 = myInt(int_fzyx, fromZ, toZ) #integrate on z
      return(val2)
    } #end x
    fromY = mxLims[2,1]
    toY = min(1-x,mxLims[2,2])
    val1 = 0
    if(fromY < toY) val1 = myInt(int_fyx, fromY  , toY)
    return(val1)
  } 
  return(myInt(int_fx,mxLims[1,1],mxLims[1,2]))
}

#3 knowns
.calcInt5p3K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fx = function(x) {  
    int_fyx = function(y) { 
      int_fzyx = function(z) { 
        int_fwzyx = function(w) { 
          return(LikFun(c(x,y,z,w,thetaOther)))
        }
        fromW = mxLims[4,1]
        toW = min( (1-x-y-z)/2,mxLims[4,2])
        val3 = 0
        if(fromW < toW) val3 = myInt(int_fwzyx, fromW, toW) #integrate on z
        return(val3)
      }
      fromZ = mxLims[3,1]
      toZ = min(1-x-y,mxLims[3,2])
      val2 = 0
      if(fromZ < toZ) val2 = myInt(int_fzyx, fromZ, toZ) #integrate on z
      return(val2)
    } #end x
    fromY = mxLims[2,1]
    toY = min(1-x,mxLims[2,2])
    val1 = 0
    if(fromY < toY) val1 = myInt(int_fyx, fromY  , toY)
    return(val1)
  } 
  return(2*myInt(int_fx,mxLims[1,1],mxLims[1,2]))
}

#2 knowns
.calcInt5p2K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fB1 = function(B1) { 
    int_fB2B1 = function(B2) {  
      int_fxB2B1 = function(x, part) { 
        int_fzxB2B1 = function(z) { 
          return(LikFun(c(B1,B2,x,z,thetaOther)))
        }
        fromZ = mxLims[4,1]
        upperZ = ifelse(part==1,x,1-B1-B2-2*x)
        toZ = min(upperZ,mxLims[4,2])
        val = 0
        if(fromZ < toZ) val = myInt(int_fzxB2B1, fromZ, toZ) #integrate on z
        return(val)
      } #end x
      xM = (1-B1-B2)/3 #Middle X depends on value of B1 and B2
      xU = (1-B1-B2)/2 #Upper X depends on value of B1 and B2
      fromX = mxLims[3,1]
      toX = min(xU,mxLims[3,2])
      val1 <- val2 <- 0
      if(fromX < xM) val1 = myInt(int_fxB2B1, fromX  , xM, part=1)
      if(xM < toX) val2 = myInt(int_fxB2B1, xM  , toX, part=2)
      return(val1 + val2)
    }
    fromB2 = mxLims[2,1]
    toB2 = min( (1-B1) ,mxLims[2,2]) #depends on value of B1
    return(myInt(int_fB2B1, fromB2, toB2))
  } 
  return(6*myInt(int_fB1,mxLims[1,1],mxLims[1,2]))
}

#1 known
.calcInt5p1K = function(LikFun,myInt,mxLims,thetaOther,...) {
  int_fB1 = function(B1) { 
    int_fyB1 = function(y) {  
      int_fxyB1 = function(x, part) { 
        int_fzxyB1 = function(z) { 
          return(LikFun(c(B1,x,z,y,thetaOther)))
        }
        fromZ = max(y, mxLims[3,1])
        upperZ = ifelse(part==1,x,1-B1-y-2*x)
        toZ = min(upperZ,mxLims[3,2])
        val = 0
        if(fromZ < toZ) val = myInt(int_fzxyB1, fromZ, toZ) #integrate on z
        return(val)
      } #end x
      
      xM = (1-B1-y)/3 #Middle X depends on value of B1 and y
      xU = (1-B1-2*y)/2 #Upper X depends on value of B1 and y
      fromX = max(y, mxLims[2,1])
      toX = min(xU,mxLims[2,2])
      val1 <- val2 <- 0
      if(fromX < xM) val1 = myInt(int_fxyB1, fromX, xM, part=1)
      if(xM < toX) val2 = myInt(int_fxyB1,   xM   , toX, part=2)
      return(val1 + val2)
    }
    fromY = mxLims[4,1]
    toY = min( (1-B1)/4 ,mxLims[4,2]) #depends on value of B1
    val0 = 0
    if(fromY < toY) val0 = myInt(int_fyB1, fromY, toY)
    return(val0)
  } 
  return(24*myInt(int_fB1,mxLims[1,1],mxLims[1,2]))
}
