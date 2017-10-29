
function (data, k) 
{
  data <- as.matrix(data)
  distdata <- dist.to.knn(data, k)
  p <- dim(distdata)[2L]
  lrddata <- reachability(distdata, k)
  lof <- rep(0, p)
  for (i in 1:p) {
    nneigh <- distdata[2, i] - distdata[1, i] + 1
    j <- seq(0, (nneigh - 1))
    local.factor <- sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
    lof[i] <- local.factor
  }
  lof
}



dist.to.knn <-
  function (dataset, neighbors) 
  {
    numrow = dim(dataset)[1]
    knndist = rep(0, 0)
    for (i in 1:numrow) {
      neighdist = knneigh.vect(dataset[i, ], dataset, neighbors)
      if (i == 2) {
        if (length(knndist) < length(neighdist)) {
          z = length(neighdist) - length(knndist)
          zeros = rep(0, z)
          knndist = c(knndist, zeros)
        }
        else if (length(knndist) > length(neighdist)) {
          z = length(knndist) - length(neighdist)
          zeros = rep(0, z)
          neighdist = c(neighdist, zeros)
        }
      }
      else {
        if (i != 1) {
          if (dim(knndist)[1] < length(neighdist)) {
            z = (length(neighdist) - dim(knndist)[1])
            zeros = rep(0, z * dim(knndist)[2])
            zeros = matrix(zeros, z, dim(knndist)[2])
            knndist = rbind(knndist, zeros)
          }
          else if (dim(knndist)[1] > length(neighdist)) {
            z = (dim(knndist)[1] - length(neighdist))
            zeros = rep(0, z)
            neighdist = c(neighdist, zeros)
          }
        }
      }
      knndist = cbind(knndist, neighdist)
    }
    return(knndist)
  }



reachability <- function(distdata,k) {
  
  p <- dim(distdata)[2]
  lrd <- rep(0,p)
  
  for (i in 1:p) {
    j <- seq(3,3+(distdata[2,i]-distdata[1,i]))
    # compare the k-distance from each observation to its kth neighbor
    # to the actual distance between each observation and its neighbors
    numneigh <- distdata[2,i]-distdata[1,i]+1
    temp <- rbind(diag(distdata[distdata[2,distdata[j,i]],distdata[j,i]]),distdata[j+numneigh,i])
    
    #calculate reachability
    reach <- 1/(sum(apply(temp,2,max))/numneigh)
    lrd[i] <- reach
  }
  lrd
}
