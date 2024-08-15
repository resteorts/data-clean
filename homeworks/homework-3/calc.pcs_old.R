###  Calculate all pairwise comparisons of records



###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################

###  calc.pcs:  

#####  a function that calculates all choose(n,2) pairwise comparisons of records
#####  the user can enter any input dataset
#####  the user can specify the type of comparison to use on each field



###  Parameters:

#####  data:  the original records, including a column with the true labels

#####  type:
#######  a vector of length equal to the number of columns in data
#######  i.e., ncol(data) == length(type)
#######  type[i] = "j" if jarowinkler similarities are wanted
#######  type[i] = "l" if levenshtein similarities are wanted
#######  type[i] = "e" if exact matching similarities are wanted
#######  so something like type = rep("l", ncol(data)) would work
#######  you do not have to use the same "j", "l", or "e" for each field
#######  if not specified, jarowinkler will be used on each field

#####  add.comp:  boolean
#######  if T, adds ".comp" to the column names
#######  if F, returns data.frame with same column names

calc.pcs = function(data, type = NULL, add.comp = T){
  nn = nrow(data)
  n.col = ncol(data)
  if(is.null(type))  type = rep("j", n.col)
  if(length(type) != ncol(data))  stop("length(type) != ncol(data)")

  
  ###  sometimes, weird things happen with NAs and formatting, depending on the data type
  for(jj in 1:n.col)  data[,jj] = ifelse(is.na(data[,jj]), NA, paste(data[,jj]))
  
  
  ###  calculate all the unique values of each field
  ###  note:  I force NA to be one of the values
  unique.vals = vector("list", n.col)
  for(jj in 1:n.col)  unique.vals[[jj]] = unique(c(NA, sort(data[,jj], dec = T)))
  
  
  ###  function for calculating the lookup index of each observation/field in the unique.vals vector
  ###  note:  I include a special case for NAs
  calc.index = function(x, col)  return(ifelse(is.na(x), 
                                               which(is.na(unique.vals[[col]])), 
                                               which(x == unique.vals[[col]])))
  
  
  ###  now, calculate the actual lookup indices for each column of data  
  indices = matrix(NA, nrow = nn, ncol = n.col)
  for(ii in 1:n.col)  indices[,ii] = unlist(sapply(data[,ii], calc.index, ii))
  indices = as.data.frame(indices)
  colnames(indices) = colnames(data)
  
  
  ###  calculate the similarity scores
  sim.scores = vector("list", n.col)
  for(ii in 1:n.col){
    if(type[ii] == "j"){
      sim.scores[[ii]] = calc.jw(unique.vals[[ii]])
    }  else if(type[ii] == "l"){
      sim.scores[[ii]] = calc.lev(unique.vals[[ii]])
    }  else if(type[ii] == "e"){
      sim.scores[[ii]] = calc.exact(unique.vals[[ii]])
    }
  }
  
  
  ###  calculate similarity matrices
  sim.mat = vector("list", n.col)
  for(ii in 1:n.col)  sim.mat[[ii]] = sim.scores[[ii]][indices[,ii], indices[,ii]]

  
  ###  fill in the resulting data.frame
  dtf = matrix(NA, nrow = choose(nn, 2), ncol = n.col+2)
  start.row = 0
  
  for(ii in 1:(nn-1)){
    ###  i add two extra cols so that the data.frame I return has the indices of the original dataset
    rows = matrix(NA, nrow = nn - ii, ncol = n.col + 2)
    for(jj in 1:n.col)  rows[,jj] = sim.mat[[jj]][(ii+1):nn, ii]
    
    ###  this fills in the indices of the comparison
    rows[,n.col+1] = ii
    rows[,n.col+2] = (ii+1):nn

    ###  rbind is a very slow function; there are definitely faster ways to do this
    #dtf = rbind(dtf, rows)
    
    ###  i fill in a big matrix rather than rbind()ing
    dtf[(start.row+1):(start.row+nn-ii), ] = rows
    start.row = start.row + nn - ii
  }
  
  ###  format the resulting data.frame
  rownames(dtf) = NULL
  dtf[which(is.na(dtf))] = 0  ## how should NAs be handled
  dtf = as.data.frame(dtf)
  if(add.comp)  names(dtf) = c(paste(colnames(data), "comp", sep = "."), "index1", "index2")
  if(!add.comp)  names(dtf) = c(colnames(data), "index1", "index2")
  
  ###  return the resulting data.frame
  return(dtf)
}



###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################


###  Function to calculate exact matching matrices
calc.exact = function(names){
  names = as.character(names)
  nn = length(names)
  if(nn < 2)  return(as.matrix(ifelse(is.na(names), 0, 1)))
  
  ex = matrix(0, nrow = nn, ncol = nn)
  
  for(ii in 1:(nn-1))
    ex[ii,(ii+1):nn] = 1*(names[ii] == names[(ii+1):nn])
  
  ex[which(is.na(ex))] = 0
  diag(ex) = 1
  ex[lower.tri(ex)] = t(ex)[lower.tri(ex)]
  return(ex)
}



###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################


###  function to calculate jarowinkler scores
calc.jw = function(names){
  names = as.character(names)
  nn = length(names)
  if(nn < 2)  return(as.matrix(ifelse(is.na(names), 0, 1)))
  
  jw = matrix(0, nrow = nn, ncol = nn)
  
  for(ii in 1:(nn-1))  jw[ii, (ii+1):nn] = jarowinkler(names[ii], names[(ii+1):nn])
  
  ###  how should NAs be handled?
  jw[which(is.na(jw))] = 0
  
  ###  the diagonal should always be a perfect match
  diag(jw) = 1
  
  ###  the lower triangle should always be the same as the upper triangle
  jw[lower.tri(jw)] = t(jw)[lower.tri(jw)]
  
  return(jw)
}



###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################


###  new function to calculate levenshtein scores
calc.lev = function(names){
  names = as.character(names)
  nn = length(names)
  if(nn < 2)  return(as.matrix(ifelse(is.na(names), 0, 1)))
  
  lev = matrix(0, nrow = nn, ncol = nn)
  
  for(ii in 1:(nn-1))  lev[ii, (ii+1):nn] = levenshteinSim(names[ii], names[(ii+1):nn])
  
  ###  how should NAs be handled?
  lev[which(is.na(lev))] = 0
  
  ###  the diagonal should always be a perfect match
  diag(lev) = 1
  
  ###  the lower triangle should always be the same as the upper triangle
  lev[lower.tri(lev)] = t(lev)[lower.tri(lev)]
  
  return(lev)
}



