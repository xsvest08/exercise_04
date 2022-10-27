
library(Biostrings)

NWScore <- function(seq1, seq2, match, mismatch, gap){
  # Stop conditions

  
  # Initialize col and rownames for matrices
  len1 <- nchar(seq1); len2 = nchar(seq2) # Save number of chars in each sequence
  seq1 <- unlist(strsplit(seq1, split = '')) # convert seq to character vector
  seq2 <- unlist(strsplit(seq2, split = '')) # convert seq to character vector
  
  # Initialize matrix M (for scores)
  M <- matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(M) <- c("-", seq1) # assign seq chars to matrix names
  colnames(M) <- c("-", seq2) # assign seq chars to matrix names
  M[1, ] <- cumsum(c(0, rep(gap, len2))) # Fill 1st row with gap penalites
  M[, 1] <- cumsum(c(0, rep(gap, len1))) # Fill 1st col with gap penalites
  
  # Initialize matrix D (for directions)
  D <- matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(D) <- c("-", seq1) # assign seq chars to matrix names
  colnames(D) <- c("-", seq2) # assign seq chars to matrix names
  D[1, ] <- rep("hor") # Fill 1st row with "hor" for horizontal moves
  D[, 1] <- rep("ver") # Fill 1st col with "ver" for vertical moves
  type <- c("dia", "hor", "ver") # Lookup vector
  
  # Compute scores and save moves
  for (i in 2:(len1 + 1)){# for every (initially zero) row
    for (j in 2:(len2 + 1)){# for every (initially zero) col
      hor <- M[i, j - 1] + gap # horizontal move = gap for seq1
      ver <- M[i - 1, j] + gap # vertical move = gap for seq2 
      dia <- ifelse(rownames(M)[i] == colnames(M)[j], # diagonal = ifelse(chars equal, match, mismatch) 
                   M[i - 1, j - 1] + match, 
                   M[i - 1, j - 1] + mismatch)
      M[i, j] = max(dia, hor, ver) # Save current (best) score in M
      D[i, j] = type[which.max(c(dia, hor, ver))] # Save direction of move in D
    }
  } 
  align1 <- c()
  align2 <- c() # Note: length of final alignments is unknown at this point
  
  while(i > 1 && j > 1){
    
    if(D[i, j] == "dia") {
      align1 = c(rownames(M)[i], align1)
      align2 = c(colnames(M)[j], align2)
      j <- j - 1
      i <- i - 1  # update indices
    } else if (D[i, j] == "ver") {
      align1 <- c(rownames(M)[i], align1)
      align2 <- c("-", align2) # vertical movement = gap for seq2
      i <- i - 1 # update indices
    } else if (D[i, j] == "hor") {
      align1 <- c("-", align1) # horizontal movement = gap for seq1
      align2 <- c(colnames(M)[j], align2) 
      j <- j - 1 # update indices
    } 
  }
  #return(list(aligned_seqs = matrix(c(align1, align2), byrow = TRUE, nrow = 2),
              #score = M[nrow(M), ncol(M)], score_matrix = M))
  return(score_matrix = M)
}

NWScore('AGTACGCA', 'TATGC', 2, -1, -2)
  


#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param seq1 DNAString object representing NT or AA sequence to align
#' @param seq2 DNAString object representing NT or AA sequence to align
#' @param align A list of DNAString objects with alignment of input sequences
#' @param match An integer value of a score for matching bases
#' @param mismatch An integer value of score for mismatching bases
#' @param gap An integer value of penalty for gap insertion
#' 
#' 
hirschberg_template <- function(seq1, seq2, align, match, mismatch, gap){
    
    seq1 <- unlist(strsplit(seq1, split = '')) # convert seq to character vector
    seq2 <- unlist(strsplit(seq2, split = '')) # convert seq to character vector
    first_align_row <- align[[1]] # initialize the first row of alignment
    second_align_row <- align[[2]] # initialize the second row of alignment
  
    if ((length(seq1)==0 && length(seq2)!=0)) # length of seq1 is equal to zero
    {
        for (i in 1:length(seq2)) # for each character in seq2
        {
            first_align_row <- c(DNAString(first_align_row),DNAString('-')) # add gap
            second_align_row <- c(DNAString(second_align_row),DNAString(seq2[i]))# add character from seq2
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if (length(seq2)==0 && length(seq1)!=0)# length of seq2 is equal to zero
    {
        for (j in 1:length(seq1))# for each character in seq1
        {
            first_align_row <- c(DNAString(first_align_row),DNAString(seq1[j])) # add character from seq1
            second_align_row <- c(DNAString(second_align_row),DNAString('-'))# add gap
        }
      
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if (length(seq1)==1 && length(seq2)==1) # length of seq1 and seq2 is equal to 1
    {
        first_align_row <- c(DNAString(first_align_row),DNAString(seq1))# add character from seq1
        second_align_row <- c(DNAString(second_align_row),DNAString(seq2))# add character from seq2
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else
    {
        x_len <- length(seq1) #delka X
        x_mid <- ceiling(x_len/2)# half of the length of seq1
        y_len <- length(seq2)# length of seq2
        
       
        
        left_score <-NWScore(seq1[1:x_mid],seq2,match,mismatch,gap) # NW score for the first half of seq1 and the whole seq2
        right_score <- NWScore(seq1[x_len:(x_mid+1)],seq2[y_len:1],match,mismatch,gap)# NW score for the second half of seq1 and the whole seq2 (both are reversed)
        count <- left_score+right_score[length(right_score):1]
        max_value <- which.max(count)
        y_mid <- max_value-1 # index of division for seq2
         
        # The first half
        if (y_mid==0)# index of division for seq2 is equal to 0
        {
            align <- hirschberg_template(seq1[1:x_mid],DNAString(),align,match,mismatch,gap)# call hirschberg function for the first half of seq1 and for an empty DNAString object
        }
        else
        {
            align <- hirschberg_template(seq1[1:x_mid],seq2[1:y_mid],align,match,mismatch,gap)# call hirschberg function for the first half of seq1 and for the first part of seq2
        }
        
        # The second half
        if ((x_mid + 1) > x_len) # seq1 cannot be further divided
        {
            align <- hirschberg_template(DNAString(),seq2[(y_mid+1):y_len],align,match,mismatch,gap)# call hirschberg function for an empty DNAString object and the second half of seq2
        }
        else if ((y_mid + 1) > y_len) # seq2 cannot be further divided
        {
            align <- hirschberg(seq1[(x_mid+1):x_len],DNAString(),align,match,mismatch,gap)# call hirschberg function for the second half of seq1 and for an empty DNAString object
        }
        else 
        {
            align <- hirschberg(seq1[(x_mid+1):x_len],Y[(y_mid+1):y_len],align,match,mismatch,gap)# call hirschberg function for the second half of seq1 and the second part of seq2
        }
    }
    return(align)
}

hirschberg_template('AGTACGCA', 'TATGC', list(DNAString(),DNAString()), 2, -1, -2)
