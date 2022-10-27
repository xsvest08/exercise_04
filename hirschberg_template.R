#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param seq1 DNAString object representing NT or AA sequence to align
#' @param seq2 DNAString object representing NT or AA sequence to align
#' @param align A list of DNAString objects with alignment of input sequences
#' @param match An integer value of a score for matching bases
#' @param mismatch An integer value of score for mismatching bases
#' @param gap An integer value of penalty for gap insertion
hirschberg_template <- function(seq1, seq2, align, match, mismatch, gap){
    
    first_align_row <- align[[1]] # initialize the first row of alignment
    second_align_row <- align[[2]] # initialize the second row of alignment
  
  
    if # length of seq1 is equal to zero
    {
        for # for each character in seq2
        {
            first_align_row <- # add gap
            second_align_row <- # add character from seq2
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if # length of seq2 is equal to zero
    {
        for # for each character in seq1
        {
            first_align_row <- # add character from seq1
            second_align_row <- # add gap
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if # length of seq1 and seq2 is equal to 1
    {
        first_align_row <- # add character from seq1
        second_align_row <- # add character from seq2
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else
    {
        x_len <- # length of seq1
        x_mid <- # half of the length of seq1
        y_len <- # length of seq2
        
        left_score <- # NW score for the first half of seq1 and the whole seq2
        right_score <- # NW score for the second half of seq1 and the whole seq2 (both are reversed)
        y_mid <- # index of division for seq2

        # The first half
        if # index of division for seq2 is equal to 0
        {
            align <- # call hirschberg function for the first half of seq1 and for an empty DNAString object
        }
        else
        {
            align <- # call hirschberg function for the first half of seq1 and for the first part of seq2
        }
        
        # The second half
        if ((x_mid + 1) > x_len) # seq1 cannot be further divided
        {
            align <- # call hirschberg function for an empty DNAString object and the second half of seq2
        }
        else if ((y_mid + 1) > y_len) # seq2 cannot be further divided
        {
            align <- # call hirschberg function for the second half of seq1 and for an empty DNAString object
        }
        else 
        {
            Align <- # call hirschberg function for the second half of seq1 and the second part of seq2
        }
    }
    return(align)
}