/*
Useful functions definition
*/

/* define a function for extracting the file name from a full path */
/* The full path will be the one defined by the user to indicate where the reference file is located */

def get_fullParent( the_string ) {
    /* pass the string as a file and get the full parent, then convert to string */
    file( the_string ).getParent().toString()
}