<?php
/*
 * RUNPRIMER.PHP
 *
 * written by Conrad Shyu (shyu4751@uidaho.edu)
 * May 6, 2005
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, ID 83843-1010
 *
 * All rights reserved. Copyrights (C) 2005.
*/
require "shared.web.inc";
drawHeader( "MiCA: Processing primer sequence prevalence analysis ..." );

// randomly generate a filename to avoid conflicts
$unique = "_" . substr( md5( uniqid( rand() ) ), 1, 8 );

/*
 * conver the parameters to match the function pointer array
 *
 * bool ( *SortOption[ 6 ] )( const stRECORD&, const stRECORD& ) =
 * {
 *     SortForwardPrimerA, // sort by forward primers in ascending order
 *     SortReversePrimerA, // sort by reverse primers in ascending order
 *     SortPrimerMatchA,   // sort by number of both primers matches in ascending order
 *     SortForwardPrimerD, // sort by forward primers in descending order
 *     SortReveresPrimerD, // sort by reverse primers in descending order
 *     SortPrimerMatchD    // sort by number of primer matches in descending order
 * };
*/
function SortOption(
    $field,         // forward fragments, reverse fragments, or organism names
    $order )        // ascending or descending order
{
    $option = array(
        0 => array( 0 => 0, 1 => 1, 2 => 2 ),
        1 => array( 0 => 3, 1 => 4, 2 => 5 ) );

    return( $option[ $order ][ $field ] );
}

/*
 * gather required parameters from the previous web page
 *
 * to reduce the possible vunlerability in PHP - shut off the global variable
 * declare local variables for the required parameters
*/
$forward    = array(
    strtoupper( GetPostValue( 'forward1' ) ),
    strtoupper( GetPostValue( 'forward2' ) ),
    strtoupper( GetPostValue( 'forward3' ) ) );
$reverse    = array(
    strtoupper( GetPostValue( 'reverse1' ) ),
    strtoupper( GetPostValue( 'reverse2' ) ),
    strtoupper( GetPostValue( 'reverse3' ) ) );
$database   = GetPostValue( 'database' );
$mismatch   = GetPostValue( 'mismatch' );
$max_base   = GetPostValue( 'max_base' );
$sort_field = GetPostValue( 'sort_field' );
$sort_order = GetPostValue( 'sort_order' );

$run_primer = true;         // indicate if the digest should be run

// verify the forward primer sequences
$forward = CheckPrimer( "Forward", $forward );
$run_primer = ( $run_primer && count( $forward ) > 0 ) ? true : false;

// verify the reverse primer sequences
$reverse = CheckPrimer( "Reverse", $reverse );
$run_primer = ( $run_primer && count( $reverse ) > 0 ) ? true : false;

/*
 * construct the parameter file even though the parameters may not be valid
*/
$param = array(
    0 => $unique,       // unique filename
    1 => $database,     // database filename
    2 => $forward,      // list of forward primer(s)
    3 => $reverse,      // list of reverse primer(s)
    4 => array(),       // list of restriction enzyme(s); not used here
    5 => $max_base,     // specify the range that allows mismatches
    6 => $mismatch,     // specify the number of allowed mismatches
    7 => 0,             // report all or shortest fragment; not used here
    8 => SortOption( $sort_field, $sort_order ),
    9 => 'not_used',    // trflp sample file; not used here
    10 => 'not_used',
    11 => 0,            // forward fragment discrepancy; not used here
    12 => 0 );          // reverse fragment discrepancy; not used here

WriteParam( $param );

// parameter errors, process will not continue
if ( !$run_primer )
{
    DrawError( "Please make sure that all required parameters have been submitted correctly." );
}
else
{
    $cmd_line = "./pspa $unique &";
    system( $cmd_line );        // call the search program

    $php_file = $unique . ".php";

    DrawSystemLoad(); DrawProgress( 7500 );

    echo "<p>You can view the results <a class=\"email\" href=\"$php_file\">here</a>.<br>\n";
}

DrawFooter();
?>
