<?php
/*
 * RUNDIGEST.PHP
 *
 * written by Conrad Shyu (shyu4751@uidaho.edu)
 * May 2, 2005
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, ID 83843-1010
 *
 * All rights reserved. Copyrights (C) 2005.
*/
require "shared.web.inc";
drawHeader( "MiCA: Processing virtual digest ..." );

// randomly generate a filename to avoid conflicts
$unique = "_" . substr( md5( uniqid( rand() ) ), 1, 8 );

/*
 * conver the parameters to match the function pointer array
 *
 * bool ( *SortOption[ 10 ] )( const stRECORD&, const stRECORD& ) =
 * {
 *     SortForwardFragmentA,   // sort by forward fragments in ascending order
 *     SortReverseFragmentA,   // sort by reverse fragments in ascending order
 *     SortShortestForwardA,   // sort by shortest forward fragment in ascending order
 *     SortShortestReverseA,   // sort by shortest reverse fragment in ascending order
 *     SortOrganismA,          // sort by organism name in ascending order
 *     SortForwardFragmentD,   // sort by forward fragments in descending order
 *     SortReverseFragmentD,   // sort by reverse fragments in descending order
 *     SortShortestForwardD,   // sort by shortest forward fragment in descending order
 *     SortShortestReverseD,   // sort by shortest reverse fragment in descending order
 *     SortOrganismD           // sort by organism name in descending order
 * };
*/
function SortOption(
    $report,    // only the shortest or all fragments
    $field,     // forward fragments, reverse fragments, or organism names
    $order )    // ascending or descending order
{
    $ascend = array(
        0 => array( 0 => 2, 1 => 3, 2 => 4 ),
        1 => array( 0 => 0, 1 => 1, 2 => 4 ) );
    $descend = array(
        0 => array( 0 => 7, 1 => 8, 2 => 9 ),
        1 => array( 0 => 5, 1 => 6, 2 => 9 ) );

    return( ( $order == 0 ) ? $ascend[ $report ][ $field ] : $descend[ $report ][ $field ] );
}

/*
 * gather required parameters from the previous web page
 *
 * to reduce the possible vunlerability in PHP - shut off the global variable
 * declare local variables for the required parameters
*/
$forward    = array( strtoupper( GetPostValue( 'forward' ) ) );
$reverse    = array( strtoupper( GetPostValue( 'reverse' ) ) );
$enzyme     = array(
    strtoupper( GetPostValue( 'enzyme1' ) ),
    strtoupper( GetPostValue( 'enzyme2' ) ),
    strtoupper( GetPostValue( 'enzyme3' ) ) );
$database   = GetPostValue( 'database' );
$mismatch   = GetPostValue( 'mismatch' );
$max_base   = GetPostValue( 'max_base' );
$sort_by    = GetPostValue( 'sort_by' );
$output_all = GetPostValue( 'output_all' );
$sort_order = GetPostValue( 'sort_order' );

$run_digest = true;         // indicate if the digest should be run

// verify the forward primer sequences
$forward = CheckPrimer( "Forward", $forward );
$run_digest = ( $run_digest && count( $forward ) > 0 ) ? true : false;

// verify the reverse primer sequences
$reverse = CheckPrimer( "Reverse", $reverse );
$run_digest = ( $run_digest && count( $reverse ) > 0 ) ? true : false;

// verify the restriction enzymes
$enzyme = CheckEnzyme( $enzyme );
$run_digest = ( $run_digest && count( $enzyme ) > 0 ) ? true : false;

/*
 * construct the parameter file even though the parameters may not be valid
*/
$param = array(
    0 => $unique,       // unique filename
    1 => $database,     // database filename
    2 => $forward,      // list of forward primer(s)
    3 => $reverse,      // list of reverse primer(s)
    4 => $enzyme,       // list of restriction enzyme(s)
    5 => $max_base,     // specify the range that allows mismatches
    6 => $mismatch,     // specify the number of allowed mismatches
    7 => $output_all,   // report only the shortest fragments or all fragments
    8 => SortOption( $output_all, $sort_by, $sort_order ),
    9 => 'not_used',    // trflp sample file; not used here
    10 => 'not_used',
    11 => 0,            // forward fragment discrepancy; not used here
    12 => 0 );          // reverse fragment discrepancy; not used here

WriteParam( $param );

// parameter errors, process will not continue
if ( !$run_digest )
{
    DrawError( "Please make sure that all required parameters have been submitted correctly." );
}
else
{
    $cmd_line = "./ispar $unique &";
    system( $cmd_line );        // call the search program

    $php_file = $unique . ".php?page=0";

    DrawSystemLoad(); DrawProgress( 7500 );

    echo "<p>You can view the results <a class=\"email\" href=\"$php_file\">here</a>.<br>\n";
}

DrawFooter();
?>
