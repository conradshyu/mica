<?php
/*
 * RUNDTRFLP.PHP
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
drawHeader( "MiCA: Processing T-RFLP analysis ..." );

// randomly generate a filename to avoid conflicts
$unique = "_" . substr( md5( uniqid( rand() ) ), 1, 8 );

/*
 * gather required parameters from the previous web page
 *
 * to reduce the possible vunlerability in PHP - shut off the global variable
 * declare local variables for the required parameters
*/
$forward        = array( strtoupper( GetPostValue( 'forward' ) ) );
$reverse        = array( strtoupper( GetPostValue( 'reverse' ) ) );
$enzyme         = array( strtoupper( GetPostValue( 'enzyme' ) ) );
$database       = GetPostValue( 'database' );
$mismatch       = GetPostValue( 'mismatch' );
$max_base       = GetPostValue( 'max_base' );
$sort_by        = GetPostValue( 'sort_by' );
$sort_order     = GetPostValue( 'sort_order' );
$forward_shift  = GetPostValue( 'forward_shift' );
$reverse_shift  = GetPostValue( 'reverse_shift' );

/*
 * handle the file uploaded by the user
*/
$run_trflp = true;      // indicate if the trflp should be run
$forward_upload = "/export/wwwroot/$unique" . '.forward';
$reverse_upload = "/export/wwwroot/$unique" . '.reverse';

if ( !move_uploaded_file( $_FILES[ 'forward_sample' ][ 'tmp_name' ], $forward_upload ) )
{
    DrawError( "Forward fragment file missing or invalid<br>" );
    $run_trflp = false;
}

if ( !move_uploaded_file( $_FILES[ 'reverse_sample' ][ 'tmp_name' ], $reverse_upload ) )
{
    DrawError( "Reverse fragment file missing or invalid<br>" );
    $run_trflp = false;
}

// verify the forward primer sequences
$forward = CheckPrimer( "Forward", $forward );
$run_trflp = ( $run_trflp && count( $forward ) > 0 ) ? true : false;

// verify the reverse primer sequences
$reverse = CheckPrimer( "Reverse", $reverse );
$run_trflp = ( $run_trflp && count( $reverse ) > 0 ) ? true : false;

$enzyme = CheckEnzyme( $enzyme );
$run_trflp = ( $run_trflp && count( $enzyme ) > 0 ) ? true : false;

$sort_order = ( 4 * $sort_order ) + $sort_by;   // there are eight possible options

/*
 * construct the parameter file even though the parameters may not be valid
*/
$param = array(
     0 => $unique,              // unique filename
     1 => $database,            // database filename
     2 => $forward,             // list of forward primer(s)
     3 => $reverse,             // list of reverse primer(s)
     4 => $enzyme,              // list of restriction enzyme(s)
     5 => $max_base,            // specify the range that allows mismatches
     6 => $mismatch,            // specify the number of allowed mismatches
     7 => 0,                    // report only the shortest fragments or all fragments
     8 => $sort_order,          // sort the abundance data with the specified order
     9 => $unique . '.forward', // trflp forward fragment sample file
    10 => $unique . '.reverse', // trflp reverse fragment sample file
    11 => $forward_shift,       // forward fragment discrepancy
    12 => $reverse_shift );     // reverse fragment discrepancy

WriteParam( $param );

// parameter errors, process will not continue
if ( !$run_trflp )
{
    DrawError( "Please make sure that all required parameters have been submitted correctly." );
}
else
{
    $cmd_line = "./trflp $unique &";
    system( $cmd_line );        // call the search program
    $php_file = $unique . '.php';

    DrawSystemLoad(); DrawProgress( 5000 );

    echo "<p>You can view the results <a class=\"email\" href=\"$php_file\">here</a>.<br>\n";
}

DrawFooter();
?>
