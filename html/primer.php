<?php
/*
 * PRIMER.PHP
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, ID 83843-1010
 *
 * this PHP file accepts the primers, restriction enzymes, and other parameters
 * for the virtual digest
 *
 * written by Conrad Shyu (shyu4751@uidaho.edu)
 * updated on May 6, 2005
 *
 * All rights reserved, 2005. Copyrights (C) University of Idaho.
*/
require "shared.web.inc";
DrawHeader( "MiCA: Primer Sequence Prevalence Analysis" );

/*
 * begin the HTML page construction
 * turn off PHP
*/
?>

<form name='primer' method='post' action='runprimer.php'>
  <table width='100%' border='0'>
    <tr>
      <td width='23%'>&nbsp;</td>
      <td width='1%'>&nbsp;</td>
      <td width='76%'>&nbsp;</td>
    </tr>

<?php
/*
 * return to PHP
*/

// select a forward primer or enter one
$forward_primer = DrawForwardPrimer();
DrawSelect( "Forward Primer:", "forward_list", "forward1", $forward_primer );
DrawSelect( "Forward Primer:", "forward_list", "forward2", $forward_primer );
DrawSelect( "Forward Primer:", "forward_list", "forward3", $forward_primer );

// select a reverse primer or enter one
$reverse_primer = DrawReversePrimer();
DrawSelect( "Reverse Primer:", "reverse_list", "reverse1", $reverse_primer );
DrawSelect( "Reverse Primer:", "reverse_list", "reverse2", $reverse_primer );
DrawSelect( "Reverse Primer:", "reverse_list", "reverse3", $reverse_primer );

// draw the dropdown selections for database
DrawDatabase( "Select Database:" );

// draw the search sensitivity parameters
DrawMismatch( 10, 15 );

// draw a separator
DrawSeparator();

// draw the checkbox selections for sorting field
DrawOption( "Sort By:", "sort_field",
    array( "Forward Primers", "Reverse Primers", "Primer Matches" ) );

// draw the checkbox selections for sorting order
DrawOption( "Sort Order:", "sort_order",
    array( "Ascending", "Descending" ) );

/*
 * turn off PHP
*/
?>
  </table><br><br>
  <center>
    <input type='submit' value='Submit'>&nbsp;&nbsp;&nbsp;
    <input type='reset' value='Reset'>
  </center>
</form>

<?php
/*
 * draw the footer of the page
*/

DrawFooter();
?>
