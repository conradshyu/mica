<?php
/*
 * PAT.PHP
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
 * updated on May 2, 2005
 * updated on July 12, 2007
 *
 * All rights reserved, 2005. Copyrights (C) University of Idaho.
*/
require "shared.web.inc";
DrawHeader( "MiCA: T-RFLP Analysis (PAT+)" );

/*
 * draw matching window
*/
function DrawMatch(
    $lower = 1,     // lower bound of the window size
    $upper = 5 )    // upper bound of the window size
{
// turn off php to generate the html code base
?>
    <tr>
      <td align='right' valign='baseline'>
        <div align='right'><font face='Arial, Helvetica, sans-serif'>
          Window (Bin) Size:</strong></font></div></td>
      <td>&nbsp;</td>
      <td align='left' valign='baseline'>
        <font face='Arial, Helvetica, sans-serif'> forward match: &plusmn; &nbsp;
          <select name='forward_shift'><?php echo GetNumericList( $lower, $upper ); ?>
          </select>&nbsp; bps.
        </font></td>
    </tr>
<?php
}

/*
 * output the html code for file upload
*/
function DrawUpload(
    $caption = "Upload Sample File:",
    $tag_name )
{
?>

    <tr>
      <td valign='baseline'><div align='right'><font face='Arial, Helvetica, sans-serif'>
          <a class="email" href="javascript:popUp( 'sample_file.html', 420, 400 )">
            <?php echo $caption; ?></a></font></div></td>
      <td>&nbsp;</td>
      <td align='left' valign='baseline'><font face='Arial, Helvetica, sans-serif'>
          <input type='file' name=<?php echo $tag_name; ?>>
          &nbsp; (CSV or TXT file ONLY!)</font></td>
    </tr>

<?php
}   // end of DrawUpload()

// begin the main procedure
?>

<form name='trflp' enctype='multipart/form-data' method='post' action='runpat.php'>
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
DrawSelect( "Forward Primer:", "forward_list", "forward", DrawForwardPrimer() );
DrawUpload( "Forward Fragments:", "forward_sample" );

// select a reverse primer or enter one
DrawSelect( "Reverse Primer:", "reverse_list", "reverse", DrawReversePrimer() );

// select a restriction enzymes
DrawSelect( "Restriction Enzyme:", "enzyme_list", "enzyme", DrawEndonuclease() );

// draw the dropdown selections for database
DrawDatabase( "Select Database:" );

// draw the search sensitivity parameters
DrawMismatch( 10, 15 );

// draw a separator
DrawSeparator();

// draw the migraton shift discrepancies
DrawMatch( 1, 5 );

// draw the checkbox selections for sorting fields
DrawOption( "Sort By:", "sort_by",
    array( "Forward fragments", "Abundance", "Name" ) );

// draw the checkbox selections for sorting order
DrawOption( "Sort Order:", "sort_order", array( "Acending", "Descending" ) );

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
