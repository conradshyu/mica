<?php
/*
 * SHARED.DATA.INC
 *
 * This file contains some commonly used PHP functions for the display of the
 * data on the web site.
 *
 * Written by Conrad Shyu (shyu4751@uidaho.edu)
 * May 2, 2005
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology
 * Department of Biological Sciences
 * University of Idaho, Moscow, Idaho 83844
 *
 * All rights reserved. Copyrights (C) 2005.
*/

// make sure the file is not included twice
if ( $_SHARED_DATA_INC == 1 ):
else:
    $_SHARED_DATA_INC = 1;

$max_list       = 250;
$max_page       = 20;
$banner_style   = "<font color=\"#FFFFFF\" size=\"2\" face=\"Arial, Helvetica, sans-serif\">";
$font_style     = "<font size=\"2\" face=\"Arial, Helvectica, sans-serif\">";
$align_center   = "<div align=\"center\">";
$align_right    = "<div align=\"right\">";
$align_left     = "<div align=\"left\">";

/*
 * Repeat()
 * repeat a given expression a few times
*/
function Repeat(
    $count,             // number of iterations
    $expr )             // experssion to be repeated
{
    $str = '';

    for ( $i = 0; $i < $count; ++$i )
    {
        $str .= $expr;
    }

    return( $str );
}   // end of Repeat()

/*
 * GetPageList()
 * construct the page list
*/
function GetPageList(
    $current,           // the current position
    $total,             // total number of records
    $filename )         // the name of the data file
{
    global $max_list;
    global $max_page;

    $filename .= ".php?page=";
    $page = ceil( $total / $max_list );

    /*
     * a very complicated pagnation algorithm to make the web page looks
     * really professional!
     *
     * added by Conrad Shyu, September 14, 2006
    */
    if ( $page > $max_page )
    {
        $begin_link = ( $current / $max_list );
        $begin_link = ( $begin_link < ( $max_page / 2 ) ) ?
            0 : $begin_link - ( $max_page / 2 );
        $end_link = $begin_link + $max_page;
        $end_link = ( $end_link > $page ) ? $page : $end_link;
    }
    else
    {
        $begin_link = 0; $end_link = $page;
    }

    $backward = ( $current < $max_list ) ? $current : $current - $max_list;
    $page_list = "Page: <a href=$filename$backward>&lt;</a> ";

    for ( $i = $begin_link; $i < $end_link; $i++ )
    {
        $index = $i * $max_list;
        $page_list .= "<a href=$filename$index>" . ( $i + 1 ) . "</a> ";
    }

    $forward = $current + $max_list;
    $forward = ( $forward > $total ) ? $current : $forward;
    $page_list .= "<a href=$filename$forward>&gt;</a> ";

    return( $page_list );
}   // end of GetPageList()

/*
 * DrawPSPA()
 * draw the output from the primer sequence prevalence analysis program (PSPA)
*/
function DrawPSPA(
    $filename )
{
    // load the data into memory
    $handle = fopen( $filename . '.dat', "r" );
    $heading = array(
        0 => "Forward Primer", 1 => "Forward Matches",
        2 => "Reverse Primer", 3 => "Reverse Matches",
        5 => "Combined Matches" );
    global $align_center, $banner_style;

    if ( !$handle )
    {
        return( false );
    }

    $delimit = ",\n\""; $list = array();

    for ( $i = 0; $buffer = fgets( $handle, 4196 ); ++$i )
    {
        $list[ $i ][ 0 ] = strtok( $buffer, $delimit );
        $list[ $i ][ 1 ] = strtok( $delimit );
        $list[ $i ][ 2 ] = strtok( $delimit );
        $list[ $i ][ 3 ] = strtok( $delimit );
        $list[ $i ][ 4 ] = strtok( $delimit );
    }

    fclose( $handle );

    DrawCaption(
        "MiCA: Primer Sequence Prevalence Analysis Output",
        "The primer sequence prevalence analysis has completed.",
        $filename );

    echo "<tr><td><table width=100% border=0>\n";
    echo "  <tr bgcolor=#999999>\n";
    echo "    <td width=25%>$align_center $banner_style\n";
    echo "      <b>Forward Primer Sequence</b></font></div></td>\n";
    echo "    <td width=15%>$align_center $banner_style\n";
    echo "      <b>Forward Matches</b></font></div></td>\n";
    echo "    <td width=25%>$align_center $banner_style\n";
    echo "      <b>Reverse Primer Sequence</b></font></div></td>\n";
    echo "    <td width=15%>$align_center $banner_style\n";
    echo "      <b>Reverse Matches</b></font></div></td>\n";
    echo "    <td>$align_center $banner_style\n";
    echo "      <b>Combined Matches</b></font></div></td>\n";
    echo "  </tr>\n";

    // now create the table to display the digest output
    $color_code = false;

    for ( $l = 0; $l < $i; ++$l )
    {
        echo DrawEntry( $list[ $l ], $color_code );
        $color_code = !$color_code;
    }

    echo "  <tr bgcolor=#CCCCCC>\n";
    echo Repeat( 5, "    <td>&nbsp;</td>\n" );
    echo "  </tr>\n";
    echo "</table></td></tr>\n";
    return( true );
}   // end of DrawPSPA()

/*
 * DrawERPA()
 * draw the output from the enzyme resolving power analysis program (ERPA)
*/
function DrawERPA(
    $filename )     // output filename without extension
{
    global $banner_style, $align_center;

    // load the data into memory
    $handle = fopen( $filename . '.dat', "r" );
    $heading = array(
        0 => "Restriction Enzyme", 1 => "Total Hits", 2 => "Success",
        3 => "5'Unique", 4 => "5'Average", 5 => "5'Std. Dev.",
        6 => "3'Unique", 7 => "3'Average", 8 => "3'Std. Dev." );

    if ( !$handle )
    {
        return( false );
    }

    $delimit = ",\n\""; $list = array();

    for ( $i = 0; $buffer = fgets( $handle, 4196 ); ++$i )
    {
        $list[ $i ][ 0 ] = strtok( $buffer, $delimit );    // get enzyme

        for ( $j = 1; $j < count( $heading ); ++$j )
        {
            $list[ $i ][ $j ] = strtok( $delimit );
        }
    }

    fclose( $handle );

    DrawCaption(
        "MiCA: Enzyme Resolving Power Analysis (ERPA) Output",
        "The enzyme resolivng power analysis has completed.",
        $filename );

    echo "<tr><td><table width=100% border=0>\n";
    echo "  <tr bgcolor=#999999>\n";
    echo "    <td width=20%>$align_center $banner_style\n";
    echo "      <b>Restriction Enzyme</b></font></div></td>\n";

    for ( $k = 1; $k < count( $heading ); ++$k )
    {
        echo "    <td width=10%>$align_center $banner_style\n";
        echo "      <b>" . $heading[ $k ] . "</b></font></div></td>\n";
    }

    // now create the table to display the digest output
    $color_code = false;

    for ( $l = 0; $l < $i; ++$l )
    {
        // display the entry
        echo DrawEntry( $list[ $l ], $color_code );
        $color_code = !$color_code;
    }

    echo "  <tr bgcolor=#CCCCCC>\n";
    echo Repeat( count( $heading ), "    <td>&nbsp;</td>\n" );
    echo "  </tr>\n";
    echo "</table></td></tr>\n";
    return( true );
}   // end of DrawERPA()

/*
 * DrawTRFLP()
 * draw the output from the trflp program (APLAUS)
*/
function DrawTRFLP(
    $filename )     // output filename without extension
{
    global $banner_style, $align_center;

    // load the data into memory
    $handle = fopen( $filename . '.dat', "r" );
    $heading = array( 0 => "Forward", 1 => "Reverse", 2 => "Forward", 3 => "Reverse",
        4 => "Abundance", 5 => "Name" );

    if ( !$handle )
    {
        return( false );
    }

    $delimit = ",\n\""; $list = array();

    for ( $i = 0; $buffer = fgets( $handle, 4196 ); ++$i )
    {
        $list[ $i ][ 0 ] = strtok( $buffer, $delimit ); // predicted forward fragments
        $list[ $i ][ 1 ] = strtok( $delimit );          // predicted reverse fragments
        $list[ $i ][ 2 ] = strtok( $delimit );          // observed forward fragments
        $list[ $i ][ 3 ] = strtok( $delimit );          // observed reverse fragments
        $list[ $i ][ 4 ] = strtok( $delimit );          // normalized relative abundance
        $list[ $i ][ 5 ] = strtok( $delimit );          // organism name
    }

    fclose( $handle );

    DrawCaption(
        "MiCA: T-RFLP Analysis (APLAUS+) Output",
        "The T-RFLP analysis (APLUAS+) has completed and returned "
        . $i . " record(s).",
        $filename );

    echo "<tr><td><table width=100% border=0>\n";
    echo "  <tr bgcolor=#999999>\n";

    for ( $k = 0; $k < count( $heading ) - 1; ++$k )
    {
        echo "    <td width=10%>$align_center $banner_style\n";
        echo "      <b>" . $heading[ $k ] . "</b></font></div></td>\n";
    }

    echo "    <td>$align_center $banner_style\n";
    echo "      <b>" . $heading[ $k ] . "</b></font></div></td>\n";

    // now create the table to display the digest output
    $color_code = false;

    for ( $l = 0; $l < $i; ++$l )
    {
        // now, display the entry
        echo DrawEntry( $list[ $l ], $color_code );
        $color_code = !$color_code;
    }

    echo "  <tr bgcolor=#CCCCCC>\n";
    echo Repeat( count( $heading ), "    <td>&nbsp;</td>\n" );
    echo "  </tr>\n";
    echo "</table></td></tr>\n";
    return( true );
}   // end of DrawTRFLP()

/*
 * DrawPAT()
 * draw the output from the trflp program (PAT+)
*/
function DrawPAT(
    $filename )     // output filename without extension
{
    global $banner_style, $align_center;

    // load the data into memory
    $handle = fopen( $filename . '.dat', "r" );
    $heading = array( 0 => "Preicted", 1 => "Observed", 2 => "Abundance", 
        3 => "Accession", 4 => "Name" );

    if ( !$handle )
    {
        return( false );
    }

    $delimit = ",\n\""; $list = array();

    for ( $i = 0; $buffer = fgets( $handle, 4196 ); ++$i )
    {
        $list[ $i ][ 0 ] = strtok( $buffer, $delimit ); // predicted forward fragments
        $list[ $i ][ 1 ] = strtok( $delimit );          // observed forward fragments
        $list[ $i ][ 2 ] = strtok( $delimit );          // normalized relative abundance
        $list[ $i ][ 3 ] = strtok( $delimit );          // accession
        $list[ $i ][ 4 ] = strtok( $delimit );          // organism name
    }

    fclose( $handle );

    DrawCaption(
        "MiCA: T-RFLP Analysis (PAT+) Output",
        "The T-RFLP analysis (PAT+) has completed and returned "
        . $i . " record(s).",
        $filename );

    echo "<tr><td><table width=100% border=0>\n";
    echo "  <tr bgcolor=#999999>\n";

    for ( $k = 0; $k < count( $heading ) - 1; ++$k )
    {
        echo "    <td width=10%>$align_center $banner_style\n";
        echo "      <b>" . $heading[ $k ] . "</b></font></div></td>\n";
    }

    echo "    <td>$align_center $banner_style\n";
    echo "      <b>" . $heading[ $k ] . "</b></font></div></td>\n";

    // now create the table to display the digest output
    $color_code = false;

    for ( $l = 0; $l < $i; ++$l )
    {
        // now, display the entry
        echo DrawEntry( $list[ $l ], $color_code );
        $color_code = !$color_code;
    }

    echo "  <tr bgcolor=#CCCCCC>\n";
    echo Repeat( count( $heading ), "    <td>&nbsp;</td>\n" );
    echo "  </tr>\n";
    echo "</table></td></tr>\n";
    return( true );
}   // end of DrawPAT()

/*
 * DrawISPAR()
 * draw the output from the digest program (ISPaR)
*/
function DrawISPAR(
    $column,        // number of columns in the output file
    $filename,      // output filename without extension
    $page = 0 )     // which page to display
{
    global $max_list, $align_center, $banner_style;

    // load the data into memory
    $handle = fopen( $filename . '.dat', "r" );

    if ( !$handle )
    {
        return( false );
    }

    $delimit = ",\n\""; $list = array();

    for ( $i = 0; $buffer = fgets( $handle, 4196 ); ++$i )
    {
        $list[ $i ][ 0 ] = strtok( $buffer, $delimit );

        // parse the data from each entry
        for ( $j = 1; $j < $column; ++$j )
        {
            $list[ $i ][ $j ] = strtok( $delimit );
        }
    }

    fclose( $handle ); $page = ( $page > $i ) ? 0 : $page;
    DrawISPARCaption(
        "MiCA: Virtual Digest (ISPaR) Output",
        "The virtual digest (ISPaR) has completed and returned $i record(s).",
        GetPageList( $page, $i, $filename ), $filename );

    echo "<tr><td><table width=100% border=0>\n";
    echo "  <tr bgcolor=#999999>\n";

    // print out the heading for the data
    for ( $k = 0; $k < ( ( $column - 3 ) * 0.5 ); ++$k )
    {
        echo "    <td width=8%><div align=center> $banner_style\n";
        echo "      <b>Forward</b></font></div></td>\n";
        echo "    <td width=8%><div align=center> $banner_style\n";
        echo "      <b>Reverse</b></font></div></td>\n";
    }

    echo "    <td width=10%>$align_center $banner_style\n";
    echo "      <b>Accession</b></font></div></td>\n";
    echo "    <td width=10% bgcolor=#999999>$align_center $banner_style\n";
    echo "      <b>Locus</b></font></div></td>\n";
    echo "    <td>$align_center $banner_style\n";
    echo "      <b>Name</b></font></div></td>\n";
    echo "  </tr>\n";

    // now create the table to display the digest output
    $color_code = false;

    for ( $l = 0; ( $l < $max_list ) && ( $page < $i ); ++$l, ++$page )
    {
        // now, display the entry
        echo DrawEntry( $list[ $page ], $color_code );
        $color_code = !$color_code;
    }

    echo "    <tr bgcolor=#CCCCCC>\n";
    echo Repeat( $column, "      <td>&nbsp;</td>\n" );
    echo "    </tr>\n";
    echo "</table></td></tr>\n";
    return( true );
}

/*
 * DrawISPARCaption()
 * draw the caption of the output page
*/
function DrawISPARCaption(
    $title,             // title of the output page
    $caption,           // caption of the output page
    $page_list,         // page list string
    $filename )         // output filename without extension
{
/*
 * turn off php
*/
?>
  <tr>
    <td><table width="100%" border="1">
        <tr>
          <td width="80%">
            <font color="#990000" size="2" face="Arial, Helvetica, sans-serif">
            <b><?php echo $title; ?></b></font></td>
          <td width="20%">
            <font size="2" face="Arial, Helvectica, sans-serif">
            <?php echo "<a href=" . $filename . ".csv>"; ?>
              CSV format (Spreadsheet)</a></font></td>
        </tr>
        <tr>
          <td><font size="2" face="Arial, Helvectica, sans-serif">
            <?php echo $caption; ?></font></td>
          <td><font size="2" face="Arial, Helvectica, sans-serif">
            <?php echo "<a href=" . $filename . ".pat>"; ?>
              PAT format</a></font></td>
        </tr>
        <tr>
          <td><div align="center">
            <font size="2" face="Arial, Helvectica, sans-serif">
            <?php echo $page_list; ?></font></div></td>
          <td><font size="2" face="Arial, Helvectica, sans-serif">
            <?php echo "<a href=" . $filename . ".txt>"; ?>
              Plain text format</a></font></td>
        </tr>
      </table></td>
  </tr>

<?php
}

/*
 * DrawCaption()
 * draw the caption of the output page
*/
function DrawCaption(
    $title,         // title of the output page
    $caption,       // caption of the output page
    $filename )     // output filename without extension
{
/*
 * turn off php
*/
?>
  <tr>
    <td><table width="100%" border="1">
        <tr>
          <td width="80%">
            <font color="#990000" size="2" face="Arial, Helvetica, sans-serif">
            <b><?php echo $title; ?></b></font></td>
          <td width="20%">
            <font size="2" face="Arial, Helvectica, sans-serif">
            <?php echo "<a href=" . $filename . ".csv>"; ?>
              CSV format (Spreadsheet)</a></font></td>
        </tr>
        <tr>
          <td><font size="2" face="Arial, Helvectica, sans-serif">
            <?php echo $caption; ?></font></td>
          <td><font size="2" face="Arial, Helvectica, sans-serif">
            <?php echo "<a href=" . $filename . ".txt>"; ?>
              Plain text format</a></font></td>
        </tr>
      </table></td>
  </tr>

<?php
}

/*
 * DrawEntry()
 * draw the data entry
*/
function DrawEntry(
    $entry,             // data entry (array)
    $color_code )       // color code
{
    global $align_right, $align_left, $font_style;

    $list = ( $color_code ) ? "<tr bgcolor=#CCCFF0>\n" : "<tr>\n";

    for ( $i = 0; $i < count( $entry ); ++$i )
    {
        $align = ( is_numeric( $entry[ $i ] ) ) ? $align_right : $align_left;
        $list .= "<td>$align $font_style" . $entry[ $i ] . "</font></div></td>\n";
    }

    $list .= "</tr>\n"; return( $list );
}

/*
 * NoCache()
 * outputs various headers to ensure the page is not cached in the browser.
*/
function NoCache()
{
    Header( "Expires: 0" );
    Header( "Last modified: now" );
    Header( "Pragma: no-cache" );
    Header( "Cache-Control: no-cache, must-validate" );
}   // end of NoCache()

/*
 * DrawHeader()
 * draw the common header for all output; the header should include the logo of
 * MiCA, dowload links for the plain text and CSV formats, brief report on the
 * duration of the analysis or digest
*/
function DrawHeader(
    $title )
{
    static $draw_header = 0;        // determine if the header has been drawn
    global $HTTP_COOKIE_VARS;
    global $USER_TYPE;

    if ( $draw_header == 1 )
        return;
    else
        $draw_header = 1;

    NoCache();

/*
 * turn off PHP
*/
?>
<html>
<head>
<title><?php echo $title; ?></title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
</head>

<body>
<table width="820" border="0">
  <tr>
    <td height="109"> 
      <table width="100%" border="0">
        <tr>
          <td width="50%" height="103"><a href="http://mica.ibest.uidaho.edu">
            <img src="images/data_banner.png" width="405" height="101"></a></td>
          <td width="50%" valign="middle">
            <font color="#990000" size="2" face="Arial, Helvetica, sans-serif">
              <b>MiCA 3: Microbial Community Analysis III</b><br><br>
              <font color="#000099">
              Initiative for Bioinformatics and Evolutionary Studies (IBEST)<br>
              Program of Bioinforamtics and Computational Biology (BCB)<br>
              Department of Biological Sciences<br>
              University of Idaho, Moscow, Idaho 83844</font></font></font></td>
        </tr>
      </table></td>
  </tr>
<?php
/*
 * turn PHP back on
*/
}

/*
 * DrawFooter()
 * draw the footer of the web page
*/
function DrawFooter()
{
/*
 * turn off PHP
*/
?>
  <tr>
    <td><div align="center">
        <font color="#996600" size="2" face="Arial, Helvetica, sans-serif">
        MiCA: Microbial Community Analysis 3<br>
        All rights reserved. Copyrights &copy; University of Idaho, 2004-2005.</font>
    </div></td>
  </tr>
</table></body></html>

<?php

/*
 * turn PHP back on
*/
}

endif;  // if _SHARED_DATA_INC == 1
?>
