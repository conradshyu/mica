<?php
/*
 * SHARED.WEB.INC
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, Idaho 83843-1010
 *
 * This file contains the commonly used PHP functions for the web site
 *
 * Written by Conrad Shyu (shyu4751@uidaho.edu)
 *
 * Copyrights 2004-2005. All rights reserved.
 *
 * last updated on May 6, 2005
 * revised on September 9, 2005
*/

// make sure the file is not included twice
if ( $_SHARED_WEB_INC == 1 ):
else:
    $_SHARED_WEB_INC = 1;

// set timeout
$SESSION_TIME_LIMIT = 120;

$banner_style   = "<font color=\"#FFFFFF\" size=\"2\" face=\"Arial, Helvetica, sans-serif\">";
$font_style     = "<font size=\"2\" face=\"Arial, Helvectica, sans-serif\">";
$align_center   = "<div align=\"center\">";
$align_right    = "<div align=\"right\">";
$align_left     = "<div align=\"left\">";

/*
 * get the values posted by the http protocol
*/
function GetPostValue(
    $key )          // the name of the variable
{
    return( $_POST[ $key ] );
}

/*
 * GetEndonuclease()
 * open the restriction enzyme file (enzyme.txt) and generate an array of 
 * restriction enzymes
*/
function GetEndonuclease()
{
    $handle = fopen( "conf/enzyme.txt", "r" );
    $enzyme = array();

    // read the file and parse the entires
    for ( $i = 0; !feof( $handle ); )
    {
        $buffer = trim( fgets( $handle, 1024 ) );

        // skip over the comments and empty line
        if ( $buffer[ 0 ] == '#' || !strlen( $buffer ) )
        {
            continue;
        }
        
        /*
         * get the name and sequence of the primers
        */
        $name = strtok( $buffer, "|\n" );     // ignore the name
        $enzyme[ $i++ ] = strtok( "|\n" );    // keep the enzyme
    }

    fclose( $handle ); return( $enzyme );
}

/*
 * GetNumericList()
 * generate a list of numbers within the specified ranges for the selection
*/
function GetNumericList(
    $lower,             // lower range of the numeric list
    $upper )            // upper range of the numeric list
{
    $select = '';

    for ( $i = $lower; !( $i > $upper ); ++$i )
    {
        $select .= "<option value=$i>&nbsp;$i</option>\n";
    }

    return( $select );
}

/*
 * CheckPrimer()
 * validate the primer; return an array of valid primers
*/
function CheckPrimer(
    $label = "Forward", // label of the primer
    $primer )           // the primer we want to check
{
    $valid = array();

    // make sure the sequence does not contain illegal characters
    for ( $i = 0, $j = 0; $i < count( $primer ); ++$i )
    {
        if ( !strlen( $primer[ $i ] ) )
        {
            continue;
        }

        echo "$label Primer:&nbsp;" . $primer[ $i ] . "&nbsp;";
        $buffer = $primer[ $i ];

        if ( strlen( strtok( $buffer, "AUGCMRWSYKVHDBNXT" ) ) )
        {
            DrawError( "(primer contain illegal character(s); skipped)" );
        }
        else
        {
            $valid[ $j++ ] = $primer[ $i ];
        }

        echo "<br>\n";
    }

    return( $valid );
}

/*
 * CheckEnzyme()
 * validate the restriction enzyme
*/
function CheckEnzyme(
    $enzyme )           // the restriction enzymes to be checked
{
    $valid = array();

    // make sure the sequence does not contain illegal characters
    for ( $i = 0, $j = 0; $i < count( $enzyme ); ++$i )
    {
        if ( !strlen( $enzyme[ $i ] ) )
        {
            continue;
        }

        echo "Restriction Enzyme:&nbsp;" . $enzyme[ $i ] . "&nbsp;";
        $buffer = $enzyme[ $i ];

        if ( strlen( strtok( $buffer, "^AUGCMRWSYKVHDBNXT" ) ) )
        {
            DrawError( "(enzyme contain illegal character(s); skipped)" );
        }
        else
        {
            $valid[ $j++ ] = $enzyme[ $i ];
        }

        echo "<br>\n";
    }

    return( $valid );
}

/*
 * write the parameters into a file
 *
 * the order of parameter
 *  0 => unique filename
 *  1 => database filename
 *  2 => list of forward primer(s)
 *  3 => list of reverse primer(s)
 *  4 => list of restriction enzyme(s)
 *  5 => specify the range that allows mismatches
 *  6 => specify the number of allowed mismatches
 *  7 => report only the shortest fragments or all fragments
 *  8 => sort option
 *  9 => trflp sample file
 * 10 => forward fragment discrepancy
 * 11 => reverse fragment discrepancy
*/
function WriteParam(
    $param )            // list of parameters
{
    $pfile = $param[ 0 ] . ".param";
    $handle = fopen( $pfile, "w" );

    if ( !$handle )
    {
        return( false );
    }

    // write the parameters into the file
    $today = date("F j, Y, g:i a" );
    fputs( $handle, "#####\n# Submitted on $today\n#####\n" );
    fputs( $handle, "filename = " . $param[ 0 ] . "\n" );
    fputs( $handle, "database = " . $param[ 1 ] . "\n" );

    // record the forward primers
    for ( $i = 0; $i < count( $param[ 2 ] ); ++$i )
    {
        fputs( $handle, "forward = " . $param[ 2 ][ $i ] . "\n" );
    }

    // record the reverse primers
    for ( $i = 0; $i < count( $param[ 3 ] ); ++$i )
    {
        fputs( $handle, "reverse = " . $param[ 3 ][ $i ] . "\n" );
    }

    // record the forward primers
    for ( $i = 0; $i < count( $param[ 4 ] ); ++$i )
    {
        fputs( $handle, "enzyme = " . $param[ 4 ][ $i ] . "\n" );
    }

    // make sure the mismatch is not larger than the max_base
    $param[ 6 ] = ( $param[ 6 ] > $param[ 5 ] ) ? $param[ 5 ] : $param[ 6 ];

    fputs( $handle, "max_base = " . $param[ 5 ] . "\n" );
    fputs( $handle, "mismatch = " . $param[ 6 ] . "\n" );
    fputs( $handle, "output_all = " . $param[ 7 ] . "\n" );
    fputs( $handle, "sort_option = " . $param[ 8 ] . "\n" );
    fputs( $handle, "forward_sample = " . $param[ 9 ] . "\n" );
    fputs( $handle, "reverse_sample = " . $param[ 10 ] . "\n" );
    fputs( $handle, "forward_shift = " . $param[ 11 ] . "\n" );
    fputs( $handle, "reverse_shift = " . $param[ 12 ] . "\n" );

    fclose( $handle ); return( true );
}

/*
 * DrawForwardPrimer()
 * open the forward primer list file (forward.txt) and generate a list of 
 * primers for the selection
*/
function DrawForwardPrimer()
{
    $select = "<option value='none'>Select a forward primer here ...</option>\n";
    $handle = fopen( "conf/forward.txt", "r" );

    // read the file and parse the entries
    while ( !feof( $handle ) )
    {
        $buffer = trim( fgets( $handle, 1024 ) );

        // skip over the comments
        if ( $buffer[ 0 ] == '#' || !strlen( $buffer ) )
        {
            continue;
        }

        /*
         * get the name and sequence of the primers
        */
        $name = strtok( $buffer, "|\n" );
        $site = strtok( "|\n" );

        /*
         * construct the selection list for the forward primers
        */
        $select .= "<option value=$site>$name&nbsp;&nbsp;&nbsp;&nbsp;$site</option>\n";
    }

    fclose( $handle ); return( $select );
}

/*
 * DrawReversePrimer()
 * open the reverse primer list file (forward.txt) and generate a list of 
 * primers for the selection
*/
function DrawReversePrimer()
{
    $select = "<option value='none'>Select a reverse primer here ...</option>\n";
    $handle = fopen( "conf/reverse.txt", "r" );

    // read the file and parse the entries
    while ( !feof( $handle ) )
    {
        $buffer = trim( fgets( $handle, 1024 ) );

        // skip over the comments
        if ( $buffer[ 0 ] == '#' || !strlen( $buffer ) )
        {
            continue;
        }

        /*
         * get the name and sequence of the primers
        */
        $name = strtok( $buffer, "|\n" );
        $site = strtok( "|\n" );

        /*
         * construct the selection list for the forward primers
        */
        $select .= "<option value=$site>$name&nbsp;&nbsp;&nbsp;&nbsp;$site</option>\n";
    }

    fclose( $handle ); return( $select );
}

/*
 * DrawEndonuclease()
 * open the restriction enzyme list file (enzyme.txt) and generate a list
 * of restriction enzymes for the selection
*/
function DrawEndonuclease()
{
    $select = "<option value='none'>Select a restriction enzyme here ...</option>\n";
    $handle = fopen( "conf/enzyme.txt", "r" );

    // read the file and parse the entries
    while ( !feof( $handle ) )
    {
        $buffer = trim( fgets( $handle, 1024 ) );

        // skip over the comments
        if ( $buffer[ 0 ] == '#' || !strlen( $buffer ) )
        {
            continue;
        }

        /*
         * get the name and sequence of the primers
        */
        $name = strtok( $buffer, "|\n" );
        $site = strtok( "|\n" );

        /*
         * construct the selection list for the forward primers
        */
        $select .= "<option value=$site>$name&nbsp;&nbsp;&nbsp;&nbsp;$site</option>\n";
    }

    fclose( $handle ); return( $select );
}

/*
 * DrawDatabase()
 * open the database list file (database.txt) and generate a list of database
 * for the dropdown menu
 *
 * written by Conrad Shyu, September 9, 2005
*/

function DrawDatabase(
    $caption = "Select a databse:" )    // caption of the selection
{
    $select = '';
    $handle = fopen( "conf/database.txt", "r" );

    // read the file and parse the entries
    while ( !feof( $handle ) )
    {
        $buffer = trim( fgets( $handle, 1024 ) );

        // skip over the comments
        if ( $buffer[ 0 ] == '#' || !strlen( $buffer ) )
        {
            continue;
        }

        /*
         * get the name and sequence of the primers
        */
        $name = strtok( $buffer, "|\n" );
        $site = strtok( "|\n" );

        /*
         * construct the selection list for the forward primers
        */
        $select .= "<option value=$site>$name</option>\n";
    }

    // draw the html tags on the page
    echo "<tr>\n";
    echo "  <td align='right' valign='baseline'>\n";
    echo "    <div align='right'><font face='Arial, Helvetica, sans-serif'>\n";
    echo "      $caption</font></div></td>\n";
    echo "  <td>&nbsp;</td>\n";
    echo "  <td align='left' valign='baseline'>\n";
    echo "    <font face='Arial, Helvetica, sans-serif'>\n";
    echo "      <select name='database'>$select</select></font></td>\n";
    echo "</tr>\n";

    fclose( $handle );
}

/*
 * list the sponsor of the MiCA website; the file is located at conf/sponsor.list
 * format: title|website link|name
*/
function DrawRelated()
{
    $handle = fopen( "conf/related.txt", "r" );

    // read the file and parse the entries
    while ( !feof( $handle ) )
    {
        $buffer = fgets( $handle, 8192 );

        // skip over the comments
        if ( $buffer[ 0 ] == '#' || !strlen( $buffer ) )
        {
            continue;
        }

        /*
         * progressively parse the string into three parts
         * the first part the title of the link, the second part is the website link
         * itself, and the last one is the caption for the website
        */
        echo "<tr><td bgcolor=\"#999999\"><div align=\"center\">\n";
        echo "<a class=\"menu\" title=\"" . strtok( $buffer, "|\n" ) . "\"\n";
        echo "  href=\"" . strtok( "|\n" ) . "\">\n";
        echo "  " . strtok( "|\n" ) . "</a></div></td></tr>\n";
    }

    fclose( $handle );
}

/*
 * DrawImage()
 * this function will generate the HTML tags for image display. it will also determines
 * the size of the image
*/
function DrawImage(
    $title,             // caption of the image
    $image )            // filename of the image
{
    $size = getimagesize( $image );
    echo "<img title=$title src=$image width=" . $size[ 0 ] . " height=" . $size[ 1 ] . ">\n";
}

/*
 * DrawSeparator()
 * draw an empty row with a horizontal bar
*/
function DrawSeparator()
{
    echo "<tr>\n";
    echo "  <td>&nbsp;</td>\n";
    echo "  <td>&nbsp;</td>\n";
    echo "  <td align='left' valign='baseline'><hr width='50%' align='left' noshade>\n";  
    echo "  </td>\n";
    echo "</tr>\n";
}

/*
 * DrawSelect()
 * draw the HTML selection and textfield tag
 * this function features a short javascript that update the textfield based on
 * what users have selected from the list
*/
function DrawSelect(
    $select_banner,     // the banner of the selection
    $select_name,       // the name of the selection
    $text_value,        // the name of the text field
    $select_list )      // the selection list
{
    echo "<tr>\n";
    echo "  <td align='right' valign='baseline'>\n";
    echo "    <div align='right'><font face='Arial, Helvetica, sans-serif'>\n";
    echo "    $select_banner</font></div></td>\n";
    echo "  <td>&nbsp;</td>\n";
    echo "  <td align='left' valign='baseline'>\n";
    echo "    <select name=$select_name onchange='$text_value.value = ";
    echo "this.options[ this.selectedIndex ].value'>\n";
    echo "$select_list</select></td>\n";
    echo "</tr>\n";

    echo "<tr>\n";
    echo "  <td align='right' valign='baseline'>\n";
    echo "    <div align='right'><font face='Arial, Helvetica, sans-serif'>\n";
    echo "    or enter here:</font></div></td>\n";
    echo "  <td>&nbsp;</td>\n";
    echo "  <td align='left' valign='baseline'>\n";
    echo "    <input type='text' name=$text_value size='32'>\n";
    echo "</tr>\n";
}

/*
 * DrawMismatch()
 * draw the number of mismatches permitted within the specified number of
 * bases of the 5' end of primer
*/
function DrawMismatch(
    $mismatch,          // maximum number of mismatches allowed
    $max_base )         // mismatches are allowed within the specified bases
{
    echo "<tr>\n";
    echo "  <td align='right' valign='baseline'>\n";
    echo "    <div align='right'><font face='Arial, Helvetica, sans-serif'>\n";
    echo "      Sensitivity:</font></div></td>\n";
    echo "  <td>&nbsp;</td>\n";
    echo "  <td align='left' valign='baseline'>\n";
    echo "    <font face='Arial, Helvetica, sans-serif'>digest allows at most &nbsp;\n";
    echo "      <select name='mismatch'>" . GetNumericList( 0, $mismatch );
    echo "        </select>&nbsp; mismatches within &nbsp;\n";
    echo "      <select name='max_base'>" . GetNumericList( 0, $max_base );
    echo "        </select>&nbsp; bases from 5' end of primer</font></td>\n";
    echo "</tr>\n";
}

/*
 * DrawOption()
 * draw the checkbox form object on the web page
*/
function DrawOption(
    $box_banner,        // the banner for the checkbox
    $box_name,          // the name of the checkbox
    $box_select )       // the banners of all possible checkbox selections
{
    echo "<tr>\n";
    echo "  <td align='right' valign='baseline'>\n";
    echo "    <div align='right'><font face='Arial, Helvetica, sans-serif'>\n";
    echo "      $box_banner</font></div></td>\n";
    echo "  <td>&nbsp;</td>\n";
    echo "  <td align='left' valign='top'>\n";
    echo "    <font face='Arial, Helvetica, sans-serif'>\n";

    // the default selection is always the first one
    $element = $box_select[ 0 ];
    echo "<input type='radio' name=$box_name value=0 checked>$element\n";

    // generate the list for all possible selections
    for ( $i = 1; $i < count( $box_select ); ++$i )
    {
        $element = $box_select[ $i ];
        echo "<input type='radio' name=$box_name value=$i>$element\n";
    }

    echo "</font></td>\n</tr>\n";
}

/*
 * draw the tips
*/
function DrawTip()
{
    $tips = array();
    $handle = fopen( "conf/tips.txt", "r" );

    // read the file and parse the entries
    for ( $i = 0; !feof( $handle ); $i++ )
    {
        $tips[ $i ] = trim( fgets( $handle, 8192 ) );
    }

    fclose( $handle );

    /*
     * an empty line will be included in the array and will be reported by 
     * the count function.
    */
    return( $tips[ rand( 2, count( $tips ) ) - 2 ] );
}

/*
 * DrawError()
 * draw the parameter error message
*/
function DrawError(
    $message )          // the error message we want to display
{
    echo "<font color='#FF0000'>$message</font>\n";
}

/*
 * DrawSystemLoad()
 * display a short message notifies the user that the process may take longer
 * to complete
*/
function DrawSystemLoad()
{
    echo "<br><b>NOTE</b>: Depending on the current system load, it may ";
    echo "take a while for the query to complete. Please wait until the ";
    echo "progress indicator completes before clicking the link below. ";
    echo "The output file will be written in the plain text and CSV";
    echo "formats. <b>The progress bar does not always give an accurate ";
    echo "indication of the current status. Please wait for a few more ";
    echo "seconds before click the link.</b><br><br>\n";

    $count = ( int ) exec( 'ps -A | grep ispar | wc -l' );

    echo "<div align='center'><strong>System Utilization Level: ";

    if ( $count < 3 )
    {
        echo "<font color='#336600'>Low</strong></font></div><br>\n";
    }
    elseif ( $count < 7 )
    {
        echo "<font color='#FF9900'>Medium</strong></font></div><br>\n";
    }
    else
    {
        echo "<font color='#FF0000'>High</strong></font></div><br>\n";
    }
}

/*
 * DrawProgress()
 * display the progress bar on the web page
*/
function DrawProgress(
    $timer )            // how long does the process take
{
    echo "<table align='center'><tr><td>\n";
    echo "<div style='font-size:10pt; padding:2px; border:solid black 1px'>\n";
    echo "  <span id='progress0'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress1'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress2'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress3'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress4'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress5'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress6'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress7'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress8'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "  <span id='progress9'> &nbsp; &nbsp; &nbsp; &nbsp; </span>\n";
    echo "</div></td></tr></table>\n";
    echo "<script language='javascript'>\n";
    echo "  setPeriod( $timer );\n";
    echo "  progressUpdate();\n";
    echo "</script>\n";
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
}

/*
 * DrawHeader()
 * this function outputs the beginnng HTML content for each page.
 * useful to keep the site layout consistent as well as easy to create new
 * pages and maintain the overall site look in one palce.
*/
function DrawHeader( $heading = "MiCA: Microbial Community Analysis" )
{
    static $draw_header = 0;        // determine if the header has been drawn
    global $HTTP_COOKIE_VARS;
    global $USER_TYPE;

    if ( $draw_header == 1 )
        return;
    else
        $draw_header = 1;

    noCache();

// output the HTML header, turn off PHP
// ************************************
?>

<html>
<head>
<title><?php echo $heading; ?></title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
<script language="JavaScript" src="shared.web.js" type="text/javascript"></script>
<style type="text/css">
<!--
h1 {
    color: #000066;
    background: #CCCCCC;
    font-size: 12pt;
    word-spacing: normal;
    border: medium double #CCCCCC;
    font-style: normal;
    text-align: center;
}
a {
    color: #ffffff;
    vertical-align: baseline;
    text-decoration: none;
}
a.email {
    color: #0000ff;
    vertical-align: baseline;
    text-decoration: none;
}
a.menu {
    color: #ffffff;
    font-family: Arial, Helvetica, sans-serif;
    background-color: #999999;
    border: thin solid #999999;
    font-size: 16;
    text-align: center;
}
a.menu:hover {
    color: #ff0000;
    font-family: Arial, Helvetica, sans-serif;
    font-size: 16;
    text-align: center;
    background-color: #999999;
}

dt {
    color: #990000;
    font-family: Arial, Helvetica, sans-serif;
    text-align: left;
    vertical-align: baseline;
    text-decoration: none;
    background-image: url(images/right_banner.png);
}
-->
</style>
</head>

<body>
<table width="827" border="0">
  <tr><td width="820"><div align="center"> 
    <table width="100%"  border="0">
      <tr><td width="82%"><a href="http://mica.ibest.uidaho.edu">
          <img src="images/mica3_logo.png" width="569" height="109"></a></td>
          <td width="18%"><a href="http://www.uidaho.edu">
          <img src="images/uofi.png" width="241" height="109"></a></td></tr>
    </table></div></td></tr>
  <tr><td bgcolor="#999999"><div align="center">
    <a class="menu" href="digest.php">Virtual Diegest (ISPaR)</a> &middot;
    <a class="menu" href="trflp.php">T-RFLP Analysis (APLAUS+)</a> &middot;
    <a class="menu" href="pat.php">T-RFLP Analysis (PAT+)</a> &middot;
    <a class="menu" href="primer.php">P.S.P.A.</a> &middot;
    <a class="menu" href="enzyme.php">E.R.P.A.</a></div></td></tr>
  <tr><td><table width="100%"  border="0">
    <tr><td bgcolor="#006699"><div align="center">
        <font color="#FFFFFF" face="Arial, Helvetica, sans-serif">
            <?php echo date("F j, Y"); ?></font></div></td>
        <td><dt>&nbsp;&nbsp;<?php echo "<strong>$heading</strong>"; ?></dt></td></tr>
      <tr><td width="20%" valign="top"><table width="100%"  border="0">
        <tr><td bgcolor="#FFFFCC"><font size="2" face="Arial, Helvetica, sans-serif">
          Do you know ...<hr align="center" width=90% size="1" noshade></hr>
            <?php echo DrawTip(); ?></font></td></tr>
          <?php DrawRelated(); ?>
  </table></td>
<!--
main contents start from here (in a table)
-->
          <td width="80%" valign="top"><font face="Arial, Helvetica, sans-serif">

<?php
// *****************************************
// finish outputing HTML, resume PHP session

}   // end of DrawHeader()

/*
 * DrawFooter()
 * output the HTML footer
*/
function DrawFooter( $update = NULL )
{
    global $PHP_SELF;

    // display the last modification date
    $update = is_null( $update ) ? date( "F d Y H:i:s.", getlastmod() ) : $update;

// output the HTML footer, turn off PHP
// ************************************
?>
          </font></td></tr>
        <tr>
          <td><div align="center">
            <img src="/cgi-bin/Count.cgi?df=mica.dat&dd=B&frgb=Gray"></div></td>
          <td><dt>&nbsp;&nbsp;m.i.c.r.o.b.i.a.l&nbsp;&nbsp;c.o.m.m.u.n.i.t.y
            &nbsp;&nbsp;a.n.a.l.y.s.i.s&nbsp;&nbsp;i.i.i</dt></td>
        </tr>
      </table></td>
  </tr>
  <tr><td bgcolor="#999999"><div align="center">
    <a class="menu" title="Return to MiCA homepage" href="default.php">
        MiCA Home</a> &middot;
    <a class="menu" title="What's new at MiCA" href="news.php">
        What's New</a> &middot;
    <a class="menu" title="Frequently asked question" href="faq.php">
        F.A.Q.</a> &middot;
    <a class="menu" title="About MiCA website" href="about.php">
        About MiCA</a> &middot;
    <a class="menu" title="MiCA online helps and tutorials"
        href="javascript:popUp( 'online_help.html', 808, 450 )">
        Online Help</a> &middot;
    <a class="menu" title="Email us your comments and suggestions"
        href="mailto:foster@uidaho.edu">Contact us</a></div></td></tr>
  <tr><td><div align="center">
    <font color="#996600" size="2" face="Arial, Helvetica, sans-serif">
      All rights reserved. Copyright &copy;, University of Idaho, 2007.<br>
      <em><font color="#990000">Last updated on <?php echo $update; ?></font></em></font></div></td>
  </tr>
</table>
</body>
</html>

<?php
// ***********************************************
// finish outputing the footer, resume PHP session
}

endif;  // if _PHP_SHARED_INC == 1
?>
