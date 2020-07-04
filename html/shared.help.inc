<?php
/*
 * SHARED.HELP.INC
 *
 * written by Conrad Shyu (shyu4751@uidaho.edu)
 *
 * Initiative for Bioinforamtics and Evolutionary Studies (IBEST)
 * Program of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, ID 83844
 *
 * last updated on May 7, 2005
 *
 * All rights reserved. Copyrights (C) 2004-2005.
*/

// make sure the file is not included twice
if ( $_SHARED_HELP_INC == 1 ):
else:
    $_SHARED_HELP_INC = 1;

/*
 * DrawFigure()
 * this function will generate the HTML tags for image display. it will also determines
 * the size of the image
*/
function DrawFigure(
    $title = "",    // caption of the image
    $image )        // filename of the image
{
    $size = getimagesize( $image );
    echo "<div align=\"center\">\n";
    echo "  <img title=\"$title\" src=\"$image\" width=" . $size[ 0 ] . " height=" . $size[ 1 ] . ">\n";
    echo "</div>\n";
}

/*
 * DrawImage()
 * this function will generate the HTML tags for image display. it will also determines
 * the size of the image
*/
function DrawImage(
    $title = "",    // caption of the image
    $image )        // filename of the image
{
    $size = getimagesize( $image );
    echo "<img title=$title src=$image width=" . $size[ 0 ] . " height=" . $size[ 1 ] . ">\n";
}

/*
 * noCache()
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
    $title = "MiCA: Online Help" )
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
<style type="text/css">
<!--
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
    font-size: 16;
    text-align: center;
}
a.menu:hover {
    color: #FF0000;
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
<table width="573"  border="0">
  <tr>
    <td><?php DrawImage( "MiCA Online Help", "images/mica3_logo.png" ); ?></td>
  </tr>
  <tr>
    <td><dt><b>&nbsp;&nbsp;<?php echo $title; ?></b></dt></td>
  </tr>
  <tr>
    <td><font size="3" face="Arial, Helvetica, sans-serif">
<!--
main contents start from here (inside a table)
-->

<?php
// *****************************************
// finish outputing HTML, resume PHP session

}   // end of DrawHeader()

/*
 * DrawFooter()
 * draw the html tags for copyrights information and other legal mumble jumbles
*/
function DrawFooter()
{
    global $PHP_SELF;
?>
    </font></td>
  </tr>
  <tr>
    <td><div align="center">
        <hr align="center" width="90%" noshade>
        <font color="#996600" size="2" face="Arial, Helvetica, sans-serif">
        MiCA: Microbial Community Analysis 3<br>
        All rights reserved. Copyrights &copy; University of Idaho, 2004-2005.</font><br>
        <em><font color="#990000">Last updated on <?php echo date( "F d Y H:i:s.", getlastmod() ); ?>
        </font></em>
    </div></td>
  </tr>
</table>
</body>
</html>

<?php
}   // end of DrawFooter()

endif;  // if _SHARED_HELP_INC = 1
?>
