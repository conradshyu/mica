<?php
require "shared.web.inc";
drawHeader( "MiCA: What's New!" );
?>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>July 12, 2007</strong>: A new feature has been added to MiCA. If 
you only use one fluorescently labeled primer (forward primer is assumed), 
please use PAT+ for analysis. PAT+ only requires the profile generated 
using the forward primer. For <i>in silico</i> amplification to perform 
properly, the reverse primer, however, is still required.</p> 

<p><img src="images/bglime.gif" width="13" height="13">
<strong>June 29, 2007</strong>: The database has been updated to RDP9.1, 
update 51. The database contains 380,520 16S rRNA sequences. In addition, 
the default database has been set to the one that only contains sequences 
with good quality (> 1,200bps). Please be aware that the time required to 
complete a query increases significantly, as the number of sequences 
increases rapidly.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>July 26, 2006</strong>: A new database has been added to MiCA. The 
database contains 95,140 16S rRNA gene sequence download from <a 
class="email" href="http://greengenes.lbl.gov">http://greengenes.lbl.gov</a>.
The greengenes database project was created by T. Z. DeSantis, I. 
Dubosarskiy, S. R. Murray, and G. L. Andersen.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>July 24, 2006</strong>: I (Conrad Shyu) am no longer the system 
administrator of MiCA. Please send all questions and comments to <a 
class="email" href="mailto:foster@uidaho.edu">Dr. James A. Foster</a>. If 
you have any technical questions about MiCA, you can still <a 
class="email" href="mailto:shyu4751@yahoo.com">email me</a>. I will try to 
answer them as soon as possible.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>July 14, 2006</strong>: I have completed a major updated on the 
APLAUS+ web interface and algorithm. The new algorithm can better estimate 
the relative abundance based on the submitted data. The output is also 
much easier to understand. The new algorithm requires separate files for 
the forward and reverse fragments. Please check out the help file for more 
information.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>May 10, 2006</strong>: The database has been updated to RDP9, 
update 38. Currently, there are 210,005 16S rRNA sequences in the 
database. Several other customized databses have also been added to 
MiCA.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>September 9, 2005</strong>: Some minor changes have been made to 
the web interface in order to support ever increasing numbers of databases 
on MiCA. Two additional databases have been added; (1) full length of high 
quality 16S sequences from the colon of humans, pigs, and the rumen, (2) 
representative sequences of the major lineages of bacteria in the colon of 
mammals.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>October 29, 2004</strong>: The MiCA database has been updated to 
Ribosomal Database Project-II Release 9. The new database contains 108,781 
unaligned and annotated Bacterial small-subunit rRNA sequences.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>March 7, 2004</strong>: The MiCA query program has been updated. 
Ambiguities are permitted in the restriction enzymes. Primer specificity 
and enzyme fidelity analysis pages are now fully functional. The community 
profile analysis pages are still under construction. We are currently 
updating the online help pages.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>September 4, 2006</strong>: MiCA is currently running on the new 
server. Everything should be fully functional. If you experience any 
problems, please contact the system administrator. The old server will be 
shut down shortly.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>September 4, 2003</strong>: We have updated the database to the 
newest version of RDP. The new database contains 62,762 bacterial and 
1,172 archaeal sequences. Therefore, queries will take much longer to run. 
In addition, we are migrating to entire web site to a new and powerful 
server. The new server will have a new domainname, 
http://mica.ibest.uidaho.edu/, and the old one will phase out shortly.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>February 9, 2003</strong>: The MiCA database has been restored. We 
are still working on migrating the updated RDP data into the MiCA 
database. If you experience any problems, please notify us.</p>

<p><img src="images/bglime.gif" width="13" height="13">
<strong>July 31, 2002</strong>: Multiple digest with at most three enzymes 
are now possible. Query can return all fragment lengths or only the 
shortest ones. The output is available in CSV (Comma Separated Values) and 
plain-text format. CSV file format allows data to be easily retrieved into 
spreadsheet or database programs. The majority of spreadsheet 
applications, such as Linux Gnumeric or Microsoft Excel, can read and 
write CSV files. Internet Explorer will invoke Excel within the browser. 
For Netscape users, please save the file first and load it with any 
spreadsheet program.</p>

<?php
DrawFooter();
?>
