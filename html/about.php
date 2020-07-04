<?php
  require "shared.web.inc";    
  drawHeader( "MiCA: About MiCA" );
?>

<p>This web site and the tools you find here were developed by students and 
faculty of the Initiative for Bioinformatics and Evolutionary Studies (IBEST) 
at the University of Idaho who are members of the Departments of Computer 
Science and Biological Sciences. The research was funded by the P&amp;G
Company.</p>

<p>The aim of this research was to develop a suite of web-based tools 
that would enable researchers to perform analyses of microbial community 
structure based on terminal-restriction fragment length polymorphisms 
(T-RFLP). MiCA enables researchers to perform the following tasks:</p>

<p>(a). <em>in silico</em> PCR amplification and restriction of 16S rRNA gene 
sequences found in public database;</p>

<p>(b). automatic retrieval of data and archival storage in a relational 
database;</p>

<p>(c). comparison of multiple T-RFLP profiles obtained from a single sample 
using different primer-enzyme combinations;</p>

<p>(d). statistical analysis of T-RFLP data and clustering of samples based 
on similarities and differences.</p>

<p>MiCA currently operates on a server computer with 1.0GB of RAM and four 
Intel Xeon 2.8GHz processor running Linux. The web server runs on Apache, 
which has been fully configured and optimized for better performance and 
security. We have devleoped several PHP scripts to facilitate the virtual 
digest. The primary database contains a large number of 16S rRNA genes 
retrieved from the Michigan State University RDP II database. We have 
gathered the most commonly used primers and restriction enzymes in the 
database. There are 19 forward and 21 reverse primers, and 53 restriction 
enzymes. A digest can incorporate at most three restriction enzymes. The 
output is available in CSV (Common Separated Values) and plain-text 
formats. CSV file format allow data to be easily retrieved into 
spreadsheet or database programs. Majority of spreadsheet applications, 
such as Linux Gnumeric or Microsoft Excel, can read and write CSV files. 
Most web browsers will invoke an appropriate application if one has been 
installed. Users may also save the file first then convert it to a desired 
format.</p>

<p>We hope you will find MiCA useful. Please <a class="email" 
href="mailto:foster@uidaho.edu">email us</a> any comments, suggestions 
or questions.</p>

<?php DrawFooter(); ?>
