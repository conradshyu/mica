<?php
  require "shared.help.inc";
  DrawHeader( "MiCA: Introduction" );
?>

<p>Microbial Community Analysis (MiCA), is a web-based suite of analytical tools comprising a database
of 16S rRNA sequences derived from the Ribosomal Database Project (RDP-II), and two computational tools,
ISPaR and APLAUS+, were developed to facilitate the analyses of microbial community structure from T-RFLP
data. MiCA requires no additional software to be installed on the client machines except a web browser
with which to access the intuitive, system independent user interface. The web interface permits the
selection of a pair of forward and reverse primers, and up to three restriction enzymes from the
drop-down menus. Alternatively, users may also enter their own primer and enzyme sequences, and stipulate
the location and number of allowed mismatches between a primer and target sequence. The search algorithm,
ISPaR, performs <em>in silico</em> amplification and restriction of 16S and 18S rRNA sequences. It
utilizes a binary-encoding scheme for DNA sequences and then search for sequences within genes that
anneal to the specified primers. Restriction sites are identified within the gene fragment delimited by
the primers. The number of base pairs between the 5' and 3' primers and the proximal restriction sites
are reported. The APLAUS+ program, on the other hand, infers plausible microbial communities from T-RFLP
data. It compares the data from one or more T-RFLP profiles to the outcomes of <em>in silico</em>
analyses of sequences in the database done using the same primers and enzymes. The outcome consists of
one or more plausible community structures, which, in theory, they would produce the T-RFLP profiles that
have been observed. For a particular combination of enzymes and primers two or more phylotypes might
theoretically produce exactly the same size fragments and cannot be distinguished. In these instances,
APLAUS+ combines the phylotypes into a single group.</p>

<p>MiCA currently includes 12 forward and 14 reverse primers, and 46 restriction enzymes. In addition, it
also features a database of 210,005 16S ribosomal RNA sequences that were obtained from the hierarchy
browser of the May release version (update 39) of the <a class="email" href="http://rdp.cme.msu.edu/">RDP</a>.
Several PHP scripts were developed to extract the identities of the organisms and validate the sequences by
replace any non-nucleotide character with the character N or X. The PHP scripts then store the extracted
and verified data into a plain-text format. The plain-text format permits easy modification and maintenance,
and quick retrieval of the sequence data. MiCA operates on a server computer with 1GB of memory and four
Intel Xeon 2.8GHz processors running Linux. The system utilizes on a customized Apache web server, which
has been optimized for better performance and security. In addition, several PHP scripts were developed to
establish an interface between the web pages and search programs. The PHP scripts are used to forward the
search parameters such as the primers and restriction enzymes to the search programs. The search program
then retrieves all sequences successively and converts sequence characters into binary streams as the search
progresses. Moreover, search programs were written entirely in C++ for fast execution and utilization of
modern processor technologies. The output from the search algorithm is sorted in CSV (Comma Separated Values)
and plain-text formats, which are easily assessable to spreadsheet or database programs. To effectively
search a large number of sequences, MiCA features an algorithm that utilizes a binary encoding scheme that
transforms nucleic acid characters into four parallel bit streams. The binary encoding scheme enables the
search to be performed with a series of binary operations, which effectively reduces the search complexity,
and increases performance. This compact representation also permits an accurate search, which is not
possible with the character-based approaches. In addition, MiCA has an online help system that provides
detailed information and instructions on the use of the computational tools available through MiCA as well
as brief tutorials to guide the selection of primers and enzymes. The MiCA online help system runs on a
separate window, and does not interfere with any tasks running through the web site.</p>

<p>We hope you will find MiCA useful for you studies in microbial diversity. Please
<a class="email" href="mailto:foster@uidaho.edu">email us</a> your comments or suggestions. Thank you.
</p>

<?php DrawFooter(); ?>
