<?php
  require "shared.help.inc";    
  DrawHeader( "MiCA: T-RFLP Analysis (APLAUS+)" );
?>

<p>Since the ability to resolve phylogenetically distinct populations 
depends on the existence of restriction site polymorphisms among the 
populations, the choice of restriction enzyme(s) used is critical. When an 
investigator can determine the optimal enzyme(s) based on empirical data, 
it is also possible to make rational choices based on prior knowledge that 
permits one to predict the phylgenetic groups likely to be present in the 
habitat. The interpretation and comparison of profiles obtained from 
multiple samples is a difficult process since single fragments are rarely 
unique for a single phylotype. Ideally, investigators would like to 
determine what known and unknown organisms could produce the observed 
community profile. To do so, one would need to search a database of DNA 
sequences from known organisms and perform <em>in silico</em> PCR amplification and 
restriction digest of the genes to determine which ones could produce the 
fragments seen in the profile. The algorithm used for this purpose must 
also account for the possible presence of phylotypes in the sample that 
are not presented in the database. For this analysis a list of "plausible 
communities" could be generated.</p>

<?php DrawFigure( "The T-RFLP analysis (APLAUS+) homepage",
    "images/mica_aplaus.png" ); ?>

<p><b>Figure 1: </b>This figure shows the web interface for the T-RFLP 
analysis using APLAUS+. The interface is very similar to that of ISPaR but 
requries an additional sample file for the analysis. In addition, two more 
options are presented to correct the size calling discrepancies (migration 
shifts on the forward and reverse fragments) on most CE systems.</p>

<p>APLAUS+ uitilizes the CSV (Comma Separated Values) format. The CSV file 
format is often used to exchange data between disparate applications. The 
file format, as it is used in Microsoft Excel, has become a widely 
accepted standard in the industry. A CSV file is a specially formatted 
plain text file, which stores spreadsheet or basic database-style 
information in a very simple format, with one record on each line, and 
each field within that record separated by a comma. CSV files are often 
used as a simple way to transfer a large volume of spreadsheet or database 
information between programs, without concerning about special file types. 
The following figure show how to create a CSV formatted sample file for 
APLAUS+.</p>

<?php DrawFigure( "The sample file must be written in the CSV format",
    "images/csv_format.png" ); ?>

<p><b>Figure 2: </b>The CSV format is often used to exchange data among 
different programs. The first column must be the abudnance of the 
fragments and second the fragment sizes. Forward and reverse fragments 
must be saved in separate files.</p>

<p>The APLAUS+ algorithm compares the data from one ore more T-RFLP profiles 
to the outcomes(s) of <em>in silico</em> analyses of sequences in the 
database done using the same primers and restriction enzymes. The outcome 
is one or more plausible community structures, one(s) that in theory would 
produce the T-RLFP profiles that have been observed. For a particular 
combiantion of enzymes and primers two or more phylotypes might 
theoretically produce exactly the same size fragments and cannot be 
distinguished. In these instances, APLAUS+ combines the phylotypes into a 
single group.</p>

<?php DrawFigure( "The APLAUS+ output page lists plausible microbial structure",
    "images/result_aplaus.png" ); ?>

<p><b>Figure 3: </b>The APLAUS+ lists one or more plausible community 
structures that in theory would produce the T-RFLP profiels that have been 
observed. The output can be sorted based on different criteria.</p>

<?php DrawFooter(); ?>
