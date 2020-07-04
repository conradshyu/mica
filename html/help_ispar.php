<?php
  require "shared.help.inc";
  DrawHeader( "MiCA: Virtual Digest (ISPaR)" );
?>

<p>Culture-based assessments of bacterial community composition are heavily biased towards readily
cultivable and fast-growing organisms and therefore only useful for limited kinds of communities.
Therefore, culture-independent methods based on analyses of 16S and 18S rRNA gene sequences have been
used in an effort to overcome these limitations, and have been widely employed to explore microbial
diversity and understand community dynamics. Culture-independent methods commonly utilize the genome
structure and genetic sequence as the major distinguishing characteristics of the family, type, and
strain of microorganism. Specific strains of microorganisms can be distinguished on the basis of their
DNA or RNA or by the DNA fragments produced when the DNA is cleaved by specific restriction enzymes.
The DNA sites recognized by restriction enzymes differ in their sequence, length, and frequency of
occurrence. As a result, different restriction enzymes cleave the DNA of a sample in different places,
yield fragments of different lengths. The cleavage of different DNA samples with one restriction enzyme
can also yield fragments of many different lengths. The differences strains of a specific organism
produced on cleavage with one or more restriction enzymes are called restriction fragment length
polymorphisms (RFLP).</p>

<p>T-RFLP uses properties of the sequences of the gene being amplified as a proxy for actual sequence
information, and provides a community profile (or fingerprint) that reflects the composition of the
numerically dominant populations in a sample. Information about specific groups of organisms can be
obtained by using group-specific PCR primers. T-RFLP takes advantage of restriction site polymorphisms
in PCR amplified rRNA genes to produce fluorescently labeled terminal restriction fragments that differ
in length. These are then resolved by either using acrylamide gel electrophoresis or by capillary
electrophoresis. DNA fragments that differ in size correspond to phylogenetically distinct populations
in the community. Since differences in the sizes of the fragments produced are solely the result of
difference in the primary sequences of rRNA genes, T-RFLP data can be interpreted using a database of
known sequences. To explore and simulate experiments with all possible recognition sites using a database
containing large numbers of sequences, a fast search algorithm that can accommodate primers with multiple
redundant positions is needed. ISPaR has been developed to meet such a challenge.</p>

<p>ISPaR performs <em>in silico</em> amplification and restriction of 16S and 18S rRNA sequences. The
query page lists all available primers and enzymes in drop-down menus. The location and number of allowed
mismatches between primers and target sequences are then stipulated to modulate the breadth and sensitivity
of the search. The output from ISPaR can be sorted in various ways. In addition, the data from ISPaR can
be exported in CSV (Comma Separated Values) and plain-text formats to other spreadsheet or database
programs. ISPaR utilizes the binary encoding scheme for the primer sequences. The algorithm first delimits
or double-digests the template sequence by the forward and reverse primers. Once the primer binding sites
are located, ISPaR then searches for the restriction sites with the same binary encoding scheme. ISPaR
permits up to three restriction enzymes. The most significant feature on the algorithm is the permission
of mismatches on the primer sequences. Ambiguities on restriction enzymes are also supported. To increase
the search sensitivity, restriction enzymes can contain any number of ambiguities. It is important to note
that the increase number of ambiguity in the restriction enzyme promotes the probability of false positive
matches.</p>

<?php DrawFigure( "The MiCA virtual digest (ISPaR) homepage", "images/mica_ispar.png" ); ?>

<p><b>Figure 1: </b>This figure shows the web interface for the virtual digest (ISPaR). The most
significant feature of the algorithm is the permission of mismatches within the primer sequences, as well
as the entire database.</p>

<?php DrawFigure( "The virtual digest (ISPaR) output page", "images/result_ispar.png" ); ?>

<p><b>Figure 2: </b>The ISPaR output page shows the forward and reverse fragment sizes in the order that
the restriction enzymes were submitted. Additional information about the organisms such as the GenBank
accession numbers, locuses, and names are also available. The output is written in three different formats
and avaiable for download for further analysis.</p>

<?php DrawFooter(); ?>
