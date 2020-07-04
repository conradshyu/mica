<?php
  require "shared.help.inc";    
  DrawHeader( "MiCA: Primer Sequence Prevalence Analysis" );
?>

<p>The ability to reproduce a target section of a DNA sequence through the use of PCR has facilitated a
wide array of amplification techniques. The success of the PCR strategy heavily depends on the small
synthetic oligonucleotide that hybridize to the complementary DNA sequence. An oligonucleotide primer is
a short segment of nucleotides that is used to amplify a section of DNA template during the PCR reaction.
Optimal primer sequences and appropriate concentration are essential for maximal specificity and
efficiency in PCR. Primers are annealed to the denatured DNA template to provide an initiation site for
the elongation of the new DNA molecule. Primers can either be specific to a particular DNA nucleotide
sequences that are very common in a particular set of DNA molecules. PCR primer design is usually based
on reverse translation of multiply aligned sequences across the conserved regions of proteins. Various
techniques have been applied to this problem, but frequent failure to amplify a desired target sequences
are often attributable to inadequate design. Primer design can be very difficult because of codon
dependency and the additional degeneracy needed to represent multiple codons at a position in the 
alignment. These degeneracy lead to complications in trying to find suitable annealing temperature and
primer lengths.</p>

<?php DrawFigure( "PSPA web interface", "images/mica_pspa.png" ); ?>

<p><b>Figure 1:</b> PSPA utilizes a simplified version of ISPaR search engine, and permits up to three
forward and three reverse primers. All possible combinations of forward and reverse primers will be used
to search the database, and generate the analysis output. The analysis of primer sequence prevalence is
considerably more computationally intensive. An analysis can take as long as 25sec, for example.</p>

<?php DrawFigure( "Output page for PSPA", "images/result_pspa.png" ); ?>

<p><b>Figure 2:</b> All possible combinations of forward and reverse primers will be used for the analysis.
PSPA reports the number of matches for given pair of forward and reverse primers. It also summarizes the
average number of mismatches for each pair of primers.</p>

<p>The development of PSPA (primer sequence prevalence analysis) (see Fig. 1 for an illustration) aims to
meet the challenge of effective primer selection and provide high-throughput <em>in silico</em> analysis
of primer combinations. PSPA utilizes a simplified version of ISPaR search algorithm, which only 
identifies the forward and reverse primer binding sites. It will accept as many as three forward and three
reverse primers. All possible combinations of given forward and reverse primers will be used for the
analysis, and the same parameter settings for mismatches will be applied uniformaly to each combination.
For example, if two forward and three reverse primers were chosen, then a total of six possible 
combinations will be used for the analysis. The analysis of primer sequence prevalence is considerably more
computationally intensive. AS an example, an analysis of three forward and three reverse primers can take
as longas 25sec. PSPA reports the number of successful amplifications for each primer individually (forward
and reverse) and for each primer pair collectively (see Fig. 2 for an illustration). It is worth noting
that the number of successful amplifications on a given forward primer will not always be the same if it is
paired with different reverse primers. This phenomenon is largely due to the nature of the search algorithm
and locations of primer binding sites. It is assumed that primer binding sites must be sufficiently distant
in order to maintain biological relevance. PSPA currently does not report the phylogenetic information on
the amplified sequences.</p>

<?php DrawFigure( "Top 10 primer pairs with largest numbers of amplicons", "images/top_primer.png" ); ?>

<p><b>Figure 3:</b> No mismatches were allowed during the searches. Therefore, the results should be 
considered as a conservative estimate of the primer sequence prevalence. The database contains a toal of
210,976 bacterial sequences, which were retrieved from the hierarchy broswer in the April 2006, update 39,
release of the RDP database.</p>

<?php DrawFigure( "Identical forward primer with different numbers of amplicons",
"images/forward_519f.png" ); ?>

<p><b>Figure 4:</b> Identical forward primers can produce different number of amplicons if paired with
different reverse primers. This phenomenon also holds for the reverse primers. It is due to the fact that
both primers are searched simultaneously, and primer binding sites must maintain a certain distance in
order to be biologically significant.</p>

<p>Selection of PCR primers is a crucial step to the production of effective restriction fragments for the
analysis of microbial diversity and species composition. For a successful amplification, there must be at
least two conserved regions in the gene sequence of interest to provide priming sites. In addition, the
primers must be far enough apart for sufficient sequence divergence to exist between them since amplicons
that are too short will result in patterns that do not reflect the true diversity of the sample because
the majority of amplicons simply will not contain a restriction site. With the avaiable analytical tool
for primer sequence, a comprehensive study on the prevalence has been performed. A total of 168 possible
pairs of forward and reverse primers (12 forward and 14 reverse primers) were used for the study. The
search for amplification sites did not permit any mismatches within both primers. Therefore, the results
should be considered as a conservative estimate of the primer sequence prevalence. If mismatches were
permitted, then the number of amplicons will increase dramatically. Fig. 3 shows the top 10 primer pairs
that give the most coverage on the sequences in the database. The primer pair, Eub-341f and Eub-536f,
amplifies nearly a half of all sequences in the database and gives a total of 54,886 amplicons. It is
worth noting that the majority of primer pairs give similar numbers of amplicons, approximately 33,000.
The next figure, Fig. 4, shows that identical forward primers may give a different number of amplicons if
paired with different reverse primers. This is due to the fact that both primers are searched
simultaneously. This approach ensures that primers are sufficiently distant.</p>

<?php DrawFooter(); ?>
