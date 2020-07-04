<?php
  require "shared.help.inc";    
  DrawHeader( "MiCA: Enzyme Resolving Power Analysis" );
?>

<p>A restriction enzyme acts like a molecular scissor that cleaves double-stranded DNA at all sites
in the molecule where the base sequence matches a particular short nucleotide sequence called the 
restriction site. By and large, the recognition site for a given restriction enzyme is unique to 
that particular enzyme. However, a few DNA sequences are recognized by different restriction enzymes
(called isoschizomers) isolated from different organisms. As an example, <em>Sau</em>3A (from 
<em>Staphylococcus aureus</em>) and <em>Mbo</em>I (from <em>Moraxella bovis</em>) both cleave the 
same sequence, but vary in their ability to cut methylated DNA. In addition to isoschizomers, the 
recognition site of one restriction enzyme may contain the recognition sequence of another. The 
discovery of restriction enzymes came from the study of host restriction and modification of phages. 
More than 1,000 different restriction enzymes have been purified from hundreds of difference 
microorganisms. In nature, restriction enzymes function as microbial immune systems that protect 
bacteria from infection by viruses. The presence of additional enzymes generally shows a 
multiplicative effect. That is, a cell with multiple independent restriction enzymes could be
virtually impregnable. The strict nucleotide sequence specificity dictates the utilization of 
restriction enzymes. In evolutionary studies, for example, restriction sites are commonly used to 
determine if DNA sequences vary between isolates. The difference in restriction fragment lengths
suggests the presence of modifications in the DNA molecules. The cut made by restriction enzymes is 
usually staggered such that the two strands of the double helix are cut a few bases apart. This 
creates single-stranded overhangs called sticky-ends at either end of the digested DNA molecule (see 
Fig. 1 for an illustration). The overhangs are said to be sticky because they are able to bind to a 
complementary single-stranded region. The sticky-ends that remain after restriction enzyme cleavage 
are actually a key aspect of recombinant DNA technology. Cleavage of DNA from any organisms with a
particular restriction enzyme always generates the same complementary ends on either side of the cut. 
The sticky-ends of a fragment from one source of DNA can then be joined to complementary ends of DNA 
cut from other source. Some restriction enzymes are also capable of producing fragments without 
single-stranded overhangs. This commonly refers to as a blunt-end (see Fig. 2 for an illustration).
As opposed to sticky cutters, no complementarities are required.</p>

<?php DrawFigure( "Sticky ends of restriction sites", "images/sticky_end.png" ); ?>

<p><b>Figure 1:</b> The restriction enzyme, <em>Eco</em>RI (5'-G^AATTC-3'), causes breaks in the
backbone of the DNA sequence and leaves a staggered cut because it results in fragments with 
single-stranded ends. Sticky ended fragments can be easily ligated to other sticky ended fragments 
with compatible single-stranded overhangs, resulting in efficient cloning.</p>

<?php DrawFigure( "Blunt ends of restriction sites", "images/blunt_end.png" ); ?>

<p><b>Figure 2:</b> The restriction enzyme, <em>Eco</em>RV (5'-GAT^ATC-3'), produces fragments without
single-stranded overhangs. This commonly refers to as a blunt end. Blunt ended fragments usually ligate
much less efficiently, making cloning more difficult. However, any blunt ended fragment can be ligated
to any other. So, blunt cutting enzymes are used when compatible sticky ended fragments cannot be
generated.</p>

<p>Restriction endonucleases that occur naturally in bacteria as host defense systems have been
widely used in laboratory applications from recombinant DNA technology to polymorphism detection for
diagnostics. The specificity of a restriction enzyme is made possible by its molecular structures.
Since the ability to resolve phylogenetically distinct populations depends on the existence of 
restriction site polymorphisms among the populations, the choice of restriction enzymes used is 
therefore critical. While it is possible to determine the optimal enzyme(s) based on empirical data,
studies have also suggested that it is also possible to make rational choices based on prior 
knowledge that permits one to predict the phylogenetic groups likely to be present in the habitat. 
Brunk et al. (1996) looked at the distribution of predicted fragments from 16S rRNA genes in a 
sequence database and recommended the use of <em>Hha</em>I, <em>Msp</em>I, <em>Rsa</em>I, and a 
combined digest of both <em>Rsa</em>I and <em>Hha</em>I. Dunbar et al. (2001) analyzed the 
phylogenetic resolution of restriction fragments from a range of enzymes and enzyme combinations.
They found that 68 percent of the <em>Rsa</em>I generated restriction fragments were specific for
less than four species of the same genus. They also concluded that the phylogenetic specific of any
restriction fragments would be greatly enhanced by the use of group specific primers. To explore the 
specificities of restriction enzymes <em>in silico</em> with a database containing a large number of 
sequences, a fast search algorithm that can accommodate primers and multiple restriction enzymes has
been developed. The enzyme resolving power analysis (ERPA) aims to provide such a tool on MiCA.</p>

<?php DrawFigure( "ERPA web interface", "images/mica_erpa.png" ); ?>

<p><b>Figure 3:</b> ERPA only accepts one pair of forward and reverse primers. The search algorithm
will iterate through the entire list of restriction enzymes in the database and reports the numbers
of unique forward and reverse restriction fragments. The parameter settings for mismatches are applied
to all possible combinations.</p>

<?php DrawFigure( "Output page for ERPA", "images/result_erpa.png" ); ?>

<p><b>Figure 4:</b> The analysis reports the number of unique forward and reverse fragments. Several
descriptive statistics, such as means and standard deviations, are also calculated. The outputs can be
sorted in various ways and are available in both the plain-text and CSV formats.</p>

<p>ERPA (enzyme resolving power analysis) provides a high performance analytical tool for the 
selection and analysis of restriction enzymes (see Fig. 3 for an illustration). It requires a pair
of forward and reverse primers in order to produce the amplicons, and then the search algorithm 
iterates through the entire lest of restriction enzymes available in MiCA. Currently there are 46 
commonly used restriction enzymes in the database. The search algorithm is the same as that of ISPaR.
ERPA reports the number of unique forward and reverse restrition fragments digested with each enzyme.
Several descriptive statistics, such as means and standard deviations, are also reported. Detailed
fragment lengths, however, are not shown on the output pages. Ideally, with a properly chosen pair of
forward and reverse primers, the restriction enzyme that gives the maximum number of unique fragments
is able to better resolve phylogenetically distinct species. Molecular markers are generally considered
constant landmarks in the genome, and are found at specific locations of the sequences. There are many
molecular techniques that rely on the restriction fragment polymorphism for the identification of 
species. Restriction fragment length polymorphism (RFLP), for example, is based on detection of 
variations in restriction fragment length of specific DNA sections among individual, as determined by
the presence or absence of restriction sites. It is then possible to locate a gene by looking for RFLP
that are almost always inherited with it. As another example, random amplification of polymorphic DNA
(RAPD) are DNA fragments, which are polymorphic in size, generated by PCR using one or two randomly 
selected primers. They are dominant markers and RAPD patterns can be used in strain identification (one
of the methods of DNA fingerprinting). In AFLP, DNA treated with restriction enzymes is amplified using
PCR. The method allows selective amplification of restriction fragments, giving rise to large amount
of useful markers that can be located on the genome relatively quickly and reliably.</p>

<?php DrawFigure( "Histogram of pooled termini digested with all restrition enzymes",
"images/histogram.png" ); ?>

<p><b>Figure 5:</b> The histogram clearly shows that the majority of fragments are less than 50bp. A
large number of forward and reverse fragments are densely clustered into few groups. It is interesting
to note that the forward and reverse terminal fragments show relatively similar distribution.</p>

<p>A comprehensive analysis on the primer sequence prevalence showed that the forward primer, Eub-341f
(5'-CCTACGGGAGGCAGCAG-3'), and the reverse primer, Eub-536f (5'-GWATTACCGCGGCKGCTG-3'), give the best
coverage on the bacterial sequences in the database. No mismatches were allowed during the searches,
and a total of 54,886 successful amplifications were reported. The analysis was performed without
incorporating any phylogenetic information from the sequences. Fig. 5 shows the histogram of the pooled
forward and reverse terminal fragments, between 50 and 200bp, digested with all available restriction
enzymes. Isoschizomers were grouped together and reported only once. The histogram shows that a large
number of forward and reverse fragmens are less than 50bp, and the majority is densely clustered into
few groups. The analysis also showed that fragments larger than 200bp are rare. Specifically 211,508
forward and 274,870 reverse fragments that are less than 50bp were reporeted, and only 719 forward and
700 reverse fragments are larger than 200bp. Both forward and reverse fragments show similar
distribution patterns.</p>

<p>The distribution of unique fragments dictates the resolving power of the restriction enzymes. 
Theoretically the restriction enzymes that give the largest number of unique fragments are more likely
to resolve phylogenetically distinct species. The resolution of a restriction enzyme depends on its
potential to differentiate species from a set of sequence variants based on the distribution of terminal
fragment sizes. Although ERPA provides a useful tool for high-throughput analysis of enzyme resolving
power, however, the specific restriction enzymes used may have to be empirically determined, especially
in specialized habitat represented by limited numbers of phylotypes. The techniques involve restriction
fragments are likely to be a very valuable screening tool when spatiotemporal changes in natural
communities with relatively low to intermediate species richness are studied. However, these techniques
still may not be sufficiently accurate for characterizing complex microbial populations.</p>

<?php DrawFooter(); ?>
