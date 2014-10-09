#!/usr/local/bin/perl

package main;
our $SEE;

package Chado_utils;

use strict;
use warnings;
## Database connection modules
use DBI;
use DBD::Pg;
## PASA Gene_obj modules
use FindBin;
use lib qw($FindBin::Bin/../PerlLib);
use Gene_obj;
use Gene_obj_indexer;
use Carp;
use URI::Escape;
use Data::Dumper;

####
sub ConnectToDb {
    my ($server, $dbtype, $user, $passwd, $db) = @_;
    my $dbproc;

    my ($schema, $port) = ("public", "5432");

    # parse $schema from $db name, if formatted like DBNAME:SCHEMA
    ($db, $schema) = split /:/, $db if ($db =~ /:/,);

    print STDERR "server: $server dbtype: $dbtype user: $user passwd: $passwd db: $db\n"
      if $SEE;

    my $dsn .= "dbi:Pg:dbname=$db;host=$server";
    $dsn .= "\;port=$port" if ($port);

    $dbproc = DBI->connect($dsn, $user, $passwd, { AutoCommit => 0 })
      or die $DBI::errstr;

    if (!defined $dbproc) {
        croak "Cannot connect to PostgreSQL server: $DBI::errstr";
    }
    $dbproc->{RaiseError} = 1;    #turn on raise error.  Must use exception handling now.

    ## set schema
    $dbproc->do("SET search_path TO $schema");

    return ($dbproc);
}

sub index_Chado_gene_objs {

    my ($dbproc, $gene_obj_indexer, $organism, $seq_id) = @_;

    # organism is optional.
    # seq_id is optional.

    $organism = "A.thaliana" if (not defined $organism);
    my ($organism_id) =
      $dbproc->selectrow_array("SELECT organism_id FROM organism WHERE abbreviation='$organism'");
    my $contigs =
      $dbproc->selectall_arrayref("SELECT f.feature_id, f.uniquename "
          . "FROM feature f "
          . "INNER JOIN cvterm c ON f.type_id=c.cvterm_id "
          . "WHERE f.organism_id=$organism_id AND (c.name='contig' OR c.name='chromosome')");

    my $hash_mode = 0;
    if (ref $gene_obj_indexer eq 'HASH') {
        $hash_mode = 1;
    }

    ## note can use either a gene_obj_indexer or a hash reference.

    my %gene_coords;
    my %asmbl_id_to_gene_id_list;
    my %transcript_to_gene;
    my %cds_phases;

    my %gene_names;
    my %loci;

    my %gene_id_to_source_type;

    my %source_tracker;

    my $counter = 0;

    foreach my $contig (@$contigs) {
        my ($contig_id, $seqid) = @$contig;
        next if ($seq_id and $seq_id ne $seqid);

        my $features = get_features($dbproc, $contig_id, $seqid);
        my $aliases = get_aliases($dbproc, $contig_id);

        foreach my $feature (@$features) {
            my (
                $feature_id, $source,  $type,   $start,
                $end,        $score,   $strand, $phase,
                $ID,         $Name,    $Parent, $Derives_from,
                $Target,     $t_start, $t_end,  $t_strand
            ) = @$feature;

            next unless ($type =~ /^(gene|mRNA|CDS|exon)$/);

            #fix semantic differences
            $start++;

            #$end++; don't iterate end since it is already correct (coordinates are space based)

            die "ERROR: No strand " . join("\t", @$feature) . "\n"
              if (!$strand);

            if (exists $source_tracker{$ID}
                && $source_tracker{$ID} ne $source)
            {
                confess
"Error, gene ID $ID is given source $source when previously encountered with source $source_tracker{$ID} ";
            }

            if ($type eq "gene") {
                $gene_names{$ID} = $Name;
                $loci{$ID}       = $aliases->{$feature_id};
                next;
            }

            if ($type eq "mRNA") {
                $transcript_to_gene{$ID} = $Parent;
                next;
            }

            my $transcript_id = $Parent;
            my $gene_id       = $transcript_to_gene{$transcript_id};

            unless (defined $gene_id) {
                print STDERR
                  "Error, no gene feature found for $transcript_id.... ignoring feature.\n"
                  if $SEE;
                next;
            }

            $gene_id_to_source_type{$gene_id} = $source;

            my ($end5, $end3) =
              ($strand == 1) ? ($start, $end) : ($end, $start);

            $gene_coords{$seqid}->{$gene_id}->{$transcript_id}->{$type}->{$end5} =
              $end3;

            if ($type eq "CDS") {
                $phase = 0 if (!$phase);
                $cds_phases{$gene_id}->{$transcript_id}->{$end5} = int($phase)
                  if ($phase =~ /^\d+$/);
            }
        }

    }

    print STDERR "\n-caching genes.\n" if $SEE;
    foreach my $asmbl_id (sort keys %gene_coords) {
        my $genes_href = $gene_coords{$asmbl_id};

        foreach my $gene_id (keys %$genes_href) {
            print STDERR "\r-indexing [$gene_id]  " if $SEE;
            my $transcripts_href = $genes_href->{$gene_id};

            my @gene_objs;

            foreach my $transcript_id (keys %$transcripts_href) {

                my $cds_coords_href = $transcripts_href->{$transcript_id}->{CDS}
                  || {};    # could be a noncoding transcript w/ no CDS
                my $exon_coords_href = $transcripts_href->{$transcript_id}->{exon};

                unless (ref $exon_coords_href) {
                    print STDERR Dumper($transcripts_href);
                    die "Error, missing exon coords for $transcript_id, $gene_id\n";
                }

                my $gene_obj = new Gene_obj();

                if (scalar(keys %$cds_coords_href) == 1) {

                    ## could be that only the cds span was provided.
                    ## break it up across the exon segments

                    my ($cds_lend, $cds_rend) =
                      sort { $a <=> $b } %$cds_coords_href;
                    my @exon_coords;
                    my $orient;
                    foreach my $exon_end5 (keys %$exon_coords_href) {
                        my $exon_end3 = $exon_coords_href->{$exon_end5};
                        push(@exon_coords, [ $exon_end5, $exon_end3 ]);
                        if ($exon_end5 < $exon_end3) {
                            $orient = '+';
                        }
                        elsif ($exon_end5 > $exon_end3) {
                            $orient = '-';
                        }
                    }

                    $gene_obj->build_gene_obj_exons_n_cds_range(\@exon_coords,
                        $cds_lend, $cds_rend, $orient);
                }
                else {

                    ## cds and exons specified separately

                    $gene_obj->populate_gene_obj($cds_coords_href, $exon_coords_href);
                }

                $gene_obj->{Model_feat_name} = $transcript_id;
                $gene_obj->{TU_feat_name}    = $gene_id;
                $gene_obj->{asmbl_id}        = $asmbl_id;

                if (my $gene_locus = $loci{$gene_id}) {
                    $gene_obj->{pub_locus} = $gene_locus;
                }
                if (my $transcript_locus = $loci{$transcript_id}) {
                    $gene_obj->{model_pub_locus} = $transcript_locus;
                }

                $gene_obj->{com_name} = $gene_names{$gene_id} || $transcript_id;

                $gene_obj->{source} = $gene_id_to_source_type{$gene_id};

                ## set CDS phase info if available from the gff
                my $cds_phases_href = $cds_phases{$gene_id}->{$transcript_id};
                if (ref $cds_phases_href) {
                    ## set the cds phases
                    my @exons = $gene_obj->get_exons();
                    foreach my $exon (@exons) {
                        if (my $cds = $exon->get_CDS_obj()) {
                            my ($end5, $end3) = $cds->get_coords();
                            my $phase = int($cds_phases_href->{$end5});
                            unless ($phase == 0 || $phase == 1 || $phase == 2) {
                                confess
"Error, should have phase set for cds $gene_id $transcript_id $end5, but I do not know what $phase is. ";
                            }
                            $cds->set_phase($phase);
                        }
                    }
                }

                push(@gene_objs, $gene_obj);
            }

            ## want single gene that includes all alt splice variants here
            my $template_gene_obj = shift @gene_objs;
            foreach my $other_gene_obj (@gene_objs) {
                $template_gene_obj->add_isoform($other_gene_obj);
            }

            $template_gene_obj->refine_gene_object();

            if ($hash_mode) {
                $gene_obj_indexer->{$gene_id} = $template_gene_obj;
            }
            else {
                $gene_obj_indexer->store_gene($gene_id, $template_gene_obj);
            }

            print STDERR "Chado_utils: stored $gene_id\n" if $SEE;

            # add to gene list for asmbl_id
            my $gene_list_aref = $asmbl_id_to_gene_id_list{$asmbl_id};
            unless (ref $gene_list_aref) {
                $gene_list_aref = $asmbl_id_to_gene_id_list{$asmbl_id} = [];
            }
            push(@$gene_list_aref, $gene_id);
        }
    }
    print STDERR "\n" if $SEE;
    return (\%asmbl_id_to_gene_id_list);
}

sub get_features {
    my $dbproc    = shift;    # db handle
    my $contig_id = shift;    # contig feature_id
    my $seqid     = shift;    # sequence id

    #get all features on contig
    my $features =
      $dbproc->selectall_arrayref("SELECT f.feature_id, r.accession, c.name, l.fmin, l.fmax, "
          . "a.significance, l.strand, l.phase, f.uniquename, f.name, f2.uniquename, "
          . "f3.uniquename, f4.uniquename, l2.fmin, l2.fmax, l2.strand "
          . "FROM feature f "
          . "INNER JOIN featureloc l ON (l.srcfeature_id = $contig_id AND f.feature_id = l.feature_id) "
          . "LEFT OUTER  JOIN feature_dbxref d ON d.feature_id = l.feature_id "
          . "INNER JOIN dbxref r ON r.dbxref_id = d.dbxref_id "
          . "INNER JOIN db b ON (b.name = 'GFF_source' AND r.db_id = b.db_id) "
          . "LEFT OUTER JOIN cvterm c ON f.type_id = c.cvterm_id "
          . "LEFT OUTER JOIN analysisfeature a ON a.feature_id = f.feature_id "
          . "LEFT OUTER JOIN feature_relationship pt ON pt.subject_id = f.feature_id "
          . "LEFT OUTER JOIN cvterm c2 ON pt.type_id = c2.cvterm_id "
          . "LEFT OUTER JOIN dbxref d2 ON d2.dbxref_id = c2.dbxref_id "
          . "LEFT OUTER JOIN feature f2 ON (c2.name = 'part_of' AND f2.feature_id = pt.object_id) "
          . "LEFT OUTER JOIN feature f3 ON (c2.name = 'derives_from' AND f3.feature_id = pt.object_id) "
          . "LEFT OUTER JOIN featureloc l2 ON (f.feature_id = l2.feature_id AND l2.srcfeature_id !=  l.srcfeature_id) "
          . "LEFT OUTER JOIN feature f4 ON l2.srcfeature_id = f4.feature_id "
          . "ORDER BY f.feature_id;");

    return $features;
}

sub get_aliases {
    my ($dbproc, $contig_id) = shift;

    return {} if (!$contig_id);

    #get Aliases
    my $aliases =
      $dbproc->selectall_arrayref("SELECT f.feature_id, s.name "
          . "FROM feature f "
          . "INNER JOIN featureloc l ON (l.srcfeature_id = $contig_id AND f.feature_id = l.feature_id) "
          . "INNER JOIN feature_synonym fs on f.feature_id = fs.feature_id "
          . "INNER JOIN synonym s ON fs.synonym_id = s.synonym_id "
          . "WHERE f.uniquename != s.name "
          . "AND f.name != s.name "
          . "AND f.name != (s.name||':'||f.uniquename) "
          . "AND f.uniquename != (s.name||'-'||f.feature_id) "
          . "AND f.uniquename != (s.name||'-'||f.feature_id) "
          . "AND f.uniquename != (s.name||'-'||f.feature_id) ");

    my %results;
    foreach my $alias (@$aliases) {
        my ($id, $nam) = @$alias;
        @{ $results{$id} } = () if (not defined $results{$id});
        push @{ $results{$id} }, $nam;
    }

    return \%results;
}

1;    #EOM
