package Chado_annot_retriever;

use strict;
use warnings;
use Carp;
use FindBin;
use lib qw($FindBin::Bin);
use lib qw($FindBin::Bin/../PerlLib);
use base qw(PASA_UPDATES::Pasa_latest_annot_retrieval_adapter);
use Chado_utils;
use Gene_obj_indexer;

## This sample annotation retriever is responsible for
## parsing gene structures from a chado database.

#### static method:
sub get_annot_retriever {
    ## the hook is always called with a hashref for the configuration values, followed by the custom parameter.

    my ($config_href, $chado_db_info) = @_;

    unless ($chado_db_info) {
        confess "Error, need chado database info: $chado_db_info";
    }

    my ($server, $type, $user, $password, $chado_db, $pred_type) =
      split(/,/, $chado_db_info);
    my ($dbproc) = Chado_utils::ConnectToDb($server, $type, $user, $password, $chado_db); #read only

    my $gene_obj_indexer = new Gene_obj_indexer({ create => "$chado_db.inx" });
    &Chado_utils::index_Chado_gene_objs($dbproc, $gene_obj_indexer);

    $dbproc->disconnect();

    my $annot_retriever_obj = new Chado_annot_retriever($gene_obj_indexer);

    return ($annot_retriever_obj);

}

## object methods below

####
sub new {
    my $packagename = shift;

    my ($gene_obj_indexer) = @_;

    my $self = $packagename->SUPER::new();
    $self->{gene_obj_indexer} = $gene_obj_indexer;

    ## get contig to gene list
    my $contig_to_genelist_aref = $self->{contig_to_genelist} = {};

    my @gene_ids = $gene_obj_indexer->get_keys();

    foreach my $gene_id (@gene_ids) {

        my $gene_obj = $gene_obj_indexer->get_gene($gene_id);

        my $contig_id = $gene_obj->{asmbl_id}
          or die "Error, gene_obj lacks asmbl_id value: " . $gene_obj->toString();

        my $gene_list_aref = $contig_to_genelist_aref->{$contig_id};
        unless (ref $gene_list_aref) {
            $gene_list_aref = $contig_to_genelist_aref->{$contig_id} = [];
        }

        push(@$gene_list_aref, $gene_id);
    }

    return ($self);
}

####
sub retrieve_gene_models_for_contig {
    my $self = shift;
    my ($contig_id) = @_;

    my @gene_objs;

    my $gene_list_aref = $self->{contig_to_genelist}->{$contig_id};
    if (ref $gene_list_aref) {
        foreach my $gene_id (@$gene_list_aref) {
            my $gene_obj = $self->{gene_obj_indexer}->get_gene($gene_id);
            unless (ref $gene_obj) {
                confess "Error, couldn't retrieve gene object for gene_id: $gene_id";
            }
            push(@gene_objs, $gene_obj);
        }
    }
    else {
        print STDERR "Warning, no genes retrieved for contig_id: $contig_id\n";
    }

    return (@gene_objs);
}

1;    #EOM

