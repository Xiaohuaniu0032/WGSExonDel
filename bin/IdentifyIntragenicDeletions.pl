#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

my ($input, $out, $help);
GetOptions(
    'Input=s' => \$input,
    'Out=s'   => \$out,
    'Help'    => \$help,
) or die usage();
if ($help) { print usage(); exit 0; }
die usage() unless defined $input && defined $out;
$input = abs_path($input) if -e $input;
die "[ERROR] Input file not found: $input\n" unless defined $input && -s $input;
die "[ERROR] Input path must be absolute\n" unless $input =~ m{^/};
die "[ERROR] Output path must be absolute\n" unless $out =~ m{^/};
make_path(dirname($out)) unless -d dirname($out);

open my $IN, '<', $input or die "[ERROR] Cannot open $input: $!\n";
my $header_line = <$IN>;
die "[ERROR] Empty input: $input\n" unless defined $header_line;
chomp $header_line; $header_line =~ s/\r$//;
my @header = split /\t/, $header_line, -1;
my %idx = map { $header[$_] => $_ } 0..$#header;
my @required = qw(SampleID Read_Name Flag Primary_Chrom Primary_Start Primary_End Primary_CIGAR Primary_MAPQ Primary_Strand SA_Chrom SA_Start SA_End SA_CIGAR SA_MAPQ SA_Strand SA_NM);
for my $col (@required) { die "[ERROR] Missing input column: $col\n" unless exists $idx{$col}; }

open my $OUT, '>', $out or die "[ERROR] Cannot write $out: $!\n";
my @out_header = qw(
    Event_ID SampleID Read_Name Chrom Start End Size_bp Strand
    Primary_Chrom Primary_Start Primary_End Primary_CIGAR Primary_MAPQ Primary_Strand
    SA_Chrom SA_Start SA_End SA_CIGAR SA_MAPQ SA_Strand SA_NM
);
print {$OUT} join("\t", @out_header), "\n";

my %stats = (input_rows=>0, different_chromosome=>0, different_strand=>0,
             unsupported_primary_cigar=>0, invalid_interval=>0, events_written=>0);
while (my $line = <$IN>) {
    chomp $line; $line =~ s/\r$//; next if $line eq '';
    $stats{input_rows}++;
    my @v = split /\t/, $line, -1;
    my %r = map { $_ => defined $v[$idx{$_}] && $v[$idx{$_}] ne '' ? $v[$idx{$_}] : 'NA' } @required;

    if ($r{Primary_Chrom} ne $r{SA_Chrom}) { $stats{different_chromosome}++; next; }
    if ($r{Primary_Strand} ne $r{SA_Strand}) { $stats{different_strand}++; next; }
    unless (is_positive_integer($r{Primary_Start}) && is_positive_integer($r{Primary_End})
            && is_positive_integer($r{SA_Start}) && is_positive_integer($r{SA_End})
            && valid_cigar($r{Primary_CIGAR}) && valid_cigar($r{SA_CIGAR})) {
        $stats{unsupported_primary_cigar}++; next;
    }

    my ($del_start, $del_end);
    if (cigar_has_terminal_softclip($r{Primary_CIGAR})) {
        # A terminal soft clip means the primary alignment is on the left.
        $del_start = $r{Primary_End} + 1;
        $del_end   = $r{SA_Start} - 1;
    }
    elsif (cigar_has_leading_softclip($r{Primary_CIGAR})) {
        # A leading soft clip means the SA alignment is on the left.
        $del_start = $r{SA_End} + 1;
        $del_end   = $r{Primary_Start} - 1;
    }
    else {
        $stats{unsupported_primary_cigar}++; next;
    }

    unless ($del_start >= 1 && $del_end >= $del_start) {
        $stats{invalid_interval}++; next;
    }

    my $event_id = sprintf('SR_EVENT_%06d', $stats{events_written} + 1);
    my @values = (
        $event_id, $r{SampleID}, $r{Read_Name}, $r{Primary_Chrom},
        $del_start, $del_end, $del_end - $del_start + 1, $r{Primary_Strand},
        @r{qw/Primary_Chrom Primary_Start Primary_End Primary_CIGAR Primary_MAPQ Primary_Strand
              SA_Chrom SA_Start SA_End SA_CIGAR SA_MAPQ SA_Strand SA_NM/}
    );
    print {$OUT} join("\t", map { defined $_ && $_ ne '' ? $_ : 'NA' } @values), "\n";
    $stats{events_written}++;
}
close $IN;
close $OUT or die "[ERROR] Cannot close $out: $!\n";
for my $key (qw/input_rows different_chromosome different_strand unsupported_primary_cigar invalid_interval events_written/) {
    print STDERR "[INFO] $key: $stats{$key}\n";
}
print STDERR "[INFO] Deletion events: $out\n";
exit 0;

sub cigar_ops {
    my ($cigar) = @_;
    return unless defined $cigar && $cigar ne '' && $cigar ne '*';
    my @ops;
    my $pos = 0;
    while ($cigar =~ /\G(\d+)([MIDNSHP=X])/gc) {
        push @ops, [$1, $2];
        $pos = pos($cigar);
    }
    return unless $pos == length($cigar) && @ops;
    return @ops;
}
sub valid_cigar { my @ops = cigar_ops($_[0]); return @ops ? 1 : 0; }
sub cigar_has_terminal_softclip {
    my @ops = cigar_ops($_[0]);
    return 0 unless @ops;
    return $ops[-1][1] eq 'S' && $ops[0][1] ne 'S' ? 1 : 0;
}
sub cigar_has_leading_softclip {
    my @ops = cigar_ops($_[0]);
    return 0 unless @ops;
    return $ops[0][1] eq 'S' && $ops[-1][1] ne 'S' ? 1 : 0;
}
sub is_positive_integer { return defined $_[0] && $_[0] =~ /^\d+$/ && $_[0] > 0; }
sub usage { return "Usage: perl IdentifyIntragenicDeletions.pl --Input split_reads.tsv --Out split_deletion_events.tsv\n"; }
