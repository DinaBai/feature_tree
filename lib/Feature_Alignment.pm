package Feature_Alignment;
use strict;
use warnings;
use List::Util qw(max);

sub new {
    my ($class, $fh) = @_;
    print(STDERR "in Feature_Alignment::new\n");
    print(STDERR " class = $class\n");
    print(STDERR " args = ", ", ".join(@_), ".\n");
    my $self = {};
    bless $self, $class;
    $self->{_seqs} = {};
    $self->{_annot} = {};
    $self->{_ids} = [];
    $self->{_is_alinged} = 0;
    $self->{_length} = 0;
    $self->{_format} = '';
    if (defined $fh) {
        print STDERR "File handle = $fh\n";
        $self->read_file($fh);
    }
    return $self;
}

sub get_ntaxa { my $self = shift; return scalar(@{$self->{_ids}})}
sub get_length { my $self = shift; return $self->{_length}}
sub is_aligned { my $self = shift; return $self->{_is_aligned}}

sub detect_format {
    my $class = shift;
    my $fh = shift;
    print STDERR "in detect_format, class=$class, fh=$fh\n";

    $_ = readline $fh;
    seek $fh, 0, 0; # reset file to beginning

    print STDERR "first line of file is :\n", $_, "\n";
    my $format = 'unknown';
    $format = 'clustal' if (/^CLUSTAL/ || /^MUSCLE/);
    $format = 'fasta' if (/^>/);
    $format = 'phylip' if (/^(\d+)\s+(\d+)\s*$/);
    $format = 'nexus' if (/\#NEXUS/);
    print STDERR "input format detected as $format\n";
    return $format
}

sub read_file {
    my $self = shift;
    my $fh = shift;
    print STDERR "in detect_format, self=$self, fh=$fh\n";
    my $format = shift;
    if (! defined $format)
    {
        $format = $self->detect_format($fh)
    }
    $self->{_format} = $format;
    if ($format eq 'unknown') {
        return undef}

    if ($format eq 'clustal') {
        my $found = 0;
        while (<$fh>) {
            if (/^CLUSTAL/ || /^MUSCLE/) {
                $found = 1;
                last;
            }
        }
        die "Format seems to be wrong, not Clustal.\n" if (!$found); 
        while (<$fh>) {
            if (/^(\S+)\s+(\S+)/) {
                push(@{$self->{_ids}}, $1) unless exists($self->{_seqs}{$1});

                $self->{_seqs}{$1} .= $2;
            }
        }
    }
    elsif ($format eq 'phylip') {
        $_ = <$fh>;
        die "Format does not seem to be phylip\n" if (!/^\s*(\d+)\s+(\d+)\s*$/);
        my $ntaxa = $1;
        my $nchar = $2;
        for my $i (1..$ntaxa) {
            $_ = <$fh>;
            /(\S+)\s+(\S.*\S)/ or die $_;
            my $id = $1;
            my $seq = $2;
            $seq =~ s/\s//g;
            $self->{_seqs}{$id} = $seq;
            push @{$self->{_ids}}, $id;
        }
        # now if there are more lines, read in same order as first set, but without identifiers
        my $index = 0;
        while (<$fh>) {
            my $seq = $_;
            $seq =~ s/\s//g;
            if ($seq) {
                my $id = $self->{_ids}[$index % $ntaxa];
                $self->{_seqs}{$id} .= $seq;
                $index += 1
            }
        }
        foreach my $id (@{$self->{_ids}}) { # phylip uses '.' as insert (unknown) character
            $self->{_seqs}{$id} =~ s/\./\-/g; # replace dot as gap char with '-'
        }
    }
    elsif ($format eq 'fasta') {
        my $id;
        while (<$fh>) {
            chomp;
            if (/^>(\S+)\s*(\S.*)?/) {
                my $id = $1;
                push @{$self->{_ids}}, $id;
                $self->{_annot}{$id} = $2 if $2;
                $self->{_seqs}{$id} = '';
            }
            else  {
                $_ =~ tr/\\s//d;
                $self->{_seqs}{$id} .= $_;
            }
        }
    }
    elsif ($format eq 'nexus')
    {
        $_ = <$fh>;
        die "Format does not seem to be NEXUS" if (!/\#NEXUS/);
        my $data;
        my $ntax;
        my $nchar;
        my $matrix;
        while (<$fh>) {
            chomp;
            s/\[[^\]]*\]//g;
            $data = 1 if (/^begin data/i);
            if ($data and !$matrix and /^dimensions/i) {
                $ntax = $1 if (/ntax=(\d+)/i);
                $nchar = $1 if (/nchar=(\d+)/i);
            }
            $matrix = 1 if ($data and /^matrix/i);
            if ($matrix) {	    
                if (/^(\S+)\s+(\S+)/) {
                    push @{$self->{_ids}}, $1 unless $self->{_seqs}{$1};
                    $self->{_seqs}{$1} .= $2;
                }
                last if (/;/);
            }
        }
    }
    my $first_id = $self->{_ids}[0];
    $self->{_length} = length($self->{_seqs}{$first_id});
    $self->{_is_aligned} = 1;
    print STDERR "now review sequences, length = $self->{_length}:\n";
    for my $id (@{$self->{_ids}}) {
        print STDERR "id $id ; len ", length($self->{_seqs}{$id}), " ; is_al=$self->{_is_aligned} \n";
        if (length($self->{_seqs}{$id}) != $self->{_length}) {
            $self->{is_aligned} = 0;
            $self->{_length} = max($self->{_length}, length($self->{_seqs}{$id}))
        }
    }
}


return 1
