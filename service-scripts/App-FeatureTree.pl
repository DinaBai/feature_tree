#
# The FeatureTree application.
#

use strict;
use Carp;
use Data::Dumper;
use File::Temp;
use File::Slurp;
use File::Basename;
use IPC::Run 'run';
use JSON;
use File::Copy 'copy';
#use Bio::KBase::AppService::AppConfig;
#use Bio::KBase::AppService::AppScript;
use Cwd;
use Feature_Alignment; # should be in lib directory

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;

#my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
# my $data_url = "http://www.alpha.patricbrc.org/api";

my $testing = 1;
print "args = ", join("\n", @ARGV), "\n";

if ($testing) {
    my $temp_params = JSON::decode_json(`cat /homes/allan/git/dev_container/modules/feature_tree/app_specs/instantiated_FeatureTree_1.json`);
    my $rc = build_tree('FeatureTree', undef, undef, $temp_params);
    exit $rc;
}


my $script = Bio::KBase::AppService::AppScript->new(\&build_tree, \&preflight);
my $rc = $script->run(\@ARGV);


sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    my $pf = {
	cpu => 8,
	memory => "128G",
	runtime => 0,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}



sub build_tree {
    my ($app, $app_def, $raw_params, $params) = @_;

    print "Proc FeatureTree build_tree ", Dumper($app_def, $raw_params, $params);
    my $time1 = `date`;

    #$global_token = $app->token();
    #$global_ws = $app->workspace;
    
    my $recipe = $params->{parameters}{recipe};
    
    my $tmpdir = File::Temp->newdir( CLEANUP => 0 );
    system("chmod", "755", "$tmpdir");
    print STDERR "$tmpdir\n";
    #$params = localize_params($tmpdir, $params);
    #print "after localize_params:\n", Dumper($params);
    #
    # Write job description.
    my $json = JSON::XS->new->pretty(1);
    my $jdesc = "$tmpdir/jobdesc.json";
    write_file($jdesc, $json->encode($params));
    

    my $seq_file_name = '';
    if ($params->{parameters}{sequences_from_local_file}) {
        $seq_file_name = basename($params->{parameters}{sequences_from_local_file});
        copy($params->{parameters}{sequences_from_local_file}, $seq_file_name) or die ("could not copy $params->{parameters}{sequences_from_local_file} to $tmpdir/$seq_file_name");
    }
    my $model = "GTR"; # default for DNA
    if ($params->{parameters}{alphabet} eq 'Protein') {
        $model = $params->{parameters}{protein_model}
    }
    my $output_name = $params->{parameters}{output_file};
    
    my @outputs;
    if ($recipe eq 'RAxML') {
        @outputs = run_raxml($seq_file_name, $model, $output_name, $tmpdir);
    } elsif ($recipe eq 'PhyML') {
        @outputs = run_phyml($params, $tmpdir);
    } else {
        die "Unrecognized recipe: $recipe \n";
    }
    
    print STDERR '\@outputs = '. Dumper(\@outputs);
    
    my $output_folder = $app->result_folder();
    # my $output_base   = $params->{output_file};
    
    #
    # Create folders first.
    #
    for my $fent (grep { $_->[1] eq 'folder' } @outputs)
    {
	my $folder = $fent->[0];
	my $file = basename($folder);
	my $path = "$output_folder/$file";
	eval {
	    $app->workspace->create( { objects => [[$path, 'folder']] } );
	};
	if ($@)
	{
	    warn "error creating $path: $@";
	}
	else
	{
	    my $type ="txt";
	    if (opendir(my $dh, $folder))
	    {
		while (my $filename = readdir($dh))
		{
		    if ($filename =~ /\.json$/)
		    {
			my $ofile = "$folder/$filename";
			my $dest = "$path/$filename";
			print STDERR "Output folder = $folder\n";
			print STDERR "Saving $ofile => $dest ...\n";
			$app->workspace->save_file_to_file($ofile, {}, $dest, $type, 1,
							   (-s "$ofile" > $shock_cutoff ? 1 : 0), # use shock for larger files
							   $global_token);
		    }
		}
	    }
	    else
	    {
		warn "Cannot open output folder $folder: $!";
	    }
	}
    }
    for my $output (@outputs)
    {
	my($ofile, $type) = @$output;
	next if $type eq 'folder';
	
	if (! -f $ofile)
	{
	    warn "Output file '$ofile' of type '$type' does not exist\n";
	    next;
	}
	
	if ($type eq 'job_result')
	{
            my $filename = basename($ofile);
            print STDERR "Output folder = $output_folder\n";
            print STDERR "Saving $ofile => $output_folder/$filename ...\n";
            $app->workspace->save_file_to_file("$ofile", {},"$output_folder/$filename", $type, 1);
	}
	else
	{
	    my $filename = basename($ofile);
	    print STDERR "Output folder = $output_folder\n";
	    print STDERR "Saving $ofile => $output_folder/$filename ...\n";
	    $app->workspace->save_file_to_file("$ofile", {}, "$output_folder/$filename", $type, 1,
					       (-s "$ofile" > $shock_cutoff ? 1 : 0), # use shock for larger files
					       $global_token);
	}
    }
    my $time2 = `date`;
    write_output("Start: $time1"."End:   $time2", "$tmpdir/DONE");
}

sub run_raxml {
    my ($alignment_file, $model, $output_name, $tmpdir) = @_;

    my $parallel = $ENV{P3_ALLOCATED_CPU};
    $parallel = 2 if $parallel < 2;
    
    my $cwd = getcwd();
    
    #
    #my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
    #my $dat = { data_api => "$data_api/genome_feature" };
    # no pretty, ensure it's on one line
    #my $sstring = encode_json($dat);

    if ($model eq 'GTR') {
        $model = 'GTRGAMMA'
    }
    else {
        $model = "PROTCAT".$model
    }

    my @cmd = ("raxmlHPC-PTHREADS-SSE3");
    push @cmd, ("-T", $parallel);
    push @cmd, ("-p", "12345");
    push @cmd, ("-m", $model);
    push @cmd, ("-s", $alignment_file);
    push @cmd, ("-n", $output_name);
    
    print STDERR "cmd = ", join(" ", @cmd) . "\n\n";
   
    chdir($tmpdir); 
    my ($rc, $out, $err) = run_cmd(\@cmd);
    chdir($cwd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    
    run("echo $tmpdir && ls -ltr $tmpdir");
    
    my @outputs;
    my @files = glob("$tmpdir/*.txt");
    @outputs = map { [ $_, 'txt' ] } @files;
    
    return @outputs;
}


sub curl_file {
    my ($url, $outfile) = @_;
    my @cmd = ("curl", curl_options(), "-o", $outfile, $url);
    print STDERR join(" ", @cmd)."\n";
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_ftp {
    my ($url, $outfile) = @_;
    my @cmd = ("curl", "-o", $outfile, $url);
    print STDERR join(" ", @cmd)."\n";
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_json {
    my ($url) = @_;
    my $out = curl_text($url);
    my $hash = JSON::decode_json($out);
    return $hash;
}

sub curl_options {
    my @opts;
    my $token = get_token()->token;
    push(@opts, "-H", "Authorization: $token");
    push(@opts, "-H", "Content-Type: multipart/form-data");
    return @opts;
}

sub run_cmd {
    my ($cmd) = @_;
    my ($out, $err);
    run($cmd, '>', \$out, '2>', \$err)
        or die "Error running cmd=@$cmd, stdout:\n$out\nstderr:\n$err\n";
    # print STDERR "STDOUT:\n$out\n";
    # print STDERR "STDERR:\n$err\n";
    return ($out, $err);
}

sub params_to_exps {
    my ($params) = @_;
    my @exps;
    for (@{$params->{paired_end_libs}}) {
        my $index = $_->{condition} - 1;
        $index = 0 if $index < 0;
        push @{$exps[$index]}, [ $_->{read1}, $_->{read2} ];
    }
    for (@{$params->{single_end_libs}}) {
        my $index = $_->{condition} - 1;
        $index = 0 if $index < 0;
        push @{$exps[$index]}, [ $_->{read} ];
    }
    return \@exps;
}

sub localize_params {
    my ($tmpdir, $params) = @_;
    for (@{$params->{paired_end_libs}}) {
        $_->{read1} = get_ws_file($tmpdir, $_->{read1}) if $_->{read1};
        $_->{read2} = get_ws_file($tmpdir, $_->{read2}) if $_->{read2};
    }
    for (@{$params->{single_end_libs}}) {
        $_->{read} = get_ws_file($tmpdir, $_->{read}) if $_->{read};
    }
    return $params;
}

sub get_ws {
    return $global_ws;
}

sub get_token {
    return $global_token;
}

sub get_ws_file {
    my ($tmpdir, $id) = @_;
    # return $id; # DEBUG
    my $ws = get_ws();
    my $token = get_token();
    
    my $base = basename($id);
    my $file = "$tmpdir/$base";
    # return $file; # DEBUG
    
    my $fh;
    open($fh, ">", $file) or die "Cannot open $file for writing: $!";

    print STDERR "GET WS => $tmpdir $base $id\n";
    system("ls -la $tmpdir");

    eval {
	$ws->copy_files_to_handles(1, $token, [[$id, $fh]]);
    };
    if ($@)
    {
	die "ERROR getting file $id\n$@\n";
    }
    close($fh);
    print "$id $file:\n";
    system("ls -la $tmpdir");

    return $file;
}

sub write_output {
    my ($string, $ofile) = @_;
    open(F, ">$ofile") or die "Could not open $ofile";
    print F $string;
    close(F);
}

sub break_fasta_lines {
    my ($fasta) = @_;
    my @lines = split(/\n/, $fasta);
    my @fa;
    for (@lines) {
        if (/^>/) {
            push @fa, $_;
        } else {
            push @fa, /.{1,60}/g;
        }
    }
    return join("\n", @fa);
}

sub verify_cmd {
    my ($cmd) = @_;
    system("which $cmd >/dev/null") == 0 or die "Command not found: $cmd\n";
}
