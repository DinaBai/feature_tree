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
use File::Copy ('copy', 'move');
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Cwd;
use Feature_Alignment; # should be in lib directory

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;

my $testing = 1;
print "args = ", join("\n", @ARGV), "\n";

if ( -e $ARGV[0]) { #testing
    #my $temp_params = JSON::decode_json(`cat /homes/allan/git/dev_container/modules/feature_tree/app_specs/instantiated_FeatureTree_1.json`);
    print "Testing with json file $ARGV[0]\n";
    my $json_data = File::Slurp::read_file($ARGV[0]);
    my $temp_params = JSON::decode_json($json_data);
    print " now call build_tree with the json object\n";
    my $rc = build_tree('FeatureTree', undef, undef, $temp_params);
    exit $rc;
}

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
#my $data_url = "http://www.alpha.patricbrc.org/api";

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

    my $recipe = $params->{recipe};
    
    my $tmpdir = File::Temp->newdir( "/tmp/FeatureTree_XXXXX", CLEANUP => 0 );
    system("chmod", "755", "$tmpdir");
    print STDERR "$tmpdir\n";
    #$params = localize_params($tmpdir, $params);
    #print "after localize_params:\n", Dumper($params);
    #
    # Write job description.
    my $json = JSON::XS->new->pretty(1);
    my $jdesc = "$tmpdir/jobdesc.json";
    write_file($jdesc, $json->encode($params));
    

    print "copy data to temp dir\n";
    my $seq_file_name = '';
    if ($params->{sequence_source} eq 'local_file') {
        $seq_file_name = basename($params->{sequences});
        print STDERR "input data file is $seq_file_name\n";
        copy($params->{sequences}, "$tmpdir/$seq_file_name") or die ("could not copy $params->{sequences} to $tmpdir/$seq_file_name");
    }
    elsif ($params->{sequence_source} eq 'ws') {
        $seq_file_name = basename($params->{sequences});
        print STDERR "input data file is $seq_file_name\n";
        $app->workspace->download_file($params->{sequences}, "$tmpdir/$seq_file_name", 1, $global_token);
    }
    elsif ($params->{sequence_source} eq 'upload') {
        $seq_file_name = $params->{output_file}."_sequence_input.txt";
        print STDERR "input data file is $seq_file_name\n";
        File::Slurp::write_file("$tmpdir/$seq_file_name", $params->{sequences})
    }
    run("echo $tmpdir && ls -ltr $tmpdir");

    my $model = "AUTO"; # default for protein
    if ($params->{alphabet} eq 'Protein' and $params->{protein_model}) {
        $model = $params->{protein_model}
    }
    elsif ($params->{alphabet} eq 'DNA') {
        $model = "GTR"
    }
    my $output_name = $params->{output_file};
    my $alphabet = $params->{alphabet};

    my @outputs;
    if (lc($recipe) eq 'raxml') {
        @outputs = run_raxml($seq_file_name, $alphabet, $model, $output_name, $tmpdir);
    } elsif (lc($recipe) eq 'phyml') {
        @outputs = run_phyml($seq_file_name, $alphabet, $model, $output_name, $tmpdir);
    } else {
        die "Unrecognized recipe: $recipe \n";
    }
    
    print STDERR '\@outputs = '. Dumper(\@outputs);
    
    my $output_folder = $app->result_folder();
    # my $output_base   = $params->{output_file};
    #$app->workspace->create( { objects => [[$path, 'folder']] } );
    
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        next if $type eq 'folder';
        
        if (! -f $ofile) {
            warn "Output file '$ofile' of type '$type' does not exist\n";
            next;
        }
        
        if ($type eq 'job_result') {
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
    my ($alignment_file, $alphabet, $model, $output_name, $tmpdir) = @_;

    my $parallel = $ENV{P3_ALLOCATED_CPU};
    $parallel = 2 if $parallel < 2;
    
    my $cwd = getcwd();
    
    if ($alphabet eq 'DNA') {
        $model = 'GTRGAMMA'
    }
    else {
        $model = "PROTCAT". uc($model)
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
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    
    my @outputs;
    my $bestTreeFile = $output_name . "_RAxML_bestTree.nwk";
    move("RAxML_bestTree.".$output_name, $bestTreeFile);
    push @outputs, ["$tmpdir/$bestTreeFile", 'nwk'];
    push @outputs, ["$tmpdir/RAxML_info.".$output_name, 'txt'];
    push @outputs, ["$tmpdir/jobdesc.json", 'json'];

    chdir($cwd);
    run("echo $tmpdir && ls -ltr $tmpdir");

    return @outputs;
}

sub run_phyml {
    my ($alignment_file, $alphabet, $model, $output_name, $tmpdir) = @_;

    my $cwd = getcwd();
    
    if ($alphabet eq 'DNA') {
        $model = 'GTR'
    }

    my @cmd = ("phyml");
    push @cmd, ("-m", $model);
    push @cmd, ("-i", $alignment_file);
    
    print STDERR "cmd = ", join(" ", @cmd) . "\n\n";
   
    chdir($tmpdir); 
    my ($rc, $out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    
    my @outputs;
    my $treeFile = $alignment_file."_phyml_tree.nwk";
    move($alignment_file."_phyml_tree.txt", $treeFile);
    my $statsFile = $alignment_file."_phyml_stats.txt";
    push @outputs, ["$tmpdir/$treeFile", 'nwk'];
    push @outputs, ["$tmpdir/$statsFile", 'txt'];
    push @outputs, ["$tmpdir/jobdesc.json", 'json'];

    chdir($cwd);
    run("echo $tmpdir && ls -ltr $tmpdir");

    return @outputs;
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

sub verify_cmd {
    my ($cmd) = @_;
    system("which $cmd >/dev/null") == 0 or die "Command not found: $cmd\n";
}
