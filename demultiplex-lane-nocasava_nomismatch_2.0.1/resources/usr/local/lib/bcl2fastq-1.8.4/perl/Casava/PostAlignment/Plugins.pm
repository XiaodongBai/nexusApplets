=head1 LICENSE

Copyright (c) 2007-2009 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=head1 NAME

Casava::PostAlignment::Plugins - Library for accessing CASAVA plugins.

=head1 DESCRIPTION

This module provides support for CASAVA plugins framework.

=head1 SUBROUTINES

=cut

package Casava::PostAlignment::Plugins;

use strict;
use warnings "all";
use Exporter 'import';
our @EXPORT_OK = qw(target2Plugin targetExecuteOrGenerateWorkflow listTargets
    targetHelpUsage getOptionsMapping getPluginScriptsDir);

use Casava::Common::Log qw(logInfo logWarning errorExit);
use Casava::PostAlignment::Sequencing::Config qw(%CONF_PROJ);
use Casava::TaskManager qw(closeTasksFile);

my %plugins = ();

my @ignoreModules = (__PACKAGE__);

sub checkInitialize() {
    if (!keys %plugins)
    {
        my $CASAVA_FULL_PERL_LIBDIR = '/usr/local/lib/bcl2fastq-1.8.4/perl'; #substituted by cmake
        errorExit "\$CASAVA_FULL_PERL_LIBDIR unconfigured. CASAVA is improperly installed"
            unless ("\@CASAVA_FULL_PERL_LIBDIR\@" ne $CASAVA_FULL_PERL_LIBDIR);
        
        my $CASAVA_PLUGINS_DIR = File::Spec->catfile(${CASAVA_FULL_PERL_LIBDIR}, 'Casava', 'PostAlignment', 'Plugins');

        logInfo("Searching for plugins in $CASAVA_PLUGINS_DIR", 5);
        foreach my $pluginFile (glob (File::Spec->catfile($CASAVA_PLUGINS_DIR, '*.pm')))
        {
            my ($volume, $path, $file) = File::Spec->splitpath($pluginFile);
            my ($moduleName) = $file =~ /([^\.]*)\.pm/;
    
            my $moduleTarget = undef;
            my $moduleQualName = "Casava::PostAlignment::Plugins::$moduleName";
            if (!grep(/^$moduleQualName$/, @ignoreModules))
            {
                logInfo("Checking package $moduleQualName", 5);
                eval "use $moduleQualName; \$moduleTarget = $moduleQualName\::getTarget()";
                logWarning ("$@") if $@;
                if (defined $moduleTarget)
                {
                    my %plugin = (pluginFile=>$pluginFile, pluginQualName=>$moduleQualName);
                    $plugins{plugins}{$moduleTarget} = \%plugin;
                    
                    logInfo("Detected plugin $moduleTarget", 0);
                }
            }
        }
    }
    $plugins{initialized} = 1;
}


sub target2Plugin($) {
    my ( $target ) = @_;
    checkInitialize();

    return $plugins{plugins}{$target};
}


=pod

=head2 targetExecuteOrGenerateWorkflow($target)

Calls executeOrGenerateWorkflow on a corresponding plugin given the target name.

=over 4

=item *

Parameters:

  $target    - name of the target

=item *

Returns:

    'unknown' in case plugin not found. Ohterwise plugin-generated user-readable target status string
    or undef if no special message is provided by plugin

=item *

Exceptions

  Plugin unknown for specified target

=back

=cut
sub targetExecuteOrGenerateWorkflow($$$\%) {
    my ( $target, $prefTarget, $projectDir, $buildChromsBinSizesRef ) = @_;
    checkInitialize();

    my $plugin = target2Plugin($target);

    if (defined $plugin)
    {
        my $pluginRet;
        eval "use $plugin->{pluginQualName}; \$pluginRet = $plugin->{pluginQualName}::executeOrGenerateWorkflow(\$prefTarget, \$projectDir, \%CONF_PROJ, \%\$buildChromsBinSizesRef)";
        errorExit "$@" if $@;
        closeTasksFile();
        return $pluginRet;
    }
    return 'unknown';
}

=pod


=head2 listTargets

Returns list of available targets.

I<Parameters:>

=over 4

=item *

None

=back

I<Returns:>

Array of target names available via installed plugins

I<Exceptions:>

=over 4

=item *

None

=back

=cut


sub listTargets() {
    checkInitialize();

    return (keys %{$plugins{plugins}});
}

=pod

=head2 targetHelpUsage($target)

Calls Pod::Usage::pod2usage on a corresponding plugin file given the target name.

=over 4

=item *

Parameters:

  $target    - name of the target

=item *

Returns:

    nothing

=item *

Exceptions

  Plugin unknown for specified target

=back

=cut
sub targetHelpUsage($) {
    my ( $target ) = @_;
    checkInitialize();

    my $plugin = target2Plugin($target);

    if (defined $plugin)
    {
        my ($volume, $path, $file) = File::Spec->splitpath($plugin->{pluginFile});

        Cwd::chdir(File::Spec->catfile($volume, $path));
        Pod::Usage::pod2usage
        (
            -exitstatus => 0,
            -verbose => 2,
            -input => $file,
        );#that should exit the process....
    }
    errorExit ("Plugin for target $target is unknown");
}

=pod

=head2 getOptionsMapping

Returns array of command line argument mappings for all loaded plugins.

=over 4

=item *

Parameters:

  PARAMS    - reference to target hash

=item *

Returns:

  mappings

=item *

Exceptions

  Plugin unknown for specified target

=back

=cut
sub getOptionsMapping(\%) {
    my ( $PARAMS ) = @_;
    checkInitialize();

    my @mappings;

    foreach my $target (keys %{$plugins{plugins}})
    {
        my $plugin = $plugins{plugins}{$target};
        my $qualName = $plugin->{pluginQualName};
        my @mapping;
        eval "use $plugin->{pluginQualName}; \@mapping = $plugin->{pluginQualName}::getOptionsMapping(\%\$PARAMS)";
        errorExit "$@" if $@;
        push @mappings, @mapping;
    }
    return @mappings;
}

=pod

=head2 getPluginScriptsDir

Returns path to folder containing the executables installed for the plugin

=over 4

=item *

Parameters:

  $pluginPackage    - qualified package name of the plugin

=item *

Returns:

  directory path

=item *

Exceptions

  unparsable package name

=back

=cut
sub getPluginScriptsDir($){
    my ($pluginPackage) = @_;

    errorExit "ERROR: Could not parse plugin module name from: $pluginPackage"
        unless ($pluginPackage =~ m/::([^:]*)$/);
        
    my $pluginPrivateScriptsDir = File::Spec->catfile('/usr/local/libexec/bcl2fastq-1.8.4', 'PostAlignment', 'Plugins', $1);
    return $pluginPrivateScriptsDir;
}

1;
__END__

=pod

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over 4

=item Standard perl modules

strict, warnings, GetOpt::Long, Pod::Usage, File::Spec

=item External perl modules

=item Casava perl modules

=back

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Roman Petrovski

=cut
