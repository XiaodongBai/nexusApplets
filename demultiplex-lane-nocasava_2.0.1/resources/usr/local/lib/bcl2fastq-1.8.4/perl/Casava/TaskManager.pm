package Casava::TaskManager;

# PROJECT: CASAVA
# MODULE:  $RCSfile: TaskManager.pm,v $
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#

=pod

=head1 NAME

Casava::TaskManager.pm - Perl libary for TaskManager(Firestarter)


=head1 SYNOPSIS

# include what functions you need... 
use Casava::TaskManager qw();  

=head1 DESCRIPTION
The library provides set of procedures and variables required to run TaskManager.

Exports:    
    
# Global variable

=head1 AUTHORSHIP

Copyright (c) 2008 Illumina

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

Created by Lukasz Szajkowski <lszajkowski@illumina.com>

=head1 SOURCE CODE

The most current release of the Perl source code 
for this module is available through CVS at genie01 a.k.a. 10.44.0.81 
cvs co BullFrog

=cut

#
# Place functions/variables you want to *export*, ie be visible
# from the caller package into @EXPORT_OK
#

use strict;
use warnings 'all';
use Exporter 'import';
our @EXPORT_OK = qw(configureWorkflow executeOrSchedule
      %taskFields checkPointInit addEverything2CheckPoint 
      executeSingleTask executeArrayTask executeArrayTaskEx
      executeCheckPoint addTask2checkPoint tryAddTask2CheckPoint closeTasksFile);

use POSIX qw(strftime);

use IO::File;
use Sys::Hostname;
use Carp;

use Casava::Common::Log;
use Casava::Common::IOLib
  qw(bufferedPrint executeCmd testIfSupports getProcessorsCount);
sub executeOrSchedule ($$$$$$$;$);
sub configureWorkflow ($$);
sub checkPointInit ($$);
sub executeSingleTask ($$$;$);
sub executeArrayTask ($$$\% );
sub executeArrayTaskEx ($$$$\% );
sub executeCheckPoint (\%);
sub addTask2checkPoint ($\%);
sub tryAddTask2CheckPoint ($\%);
sub addEverything2CheckPoint (\%);

my %taskManagerConf = ();
my %numberOfTasks   = ();
$taskManagerConf{workflowFileName} = "tasks.txt";
$taskManagerConf{cmdExecuteDef}    = 1;
$taskManagerConf{NA}               = 'N/A       ';
$taskManagerConf{defaultStatus}    = 'waiting';

our %taskFields = (
    queueId   => 0,
    taskId    => 1,
    taskType  => 2,
    appTag    => 3,
    status    => 4,
    Reuse     => 1,
    ReuseAddr => 1,
    waitFor   => 5,
    cmd       => 6,
    submitTime   => 7,
    startTime   => 8,
    host  => 9,
    finishTime  => 10,
);

my $cmdFile = undef;
my $NA    = 'N/A';

=pod

=head1 The procedure configures the workflow.

=over 4

=item configureWorkflow($isTaskManager, $taskDescriptionFile)

The procedure determines if the executeOrSchedule shuld execute the command or
write the task definition to the workflow file.


Parameters:
    $isWorkflow            - is TaskManager on
    $taskDescriptionFile   - path to athe file with task descriptions

Returns:
    nothing
        
=back

=cut

sub configureWorkflow ($$) {
    my ( $isTaskManager, $taskDescriptionFile ) = @_;
    if ( defined $isTaskManager && $isTaskManager == 1 ) {
        $taskManagerConf{cmdExecuteDef}    = 0;
        $taskManagerConf{workflowFileName} = $taskDescriptionFile;
    }    # if
    else {
        $taskManagerConf{cmdExecuteDef} = 1;
    }    # else
}

=pod

=head1 The procedure initialize checkPoint

=over 4

=item checkPointInit($checkPointId, $checkPointName)

The procedure initialize checkPoint.


Parameters:
    $checkPointId            - checkpoint id (use it to depend other check points or tasks)
    $checkPointName          - Human-readable description of the checkpoint

Returns:
    HASH MAP REF to new checkPoint or undef if the checkpoint already exists.
        
=back

=cut

sub checkPointInit ($$) {
    my ( $checkPointId, $checkPointName ) = @_;
    if (!defined $numberOfTasks{$checkPointId})
    {
        my %checkPoint = ();
        $checkPoint{checkPointId}   = $checkPointId;
        $checkPoint{checkPointName} = $checkPointName;
        return \%checkPoint;
    }
    return undef;
}

=pod

=head1 The procedure adds single task to the checkpoint
if the task exists in the workflow

=over 4

=item tryAddTask2CheckPoint($queueId, $checkPointRef)

The procedure initialize checkPoint.


Parameters:
    $queueId            - task's queue id 
    $checkPointRef      - HASH MAP REF to checkPoint 

Returns:
    undef if the task does not exist
        
=back

=cut

sub tryAddTask2CheckPoint ($\%) {
    my ( $queueId, $checkPointRef ) = @_;
    createOrOpenTaskFile ($queueId) if ( ! $cmdFile ); 
    return undef if ( !$numberOfTasks{$queueId} );

    my $task .= "${queueId}[" . ( $numberOfTasks{$queueId} - 1 ) . "],";
    push @{ $checkPointRef->{tasks} }, $task;
    return 1;
}


=pod

=head1 The procedure adds single task to the checkpoint

=over 4

=item addTask2checkPoint($checkPointId, $checkPointRef)

The procedure initialize checkPoint.


Parameters:
    $queueId            - task's queue id 
    $checkPointRef      - HASH MAP REF to checkPoint 

Returns:
    nothing
        
=back

=cut

sub addTask2checkPoint ($\%) {
    my ( $queueId, $checkPointRef ) = @_;
    errorExit "Unable to add undefined queueId $queueId to checkpoint $checkPointRef->{checkPointId}"
        if ( !tryAddTask2CheckPoint($queueId, %$checkPointRef) );
}

=pod

=head2 closeTasksFile

Flushes and closes cached workflow file handle

=over 4

=item *

Parameters:

  none

=item *

Returns:

  nothing

=item *

Exceptions

  none

=back

=cut

sub closeTasksFile() {
    if ($cmdFile) 
    {
        $cmdFile->flush(); # There might be some pending data. Flush file before reading it.
        $cmdFile->close();
        undef $cmdFile;
    }
}

=pod

=head1 The procedure adds every leaf (task or checkpoint) to the checkpoint

=over 4

=item addEverything2CheckPoint($checkPointRef)

The procedure initialize checkPoint.


Parameters:
    $checkPointRef      - HASH MAP REF to checkPoint that will receive the new
                          dependencies

Returns:
    nothing
        
=back

=cut

sub addEverything2CheckPoint (\%) {
    my ( $checkPointRef ) = @_;

    # the file is going to be re-parsed. Close cached handle if we have one
    closeTasksFile();
    
    open( FILE, "<$taskManagerConf{workflowFileName}" ) 
    || errorExit "Couldn't open file handle for $taskManagerConf{workflowFileName} $!\n";

    my %queueEndLeafCandidateTaskNumbers = ();# collects the queue ends
    my @checkPointDependencies = ();# collects the dependencies of checkpoints
    
    while (<FILE>)
    {
        chomp;
        my @taskDesc = split '\t', $_;
        my $taskQueueId = $taskDesc[$taskFields{queueId}] ;
        
        if( !defined $queueEndLeafCandidateTaskNumbers{ $taskQueueId } )
        {
            $queueEndLeafCandidateTaskNumbers{ $taskQueueId } = 0;
        }
        elsif (-1 != $queueEndLeafCandidateTaskNumbers{ $taskQueueId })
        {
            # unless we've already detected that this queue merges with another one...
            ++$queueEndLeafCandidateTaskNumbers{ $taskQueueId };
        }
        
        my $waitForQueueId = $taskDesc[ $taskFields{waitFor} ];
        if ( $waitForQueueId !~ /$NA/ )
        {
            # signal the fact that the queue is merged with another queue 
            # so it can't have a leaf
            $queueEndLeafCandidateTaskNumbers{ $waitForQueueId } = -1;
        }

        my $taskType = $taskDesc[$taskFields{taskType}];

        if ('synch' eq $taskType)
        {
            my $cmd = $taskDesc[ $taskFields{cmd} ];
            my ($dependencies) = $cmd =~ /TASKS\s+\((\S+)\)/;
            push @checkPointDependencies, $dependencies;  
        }
    }
    close FILE;

    foreach my $checkpointTaskDependency (@checkPointDependencies)
    {
        if (defined $checkpointTaskDependency)
        {
            foreach my $waitFor (split ',', $checkpointTaskDependency)
            {
                my ($waitForQueueId, $waitForTaskId) = $waitFor =~ /(\S+)\[(\S+)?\]/;
    
                if ($waitForTaskId == $queueEndLeafCandidateTaskNumbers{ $waitForQueueId })
                {
                    # signal the fact that this queue is merged with checkpoint
                    $queueEndLeafCandidateTaskNumbers{ $waitForQueueId } = -1;
                }
            }
        }
    }
    
    # by now all queue ends that are not -1 are leafs. Add them as dependencies
    foreach my $taskQueueName (keys %queueEndLeafCandidateTaskNumbers)
    {
        my $taskId = $queueEndLeafCandidateTaskNumbers{$taskQueueName};
        if (-1 != $taskId)
        {
            my $task = "$taskQueueName\[" . $taskId . "\],";
            push @{ $checkPointRef->{tasks} }, $task;
        }                
    }
}



=pod

=head1 The procedure executes or schedules the task. The task also gets added to checkPoint.

=over 4

=item executeSingleTask( $command, $queueId, $queueName; $waitFor )

The procedure executes or schedules the task.


Parameters:
    $command            - command to be executed
    $queueId            - queue Id where command should be scheduled
    $queueName          - human-readable decription of task. tabs and new-lines are not allowed
    $waitFor            - the command has to wait for $waitFor = CheckPointId or $waitFor = $queueId or N/A or undef

Returns:
    nothing
        
=back

=cut

sub executeSingleTask ($$$;$ ) {
    my ( $command, $queueId, $queueName, $waitFor ) = @_;
    $waitFor = $NA unless defined $waitFor;
    
    executeOrSchedule( $command, 0, $queueName, 0, 1, $queueId, $waitFor );
}

=pod

=head1 The procedure executes or schedules the task. The task also gets added to checkPoint.

=over 4

=item executeArrayTaskEx( $command, $queueId, $queueName, $waitForCheckPoint, $checkPointRef )

The procedure executes or schedules the task. The task also gets added to checkPoint.


Parameters:
    $command            - command to be executed
    $queueId            - queue Id where command should be scheduled
    $queuName           - human-readable decription of task. tabs and new-lines are not allowed
    $waitForCheckPoint  - the command has to wait for  $waitForCheckPoint
    $checkPointRef      - the command whould be added to $checkPointRef

Returns:
    nothing
        
=back

=cut

sub executeArrayTaskEx ($$$$\% ) {
    my ( $command, $queueId, $queueName, $waitForCheckPoint, $checkPointRef ) = @_;
    $waitForCheckPoint = $NA unless defined $waitForCheckPoint;

    if ( !$queueName || !defined $checkPointRef->{checkPointId} )
    {
        errorExit "ERROR: executeArrayTask - task name or checkPointId not defined\n";
    }    # if
    executeOrSchedule( $command, 0, $queueName, 0, 1, $queueId, $waitForCheckPoint );
    addTask2checkPoint( $queueId, %{$checkPointRef} );
}

=pod

=head1 The procedure executes or schedules the task. The task also gets added to checkPoint.

=over 4

=item executeArrayTask( $command, $queueId, $waitForCheckPoint, $checkPointRef )

The procedure executes or schedules the task. The task also gets added to checkPoint.


Parameters:
    $command            - command to be executed
    $queueId            - queue Id where command should be scheduled
    $waitForCheckPoint  - the command has to wait for  $waitForCheckPoint
    $checkPointRef      - the command whould be added to $checkPointRef

Returns:
    nothing
        
=back

=cut

sub executeArrayTask ($$$\% ) {
    my ( $command, $queueId, $waitForCheckPoint, $checkPointRef ) = @_;

    if ( !defined $checkPointRef->{checkPointName} )
    {
        errorExit "ERROR: executeArrayTask - checkPointName or checkPointId not defined\n";
    }    # if
    executeArrayTaskEx ( $command, $queueId, $checkPointRef->{checkPointName}, $waitForCheckPoint, %{$checkPointRef} ); 
}

=pod

=head1 The procedure executes or schedules the checkPoint. 

=over 4

=item executeCheckPoint($checkPointRef)

The procedure executes or schedules the checkPoint. 


Parameters:
    $checkPointRef      - HASH MAP REF to checkPoint 

Returns:
    nothing
        
=back

=cut

sub executeCheckPoint (\%) {
    my ($checkPointRef) = @_;
#    use Data::Dumper;
#    print Dumper($checkPointRef);
    my $task            = "TASKS (";
    my $taskType        = "printNumbers";
    my @taskIds         = defined $checkPointRef->{tasks} ? @{ $checkPointRef->{tasks} } : ();
    foreach my $taskId (@taskIds) {
        $task .= $taskId;
    }    # foreach
    $task .= ")";
    my $status = (@taskIds) ? undef : 'finished';
    executeOrSchedule( $task, 0, $checkPointRef->{checkPointName},
        0, 2, $checkPointRef->{checkPointId}, "N/A", $status );
}    # executeCheckPoint

=pod

=head1 The procedure executes the command or writes the command to workflow file.

=over 4

=item executeOrSchedule($command, $verbose, $appTag, $cmdExecute, 
        $verboseCmd, $queueId, $waitFor, [$status])

The procedure executes the command or writes the command to workflow file.
Calling procedure has to manage the dependencies.


Parameters:
    $command    - command to be process 
    $verbose    - verbose level
    $appTag     - application tag (e.g.: application tag id or stage id)
    $cmdExecute - if 1 then command will be executed
    $verboseCmd - ?
    $queueId    - name of queue where the command should be scheduled 
    $waitFor    - the command has to wait for $waitFor SYNCH to finish
    $status     - 'waiting' unless specified

Returns:
    nothing
        
=back

=cut

sub executeOrSchedule ($$$$$$$;$) {
    my ( $command, $verbose, $appTag, $cmdExecute, $verboseCmd, $queueId, $waitFor, $status ) = @_;
    if ( !defined $waitFor || $waitFor eq 'N/A' ) {
        $waitFor = $taskManagerConf{NA};
    }    # if
    $cmdExecute = $taskManagerConf{cmdExecuteDef};
    $status = $taskManagerConf{defaultStatus} unless $status;

    if ( !$cmdFile && $verboseCmd > 0 ) {
        createOrOpenTaskFile ($queueId);
    }    # if
    if ( !defined( $numberOfTasks{$queueId} ) ) {
        $numberOfTasks{$queueId} = 0;
    }    # if
    if ( $verboseCmd == 1 && defined($appTag) ) {
        $appTag = sprintf( "\%-10s", $appTag );
        print $cmdFile ( $queueId . "\t"
              . $numberOfTasks{$queueId} . "\t"
              . "task\t"
              . "$appTag\t"
              . "$status\t"
              . "$waitFor\t"
              . $command
              . "\n" );
        $numberOfTasks{$queueId}++;
    }    # if
    if ( $verboseCmd == 2 && defined($appTag) ) {
        $appTag = sprintf( "\%-10s", $appTag );
        print $cmdFile ( $queueId . "\t"
              . $numberOfTasks{$queueId} . "\t"
              . "synch\t"
              . "$appTag\t"
              . "$status\t"
              . "$taskManagerConf{NA}\t"
              . $command
              . "\n" );
        $numberOfTasks{$queueId}++;
    }    # if
    if ( $cmdExecute > 0 && $verboseCmd != 2 ) {
        executeCmd( $command, $verbose );
    }    # if
}     # executeOrSchedule

sub createOrOpenTaskFile
{
    my $queueId       = shift;

    %numberOfTasks = ();

    #my $header = "CASAVAWorkflow\t$numberOfTask\n";
    my $opt = ">";
    if ( -e $taskManagerConf{workflowFileName} )
    {
        open( FILE, "<$taskManagerConf{workflowFileName}" )
          || errorExit "Couldn't open file handle for $opt$taskManagerConf{workflowFileName} $!\n";
        while (<FILE>)
        {
            my $line = $_;
            chomp($line);
            my @taskDesc = split '\t', $line;
            $numberOfTasks{ $taskDesc[0] }++;
        }    # while
        close FILE;
        $opt = ">>";
    }
    else
    {
        $numberOfTasks{$queueId} = 0;
    }    # else

    $cmdFile = IO::File->new("$opt$taskManagerConf{workflowFileName}")
        or errorExit "$0:ERROR: Couldn't open $taskManagerConf{workflowFileName} $!\n";
}


1;    # says use was ok
__END__
