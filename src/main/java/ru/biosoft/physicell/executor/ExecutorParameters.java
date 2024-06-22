package ru.biosoft.physicell.executor;

import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellContainerExperimental;
import ru.biosoft.physicell.core.CellContainerParallel;

public class ExecutorParameters
{
    String resultPath;
    boolean localResult = false;
    boolean parallelDiffusion;
    String cellType = CellContainer.DEFAULT_NAME;
    boolean experimentalCell;
    boolean quiet = false;
    boolean trackTime = false;
    String project = null;
    long seed = (long) ( Math.random() * 1E12 );
    String settingdPath;
    int finalTime = -1;
    boolean printInfo;
    boolean simulate;
    boolean showHelp;
    boolean listProjects;

    public static String PROJECTS_CMD = "--projects";
    public static String HELP_CMD = "--help";
    public static String SHORT_HELP_CMD = "-h";

    public ExecutorParameters(String ... args)
    {
        for( int i = 0; i < args.length; i++ )
        {
            String arg = args[i];

            if( arg.equals( PROJECTS_CMD ) )
            {
                this.listProjects = true;
                return;
            }
            if( HELP_CMD.equals( arg ) || SHORT_HELP_CMD.equals( arg ) )
            {
                this.showHelp = true;
                return;
            }

            if( arg.equals( "--output-dir" ) || arg.equals( "-o" ) )
            {
                String result = args[++i];
                if( !result.contains( "/" ) && !result.contains( "\\" ) )
                    localResult = true;
                resultPath = result;
            }
            else if( arg.equals( "--quiet" ) )
            {
                quiet = true;
            }
            else if( arg.equals( "--track-time" ) )
            {
                trackTime = true;
            }
            else if( arg.equals( "--parallel" ) )
            {
                cellType = CellContainerParallel.PARALLEL_CONTAINER_NAME;
            }
            else if( arg.equals( "--experimental" ) )
            {
                cellType = CellContainerExperimental.EXPERIMENTAL_CONTAINER_NAME;
            }
            else if( arg.equals( "--diff-parallel" ) )
            {
                parallelDiffusion = true;
            }
            else if( arg.equals( "--run" ) )
            {
                simulate = true;
            }
            else if( arg.equals( "--name" ) || arg.equals( "-n" ) )
            {
                project = args[++i];
            }
            else if( arg.equals( "--info" ) )
            {
                printInfo = true;
            }
            else if( arg.equals( "-s" ) || arg.equals( "--seed" ) )
            {
                seed = Long.parseLong( args[++i] );
            }
            else if( arg.equals( "--settings" ) )
            {
                this.settingdPath = args[++i];
            }
            else if( arg.equals( "--tf" ) )
            {
                this.finalTime = Integer.parseInt( args[++i] );
            }
            else if( arg.startsWith( "-" ) && !arg.subSequence( 0, 2 ).equals( "--" ) )
            {
                if( arg.contains( "p" ) )
                {
                    cellType = CellContainerParallel.PARALLEL_CONTAINER_NAME;
                }
                if( arg.contains( "x" ) )
                {
                    cellType = CellContainerExperimental.EXPERIMENTAL_CONTAINER_NAME;
                }
                if( arg.contains( "d" ) )
                {
                    parallelDiffusion = true;
                }
                if( arg.contains( "q" ) )
                {
                    quiet = true;
                }
                if( arg.contains( "t" ) )
                {
                    trackTime = true;
                }
                if( arg.contains( "r" ) )
                {
                    this.simulate = true;
                }
                if( arg.contains( "i" ) )
                {
                    this.printInfo = true;
                }
            }
        }
    }
}