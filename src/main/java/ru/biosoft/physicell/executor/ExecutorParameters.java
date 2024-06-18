package ru.biosoft.physicell.executor;

import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellContainerExperimental;
import ru.biosoft.physicell.core.CellContainerParallel;

public class ExecutorParameters
{
    String resultPath;
    boolean parallelDiffusion;
    String cellType = CellContainer.DEFAULT_NAME;
    boolean experimentalCell;
    boolean quiet = false;
    boolean trackTime = false;
    String project = null;
    long seed = (long) ( Math.random() * 1E12 );

    public ExecutorParameters(String ... args)
    {
        for( int i = 0; i < args.length; i++ )
        {
            String arg = args[i];
            if( arg.equals( "--output-dir" ) )
            {
                resultPath = args[++i];
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
            else if( arg.equals( "--run" ) || arg.equals( "-r" ) )
            {
                project = args[++i];
            }
            else if( arg.equals( "-s" ) || arg.equals( "--seed" ) )
            {
                seed = Long.parseLong( args[++i] );
            }
            else if( arg.startsWith( "-" ) )
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
            }
        }
    }
}