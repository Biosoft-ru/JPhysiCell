package ru.biosoft.physicell.executor;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import ru.biosoft.physicell.biofvm.ConstantCoefficientsLOD3D;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.xml.ModelReader;

public class Executor
{
    private static String PROJECTS_FILE = "projects";
    public static String PROJECTS_CMD = "-projects";
    public static String HELP_CMD = "-help";
    public static String PATH_TO_SETTINGS = "config/PhysiCell_settings.xml";
    private static String HELP_FILE = "help";
    public static String help;
    public static Map<String, String> projects;

    public static void main(String ... args)
    {
        if( projects == null )
            findProjects();

        if( help == null )
            findHelp();

        if( PROJECTS_CMD.equals( args[0] ) )
        {
            System.out.println( "There are " + projects.size() + " projects: " + String.join( ", ", projects.keySet() ) );
            return;
        }
        else if( HELP_CMD.equals( args[0] ) )
        {
            System.out.println( help );
            return;
        }

        ExecutorParameters parameters = new ExecutorParameters( args );
        String project = parameters.project;
        if( project == null )
        {
            System.out.println( "Unknown command, please specify project. Use --help for more information." );
            return;
        }

        if( !projects.keySet().contains( project ) )
        {
            System.out.println( "Unknown project " + project + ". Use --projects to see available projects." );
            return;
        }

        if( projects.keySet().contains( project ) )
        {
            System.out.println( PhysiCellUtilities.getCurrentTime() + "Running " + project + " project..." );
            String path = projects.get( project );
            try
            {
                Class c = Class.forName( path );
                InputStream stream = c.getResourceAsStream( PATH_TO_SETTINGS );
                Model model = new ModelReader().read( stream, c );
                runProject( model, parameters );
            }
            catch( Exception ex )
            {
                ex.printStackTrace();
            }
        }
    }

    public static void runProject(Model model, ExecutorParameters params)
    {
        try
        {
            double startTime = System.nanoTime();
            PhysiCellUtilities.setSeed( params.seed );
            ( (ConstantCoefficientsLOD3D)model.getMicroenvironment().getSolver() ).setPrallel( params.parallelDiffusion );
            model.createContainer( 30, params.cellType );

            if( params.resultPath != null )
            {
                model.setResultFolder( params.resultPath );
                model.setSaveFull( true );
                model.setSaveImg( true );
                model.setWriteDensity( true );
                for( int i = 0; i < model.getMicroenvironment().densityNames.length; i++ )
                    model.addGIFVisualizer( 0, model.getMicroenvironment().densityNames[i] ).setStubstrateIndex( i ).setMaxDensity( 0.5 );
            }
            else
            {
                model.setSaveFull( false );
                model.setSaveImg( false );
                model.setWriteDensity( false );
            }
            if( params.quiet )
                model.setVerbose( false );

            model.init();
            double tStart = System.nanoTime();
            model.simulate();
            tStart = System.nanoTime() - tStart;

            System.out.println(
                    PhysiCellUtilities.getCurrentTime() + "Completed, elapsed time: " + ( System.nanoTime() - startTime ) / 1E9 + " s." );

            //            System.out.println( "Total " + tStart / 1E9 );
            //            System.out.println( "Diffusion " + Model.tDiffusion / 1E9 );
            //            System.out.println( "Secretion " + CellContainer.tSecretion / 1E9 );
            //            System.out.println( "Phenotype " + CellContainer.tPhenotype / 1E9 );
            //            System.out.println( "Velocity " + CellContainer.tVelocity / 1E9 );
            //            System.out.println( "Interaction " + CellContainer.tInteraction / 1E9 );
            //            System.out.println( "Contact " + CellContainer.tContact / 1E9 );
            //            System.out.println( "Custom " + CellContainer.tCustom / 1E9 );
            //            System.out.println( "Attachment " + CellContainer.tAttachment / 1E9 );
            //            System.out.println( "Divide " + CellContainer.tDivide / 1E9 );
            //            System.out.println( "Rest " + CellContainer.tRestAll / 1E9 );
            //            System.out.println( "Gradient " + CellContainer.tGradient / 1E9 );
            //            System.out.println( "Total cell " + CellContainer.tTotal / 1E9 );
        }
        catch( Exception ex )
        {
            ex.printStackTrace();
        }
    }

    public static void findHelp()
    {
        InputStream stream = Executor.class.getResourceAsStream( HELP_FILE );
        BufferedReader reader = new BufferedReader( new InputStreamReader( stream ) );
        help = reader.lines().collect( Collectors.joining( "\n" ) );
    }

    public static void findProjects()
    {
        projects = new TreeMap<>();
        InputStream stream = Executor.class.getResourceAsStream( PROJECTS_FILE );
        BufferedReader reader = new BufferedReader( new InputStreamReader( stream ) );
        reader.lines().map( s -> s.split( "\t" ) ).forEach( s -> projects.put( s[0], s[1] ) );
    }
}