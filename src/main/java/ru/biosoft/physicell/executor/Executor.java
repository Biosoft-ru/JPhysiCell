package ru.biosoft.physicell.executor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import ru.biosoft.physicell.biofvm.ConstantCoefficientsLOD3D;
import ru.biosoft.physicell.biouml.BioUMLIntraReader;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.xml.ModelReader;

public class Executor
{
    private static String PROJECTS_FILE = "projects";
    private static String PATH_TO_SETTINGS = "config/PhysiCell_settings.xml";
    private static String HELP_FILE = "help";
    private static String VERSION_FILE = "version";
    private static String LICENSE_FILE = "license";
    private static Map<String, String> projects;

    public static void main(String ... args)
    {
        ExecutorParameters parameters = new ExecutorParameters( args );

        if( projects == null )
            findProjects();

        if( parameters.listProjects )
        {
            if( projects == null )
                findProjects();
            log( "There are " + projects.size() + " projects: " + String.join( ", ", projects.keySet() ) );
            return;
        }

        if( parameters.showHelp )
        {
            log( getResource( HELP_FILE ) );
            return;
        }

        if( parameters.showVersion )
        {
            log( getResource( VERSION_FILE ) );
            return;
        }

        if( parameters.showLicense )
        {
            log( getResource( LICENSE_FILE ) );
            return;
        }

        String project = parameters.project;
        if( project == null )
        {
            log( "Unknown command, please specify project. Use --help for more information." );
            return;
        }

        if( !projects.keySet().contains( project ) )
        {
            log( "Unknown project " + project + ". Use --projects to see available projects." );
            return;
        }

        if( projects.keySet().contains( project ) )
        {
            log( "Selected project: " + project + "." );
            String path = projects.get( project );
            try
            {
                Class<?> c = Class.forName( path );
                InputStream stream = null;
                if( parameters.settingdPath != null )
                {
                    stream = new FileInputStream( new File( parameters.settingdPath ) );
                }
                else
                {
                    stream = c.getResourceAsStream( PATH_TO_SETTINGS );
                }

                ModelReader reader = new ModelReader();
                reader.setIntracellularReader( new BioUMLIntraReader() );
                Model model = reader.read( stream, c );
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

            ( (ConstantCoefficientsLOD3D)model.getMicroenvironment().getSolver() ).setPrallel( params.parallelDiffusion );
            model.createContainer( 30, params.cellType );

            if( params.resultPath != null )
            {
                String path = params.resultPath;
                //                if( params.localResult )
                //                {
                //                    String jarPath = Executor.class.getProtectionDomain().getCodeSource().getLocation().getFile();
                //                    String parent = new File( jarPath ).getParent();
                //                    path = parent + "/" + path;
                //                }
                model.setResultFolder( path );
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

            if( params.finalTime > 0 )
                model.setTMax( params.finalTime );

            model.init();
            model.setSeed( params.seed );
            if( params.printInfo )
            {
                log( model.display() );
            }
            else if( params.simulate )
            {
                double tStart = System.nanoTime();
                model.simulate();
                tStart = System.nanoTime() - tStart;

                log( "Completed, elapsed time: " + ( System.nanoTime() - startTime ) / 1E9 + " s." );
            }
            else
            {
                log( "Nothing to do with selected project. Use -r to simulate or -i to print info." );
            }
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

    private static void log(String s)
    {
        System.out.println( PhysiCellUtilities.getCurrentTime() + s );
    }

    public static String getResource(String name)
    {
        InputStream stream = Executor.class.getResourceAsStream( name );
        if( stream == null )
            return "Can not find resource :" + name;
        BufferedReader reader = new BufferedReader( new InputStreamReader( stream ) );
        return reader.lines().collect( Collectors.joining( "\n" ) );
    }

    public static void findProjects()
    {
        projects = new TreeMap<>();
        InputStream stream = Executor.class.getResourceAsStream( PROJECTS_FILE );
        BufferedReader reader = new BufferedReader( new InputStreamReader( stream ) );
        reader.lines().map( s -> s.split( "\t" ) ).forEach( s -> projects.put( s[0], s[1] ) );
    }
}