package ru.biosoft.physicell.sample_projects.biorobots;

import java.awt.Color;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.standard.standard_elastic_contact_function;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.xml.ModelReader;

public class CustomBiorobots
{
    public static void setColors(Model model) throws Exception
    {
        Color cargoColor = ModelReader.readColor( model.getParameter( "cargo_color" ) );
        Color workerColor = ModelReader.readColor( model.getParameter( "worker_color" ) );
        Color directorColor = ModelReader.readColor( model.getParameter( "director_color" ) );
        for( Visualizer visualizer : model.getVisualizers() )
        {
            visualizer.setColorType( CellDefinition.getCellDefinition( "cargo cell" ).type, cargoColor );
            visualizer.setColorType( CellDefinition.getCellDefinition( "worker cell" ).type, workerColor );
            visualizer.setColorType( CellDefinition.getCellDefinition( "director cell" ).type, directorColor );
        }
    }

    public static void init(Model m) throws Exception
    {
        PhysiCellUtilities.setSeed( m.getParameterInt( "random_seed" ) );
        createCellTypes( m );
        setupTissue( m );
        setColors( m );
    }

    static void createCellTypes(Model m) throws Exception
    {
        SignalBehavior.setup_signal_behavior_dictionaries( m.getMicroenvironment() );
        CellDefinition pCD = CellDefinition.getCellDefinition( "director cell" );
        pCD.functions.updatePhenotype = new DirectorCellRule();

        pCD = CellDefinition.getCellDefinition( "cargo cell" );
        pCD.functions.updatePhenotype = new CargoCellRule( m );
        pCD.functions.contact_function = new standard_elastic_contact_function();

        pCD = CellDefinition.getCellDefinition( "worker cell" );
        pCD.functions.updatePhenotype = new WorkerCellRule( m );
        pCD.functions.contact_function = new standard_elastic_contact_function();
    }

    static void setupTissue(Model model) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Zmin = m.mesh.boundingBox[2];
        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];
        double Zmax = m.mesh.boundingBox[5];
        if( m.options.simulate2D )
        {
            Zmin = 0.0;
            Zmax = 0.0;
        }
        double Xrange = Xmax - Xmin;
        double Yrange = Ymax - Ymin;

        // create some of each type of cell 
        for( int k = 0; k < CellDefinition.getDefinitionsCount(); k++ )
        {
            CellDefinition pCD = CellDefinition.getCellDefinitionByIndex( k );
            System.out.println( "Placing cells of type " + pCD.name + " ... " );
            int numberCells = model.getParameterInt( "number_of_cells" );
            for( int n = 0; n < numberCells; n++ )
            {
                double[] position = new double[] {PhysiCellUtilities.UniformRandom( Xmin, Xmax ),
                        PhysiCellUtilities.UniformRandom( Ymin, Ymax ), PhysiCellUtilities.UniformRandom( Zmin, Zmax )};
                Cell.createCell( pCD, m, position );
            }
        }
        /* custom loading here */
        int directorsNumber = model.getParameterInt( "number_of_directors" ); // 15;  
        int cargoClustersNumber = model.getParameterInt( "number_of_cargo_clusters" ); // 100;  
        int workersNumber = model.getParameterInt( "number_of_workers" ); // 50;  
        CellDefinition pCargoDef = CellDefinition.getCellDefinition( "cargo cell" );
        CellDefinition pDirectorDef = CellDefinition.getCellDefinition( "director cell" );
        CellDefinition pWorkerDef = CellDefinition.getCellDefinition( "worker cell" );

        System.out.println( "Placing cells ... " );
        // randomly place seed cells 
        double[] position = new double[3];
        double relative_margin = 0.2;
        double relative_outer_margin = 0.02;

        System.out.println( "\tPlacing " + directorsNumber + " director cells ... " );
        for( int i = 0; i < directorsNumber; i++ )
        {
            position[0] = m.options.X_range[0]
                    + Xrange * ( relative_margin + ( 1.0 - 2 * relative_margin ) * PhysiCellUtilities.UniformRandom() );

            position[1] = m.options.Y_range[0]
                    + Yrange * ( relative_outer_margin + ( 1.0 - 2 * relative_outer_margin ) * PhysiCellUtilities.UniformRandom() );
            Cell cell = Cell.createCell( pDirectorDef, m, position );
            SignalBehavior.setSingleBehavior( cell, "movable", 0 );
        }

        // place cargo clusters on the fringes 
        System.out.println( "\tPlacing cargo cells ... " );
        for( int i = 0; i < cargoClustersNumber; i++ )
        {
            position[0] = m.options.X_range[0]
                    + Xrange * ( relative_outer_margin + ( 1 - 2.0 * relative_outer_margin ) * PhysiCellUtilities.UniformRandom() );

            position[1] = m.options.Y_range[0]
                    + Yrange * ( relative_outer_margin + ( 1 - 2.0 * relative_outer_margin ) * PhysiCellUtilities.UniformRandom() );

            if( PhysiCellUtilities.UniformRandom() < 0.5 )
                Cell.createCell( pCargoDef, m, position );
            else
                createCargoCluster7( position, m );
        }

        // place workers
        System.out.println( "\tPlacing worker cells ... " );
        for( int i = 0; i < workersNumber; i++ )//number_of_workers; i++ )
        {
            position[0] = m.options.X_range[0]
                    + Xrange * ( relative_margin + ( 1.0 - 2 * relative_margin ) * PhysiCellUtilities.UniformRandom() );

            position[1] = m.options.Y_range[0]
                    + Yrange * ( relative_outer_margin + ( 1.0 - 2 * relative_outer_margin ) * PhysiCellUtilities.UniformRandom() );
            Cell.createCell( pWorkerDef, m, position );
        }
        System.out.println( "Done!" );
    }

    /** 
     * Create a hollow cluster at position, with random orientation 
     */
    static void createCargoCluster6(double[] center, Microenvironment m)
    {
        CellDefinition pCargoDef = CellDefinition.getCellDefinition( "cargo cell" );
        double spacing = 0.95 * pCargoDef.phenotype.geometry.radius * 2.0;
        double dTheta = 1.047197551196598; // 2*pi / 6.0 

        double theta = 6.283185307179586 * PhysiCellUtilities.UniformRandom();
        double[] position = new double[3];
        for( int i = 0; i < 6; i++ )
        {
            position[0] = center[0] + spacing * Math.cos( theta );
            position[1] = center[1] + spacing * Math.sin( theta );
            Cell.createCell( pCargoDef, m, position );
            theta += dTheta;
        }
    }

    /**
     * Creates a filled cluster at position, with random orientation 
     */
    static void createCargoCluster7(double[] center, Microenvironment m)
    {
        CellDefinition pCargoDef = CellDefinition.getCellDefinition( "cargo cell" );
        createCargoCluster6( center, m );
        Cell.createCell( pCargoDef, m, center );
    }

    /**
     *  Creates a small cluster at position, with random orientation 
     */
    static void createCargoCluster3(double[] center, Microenvironment m)
    {
        CellDefinition pCargoDef = CellDefinition.getCellDefinition( "cargo cell" );
        double spacing = 0.95 * pCargoDef.phenotype.geometry.radius * 1.0;
        double d_Theta = 2.094395102393195; // 2*pi / 3.0 
        double theta = 6.283185307179586 * PhysiCellUtilities.UniformRandom();
        double[] position = new double[3];
        for( int i = 0; i < 3; i++ )
        {
            position[0] = center[0] + spacing * Math.cos( theta );
            position[1] = center[1] + spacing * Math.sin( theta );
            Cell.createCell( pCargoDef, m, position );
            theta += d_Theta;
        }
    }
}
