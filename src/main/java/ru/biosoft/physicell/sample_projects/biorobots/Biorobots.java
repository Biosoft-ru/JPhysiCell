package ru.biosoft.physicell.sample_projects.biorobots;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.standard.StandardElasticContact;

public class Biorobots extends Model
{
    @Override
    public void init() throws Exception
    {
        super.init();
        //        PhysiCellUtilities.setSeed( (int) ( Math.random() * 1000 ) );
        //        PhysiCellUtilities.setSeed( getParameterInt( "random_seed" ) );
        createCellTypes();
        setupTissue();
        getVisualizers().forEach( v -> v.setAgentVisualizer( new BiorobotsVisualizer( this ) ) );
    }

    void createCellTypes() throws Exception
    {
        signals.setupDictionaries( this );
        CellDefinition pCD = getCellDefinition( "director cell" );
        pCD.functions.updatePhenotype = new DirectorCellRule();

        pCD = getCellDefinition( "cargo cell" );
        pCD.functions.updatePhenotype = new CargoCellRule( this );
        pCD.functions.contact = new StandardElasticContact();

        pCD = getCellDefinition( "worker cell" );
        pCD.functions.updatePhenotype = new WorkerCellRule( this );
        pCD.functions.contact = new StandardElasticContact();
    }

    void setupTissue() throws Exception
    {
        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];

        double Xrange = Xmax - Xmin;
        double Yrange = Ymax - Ymin;

        int directorsNumber = getParameterInt( "number_of_directors" ); // 15;  
        int cargoClustersNumber = getParameterInt( "number_of_cargo_clusters" ); // 100;  
        int workersNumber = getParameterInt( "number_of_workers" ); // 50;  
        CellDefinition pCargoDef = getCellDefinition( "cargo cell" );
        CellDefinition pDirectorDef = getCellDefinition( "director cell" );
        CellDefinition pWorkerDef = getCellDefinition( "worker cell" );

        // randomly place seed cells 
        double[] position = new double[3];
        double relative_margin = 0.2;
        double relative_outer_margin = 0.02;

        for( int i = 0; i < directorsNumber; i++ )
        {
            position[0] = m.options.X_range[0]
                    + Xrange * ( relative_margin + ( 1.0 - 2 * relative_margin ) * rng.UniformRandom() );

            position[1] = m.options.Y_range[0]
                    + Yrange * ( relative_outer_margin + ( 1.0 - 2 * relative_outer_margin ) * rng.UniformRandom() );
            Cell cell = Cell.createCell( pDirectorDef, this, position );
            getSignals().setSingleBehavior( cell, "movable", 0 );
        }

        // place cargo clusters on the fringes 
        for( int i = 0; i < cargoClustersNumber; i++ )
        {
            position[0] = m.options.X_range[0]
                    + Xrange * ( relative_outer_margin + ( 1 - 2.0 * relative_outer_margin ) * rng.UniformRandom() );

            position[1] = m.options.Y_range[0]
                    + Yrange * ( relative_outer_margin + ( 1 - 2.0 * relative_outer_margin ) * rng.UniformRandom() );

            if( rng.UniformRandom() < 0.5 )
                Cell.createCell( pCargoDef, this, position );
            else
                createCargoCluster7( position, m );
        }

        // place workers
        for( int i = 0; i < workersNumber; i++ )
        {
            position[0] = m.options.X_range[0]
                    + Xrange * ( relative_margin + ( 1.0 - 2 * relative_margin ) * rng.UniformRandom() );

            position[1] = m.options.Y_range[0]
                    + Yrange * ( relative_outer_margin + ( 1.0 - 2 * relative_outer_margin ) * rng.UniformRandom() );
            Cell.createCell( pWorkerDef, this, position );
        }
    }

    /** 
     * Create a hollow cluster at position, with random orientation 
     */
    void createCargoCluster6(double[] center, Microenvironment m)
    {
        CellDefinition pCargoDef = getCellDefinition( "cargo cell" );
        double spacing = 0.95 * pCargoDef.phenotype.geometry.radius * 2.0;
        double dTheta = 1.047197551196598; // 2*pi / 6.0 

        double theta = 6.283185307179586 * rng.UniformRandom();
        double[] position = new double[3];
        for( int i = 0; i < 6; i++ )
        {
            position[0] = center[0] + spacing * Math.cos( theta );
            position[1] = center[1] + spacing * Math.sin( theta );
            Cell.createCell( pCargoDef, this, position );
            theta += dTheta;
        }
    }

    /**
     * Creates a filled cluster at position, with random orientation 
     */
    void createCargoCluster7(double[] center, Microenvironment m)
    {
        CellDefinition pCargoDef = getCellDefinition( "cargo cell" );
        createCargoCluster6( center, m );
        Cell.createCell( pCargoDef, this, center );
    }

    /**
     *  Creates a small cluster at position, with random orientation 
     */
    void createCargoCluster3(double[] center, Model model)
    {
        CellDefinition pCargoDef = getCellDefinition( "cargo cell" );
        double spacing = 0.95 * pCargoDef.phenotype.geometry.radius * 1.0;
        double d_Theta = 2.094395102393195; // 2*pi / 3.0 
        double theta = 6.283185307179586 * rng.UniformRandom();
        double[] position = new double[3];
        for( int i = 0; i < 3; i++ )
        {
            position[0] = center[0] + spacing * Math.cos( theta );
            position[1] = center[1] + spacing * Math.sin( theta );
            Cell.createCell( pCargoDef, model, position );
            theta += d_Theta;
        }
    }
}
