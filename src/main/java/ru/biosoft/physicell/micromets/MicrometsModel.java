package ru.biosoft.physicell.micromets;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import one.util.streamex.IntStreamEx;
import one.util.streamex.StreamEx;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.standard.StandardElasticContact;
import ru.biosoft.physicell.ui.AgentColorer;

public class MicrometsModel extends Model
{
    int recruited_Tcells = 0;
    int recruited_macrophages = 0;
    int recruited_CD4Tcells = 0;
    int recruited_DCs = 0;

    double first_macrophage_recruitment_time = 9e9;
    double first_CD8_T_cell_recruitment_time = 9e9;
    double first_CD4_T_cell_recruitment_time = 9e9;
    double first_DC_recruitment_time = 9e9;

    private List<Cell> cells_to_move_from_edge = new ArrayList<>();
    private List<double[]> valid_position = new ArrayList<>();
    List<Integer> vascularized_voxel_indices = new ArrayList<>();

    LymphNode lympNode = new LymphNode();
    external_immue external_immue;
    dc_history dc_history;

    History historyTc = new History( 60 );
    History historyTh = new History( 60 );
    History history = new History( 72000 );

    public String getLogInfo() throws Exception
    {
        int[] cells = cell_count();
        return Math.round( getCurrentTime()) + " " +IntStreamEx.of( cells ).joining(" ");
//        return PhysiCellUtilities.getCurrentTime() + "\tElapsed\t" + ( System.currentTimeMillis() - startTime ) / 1000 + "\tTime:\t"
//                + (int)Math.round( curTime ) + "\tCells\t" + m.getAgentsCount() + " CD8: " + getCells( "CD8: " ) + " macrophages: "
//                + getCells( "macrophage" ) +" activated macropaghes: "+ getActivatedCells("macrophage")+" DC: " + getCells( "DC" ) + " CD4: " + getCells( "CD4 Tcell" );
    }

    private long getCells(String type)
    {
        return StreamEx.of( m.getAgents( Cell.class ) ).filter( c -> c.typeName.equals( type ) ).count();
    }
    
    private long getActivatedCells(String type)
    {
        return StreamEx.of( m.getAgents( Cell.class ) ).filter( c -> c.typeName.equals( type ) ).filter( c->c.customData.get(  "activated_immune_cell" ) > 0.5  ).count();
    }

    @Override
    public void setupInitial() throws Exception
    {
        immune_submodels_setup();
        create_cell_types();
        setup_tissue();
        
        lympNode.InitialCondition( getParameterDouble("initial_DCs_lymphNode"),10.0,1.0,1.0,0.0,0.0);
    }

    @Override
    public AgentColorer getDefaultColorer()
    {
        return new MicrometsColorer();
    }

    public MicrometsModel()
    {
        external_immue = new external_immue( this );
        dc_history = new dc_history( this );
    }

    class LymphNode
    {
        int DCAMOUNT = 0;
        double GridCOUNT = 1.0;
        double DM;
        double TC;
        double TH1;
        double TH2;
        double TCt;
        double Tht;

        void InitialCondition(double dm, double tc, double th1, double th2, double tct, double tht)
        {
            DM = dm;
            TC = tc;
            TH1 = th1;
            TH2 = th2;
            TCt = tct;
            Tht = tht;
        }
    }

    @Override
    public void stepBeforeCells() throws Exception
    {
        dc_history.DC_history_main_model( diffusionStep );

        external_immue.external_immune_model( diffusionStep );
        clear_cells_to_move_from_edge();
        include_tumor_cells();
        vaccine();
    }

    @Override
    public void stepAfterCells() throws Exception
    {
        divide_custom_data();
        process_tagged_cells_on_edge();
        immune_cell_recruitment( diffusionStep );
    }

    void create_cell_types()
    {
        // set the random seed
        //        SeedRandom( parameters.ints("random_seed") );

        /*
        Put any modifications to default cell definition here if you
        want to have "inherited" by other cell types.
        
        This is a good place to set default functions.
        */
        //        initialize_default_cell_definition();
        //        cell_defaults.functions.volume_update_function = standard_volume_update_function;
        //        cell_defaults.functions.update_velocity = standard_update_cell_velocity;
        //
        //        cell_defaults.functions.update_migration_bias = NULL;
        //        cell_defaults.functions.update_phenotype = NULL;
        //        cell_defaults.functions.custom_cell_rule = NULL;
        //
        //        cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
        //        cell_defaults.functions.calculate_distance_to_membrane = NULL;

        /*
        This parses the cell definitions in the XML config file.
        */
        // Variables of plastoelastic mechanics
        //        std::vector<double> myvec(3,0.0);
        //        double[] myvec = new double[] {3, 0, 0};
        //        cell_defaults.custom_data.add_vector_variable( "ECM_attachment_point", "micron", myvec );
        //        cell_defaults.custom_data.add_vector_variable( "mechanical_strain_displacement", "micron", myvec );
        //        // assume elastic movement on the order of 10 min at maximum 10 micron elongation
        //        cell_defaults.custom_data.add_variable( "spring_constant", "1/min", 0.05 ); // 0.05    // (1.0/10.0) * (1.0/10.0)
        //        // assume plastic movement on the order of 1 day at maximum 10 micron elongation
        //        cell_defaults.custom_data.add_variable( "mechanical_relaxation_rate", "1/min", 0.0005 ); // 0.0005  // (1.0/10.0) * (1.0/(24.0*60.0)
        //        cell_defaults.custom_data.add_variable( "mechanical_strain", "micron", 0.0 );
    }

    private void setVectorData(CellDefinition cd)
    {
        double[] myvec = new double[] {3, 0, 0};
        cd.custom_data.addVectorVariable( "ECM_attachment_point", "micron", myvec.clone() );
        cd.custom_data.addVectorVariable( "mechanical_strain_displacement", "micron", myvec.clone() );
        // assume elastic movement on the order of 10 min at maximum 10 micron elongation
        cd.custom_data.addVariable( "spring_constant", "1/min", 0.05 ); // 0.05    // (1.0/10.0) * (1.0/10.0)
        // assume plastic movement on the order of 1 day at maximum 10 micron elongation
        cd.custom_data.addVariable( "mechanical_relaxation_rate", "1/min", 0.0005 ); // 0.0005  // (1.0/10.0) * (1.0/(24.0*60.0)
        cd.custom_data.addVariable( "mechanical_strain", "micron", 0.0 );
    }

    void immune_submodels_setup()
    {

        CellDefinition pCD = getCellDefinition( "CD8 Tcell" );
        pCD.functions.updatePhenotype = new CD8_Tcell_phenotype();
        pCD.functions.customCellRule = new CD8_Tcell_mechanics();// CD8_submodel_info.mechanics_function;
        pCD.functions.contact = new StandardElasticContact();// CD8_Tcell_contact_function;
        pCD.functions.updateMigration = null;
        setVectorData( pCD );
        // set up macrophages
        //        Macrophage_submodel_info = ; // much shared information
        //        // set version info
        //        Macrophage_submodel_info.name = "macrophage model";
        //        Macrophage_submodel_info.version = immune_submodels_version;
        //        // set functions
        //        Macrophage_submodel_info.main_function = NULL;
        //        Macrophage_submodel_info.phenotype_function = macrophage_phenotype;
        //        Macrophage_submodel_info.mechanics_function = macrophage_mechanics;
        //        // what microenvironment variables do you expect?
        //        // nothing unique
        //        // what custom data do I need?CD8_submodel_info
        //        //CD8_submodel_info.cell_variables.push_back( "something" );
        //        // register the submodel
        //        Macrophage_submodel_info.register_model();
        // set functions for the corresponding cell definition
        pCD = getCellDefinition( "macrophage" );
        pCD.functions.updatePhenotype = new macrophage_phenotype();
        pCD.functions.customCellRule = new macrophage_mechanics();
        pCD.functions.updateMigration = new immune_cell_motility_direction();
        setVectorData( pCD );

        // (Adrianne) set up DC submodel info
        //        DC_submodel_info = CD8_submodel_info; // much shared information
        //        // set version info
        //        DC_submodel_info.name = "DC model";
        //        DC_submodel_info.version = immune_submodels_version;
        //        // set functions
        //        DC_submodel_info.main_function = NULL;
        //        DC_submodel_info.phenotype_function = DC_phenotype;
        //        DC_submodel_info.mechanics_function = DC_mechanics;
        // what microenvironment variables do you expect?
        // nothing unique
        // what custom data do I need?
        //CD8_submodel_info.cell_variables.push_back( "something" );
        // register the submodel
        //        DC_submodel_info.register_model();
        // set functions for the corresponding cell definition
        pCD = getCellDefinition( "DC" );
        pCD.functions.updatePhenotype = new DC_phenotype();
        pCD.functions.customCellRule = new DC_mechanics();// DC_submodel_info.mechanics_function;
        pCD.functions.contact = new StandardElasticContact();
        pCD.functions.updateMigration = new immune_cell_motility_direction();
        setVectorData( pCD );
        // (Adrianne) set up CD4 Tcells ** we don't want CD4's to do anything expect be recruited to the tissue and migrate in tissue
        // set version info
        //        CD4_submodel_info = CD8_submodel_info; // much shared information
        //        CD4_submodel_info.name = "CD4 Tcell model";
        //        CD4_submodel_info.version = immune_submodels_version;
        //
        //        // set functions
        //        CD4_submodel_info.main_function = NULL;
        //        CD4_submodel_info.phenotype_function = CD4_Tcell_phenotype;
        //        CD4_submodel_info.mechanics_function = CD4_Tcell_mechanics;
        //        // what microenvironment variables do you expect?
        //        CD4_submodel_info.microenvironment_variables.push_back( "TNF" );
        //        // what custom data do I need?
        //        //CD8_submodel_info.cell_variables.push_back( "something" );
        //        // register the submodel
        //        CD4_submodel_info.register_model();
        // set functions for the corresponding cell definition
        pCD = getCellDefinition( "CD4 Tcell" );
        pCD.functions.updatePhenotype = new CD4_Tcell_plenotype();
        pCD.functions.customCellRule = new CD4_Tcell_mechanics();
        pCD.functions.updateMigration = null;
        setVectorData( pCD );

        pCD = getCellDefinition( "cancer cell" );
        pCD.functions.updatePhenotype = new cancer_phenotype();
        pCD.functions.customCellRule = new cancer_mechanics();
        pCD.functions.contact = new StandardElasticContact();
        pCD.functions.updateMigration = null;
        setVectorData( pCD );

        pCD = getCellDefinition( "lung cell" );
        pCD.functions.updatePhenotype = new epithelium_phenotype();
        pCD.functions.customCellRule = new epithelium_mechanics();
        //        pCD.functions.updateVelocity = null;
        pCD.functions.contact = new StandardElasticContact();
        pCD.functions.updateMigration = null;
        setVectorData( pCD );
    }

    void setup_tissue() throws Exception
    {
        choose_initialized_voxels();
        // create some cells near the origin

        Cell pC;

        // hexagonal cell packing
        CellDefinition pCD = getCellDefinition( "lung cell" );

        double cell_radius = pCD.phenotype.geometry.radius;
        double spacing = 0.95 * cell_radius * 2.0;

        double x_min = m.mesh.boundingBox[0] + cell_radius;
        double x_max = m.mesh.boundingBox[3] - cell_radius;

        double y_min = m.mesh.boundingBox[1] + cell_radius;
        double y_max = m.mesh.boundingBox[4] - cell_radius;

        double x = x_min;
        double y = y_min;

        double center_x = 0.5 * ( x_min + x_max );
        double center_y = 0.5 * ( y_min + y_max );

        double triangle_stagger = Math.sqrt( 3.0 ) * spacing * 0.5;

        // find the cell nearest to the center
        double nearest_distance_squared = 9e99;
        //        Cell pNearestCell = null;

        // temp_position(3, 0.0);

        int n = 0;
        while( y < y_max )
        {
            while( x < x_max )
            {
                double[] temp_position = new double[3];
                temp_position[0] = x;
                temp_position[1] = y;
                temp_position[2] = 0.0;
                valid_position.add( temp_position );

                double dx = x - center_x;
                double dy = y - center_y;

                double temp = dx * dx + dy * dy;
                if( temp < nearest_distance_squared )
                {
                    nearest_distance_squared = temp;
                    //                    pNearestCell = pC;
                }
                x += spacing;
            }
            x = x_min;

            n++;
            y += triangle_stagger;
            // in odd rows, shift
            if( n % 2 == 1 )
            {
                x += 0.5 * spacing;
            }
        }
        int Max_number_of_cell = valid_position.size();
        // extern double GridCOUNT;
        lympNode.GridCOUNT = Max_number_of_cell;
        // place immune cells
        initial_immune_cell_placement();

        int ECM_attachment_point_index = pCD.custom_data.findVectorVariableIndex( "ECM_attachment_point" );
        // cancer cells
        if( getParameterDouble( "time_add_cancer_cell" ) == 0.0 )
        {
            for( int i = 0; i < getParameterInt( "number_of_cancer_cells" ); i++ )
            {
                int index_sample = (int) ( rng.UniformRandom() * ( valid_position.size() - 1 ) );
                pC = Cell.createCell( getCellDefinition( "cancer cell" ), this, valid_position.get( index_sample ) );

                pC.customData.vectorVariables.get( ECM_attachment_point_index ).value = pC.position.clone();
                valid_position.remove( index_sample );
                //                valid_position.erase( valid_position.begin() + index_sample );
                //std::cout << "SIZE: " << valid_position.size() << " Index: " << index_sample <<  std::endl;
            }
        }
        //        // Ephitelium cells
        for( int i = 0; i < getParameterDouble( "cell_confluence_lung_cells" ) * Max_number_of_cell; i++ )
        {
            int index_sample = (int) ( rng.UniformRandom() * ( valid_position.size() - 1 ) );
            pC = Cell.createCell( getCellDefinition( "lung cell" ), this, valid_position.get( index_sample ) );
            pC.customData.vectorVariables.get( ECM_attachment_point_index ).value = pC.position.clone();
            valid_position.remove( index_sample );
            //            valid_position.erase( valid_position.begin() + index_sample );
            if( valid_position.size() == 0 )
                break;
        }

        return;
    }

    public void addCellsToMoveFromEdge(Cell cell)
    {
        cells_to_move_from_edge.add( cell );
    }

    public void process_tagged_cells_on_edge()
    {
        int cancer_type = getCellDefinition( "cancer cell" ).type;
        for( Cell pC : cells_to_move_from_edge )
        {
            if( pC.type == cancer_type )
            {
                Cell.deleteCell( pC );
            }
            else
            {
//                System.out.println( "Nudge " + pC.typeName );
                nudge_out_of_bounds_cell( pC, 10.0 );
            }
        }
        return;
    }

    public void clear_cells_to_move_from_edge()
    {
        cells_to_move_from_edge.clear();
    }

    void divide_custom_data()
    {
        int cancer_type = getCellDefinition( "cancer cell" ).type;
        double tolerance = 0.001;

        //      #pragma omp parallel for
        for( Cell pCell : m.getAgents( Cell.class ) )
        {
            //            Cell pCell;
            //            pCell = (*all_cells)[i];

            // if cell is dead or is not cancer cell, skip it
            if( pCell.phenotype.death.dead == true || pCell.type != cancer_type )
            {
                continue;
            }

            int last_cycle_index = pCell.customData.findVariableIndex( "last_cycle_entry_time" );
            int generation_index = pCell.customData.findVariableIndex( "division_generation" );
            if( pCell.phenotype.cycle.data.elapsedTimePhase < tolerance
                    && Math.abs( getCurrentTime() - pCell.customData.get( last_cycle_index ) ) >= phenotypeStep )
            {
                // add generation by 1
                pCell.customData.set( generation_index, pCell.customData.get( generation_index ) + 1 );
                pCell.customData.set( last_cycle_index, getCurrentTime() );
            }
        }

        return;
    }

    void nudge_out_of_bounds_cell(Cell pC, double tolerance)
    {
        double[] nudge = set_nudge_from_edge( pC, tolerance );
        // remove attachments
        pC.removeAllAttachedCells();
        // set velocity away rom edge
        pC.velocity = nudge;

        // set new position
        VectorUtil.prod( nudge, tolerance );
        VectorUtil.sum( pC.position, nudge );

        // update in the data structure
        pC.updateVoxelInContainer();

        // allow that cell to move and be movable
        pC.isOutOfDomain = false;
        pC.isActive = true;
        pC.isMovable = true;

        return;
    }

    double[] set_nudge_from_edge(Cell pC, double tolerance) // return {push_x,push_y,push_z} of direction to nudge cell
    {

        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Zmin = m.mesh.boundingBox[2];

        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];
        double Zmax = m.mesh.boundingBox[5];

        boolean two_dimensions = m.options.simulate2D;

        boolean setup_done = false;
        if( two_dimensions && setup_done == false )
        {
            Zmin = 0.0;
            Zmax = 0.0;
            setup_done = true;
        }

        double[] nudge = {0, 0, 0};

        if( pC.position[0] < Xmin + tolerance )
        {
            nudge[0] += 1;
        }
        if( pC.position[0] > Xmax - tolerance )
        {
            nudge[0] -= 1;
        }

        if( pC.position[1] < Ymin + tolerance )
        {
            nudge[1] += 1;
        }
        if( pC.position[1] > Ymax - tolerance )
        {
            nudge[1] -= 1;
        }

        if( two_dimensions )
        {
            VectorUtil.normalize( nudge );
            return nudge;
        }

        if( pC.position[2] < Zmin + tolerance )
        {
            nudge[2] += 1;
        }
        if( pC.position[2] > Zmax - tolerance )
        {
            nudge[2] -= 1;
        }

        VectorUtil.normalize( nudge );
        return nudge;
    }


    void include_tumor_cells() throws Exception
    {
        if( getCurrentTime() == 0
                || Math.abs( getCurrentTime() - this.getParameterDouble( "time_add_cancer_cell" ) ) >= 0.01 * diffusionStep )
            return;
        Cell pC;
        CellDefinition pCD = getCellDefinition( "cancer cell" );
        int ECM_attachment_point_index = pCD.custom_data.findVectorVariableIndex( "ECM_attachment_point" );
        for( int i = 0; i < getParameterInt( "number_of_cancer_cells" ); i++ )
        {
            int index_sample = (int) ( rng.UniformRandom() * ( valid_position.size() - 1 ) );
            pC = Cell.createCell( getCellDefinition( "cancer cell" ), this, valid_position.get( index_sample ) );
            pC.customData.vectorVariables.get( ECM_attachment_point_index ).value = pC.position.clone();
            valid_position.remove( index_sample );
            //            valid_position.erase( valid_position.begin() + index_sample );
        }
    }

    void vaccine() throws Exception
    {
        if( getParameterInt( "number_of_cells_vaccine" ) == 0 )
            return; // No vaccine
        if( Math.abs( getCurrentTime() - ( getParameterDouble( "time_vaccine" )
                + getParameterInt( "count_vaccine" ) * getParameterDouble( "vaccine_interval" ) ) ) >= 0.01 * this.diffusionStep ) // Check right time
            return;

        getParameter( "count_vaccine" ).setValue( String.valueOf( getParameterInt( "count_vaccine" ) + 1 ) );

        System.out.println( "Take shot: " + getCurrentTime() + " min - dose: " + getParameterInt( "number_of_cells_vaccine" ) + " cells" );

        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];
        double Xrange = Xmax - Xmin;
        double Yrange = Ymax - Ymin;
        Cell pC;
        CellDefinition pCD = getCellDefinition( "cancer cell" );

        for( int i = 0; i < getParameterInt( "number_of_cells_vaccine" ); i++ )
        {
            double[] position = {0, 0, 0};
            position[0] = Xmin + rng.UniformRandom() * Xrange;
            position[1] = Ymin + rng.UniformRandom() * Yrange;
            pC = Cell.createCell( getCellDefinition( "cancer cell" ), this, position );

            int apoptosis_index = pC.phenotype.death.findDeathModelIndex( "Apoptosis" );
            int debris_index = m.findDensityIndex( "debris" );
            pC.startDeath( apoptosis_index );
            pC.phenotype.volume.total = 478;
            pC.phenotype.volume.nuclear = 47.8;
            pC.phenotype.geometry.update( pC, pC.phenotype, 0.0 );
            pC.phenotype.cycle.data.setTransitionRate( 0, 1, 0.0 );
            pC.phenotype.secretion.secretionRates[debris_index] = pC.customData.get( "debris_secretion_rate" );
            pC.functions.updateVolume = null;
            pC.functions.updatePhenotype = null;
        }
    }

    int[] cell_count()
    {
        int lung_cell_type = getCellDefinition( "lung cell" ).type;
        int cancer_cell_type = getCellDefinition( "cancer cell" ).type;
        int CD8_type = getCellDefinition( "CD8 Tcell" ).type;
        int CD4_type = getCellDefinition( "CD4 Tcell" ).type;
        int macrophage_type = getCellDefinition( "macrophage" ).type;
        int DC_type = getCellDefinition( "DC" ).type;

        int[] NumberofCells = new int[16];
        for( Cell cell : this.getMicroenvironment().getAgents( Cell.class ) )
        {
            if( cell.phenotype.death.dead == true )
            {
                if( cell.type == lung_cell_type )
                    NumberofCells[1]++; // Dead lung
                else if( cell.type == cancer_cell_type )
                    NumberofCells[2]++; // Dead cancer
                else if( cell.type == DC_type )
                    NumberofCells[3]++; // Dead DC
                else if( cell.type == macrophage_type )
                    NumberofCells[4]++; // Dead macrophage
                else if( cell.type == CD4_type )
                    NumberofCells[5]++; // Dead CD4
                else if( cell.type == CD8_type )
                    NumberofCells[6]++; // Dead CD8
            }
            else if( cell.type == lung_cell_type )
            {
                NumberofCells[0]++;
            }
            else if( cell.type == cancer_cell_type )
            {
                NumberofCells[7]++;
            }
            else if( cell.type == CD8_type )
            {
                NumberofCells[15]++;
            }
            else if( cell.type == macrophage_type )
            {
                if( cell.customData.get( "ability_to_phagocytose_cancer_cell" ) == 1 )
                {
                    NumberofCells[11]++;
                } //hyperactivated macrophage
                else
                {
                    if( cell.phenotype.volume.total > ( cell.customData.get( "threshold_macrophage_volume" ) ) )
                    {
                        NumberofCells[10]++;
                    } //exhausted macrophage
                    else
                    {
                        if( cell.customData.get( "activated_immune_cell" ) > 0.5 )
                        {
                            NumberofCells[14]++;
                        } //activated macrophage
                        else
                            NumberofCells[9]++; //inactivated macrophage
                    }
                }
            }
            else if( cell.type == DC_type )
            {
                if( cell.customData.get( "activated_immune_cell" ) > 0.5 )
                {
                    NumberofCells[13]++;
                } //activated DC
                else
                    NumberofCells[8]++; //inactivated DC
            }
            else if( cell.type == CD4_type )
            {
                NumberofCells[12]++;
            }
        }
        return NumberofCells;
    }

    void initial_immune_cell_placement() throws Exception
    {
        CellDefinition pCD8 = getCellDefinition( "CD8 Tcell" );
        CellDefinition pMF = getCellDefinition( "macrophage" );
        CellDefinition pDC = getCellDefinition( "DC" );
        CellDefinition pCD4 = getCellDefinition( "CD4 Tcell" );

        // CD8+ T cells;
        for( int n = 0; n < getParameterInt( "number_of_CD8_Tcells" ); n++ )
        {
            create_infiltrating_immune_cell( pCD8 );
        }

        // macrophages
        for( int n = 0; n < getParameterInt( "number_of_macrophages" ); n++ )
        {
            create_infiltrating_immune_cell_initial( pMF );
        }

        // DC
        for( int n = 0; n < getParameterInt( "number_of_DCs" ); n++ )
        {
            create_infiltrating_immune_cell_initial( pDC );
        }

        // CD4+ T cells
        for( int n = 0; n < getParameterInt( "number_of_CD4_Tcells" ); n++ )
        {
            create_infiltrating_immune_cell_initial( pCD4 );
        }

        return;
    }

    void choose_initialized_voxels()
    {
        // read in percentage of tissue that's vascularized
        double percentage_vascularised = getParameterDouble( "percentage_tissue_vascularized" );
        int max_voxel_index = m.mesh.voxels.length - 1;
        int number_of_vascularized_voxels = (int) ( percentage_vascularised / 100.0 * ( max_voxel_index + 1 ) );

        // Sample with no replacement
        List<Integer> voxelNoVacularized = new ArrayList<>();//int[m.mesh.voxels.length];
        for( int i = 0; i < m.mesh.voxels.length; i++ )
            voxelNoVacularized.add( i );
        //        std::iota(voxelNoVacularized.begin(), voxelNoVacularized.end(), 0); // Enumerate from 0 to max_voxel_index

        // choose which voxels are veins
        for( int n = 0; n < number_of_vascularized_voxels; n++ )
        {
            int index_NoVasc_voxel = (int) ( rng.UniformRandom() * ( voxelNoVacularized.size() - 1 ) );
            vascularized_voxel_indices.add( voxelNoVacularized.get( index_NoVasc_voxel ) );
            // Remove from No vascularized vector
            voxelNoVacularized.remove( index_NoVasc_voxel );
            //            voxelNoVacularized.erase (voxelNoVacularized.begin()+index_NoVasc_voxel);
        }
    }

    double[] choose_vascularized_position()
    {
        // Randomly select the vessel to extravasate the cell
        int my_voxel_index = (int) ( rng.UniformRandom() * ( vascularized_voxel_indices.size() - 1 ) );
        int n = vascularized_voxel_indices.get( my_voxel_index );

        return m.mesh.voxels[n].center;
    }

    void create_infiltrating_immune_cell_initial(CellDefinition pCD) throws Exception
    {
        // randomly select an place from the mesh to create cell (Initial condition)
        int index_sample = (int) ( rng.UniformRandom() * ( valid_position.size() - 1 ) );
        Cell pC = Cell.createCell( pCD, this, valid_position.get( index_sample ) );
        valid_position.remove( index_sample );
        //        valid_position.erase( valid_position.begin() + index_sample );
        return;
    }

    void create_infiltrating_immune_cell(CellDefinition pCD) throws Exception
    {
        Cell pC = Cell.createCell( pCD, this, choose_vascularized_position() );
    }

    void create_infiltrating_immune_cell(String cell_name) throws Exception
    {
        create_infiltrating_immune_cell( getCellDefinition( cell_name ) );

        return;
    }

    void create_infiltrating_Tcell() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "CD8 Tcell" );
        create_infiltrating_immune_cell( pCD );

        return;
    }

    void create_infiltrating_CD4Tcell() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "CD4 Tcell" );
        create_infiltrating_immune_cell( pCD );

        return;
    }

    void create_infiltrating_DC() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "DC" );
        create_infiltrating_immune_cell( pCD );

        return;
    }

    void create_infiltrating_macrophage() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "macrophage" );
        create_infiltrating_immune_cell( pCD );

        return;
    }

    double t_immune = 0.0;
    double t_last_immune = 0.0;
    double t_next_immune = 0.0;
    
    void immune_cell_recruitment(double dt) throws Exception
    {
        int TNF_index = m.findDensityIndex( "TNF" );

        double dt_immune = getParameterDouble( "immune_dt" );

        double tolerance = 0.1 * diffusionStep;

        // is it time for the next immune recruitment?
        if( t_immune > t_next_immune - tolerance )
        {
            double elapsed_time = ( t_immune - t_last_immune );

            // macrophage recruitment
            double macrophage_recruitment_rate = getParameterDouble( "macrophage_max_recruitment_rate" );
            double M_min_signal = getParameterDouble( "macrophage_recruitment_min_signal" );
            double M_sat_signal = getParameterDouble( "macrophage_recruitment_saturation_signal" );
            double M_max_minus_min = M_sat_signal - M_min_signal;
            // integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV
            double total_rate = 0;
            double total_scaled_signal = 0.0;
            for( int n = 0; n < m.mesh.voxels.length; n++ )
            {
                // (signal(x)-signal_min)/(signal_max/signal_min)
                double dRate = ( m.get( n )[TNF_index]  - M_min_signal );
                dRate /= M_max_minus_min;
                // crop to [0,1]
                if( dRate > 1 )
                {
                    dRate = 1;
                }
                if( dRate < 0 )
                {
                    dRate = 0;
                }
                total_rate += dRate;
            }
            total_scaled_signal = total_rate;
            // multiply by dV and rate_max
            total_rate *= m.mesh.dV;
            total_rate *= macrophage_recruitment_rate;

            // expected number of new macrophages
            double number_of_new_cells_prob = total_rate * elapsed_time;
            int number_of_new_cells_int = (int)Math.floor( number_of_new_cells_prob );
            double alpha = number_of_new_cells_prob - number_of_new_cells_int;

            //STOCHASTIC PORTION
//            System.out.println( "Time "+t_immune+" Probability "+alpha );
            if( rng.UniformRandom() < alpha )
            {
                number_of_new_cells_int++;
            }
            recruited_macrophages += number_of_new_cells_int;

            if( number_of_new_cells_int > 0 )
            {
                if( t_immune < first_macrophage_recruitment_time )
                {
                    first_macrophage_recruitment_time = t_immune;
                }
                System.out.println( "\tRecruiting " + number_of_new_cells_int + " macrophages ... " );
                for( int n = 0; n < number_of_new_cells_int; n++ )
                {
                    create_infiltrating_macrophage();
                }
            }

            // CD8 Tcell recruitment (Michael) changed to take floor of ODE value
            //            extern std::vector<int>historyTc;
            int number_of_new_cells = (int)Math.floor( lympNode.TCt );
            lympNode.TCt -= number_of_new_cells;
            //            std::rotate(historyTc.rbegin(),historyTc.rbegin()+1,historyTc.rend());
            historyTc.addFront( number_of_new_cells );
            recruited_Tcells += historyTc.getBack();

            if( historyTc.getBack() > 0 )
            {
                if( t_immune < first_CD8_T_cell_recruitment_time )
                {
                    first_CD8_T_cell_recruitment_time = t_immune;
                }
                System.out.println( "\tRecruiting " + historyTc.getBack() + " CD8 T cells ... " );
                for( int n = 0; n < historyTc.getBack(); n++ )
                {
                    create_infiltrating_Tcell();
                }
            }

            // CD4 recruitment (Michael) changed to take floor of ODE value
            //            extern std::vector<int>historyTh;
            number_of_new_cells = (int)Math.floor( lympNode.Tht );

            lympNode.Tht -= number_of_new_cells;
            //            std::rotate(historyTh.rbegin(),historyTh.rbegin()+1,historyTh.rend());
            historyTh.addFront( number_of_new_cells );
            recruited_CD4Tcells += historyTh.getBack();

            if( historyTh.getBack() > 0 )
            {
                if( t_immune < first_CD4_T_cell_recruitment_time )
                {
                    first_CD4_T_cell_recruitment_time = t_immune;
                }
                System.out.println( "\tRecruiting " + historyTh.getBack() + " CD4 T cells ... " );
                for( int n = 0; n < historyTh.getBack(); n++ )
                {
                    create_infiltrating_CD4Tcell();
                }
            }

            // (Adrianne) DC recruitment - *** This section will be changed to be Tarun's model  so I've left recruitment parameters to be mac cell parameters**
            double DC_recruitment_rate = getParameterDouble( "DC_max_recruitment_rate" );
            double DC_min_signal = getParameterDouble( "DC_recruitment_min_signal" );
            double DC_sat_signal = getParameterDouble( "DC_recruitment_saturation_signal" );
            double DC_max_minus_min = DC_sat_signal - DC_min_signal;
            // integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV
            total_rate = 0;
            total_scaled_signal = 0.0;
            for( int n = 0; n < m.mesh.voxels.length; n++ )
            {
                // (signal(x)-signal_min)/(signal_max/signal_min)
                double dRate = ( m.get( n )[TNF_index] - DC_min_signal );
                dRate /= DC_max_minus_min;
                // crop to [0,1]
                if( dRate > 1 )
                {
                    dRate = 1;
                }
                if( dRate < 0 )
                {
                    dRate = 0;
                }
                total_rate += dRate;
            }
            total_scaled_signal = total_rate;
            // multiply by dV and rate_max
            total_rate *= m.mesh.dV;
            total_rate *= DC_recruitment_rate;

            // expected number of new DCs
            number_of_new_cells_prob = total_rate * elapsed_time;
            number_of_new_cells_int = (int)Math.floor( number_of_new_cells_prob );
            alpha = number_of_new_cells_prob - number_of_new_cells_int;

            //STOCHASTIC PORTION
            if( rng.UniformRandom() < alpha )
            {
                number_of_new_cells_int++;
            }
            recruited_DCs += number_of_new_cells_int;

            if( number_of_new_cells_int > 0 )
            {
                if( t_immune < first_DC_recruitment_time )
                {
                    first_DC_recruitment_time = t_immune;
                }
                System.out.println( "\tRecruiting " + number_of_new_cells_int + " DCs ... " );
                for( int n = 0; n < number_of_new_cells_int; n++ )
                {
                    create_infiltrating_DC();
                }
            }

            t_last_immune = t_immune;
            t_next_immune = t_immune + dt_immune;

        }
        t_immune += dt;

        return;
    }
}
