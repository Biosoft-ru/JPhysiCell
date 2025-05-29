package ru.biosoft.physicell.covid;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.standard.StandardElasticContact;
import ru.biosoft.physicell.core.standard.StandardUpdateVelocity;
import ru.biosoft.physicell.core.standard.StandardVolumeUpdate;
import ru.biosoft.physicell.ui.AgentColorer;

public class ModelCovid extends Model
{
    private ExternalImmune externalImmune = new ExternalImmune();
    private ReceptorDynamics receptorsDynamics = new ReceptorDynamics();

    List<Integer> vascularized_voxel_indices = new ArrayList<>();
    public Set<Cell> cells_to_move_from_edge = new HashSet<>();
    private Set<double[]> newResiduals = new HashSet<>();

    int recruited_Tcells = 0;
    int recruited_neutrophils = 0;
    int recruited_macrophages = 0;
    int recruited_CD4Tcells = 0;
    int recruited_DCs = 0;
    int recruited_fibroblasts = 0;

    double first_macrophage_recruitment_time = 9e9;
    double first_neutrophil_recruitment_time = 9e9;
    double first_CD8_T_cell_recruitment_time = 9e9;
    double first_CD4_T_cell_recruitment_time = 9e9;
    double first_DC_recruitment_time = 9e9;
    double first_fibroblast_cell_recruitment_time = 9e9;

    int DCAMOUNT = 0;
    double DM = 0; // global ICs
    double TC = 10;
    double TH1 = 1;
    double TH2 = 1;
    double TCt = 0;
    double Tht = 0;
    double Bc = 20;
    double Ps = 0;
    double Ig = 0;
    double EPICOUNT = 1;

    History historyTc = new History( 60 );
    History historyTh = new History( 60 );
    History history = new History( 72000 );

    double dt_immune;
    double t_immune;
    double t_last_immune;
    double t_next_immune;

    public void addResidual(double[] pos)
    {
        this.newResiduals.add( pos );
    }
    
    @Override
    public void stepBeforeCells()
    {
        DC_history_main_model( diffusionStep );

        externalImmune.doStep( diffusionStep );

        // receptor dynamics    
        receptorsDynamics.doStep( diffusionStep );

        // detach dead cells 
        detach_all_dead_cells( diffusionStep );

        cells_to_move_from_edge.clear();
    }

    @Override
    public void stepAfterCells() throws Exception
    {
        /*
        Custom add-ons could potentially go here. 
        */
        
        placeResiduals();
        
        process_tagged_cells_on_edge();

        move_exported_to_viral_field();

        immune_cell_recruitment( diffusionStep );
    }


    void detach_all_dead_cells(double dt)
    {
        double dt_detach = 0.1;
        double next_time = 0.0;

        double t_detach = 0.0;
        double t_last = 0.0;
        double t_next = 0.0;

        double tolerance = 0.1 * dt;

        // is it time for the next immune recruitment? 
        if( dt_detach > next_time - tolerance )
        {
            double elapsed_time = ( t_detach - t_last );

            detach_all_dead_cells();

            t_last = t_detach;
            next_time = t_detach + dt_detach;
        }
        dt_detach += dt;
    }

    void detach_all_dead_cells()
    {
        //        Cell* pC;
        //        for( int n = 0 ; n < (*all_cells).size() ; n++ )
        for( Cell cell : getMicroenvironment().getAgents( Cell.class ) )
        {
            //            pC = (*all_cells)[n]; 
            if( cell.phenotype.death.dead )
            {
                if( cell.state.neighbors.size() > 0 )
                {
                    //                    std::cout << "remove all attachments for " << pC << " " << pC.type_name << std::endl; 
                    cell.removeAllAttachedCells();
                }
            }
        }
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
        //        nudge *= tolerance;
        VectorUtil.sum( pC.position, nudge );
        //        pC.position += nudge;

        // update in the data structure 
        pC.updateVoxelInContainer();

        // allow that cell to move and be movable 
        pC.isOutOfDomain = false;
        pC.isActive = true;
        pC.isMovable = true;
    }

    void replace_out_of_bounds_cell(Cell pC, double tolerance) throws Exception
    {
        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Zmin = m.mesh.boundingBox[2];

        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];
        double Zmax = m.mesh.boundingBox[5];

        boolean setup_done = false;
        if( setup_done == false )
        {
            Xmin += tolerance;
            Ymin += tolerance;
            Zmin += tolerance;

            Xmax -= tolerance;
            Ymax -= tolerance;
            Zmax -= tolerance;

            if( m.options.simulate2D )
            {
                Zmin = 0.0;
                Zmax = 0.0;
            }
            setup_done = true;
        }

        double Xrange = Xmax - Xmin;
        double Yrange = Ymax - Ymin;
        double Zrange = Zmax - Zmin;

        double[] position = {Xmin, Ymin, Zmin}; // 
        position[0] += Xrange * getRNG().UniformRandom();
        position[1] += Yrange * getRNG().UniformRandom();
        position[2] += Zrange * getRNG().UniformRandom() + getParameterDouble( "immune_z_offset" );

        //        #pragma omp critical
        //        {
        // std::cout << "moving cell from edge " << pC << " " << pC.type_name << std::endl; 
        // create a new cell of same type 
        Cell.createCell( getCellDefinition( pC.typeName ), this, position );
        //            pNewCell.assignPosition( position ); 
        // pNewCell.custom_data = pC.custom_data; // enable in next testing 

        // get rid of the old one 
        pC.lyseCell();
        //        }    
    }

    public void placeResiduals() throws Exception
    {
        for (double[] pos: newResiduals)
        {
            create_secreting_agentcall( this, pos );
        }
        newResiduals.clear();
    }
    
    void process_tagged_cells_on_edge()
    {
        for( Cell cell : cells_to_move_from_edge )
            nudge_out_of_bounds_cell( cell, 10.0 );
    }

    // return {push_x,push_y,push_z} of direction to nudge cell 
    double[] set_nudge_from_edge(Cell pC, double tolerance)
    {
        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Zmin = m.mesh.boundingBox[2];

        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];
        double Zmax = m.mesh.boundingBox[5];

        boolean two_dimensions = m.options.simulate2D;

        boolean setup_done = false;
        if( m.options.simulate2D && !setup_done )
        {
            Zmin = 0.0;
            Zmax = 0.0;
            setup_done = true;
        }

        double[] nudge = new double[3];

        if( pC.position[0] < Xmin + tolerance )
            nudge[0] += 1;
        if( pC.position[0] > Xmax - tolerance )
            nudge[0] -= 1;

        if( pC.position[1] < Ymin + tolerance )
            nudge[1] += 1;
        if( pC.position[1] > Ymax - tolerance )
            nudge[1] -= 1;

        if( two_dimensions )
        {
            VectorUtil.normalize( nudge );
            return nudge;
        }

        if( pC.position[2] < Zmin + tolerance )
            nudge[2] += 1;
        if( pC.position[2] > Zmax - tolerance )
            nudge[2] -= 1;

        VectorUtil.normalize( nudge );
        return nudge;
    }

    void move_exported_to_viral_field()
    {
        int nV = m.findDensityIndex( "virion" );
        int nA = m.findDensityIndex( "assembled virion" );

        //        #pragma omp parallel for 
        for( int n = 0; n < m.numberVoxels(); n++ )
        {
            m.get( n )[nV] += m.get( n )[nA];
            m.get( n )[nA] = 0;
        }
    }

    void immune_cell_recruitment(double dt) throws Exception
    {
        int proinflammatory_cytokine_index = m.findDensityIndex( "pro-inflammatory cytokine" );

        int antiinflammatory_cytokine_index = m.findDensityIndex( "anti-inflammatory cytokine" );

        //        double dt_immune = getParameterDouble( "immune_dt" );
        //        double t_immune = 0.0;
        //        double t_last_immune = 0.0;
        //        double t_next_immune = 0.0;

        double tolerance = 0.1 * diffusionStep;

        // is it time for the next immune recruitment? 
        if( t_immune > t_next_immune - tolerance )
        {
            //System.out.println( "STEP " + t_immune );
            double elapsed_time = ( t_immune - t_last_immune );

            // macrophage recruitment 

            double macrophage_recruitment_rate = getParameterDouble( "macrophage_max_recruitment_rate" );
            double M_min_signal = getParameterDouble( "macrophage_recruitment_min_signal" );
            double M_sat_signal = getParameterDouble( "macrophage_recruitment_saturation_signal" );
            double M_max_minus_min = M_sat_signal - M_min_signal;

            double total_rate = 0;
            // integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV 
            double total_scaled_signal = 0.0;
            for( int n = 0; n < m.mesh.voxels.length; n++ )
            {
                // (signal(x)-signal_min)/(signal_max/signal_min)
                double dRate = ( m.get( n )[proinflammatory_cytokine_index] - M_min_signal );
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
            // multiply by dV and rate_max 
            total_scaled_signal = total_rate;

            total_rate *= m.mesh.dV;
            total_rate *= macrophage_recruitment_rate;

            // expected number of new neutrophils 
            double number_of_new_cells_prob = total_rate * elapsed_time;
            // recruited_DCs += number_of_new_cells;        

            int number_of_new_cells_int = (int)Math.floor( number_of_new_cells_prob );
            double alpha = number_of_new_cells_prob - number_of_new_cells_int;

            //STOCHASTIC PORTION        

            if( getRNG().UniformRandom() < alpha )
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

                System.out.println( "Recruiting " + number_of_new_cells_int + " macrophages ... " );

                for( int n = 0; n < number_of_new_cells_int; n++ )
                {
                    create_infiltrating_macrophage();
                }
            }

            // neutrophil recruitment 
            double neutrophil_recruitment_rate = getParameterDouble( "neutrophil_max_recruitment_rate" );
            double NR_min_signal = getParameterDouble( "neutrophil_recruitment_min_signal" );
            double NR_sat_signal = getParameterDouble( "neutrophil_recruitment_saturation_signal" );
            double NR_max_minus_min = NR_sat_signal - NR_min_signal;

            total_rate = 0;
            // integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV 
            total_scaled_signal = 0.0;
            for( int n = 0; n < m.mesh.voxels.length; n++ )
            {
                // (signal(x)-signal_min)/(signal_max/signal_min)
                double dRate = ( m.get( n )[proinflammatory_cytokine_index] - NR_min_signal );
                dRate /= NR_max_minus_min;
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
            // multiply by dV and rate_max 
            total_scaled_signal = total_rate;

            total_rate *= m.mesh.dV;
            total_rate *= neutrophil_recruitment_rate;

            // expected number of new neutrophils 
            number_of_new_cells_prob = total_rate * elapsed_time;
            // recruited_DCs += number_of_new_cells;        

            number_of_new_cells_int = (int)Math.floor( number_of_new_cells_prob );
            alpha = number_of_new_cells_prob - number_of_new_cells_int;

            //STOCHASTIC PORTION        

            if( getRNG().UniformRandom() < alpha )
            {
                number_of_new_cells_int++;
            }
            recruited_neutrophils += number_of_new_cells_int;

            if( number_of_new_cells_int > 0 )
            {
                if( t_immune < first_neutrophil_recruitment_time )
                {
                    first_neutrophil_recruitment_time = t_immune;
                }

                System.out.println( "Recruiting " + number_of_new_cells_int + " neutrophils ... " );

                for( int n = 0; n < number_of_new_cells_int; n++ )
                {
                    create_infiltrating_neutrophil();
                }
            }

            // CD8 Tcell recruitment (Michael) changed to take floor of ODE value

            //            double TCt = 0;
            //            int[] historyTc;

            int number_of_new_cells = (int)Math.floor( TCt );
            TCt -= number_of_new_cells;

            //            std::rotate(historyTc.rbegin(),historyTc.rbegin()+1,historyTc.rend());
            //            historyTc.front() = number_of_new_cells;

            historyTc.addFront( number_of_new_cells );
            //            System.out.println( "historyTc " + number_of_new_cells );
            recruited_Tcells += historyTc.getBack();


            int nAb = m.findDensityIndex( "Ig" );
            int nV = m.findDensityIndex( "virion" );

            if( historyTc.isFull() && historyTc.getBack() != 0 )
            {
                if( t_immune < first_CD8_T_cell_recruitment_time )
                {
                    first_CD8_T_cell_recruitment_time = t_immune;
                }

                System.out.println( "Recruiting " + historyTc.getBack() + " CD8 T cells ... " );

                for( int n = 0; n < historyTc.getBack(); n++ )
                {
                    create_infiltrating_Tcell();
                }
            }


            // CD4 recruitment (Michael) changed to take floor of ODE value
            //            double Tht = 0;
            //             int[] historyTh;

            number_of_new_cells = (int)Math.floor( Tht );
            Tht -= number_of_new_cells;

            historyTh.addFront( number_of_new_cells );
            //            System.out.println( "historyTh " + number_of_new_cells );
            //            std::rotate(historyTh.rbegin(),historyTh.rbegin()+1,historyTh.rend());
            //            historyTh.front() = number_of_new_cells;

            recruited_CD4Tcells += historyTh.getBack();

            if( historyTh.isFull() && historyTh.getBack() != 0 )
            {
                if( t_immune < first_CD4_T_cell_recruitment_time )
                {
                    first_CD4_T_cell_recruitment_time = t_immune;
                }

                System.out.println( "Recruiting " + historyTh.getBack() + " CD4 T cells ... " );

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

            total_rate = 0;
            // integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV 
            total_scaled_signal = 0.0;
            for( int n = 0; n < m.mesh.voxels.length; n++ )
            {
                // (signal(x)-signal_min)/(signal_max/signal_min)
                double dRate = ( m.get( n )[proinflammatory_cytokine_index] - DC_min_signal );
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
            // multiply by dV and rate_max 
            total_scaled_signal = total_rate;

            total_rate *= m.mesh.dV;
            total_rate *= DC_recruitment_rate;

            // expected number of new neutrophils 
            number_of_new_cells_prob = total_rate * elapsed_time;
            // recruited_DCs += number_of_new_cells;        

            number_of_new_cells_int = (int)Math.floor( number_of_new_cells_prob );
            alpha = number_of_new_cells_prob - number_of_new_cells_int;

            //STOCHASTIC PORTION        

            if( getRNG().UniformRandom() < alpha )
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

                System.out.println( "Recruiting " + number_of_new_cells_int + " DCs ... " );

                for( int n = 0; n < number_of_new_cells_int; n++ )
                {
                    create_infiltrating_DC();
                }
            }

            //  fibroblast recruitment
            double fibroblast_recruitment_rate = getParameterDouble( "fibroblast_max_recruitment_rate" );
            double f_min_signal = getParameterDouble( "fibroblast_recruitment_min_signal" );
            double f_sat_signal = getParameterDouble( "fibroblast_recruitment_saturation_signal" );
            double f_max_minus_min = f_sat_signal - f_min_signal;

            total_rate = 0;
            // integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV
            total_scaled_signal = 0.0;
            for( int n = 0; n < m.mesh.voxels.length; n++ )
            {
                // (signal(x)-signal_min)/(signal_max/signal_min)
                double TGF_beta = m.get( n )[antiinflammatory_cytokine_index];
                double dRate = ( 0.0492 * Math.pow( TGF_beta, 3 ) - 0.9868 * Math.pow( TGF_beta, 2 ) + 6.5408 * TGF_beta + 7
                        - f_min_signal );
                dRate /= f_max_minus_min;
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

            // multiply by dV and rate_max
            total_scaled_signal = total_rate;

            total_rate *= m.mesh.dV;
            total_rate *= fibroblast_recruitment_rate;

            // expected number of new fibroblast
            number_of_new_cells_prob = total_rate * elapsed_time;
            number_of_new_cells_int = (int)Math.floor( number_of_new_cells_prob );
            alpha = number_of_new_cells_prob - number_of_new_cells_int;

            //STOCHASTIC PORTION        

            if( getRNG().UniformRandom() < alpha )
            {
                number_of_new_cells_int++;
            }
            recruited_fibroblasts += number_of_new_cells_int;

            if( number_of_new_cells_int > 0 )
            {
                if( t_immune < first_fibroblast_cell_recruitment_time )
                {
                    first_fibroblast_cell_recruitment_time = t_immune;
                }

                System.out.println( "Recruiting " + number_of_new_cells_int + " fibroblast cells ... " );

                for( int n = 0; n < number_of_new_cells_int; n++ )
                {
                    create_infiltrating_fibroblast();
                }
            }

            t_last_immune = t_immune;
            t_next_immune = t_immune + dt_immune;

        }


        t_immune += dt;
    }

    void create_infiltrating_immune_cell(CellDefinition pCD) throws Exception
    {
        Cell.createCell( pCD, this, choose_vascularized_position() );
    }

    double[] choose_vascularized_position()
    {
        //extern std::vector<int> vascularized_voxel_indices;
        int my_voxel_index = (int) ( getRNG().UniformRandom() * ( vascularized_voxel_indices.size() - 1 ) );
        int n = vascularized_voxel_indices.get( my_voxel_index );
        return m.mesh.voxels[n].center;
    }

    void create_infiltrating_immune_cell(String cell_name) throws Exception
    {
        create_infiltrating_immune_cell( getCellDefinition( cell_name ) );
    }

    void create_infiltrating_neutrophil() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "neutrophil" );
        create_infiltrating_immune_cell( pCD );
    }

    void create_infiltrating_Tcell() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "CD8 Tcell" );
        create_infiltrating_immune_cell( pCD );
    }

    // (Adrianne) creating infiltrating CD4 cells
    void create_infiltrating_CD4Tcell() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "CD4 Tcell" );
        create_infiltrating_immune_cell( pCD );
    }

    // (Adrianne) creating infiltrating DCs
    void create_infiltrating_DC() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "DC" );
        create_infiltrating_immune_cell( pCD );
    }
    void create_infiltrating_macrophage() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "macrophage" );
        create_infiltrating_immune_cell( pCD );
    }

    void create_infiltrating_fibroblast() throws Exception
    {
        CellDefinition pCD = getCellDefinition( "fibroblast" );
        create_infiltrating_immune_cell( pCD );
        ;
    }

    void create_secreting_agent(CellDefinition pCD, Model model, double[] pos) throws Exception
    {
        Cell pC = Cell.createCell( pCD, model, pos );
        pC.isMovable = false;
    }

    /* void create_secreting_agent( std::string cell_name )
    {
        create_secreting_agent( find_cell_definition( cell_name ) ); 
        
        return;
    } */

    void create_secreting_agentcall(Model model, double[] pos) throws Exception
    {
        CellDefinition pCD = model.getCellDefinition( "residual" );
        create_secreting_agent( pCD, model, pos );
    }

    private void setup_tissue() throws Exception
    {
        int nV = m.findDensityIndex( "virion" );

        choose_initialized_voxels();

        // create some cells near the origin

        Cell pC;

        // hexagonal cell packing 
        CellDefinition pCD = getCellDefinition( "lung epithelium" );

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
        Cell pNearestCell = null;

        int n = 0;
        while( y < y_max )
        {
            while( x < x_max )
            {
                pC = Cell.createCell( getCellDefinition( "lung epithelium" ), this, new double[] {x, y, 0.0} );
                //                pC.assign_position( x,y, 0.0 );

                double dx = x - center_x;
                double dy = y - center_y;

                double temp = dx * dx + dy * dy;
                if( temp < nearest_distance_squared )
                {
                    nearest_distance_squared = temp;
                    pNearestCell = pC;
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

        EPICOUNT = m.getAgentsCount();// (*all_cells).size();
        int number_of_virions = (int) ( getParameterDouble( "multiplicity_of_infection" ) * m.getAgentsCount() );
        double single_virion_density_change = 1.0 / m.mesh.dV;

        // infect the cell closest to the center  

        if( getParameterBoolean( "use_single_infected_cell" ) == true )
        {
            //            std::cout << "Infecting center cell with one virion ... " << std::endl; 
            pNearestCell.phenotype.molecular.internSubstrates[nV] = 1.0;
        }
        else
        {
            //            std::cout << "Placing " << number_of_virions << " virions ... " << std::endl; 
            for( int l = 0; l < number_of_virions; l++ )
            {
                // pick a random voxel 
                double[] position = {0, 0, 0};
                position[0] = x_min + ( x_max - x_min ) * getRNG().UniformRandom();
                position[1] = y_min + ( y_max - y_min ) * getRNG().UniformRandom();

                int i = m.nearestVoxelIndex( position );

                // int n = (int) ( ( microenvironment.number_of_voxels()-1.0 ) * UniformRandom() ); 
                // microenvironment(i,j)[nV] += single_virion_density_change; 
                m.get( i )[nV] += single_virion_density_change;
            }
        }

        // now place immune cells 

        initial_immune_cell_placement();
    }

    void initial_immune_cell_placement() throws Exception
    {
        CellDefinition pCD8 = getCellDefinition( "CD8 Tcell" );
        CellDefinition pMF = getCellDefinition( "macrophage" );
        CellDefinition pN = getCellDefinition( "neutrophil" );
        CellDefinition pDC = getCellDefinition( "DC" );
        CellDefinition pCD4 = getCellDefinition( "CD4 Tcell" );
        CellDefinition pF = getCellDefinition( "fibroblast" );

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

        // neutrophils  
        for( int n = 0; n < getParameterInt( "number_of_neutrophils" ); n++ )
        {
            create_infiltrating_immune_cell( pN );
        }

        // DC   
        for( int n = 0; n < getParameterInt( "number_of_DCs" ); n++ )
        {
            create_infiltrating_immune_cell_initial( pDC );
        }

        // fibroblast   
        for( int n = 0; n < getParameterInt( "number_of_fibroblast" ); n++ )
        {
            create_infiltrating_immune_cell( pF );
        }

        // CD4+ T cells 
        for( int n = 0; n < getParameterInt( "number_of_CD4_Tcells" ); n++ )
        {
            create_infiltrating_immune_cell_initial( pCD4 );
        }
        return;
    }

    void create_infiltrating_immune_cell_initial(CellDefinition pCD) throws Exception
    {
        // randomly place cell intially
        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Zmin = m.mesh.boundingBox[2];

        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];
        double Zmax = m.mesh.boundingBox[5];

        if( m.options.simulate2D == true )
        {
            Zmin = 0.0;
            Zmax = 0.0;
        }

        double Xrange = ( Xmax - Xmin );
        double Yrange = ( Ymax - Ymin );
        double Zrange = ( Zmax - Zmin );

        // keep cells away from the outer edge 

        Xmin += 0.1 * Xrange;
        Ymin += 0.1 * Yrange;
        Zmin = 0;

        Xrange *= 0.8;
        Yrange *= 0.8;
        Zrange = 0.0;

        // create some of each type of cell 

        double[] position = {0, 0, 0};
        position[0] = Xmin + getRNG().UniformRandom() * Xrange;
        position[1] = Ymin + getRNG().UniformRandom() * Yrange;

        Cell.createCell( pCD, this, position );
    }

    void choose_initialized_voxels()
    {
        // read in percentage of tissue that's vascularised
        double percentage_vascularised = getParameterDouble( "perecentage_tissue_vascularized" );
        int max_voxel_index = m.mesh.voxels.length - 1;
        int number_of_vascularized_voxels = (int) ( percentage_vascularised / 100.0 * ( max_voxel_index + 1 ) );

        // choose which voxels are veins
        for( int n = 0; n < number_of_vascularized_voxels; n++ )
        {
            int index_vascularised_voxel = (int) ( getRNG().UniformRandom() * max_voxel_index );
            vascularized_voxel_indices.add( index_vascularised_voxel );
        }
    }


    public static boolean check_for_out_of_bounds(Cell pC, double tolerance)
    {
        double Xmin = pC.getMicroenvironment().mesh.boundingBox[0];
        double Ymin = pC.getMicroenvironment().mesh.boundingBox[1];
        double Zmin = pC.getMicroenvironment().mesh.boundingBox[2];

        double Xmax = pC.getMicroenvironment().mesh.boundingBox[3];
        double Ymax = pC.getMicroenvironment().mesh.boundingBox[4];
        double Zmax = pC.getMicroenvironment().mesh.boundingBox[5];

        boolean two_dimensions = pC.getMicroenvironment().options.simulate2D;

        boolean setup_done = false;
        if( pC.getMicroenvironment().options.simulate2D && !setup_done )
        {
            Zmin = 0.0;
            Zmax = 0.0;
            setup_done = true;
        }

        if( pC.position[0] < Xmin + tolerance )
        {
            return true;
        }
        if( pC.position[0] > Xmax - tolerance )
        {
            return true;
        }

        if( pC.position[1] < Ymin + tolerance )
        {
            return true;
        }
        if( pC.position[1] > Ymax - tolerance )
        {
            return true;
        }

        if( two_dimensions )
        {
            return false;
        }

        if( pC.position[2] < Zmin + tolerance )
        {
            return true;
        }
        if( pC.position[2] > Zmax - tolerance )
        {
            return true;
        }

        return false;
    }



    void DC_history_model(Cell pCell, Phenotype phenotype, double dt)
    {
        int DC_type = getCellDefinition( "DC" ).type;

        // bookkeeping -- find microenvironment variables we need

        // bookkeeping -- find custom data we need 
        double DCprob = getParameterDouble( "DC_leave_prob" );
        //         double DCAMOUNT; //declare existance of counter
        // do nothing if dead 
        if( phenotype.death.dead )
        {
            return;
        }

        // if not DC, do nothing 
        if( pCell.type != DC_type )
        {
            return;
        }

        //  if( pCell.customData.get( "activated_immune_cell" ) > 0.5 )
        //     System.out.println( "POSSBIBLE " +  pCell.customData.get( "activated_immune_cell" ));
        // (Adrianne) if DC is already activated, then check whether it leaves the tissue
        if( pCell.customData.get( "activated_immune_cell" ) > 0.5 && getRNG().UniformRandom() < DCprob )
        {
            // (Adrianne) DC leaves the tissue and so we lyse that DC
            System.out.println( "DC leaves tissue" );
            pCell.lyseCell();
            //            #pragma omp critical 
            {
                DCAMOUNT++;
            } // add one  
            return;

        }
    }

    void DC_history_main_model(double dt)
    {
        //        extern double DCAMOUNT;
        //        extern std::vector<int>history;
        //        DCAMOUNT = 0;

        //        #pragma omp parallel for 
        //        for( int n=0; n < (*all_cells).size() ; n++ )
        //        {
        for( Cell pC : m.getAgents( Cell.class ) )
        {
            //            Cell* pC = (*all_cells)[n]; 
            if( pC.phenotype.death.dead == false )
            {
                DC_history_model( pC, pC.phenotype, dt );
            }
        }
        history.addFront( DCAMOUNT );
        //        std::rotate(history.rbegin(),history.rbegin()+1,history.rend());
        //        history.front() = DCAMOUNT;

        /* std::copy(history.begin(), history.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl; 
         */
        //        return; 
    }

    @Override
    public void setupInitial() throws Exception
    {
        receptorsDynamics.setup( this );
        externalImmune.setup( this );
        create_cell_types();
        setup_tissue();

        dt_immune = getParameterDouble( "immune_dt" );
        t_immune = 0.0;
        t_last_immune = 0.0;
        t_next_immune = 0.0;
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
        for( CellDefinition cd : this.getCellDefinitions() )
        {
            cd.functions.updateVolume = new StandardVolumeUpdate();
            cd.functions.updateVelocity = new StandardUpdateVelocity();

            cd.functions.updateMigration = null;
            cd.functions.updatePhenotype = null;
            cd.functions.customCellRule = null;

            cd.functions.membraneInteraction = null;
            cd.functions.membraneDistanceCalculator = null;
        }


        int virion_index = m.findDensityIndex( "virion" );
        int assembled_virion_index = m.findDensityIndex( "assembled virion" );

        /*
           This parses the cell definitions in the XML config file. 
        */

        //        initialize_cell_definitions_from_pugixml();

        /* 
           Put any modifications to individual cell definitions here. 
           
           This is a good place to set custom functions. 
        */

        // register the submodels 
        // (which ensures that the cells have all the internal variables they need) 

        CellDefinition pCD = getCellDefinition( "lung epithelium" );
        //        pCD.phenotype.molecular.fractionReleasedDeath[virion_index] = getParameterDouble( "virus_fraction_released_at_death" );
        //        pCD.phenotype.molecular.fractionReleasedDeath[assembled_virion_index] = getParameterDouble( "virus_fraction_released_at_death" );

        pCD.phenotype.molecular.fractionReleasedDeath[virion_index] = 0;//getParameterDouble( "virus_fraction_released_at_death" ); TODO: ?
        pCD.phenotype.molecular.fractionReleasedDeath[assembled_virion_index] = 0;//getParameterDouble( "virus_fraction_released_at_death" );

        immune_submodels_setup();

        //        submodel_registry.display( std::cout ); 

        /*
           This builds the map of cell definitions and summarizes the setup. 
        */

        //        build_cell_definitions_maps(); 
        //        display_cell_definitions( std::cout ); 

        //        return; 
    }

    void immune_submodels_setup()
    {
        //            CellDefinition pCD;
        //
        //            // 
        //            // set up CD8 Tcells
        //            // set version info 
        //            CD8_submodel_info.name = "CD8 Tcell model";
        //            CD8_submodel_info.version = immune_submodels_version;
        //            // set functions 
        //            CD8_submodel_info.main_function = null;
        //            CD8_submodel_info.phenotype_function = CD8_Tcell_phenotype;
        //            CD8_submodel_info.mechanics_function = CD8_Tcell_mechanics;
        //            // what microenvironment variables do you expect? 
        //            CD8_submodel_info.microenvironment_variables.push_back( "virion" );
        //            CD8_submodel_info.microenvironment_variables.push_back( "interferon 1" );
        //            CD8_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" );
        //            CD8_submodel_info.microenvironment_variables.push_back( "chemokine" );
        //            // what custom data do I need? 
        //            //CD8_submodel_info.cell_variables.push_back( "something" ); 
        //            // register the submodel  
        //            CD8_submodel_info.register_model();
        //            // set functions for the corresponding cell definition 
        CellDefinition pCD = getCellDefinition( "CD8 Tcell" );
        pCD.functions.updatePhenotype = new CD8_Tcell_phenotype();
        pCD.functions.customCellRule = new CD8_Tcell_mechanics();
        pCD.functions.contact = new CD8_Tcell_contact_function();
        //
        //            // set up macrophages
        //            Macrophage_submodel_info = CD8_submodel_info; // much shared information 
        //            // set version info 
        //            Macrophage_submodel_info.name = "macrophage model";
        //            Macrophage_submodel_info.version = immune_submodels_version;
        //            // set functions 
        //            Macrophage_submodel_info.main_function = null;
        //            Macrophage_submodel_info.phenotype_function = macrophage_phenotype;
        //            Macrophage_submodel_info.mechanics_function = macrophage_mechanics;
        //            // what microenvironment variables do you expect? 
        //            // nothing unique 
        //            // what custom data do I need? 
        //            //CD8_submodel_info.cell_variables.push_back( "something" ); 
        //            // register the submodel  
        //            Macrophage_submodel_info.register_model();
        //            // set functions for the corresponding cell definition 
        pCD = getCellDefinition( "macrophage" );
        pCD.functions.updatePhenotype = new macrophage_phenotype();
        pCD.functions.customCellRule = new macrophage_mechanics();
        pCD.functions.updateMigration = new immune_cell_motility_direction();
        //
        //            // set up neutrophils 
        //            // set up macrophages
        //            Neutrophil_submodel_info = CD8_submodel_info; // much shared information 
        //            // set version info 
        //            Neutrophil_submodel_info.name = "neutrophil model";
        //            Neutrophil_submodel_info.version = immune_submodels_version;
        //            // set functions 
        //            Neutrophil_submodel_info.main_function = null;
        //            Neutrophil_submodel_info.phenotype_function = neutrophil_phenotype;
        //            Neutrophil_submodel_info.mechanics_function = neutrophil_mechanics;
        //            // what microenvironment variables do you expect? 
        //            // nothing unique 
        //            // what custom data do I need? 
        //            //CD8_submodel_info.cell_variables.push_back( "something" ); 
        //            // register the submodel  
        //            Neutrophil_submodel_info.register_model();
        //            // set functions for the corresponding cell definition 
        pCD = getCellDefinition( "neutrophil" );
        pCD.functions.updatePhenotype = new neutrophil_phenotype();
        pCD.functions.customCellRule = new neutrophil_mechanics();
        pCD.functions.updateMigration = new immune_cell_motility_direction();
        //
        //            // (Adrianne) set up DC submodel info
        //            DC_submodel_info = CD8_submodel_info; // much shared information 
        //            // set version info 
        //            DC_submodel_info.name = "DC model";
        //            DC_submodel_info.version = immune_submodels_version;
        //            // set functions 
        //            DC_submodel_info.main_function = null;
        //            DC_submodel_info.phenotype_function = DC_phenotype;
        //            DC_submodel_info.mechanics_function = DC_mechanics;
        //            // what microenvironment variables do you expect? 
        //            // nothing unique 
        //            // what custom data do I need? 
        //            //CD8_submodel_info.cell_variables.push_back( "something" ); 
        //            // register the submodel  
        //            DC_submodel_info.register_model();
        //            // set functions for the corresponding cell definition 
        pCD = getCellDefinition( "DC" );
        pCD.functions.updatePhenotype = new DC_phenotype();//DC_submodel_info.phenotype_function;
        pCD.functions.customCellRule = new DC_mechanics();//DC_submodel_info.mechanics_function;
        pCD.functions.updateMigration = new immune_cell_motility_direction();
        //
        //
        //            // (Adrianne) set up CD4 Tcells ** we don't want CD4's to do anything expect be recruited to the tissue and migrate in tissue
        //            // set version info 
        //            CD4_submodel_info = CD8_submodel_info; // much shared information 
        //            CD4_submodel_info.name = "CD4 Tcell model";
        //            CD4_submodel_info.version = immune_submodels_version;
        //
        //            // set functions 
        //            CD4_submodel_info.main_function = null;
        //            CD4_submodel_info.phenotype_function = CD4_Tcell_phenotype;
        //            CD4_submodel_info.mechanics_function = CD4_Tcell_mechanics;
        //            // what microenvironment variables do you expect? 
        //            CD4_submodel_info.microenvironment_variables.push_back( "virion" );
        //            CD4_submodel_info.microenvironment_variables.push_back( "interferon 1" );
        //            CD4_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" );
        //            CD4_submodel_info.microenvironment_variables.push_back( "chemokine" );
        //            // what custom data do I need? 
        //            //CD8_submodel_info.cell_variables.push_back( "something" ); 
        //            // register the submodel  
        //            CD4_submodel_info.register_model();
        //            // set functions for the corresponding cell definition 
        pCD = getCellDefinition( "CD4 Tcell" );
        pCD.functions.updatePhenotype = new CD4_Tcell_phenotype();//CD4_submodel_info.phenotype_function;
        pCD.functions.customCellRule = new CD4_Tcell_mechanics();//CD4_submodel_info.mechanics_function;
        //
        //            // set up fibroblast
        //            fibroblast_submodel_info = CD8_submodel_info; // much shared information 
        //            fibroblast_submodel_info.name = "fibroblast model";
        //            fibroblast_submodel_info.version = immune_submodels_version;
        //
        //            fibroblast_submodel_info.main_function = null;
        //            fibroblast_submodel_info.phenotype_function = fibroblast_phenotype;
        //            fibroblast_submodel_info.mechanics_function = fibroblast_mechanics;
        //
        //            fibroblast_submodel_info.register_model();
        //            // set functions for the corresponding cell definition 
        pCD = getCellDefinition( "fibroblast" );
        pCD.functions.updatePhenotype = new fibroblast_phenotype();//fibroblast_submodel_info.phenotype_function;
        pCD.functions.customCellRule = new fibroblast_mechanics();//fibroblast_submodel_info.mechanics_function;

        pCD.functions.updateMigration = new immune_cell_motility_direction();

        pCD = getCellDefinition( "lung epithelium" );
        pCD.functions.updatePhenotype = new epithelium_phenotype();
        pCD.functions.customCellRule = new epithelium_mechanics();
        pCD.functions.contact = new StandardElasticContact();
    }
    
    @Override
    public AgentColorer getDefaultColorer()
    {
        return new CovidColorer();
    }
}