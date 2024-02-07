package ru.biosoft.physicell.sample_projects.cancer_immune;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.contact_function;
import ru.biosoft.physicell.core.CellFunctions.custom_cell_rule;
import ru.biosoft.physicell.core.CellFunctions.update_migration_bias;
import ru.biosoft.physicell.core.CellFunctions.update_phenotype;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Model.Event;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.Rules;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.ui.AgentVisualizer;
import ru.biosoft.physicell.ui.Visualizer;

/**
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

public class CancerImmune
{
    public static boolean use2D = true;

    public static void init(Model model) throws Exception
    {
        SignalBehavior.setup_signal_behavior_dictionaries( model.getMicroenvironment() );
        createCancerCell();
        create_immune_cell_type( model );
        setup_tissue( model );
        model.addEvent( new ImmunityEvent( model.getParameterDouble( "immune_activation_time" ) ) );
        for( Visualizer visualizer : model.getVisualizers() )
        {
            visualizer.setAgentVisualizer( new CancerImmunityVisualizer() );
        }
    }

    static void create_immune_cell_type(Model model)
    {
        Microenvironment microenvironment = model.getMicroenvironment();
        CellDefinition pImmuneCell = CellDefinition.getCellDefinition( "immune cell" );

        CellDefinition cancerCellCD = CellDefinition.getCellDefinition( "cancer cell" );

        int oxygen_ID = microenvironment.findDensityIndex( "oxygen" );
        int immuno_ID = microenvironment.findDensityIndex( "immunostimulatory factor" );

        // reduce o2 uptake 

        pImmuneCell.phenotype.secretion.uptakeRates[oxygen_ID] *= model.getParameterDouble( "immune_o2_relative_uptake" );

        pImmuneCell.phenotype.mechanics.cell_cell_adhesion_strength *= model.getParameterDouble( "immune_relative_adhesion" );
        pImmuneCell.phenotype.mechanics.cell_cell_repulsion_strength *= model.getParameterDouble( "immune_relative_repulsion" );

        // figure out mechanics parameters 

        pImmuneCell.phenotype.mechanics.relative_maximum_attachment_distance = cancerCellCD.custom_data.get( "max_attachment_distance" )
                / pImmuneCell.phenotype.geometry.radius;

        pImmuneCell.phenotype.mechanics.attachment_elastic_constant = cancerCellCD.custom_data.get( "elastic_coefficient" );

        pImmuneCell.phenotype.mechanics.relative_detachment_distance = cancerCellCD.custom_data.get( "max_attachment_distance" )
                / pImmuneCell.phenotype.geometry.radius;

        // set functions 

        pImmuneCell.functions.updatePhenotype = null;
        pImmuneCell.functions.custom_cell_rule = new immune_cell_rule();
        pImmuneCell.functions.update_migration_bias = new immune_cell_motility();
        pImmuneCell.functions.contact_function = new adhesion_contact_function();
        // set custom data values 
    }

    public static void createCancerCell() throws Exception
    {
        //        Microenvironment microenvironment = model.getMicroenvironment();
        // housekeeping 
        CellDefinition cell_defaults = CellDefinition.getCellDefinition( "cancer cell" );
        //        CellDefinition cell_defaults = StandardModels.getDefaultCellDefinition();
        cell_defaults.parameters.o2_proliferation_saturation = 38.0;
        cell_defaults.parameters.o2_reference = 38.0;

        //        int oxygen_ID = microenvironment.findDensityIndex( "oxygen" ); // 0 
        //        int immuno_ID = microenvironment.findDensityIndex( "immunostimulatory factor" ); // 1

        /*
           This parses the cell definitions in the XML config file. 
        */
        //        initialize_cell_definitions_from_pugixml();

        // change the max cell-cell adhesion distance 
        cell_defaults.phenotype.mechanics.relative_maximum_attachment_distance = cell_defaults.custom_data.get( "max_attachment_distance" )
                / cell_defaults.phenotype.geometry.radius;

        cell_defaults.phenotype.mechanics.relative_detachment_distance = cell_defaults.custom_data.get( "max_attachment_distance" )
                / cell_defaults.phenotype.geometry.radius;

        cell_defaults.phenotype.mechanics.attachment_elastic_constant = cell_defaults.custom_data.get( "elastic_coefficient" );

        cell_defaults.functions.updatePhenotype = new tumor_cell_phenotype_with_and_immune_stimulation();
        cell_defaults.functions.custom_cell_rule = null;
        cell_defaults.functions.contact_function = new adhesion_contact_function();
        cell_defaults.functions.update_migration_bias = null;
    }

    static void create_cell_types(Model model) throws Exception
    {
        Microenvironment microenvironment = model.getMicroenvironment();
        // use the same random seed so that future experiments have the 
        // same initial histogram of oncoprotein, even if threading means 
        // that future division and other events are still not identical 
        // for all runs 
        PhysiCellUtilities.setSeed( model.getParameterInt( "random_seed" ) );

        // create the immune cell type 
        create_immune_cell_type( model );

        //        build_cell_definitions_maps();

        /*
           This intializes cell signal and response dictionaries 
        */

        SignalBehavior.setup_signal_behavior_dictionaries( microenvironment );

        /*
           Cell rule definitions 
        */
        Rules.setup_cell_rules( model );


        //    display_cell_definitions( std::cout ); 
        //        
        //        return; 
    }

    void setup_microenvironment(Model model)
    {
        Microenvironment m = model.getMicroenvironment();
        if( m.options.simulate_2D == true )
        {
            //            std::cout << "Warning: overriding 2D setting to return to 3D" << std::endl; 
            m.options.simulate_2D = false;
        }
        Microenvironment.initialize_microenvironment( m );
        //        initialize_microenvironment();

        //        return; 
    }

    public static void introduce_immune_cells(Model model)
    {
        CellDefinition cd = CellDefinition.getCellDefinition( "immune cell" );
        Microenvironment m = model.getMicroenvironment();
        Set<Cell> cells = m.getAgents( Cell.class );
        double tumor_radius = -9e9; // 250.0; 
        double temp_radius = 0.0;

        // for the loop, deal with the (faster) norm squared 
        for( Cell cell : cells )//int i=0; i < (all_cells).size() ; i++ )
        {
            temp_radius = VectorUtil.norm_squared( cell.position );
            if( temp_radius > tumor_radius )
            {
                tumor_radius = temp_radius;
            }
        }
        // now square root to get to radius 
        tumor_radius = Math.sqrt( tumor_radius );

        // if this goes wackadoodle, choose 250 
        if( tumor_radius < 250.0 )
        {
            tumor_radius = 250.0;
        }

        System.out.println( "current tumor radius: " + tumor_radius );

        // now seed immune cells 

        int number_of_immune_cells = model.getParameterInt( "number_of_immune_cells" ); // 7500; // 100; // 40; 
        double radius_inner = tumor_radius + model.getParameterDouble( "initial_min_immune_distance_from_tumor" );// 30.0; // 75 // 50; 
        double radius_outer = radius_inner + model.getParameterDouble( "thickness_of_immune_seeding_region" ); // 75.0; // 100; // 1000 - 50.0; 

        double mean_radius = 0.5 * ( radius_inner + radius_outer );
        double std_radius = 0.33 * ( radius_outer - radius_inner ) / 2.0;

        if( use2D )
        {
            number_of_immune_cells /= 10;
            for( int i = 0; i < number_of_immune_cells; i++ )
            {
                double theta = PhysiCellUtilities.UniformRandom() * 6.283185307179586476925286766559;
                //                double phi = Math.acos( 2.0 * PhysiCellUtilities.UniformRandom() - 1.0 );

                double radius = PhysiCellUtilities.NormalRandom( mean_radius, std_radius );

                double[] position = new double[] {radius * Math.cos( theta ), radius * Math.sin( theta ), 0};
                Cell.createCell( cd, m, position );
            }
        }
        else
        {
        for( int i = 0; i < number_of_immune_cells; i++ )
        {
            double theta = PhysiCellUtilities.UniformRandom() * 6.283185307179586476925286766559;
            double phi = Math.acos( 2.0 * PhysiCellUtilities.UniformRandom() - 1.0 );

            double radius = PhysiCellUtilities.NormalRandom( mean_radius, std_radius );

            double[] position = new double[] {radius * Math.cos( theta ) * Math.sin( phi ), radius * Math.sin( theta ) * Math.sin( phi ),
                    radius * Math.cos( phi )};
            Cell.createCell( cd, m, position );
        }
    }
    }

    static List<double[]> create_cell_sphere_positions(double cell_radius, double sphere_radius)
    {
        List<double[]> cells = new ArrayList<>();
        int xc = 0, yc = 0, zc = 0;
        double x_spacing = cell_radius * Math.sqrt( 3 );
        double y_spacing = cell_radius * 2;
        double z_spacing = cell_radius * Math.sqrt( 3 );

        //(3,0.0);
        // std::vector<double> cylinder_center(3,0.0);
        if( use2D )
        {
            for( double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++ )
            {
                for( double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++ )
                {
                    double[] tempPoint = new double[3];
                    tempPoint[0] = x + ( zc % 2 ) * 0.5 * cell_radius;
                    tempPoint[1] = y + ( xc % 2 ) * cell_radius;
                    tempPoint[2] = 0;

                    if( Math.sqrt( VectorUtil.norm_squared( tempPoint ) ) < sphere_radius )
                    {
                        cells.add( tempPoint );
                    }
                }
            }
        }
        else
        {
            for( double z = -sphere_radius; z < sphere_radius; z += z_spacing, zc++ )
            {
                for( double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++ )
                {
                    for( double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++ )
                    {
                        double[] tempPoint = new double[3];
                        tempPoint[0] = x + ( zc % 2 ) * 0.5 * cell_radius;
                        tempPoint[1] = y + ( xc % 2 ) * cell_radius;
                        tempPoint[2] = z;

                        if( Math.sqrt( VectorUtil.norm_squared( tempPoint ) ) < sphere_radius )
                        {
                            cells.add( tempPoint );
                        }
                    }

                }
            }
        }
        return cells;
    }

    static void setup_tissue(Model model) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        // place a cluster of tumor cells at the center 
        CellDefinition cell_defaults = CellDefinition.getCellDefinition( "cancer cell" );// StandardModels.getDefaultCellDefinition();
        double cell_radius = cell_defaults.phenotype.geometry.radius;
        double cell_spacing = 0.95 * 2.0 * cell_radius;

        double tumor_radius = 250;//model.getParameterDouble( "tumor_radius" );// 250.0;  

        List<double[]> positions = create_cell_sphere_positions( cell_radius, tumor_radius );
        System.out.println( "creating " + positions.size() + " closely-packed tumor cells ... " );

        double imm_mean = model.getParameterDouble( "tumor_mean_immunogenicity" );
        double imm_sd = model.getParameterDouble( "tumor_immunogenicity_standard_deviation" );

        for( int i = 0; i < positions.size(); i++ )
        {
            Cell pCell = Cell.createCell( cell_defaults, m, positions.get( i ) ); // tumor cell 
            pCell.custom_data.set( "oncoprotein", PhysiCellUtilities.NormalRandom( imm_mean, imm_sd ) );
            if( pCell.custom_data.get( "oncoprotein" ) < 0.0 )
            {
                pCell.custom_data.set( "oncoprotein", 0.0 );
            }
        }
        double sum = 0.0;
        double min = 9e9;
        double max = -9e9;
        Set<Cell> cells = m.getAgents( Cell.class );
        for( Cell cell : cells )
        {
            double r = cell.custom_data.get( "oncoprotein" );
            sum += r;
            if( r < min )
            {
                min = r;
            }
            if( r > max )
            {
                max = r;
            }
        }
        double mean = sum / ( cells.size() + 1e-15 );
        // compute standard deviation 
        sum = 0.0;
        for( Cell cell : cells )
        {
            sum += ( cell.custom_data.get( "oncoprotein" ) - mean ) * ( cell.custom_data.get( "oncoprotein" ) - mean );
        }
        double standard_deviation = Math.sqrt( sum / ( cells.size() - 1.0 + 1e-15 ) );

        System.out.println( "Oncoprotein summary: " );
        System.out.println( "===================" );
        System.out.println( "mean: " + mean );
        System.out.println( "standard deviation: " + standard_deviation );
        System.out.println( "[min max]: [" + min + " " + max + "]" );

    }

    // custom cell phenotype function to scale immunostimulatory factor with hypoxia 
    public static class tumor_cell_phenotype_with_and_immune_stimulation implements update_phenotype
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            Microenvironment m = pCell.getMicroenvironment();
            int cycle_start_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
            int cycle_end_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
            int oncoprotein_i = pCell.custom_data.find_variable_index( "oncoprotein" );

            // update secretion rates based on hypoxia 
            int o2_index = m.findDensityIndex( "oxygen" );
            int immune_factor_index = m.findDensityIndex( "immunostimulatory factor" );
            double o2 = pCell.nearest_density_vector()[o2_index];

            phenotype.secretion.secretionRates[immune_factor_index] = 10.0;

            new StandardModels.update_cell_and_death_parameters_O2_based().execute( pCell, phenotype, dt );

            // if cell is dead, don't bother with future phenotype changes. 
            // set it to secrete the immunostimulatory factor 
            if( phenotype.death.dead == true )
            {
                phenotype.secretion.secretionRates[immune_factor_index] = 10;
                pCell.functions.updatePhenotype = null;
                return;
            }

            // multiply proliferation rate by the oncoprotein 
            //        double transitionRate = phenotype.cycle.data.
            //            double transitionRate = phenotype.cycle.data.getTransitionRate( cycle_start_index, cycle_end_index );
            //            transitionRate *= pCell.custom_data.get( oncoprotein_i );
            //            phenotype.cycle.data.setTransitionRate( cycle_start_index, cycle_end_index, transitionRate );// *= pCell.custom_data[oncoprotein_i];
            phenotype.cycle.data.modifyTransitionRate( cycle_start_index, cycle_end_index, pCell.custom_data.get( oncoprotein_i ) );
        }
    }

    public static class CancerImmunityVisualizer extends AgentVisualizer
    {
        @Override
        public Color findBorderColor(Cell cell)
        {
            return Color.black;
        }

        @Override
        public Color findColor(Cell cell)
        {
            int oncoprotein_i = cell.custom_data.find_variable_index( "oncoprotein" );
            Color result = Color.black; // immune are black

            if( cell.type == 1 )
            {
                return new Color( 50, 205, 50 );
            }

            // if I'm under attack, color me 
            if( cell.state.attachedCells.size() > 0 )
            {
                return new Color( 128, 0, 128 );
            }
            // live cells are green, but shaded by oncoprotein value 

            if( cell.phenotype.death.dead == false )
            {
                int oncoprotein = (int)Math.round( 0.5 * cell.custom_data.get( oncoprotein_i ) * 255.0 );
                return new Color( oncoprotein, oncoprotein, 255 - oncoprotein );
                //                        char szTempString [128];
                //                        sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
                //                        output[0].assign( szTempString );
                //                        output[1].assign( szTempString );

                //                        sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
                //                        output[2].assign( szTempString );

                //                        return output; 
            }

            // if not, dead colors 

            if( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.apoptotic ) // Apoptotic - Red
            {
                return new Color( 255, 0, 0 );
                //            output[0] = "rgb(255,0,0)";
                //            output[2] = "rgb(125,0,0)";
            }

            //        // Necrotic - Brown
            if( cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                    || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                    || cell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
            {
                return new Color( 250, 138, 138 );
                //                        output[0] = "rgb(250,138,38)";
                //                        output[2] = "rgb(139,69,19)";
            }
            //        
            return result;
        }
    }

    /*
    void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
    {
        std::vector<double> displacement = pAttachedTo.position - pActingOn.position; 
        axpy( &(pActingOn.velocity) , elastic_constant , displacement ); 
        
        return; 
    }
    
    void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
    {
        for( int i=0; i < pCell.state.attached_cells.size() ; i++ )
        {
            add_elastic_velocity( pCell, pCell.state.attached_cells[i], pCell.custom_data["elastic_coefficient"] ); 
        }
    
        return; 
    }   
    
    void attach_cells( Cell* pCell_1, Cell* pCell_2 )
    {
        #pragma omp critical
        {
            
        bool already_attached = false; 
        for( int i=0 ; i < pCell_1.state.attached_cells.size() ; i++ )
        {
            if( pCell_1.state.attached_cells[i] == pCell_2 )
            { already_attached = true; }
        }
        if( already_attached == false )
        { pCell_1.state.attached_cells.push_back( pCell_2 ); }
        
        already_attached = false; 
        for( int i=0 ; i < pCell_2.state.attached_cells.size() ; i++ )
        {
            if( pCell_2.state.attached_cells[i] == pCell_1 )
            { already_attached = true; }
        }
        if( already_attached == false )
        { pCell_2.state.attached_cells.push_back( pCell_1 ); }
    
        }
    
        return; 
    }
    
    void dettach_cells( Cell* pCell_1 , Cell* pCell_2 )
    {
        #pragma omp critical
        {
            bool found = false; 
            int i = 0; 
            while( !found && i < pCell_1.state.attached_cells.size() )
            {
                // if cell 2 is in cell 1's list, remove it
                if( pCell_1.state.attached_cells[i] == pCell_2 )
                {
                    int n = pCell_1.state.attached_cells.size(); 
                    // copy last entry to current position 
                    pCell_1.state.attached_cells[i] = pCell_1.state.attached_cells[n-1]; 
                    // shrink by one 
                    pCell_1.state.attached_cells.pop_back(); 
                    found = true; 
                }
                i++; 
            }
        
            found = false; 
            i = 0; 
            while( !found && i < pCell_2.state.attached_cells.size() )
            {
                // if cell 1 is in cell 2's list, remove it
                if( pCell_2.state.attached_cells[i] == pCell_1 )
                {
                    int n = pCell_2.state.attached_cells.size(); 
                    // copy last entry to current position 
                    pCell_2.state.attached_cells[i] = pCell_2.state.attached_cells[n-1]; 
                    // shrink by one 
                    pCell_2.state.attached_cells.pop_back(); 
                    found = true; 
                }
                i++; 
            }
    
        }
        
        return; 
    }
    */

    public static class immune_cell_motility implements update_migration_bias
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            Microenvironment microenvironment = pCell.getMicroenvironment();
            // if attached, biased motility towards director chemoattractant 
            // otherwise, biased motility towards cargo chemoattractant 

            int immune_factor_index = microenvironment.findDensityIndex( "immunostimulatory factor" );

            // if not docked, attempt biased chemotaxis 
            if( pCell.state.attachedCells.size() == 0 )
            {
                phenotype.motility.is_motile = true;

                phenotype.motility.migration_bias_direction = pCell.nearest_gradient( immune_factor_index ).clone();
                VectorUtil.normalize( ( phenotype.motility.migration_bias_direction ) );
            }
            else
            {
                phenotype.motility.is_motile = false;
            }
        }
    }

    public static Cell immune_cell_check_neighbors_for_attachment(Cell pAttacker, double dt)
    {
        //        std::vector<Cell> nearby = pAttacker.cells_in_my_container(); 
        for( Cell nearbyCell : pAttacker.cells_in_my_container() )
        {
            //        int i = 0; 
            //        while( i < nearby.size() )
            //        {
            // don't try to kill yourself 
            if( nearbyCell != pAttacker )
            {
                if( immune_cell_attempt_attachment( pAttacker, nearbyCell, dt ) )
                {
                    return nearbyCell;
                }
            }
            //            i++; 
        }

        return null;
    }

    static boolean immune_cell_attempt_attachment(Cell pAttacker, Cell pTarget, double dt)
    {
        int oncoprotein_i = pTarget.custom_data.find_variable_index( "oncoprotein" );
        int attach_rate_i = pAttacker.custom_data.find_variable_index( "attachment_rate" );

        double oncoprotein_saturation = pAttacker.custom_data.get( "oncoprotein_saturation" );
        double oncoprotein_threshold = pAttacker.custom_data.get( "oncoprotein_threshold" );
        double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;

        double max_attachment_distance = pAttacker.custom_data.get( "max_attachment_distance" );
        double min_attachment_distance = pAttacker.custom_data.get( "min_attachment_distance" );
        double attachment_difference = max_attachment_distance - min_attachment_distance;

        if( pTarget.custom_data.get( oncoprotein_i ) > oncoprotein_threshold && pTarget.phenotype.death.dead == false )
        {
            double[] displacement = VectorUtil.newDiff( pTarget.position, pAttacker.position );
            double distance_scale = VectorUtil.norm( displacement );
            if( distance_scale > max_attachment_distance )
            {
                return false;
            }

            double scale = pTarget.custom_data.get( oncoprotein_i );
            scale -= oncoprotein_threshold;
            scale /= oncoprotein_difference;
            if( scale > 1.0 )
            {
                scale = 1.0;
            }

            distance_scale *= -1.0;
            distance_scale += max_attachment_distance;
            distance_scale /= attachment_difference;
            if( distance_scale > 1.0 )
            {
                distance_scale = 1.0;
            }

            if( PhysiCellUtilities.UniformRandom() < pAttacker.custom_data.get( attach_rate_i ) * scale * dt * distance_scale )
            {
                //              std::cout << "\t attach!" << " " << pTarget.custom_data[oncoprotein_i] << std::endl; 
                Cell.attach_cells( pAttacker, pTarget );
            }

            return true;
        }

        return false;
    }

    static boolean immune_cell_attempt_apoptosis(Cell pAttacker, Cell pTarget, double dt)
    {
        int oncoprotein_i = pTarget.custom_data.find_variable_index( "oncoprotein" );
        int apoptosis_model_index = pTarget.phenotype.death.find_death_model_index( "apoptosis" );
        int kill_rate_index = pAttacker.custom_data.find_variable_index( "kill_rate" );

        double oncoprotein_saturation = pAttacker.custom_data.get( "oncoprotein_saturation" ); // 2.0; 
        double oncoprotein_threshold = pAttacker.custom_data.get( "oncoprotein_threshold" ); // 0.5; // 0.1; 
        double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;

        // new 
        if( pTarget.custom_data.get( oncoprotein_i ) < oncoprotein_threshold )
        {
            return false;
        }

        // new 
        double scale = pTarget.custom_data.get( oncoprotein_i );
        scale -= oncoprotein_threshold;
        scale /= oncoprotein_difference;
        if( scale > 1.0 )
        {
            scale = 1.0;
        }

        if( PhysiCellUtilities.UniformRandom() < pAttacker.custom_data.get( kill_rate_index ) * scale * dt )
        {
            //          std::cout << "\t\t kill!" << " " << pTarget.custom_data[oncoprotein_i] << std::endl; 
            return true;
        }
        return false;
    }

    static boolean immune_cell_trigger_apoptosis(Cell pAttacker, Cell pTarget)
    {
        int apoptosis_model_index = pTarget.phenotype.death.find_death_model_index( "apoptosis" );

        // if the Target cell is already dead, don't bother!
        if( pTarget.phenotype.death.dead == true )
        {
            return false;
        }
        //        System.out.println( "Start death " + pTarget.toString() );
        pTarget.start_death( apoptosis_model_index );
        return true;
    }

    public static class immune_cell_rule implements custom_cell_rule
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            int attach_lifetime_i = pCell.custom_data.find_variable_index( "attachment_lifetime" );

            if( phenotype.death.dead == true )
            {
                // the cell death functions don't automatically turn off custom functions, 
                // since those are part of mechanics. 

                // Let's just fully disable now. 
                pCell.functions.custom_cell_rule = null;
                return;
            }

            // if I'm docked
            if( pCell.state.numberAttachedCells() > 0 )
            {
                // attempt to kill my attached cell
                Cell attached = pCell.state.attachedCells.iterator().next();///[0];
                boolean detach_me = false;

                if( immune_cell_attempt_apoptosis( pCell, attached, dt ) )
                {
                    immune_cell_trigger_apoptosis( pCell, attached );
                    detach_me = true;
                }

                // decide whether to detach 

                if( PhysiCellUtilities.UniformRandom() < dt / ( pCell.custom_data.get( attach_lifetime_i ) + 1e-15 ) )
                {
                    detach_me = true;
                }

                // if I dettach, resume motile behavior 

                if( detach_me )
                {
                    Cell.detach_cells( pCell, attached );
                    phenotype.motility.is_motile = true;
                }
                return;
            }

            // I'm not docked, look for cells nearby and try to docked

            // if this returns non-NULL, we're now attached to a cell 
            if( immune_cell_check_neighbors_for_attachment( pCell, dt ) != null )
            {
                // set motility off 
                phenotype.motility.is_motile = false;
                return;
            }
            phenotype.motility.is_motile = true;
        }
    }

    public static class adhesion_contact_function implements contact_function
    {
        public void execute(Cell pActingOn, Phenotype pao, Cell pAttachedTo, Phenotype pat, double dt)
        {
            double[] displacement = VectorUtil.newDiff( pAttachedTo.position, pActingOn.position );

            double max_elastic_displacement = pao.geometry.radius * pao.mechanics.relative_detachment_distance;
            double max_displacement_squared = max_elastic_displacement * max_elastic_displacement;

            // detach cells if too far apart 

            if( VectorUtil.norm_squared( displacement ) > max_displacement_squared )
            {
                Cell.detach_cells( pActingOn, pAttachedTo );
                return;
            }
            VectorUtil.axpy( ( pActingOn.velocity ), pao.mechanics.attachment_elastic_constant, displacement );
        }
    }

    public static class ImmunityEvent extends Event
    {
        public ImmunityEvent(double executionTime)
        {
            super( executionTime );
        }

        @Override
        public void execute(Model model) throws Exception
        {
            System.out.println( "Therapy started!" );
            model.setSaveInterval( model.getParameterDouble( "save_interval_after_therapy_start" ) ); // 3.0; 
            introduce_immune_cells( model );
        }
    }
}