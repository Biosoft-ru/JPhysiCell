package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.CellFunctions.add_cell_basement_membrane_interactions;
import ru.biosoft.physicell.core.CellFunctions.calculate_distance_to_membrane;
import ru.biosoft.physicell.core.CellFunctions.contact_function;
import ru.biosoft.physicell.core.CellFunctions.update_migration_bias;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

/*
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
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
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
public class StandardModels
{
    static boolean PhysiCell_standard_models_initialized = false;
    static boolean PhysiCell_standard_death_models_initialized = false;
    static boolean PhysiCell_standard_cycle_models_initialized = false;

    public static CycleModel Ki67_advanced = new CycleModel();
    public static CycleModel Ki67_basic = new CycleModel();
    public static CycleModel live = new CycleModel();
    public static CycleModel apoptosis = new CycleModel();
    public static CycleModel necrosis = new CycleModel();
    public static CycleModel cycling_quiescent = new CycleModel();

    public static CellDefinition defaults;// = new CycleModel();

    static DeathParameters apoptosis_parameters = new DeathParameters();
    static DeathParameters necrosis_parameters = new DeathParameters();

    // new cycle models:
    public static CycleModel flow_cytometry_cycle_model = new CycleModel();
    public static CycleModel flow_cytometry_separated_cycle_model = new CycleModel();

    public static CellDefinition createFromDefault(String name, int type, Microenvironment m) throws Exception
    {
        if( defaults == null )
            createDefaultCellDefinition();
        return defaults.clone( name, type, m );
    }

    public static CellDefinition getDefaultCellDefinition() throws Exception
    {
        if( defaults == null )
            createDefaultCellDefinition();
        return defaults;
    }

    public static CellDefinition createDefaultCellDefinition() throws Exception
    {
        defaults = new CellDefinition();
        create_standard_cycle_and_death_models(); // If the standard models have not yet been created, do so now. 

        //        result.parameters.pReference_live_phenotype = result.phenotype;

        // set up the default custom data 
        // the default Custom_Cell_Data constructor should take care of this

        // set up the default functions 
        defaults.phenotype.cycle = Ki67_advanced;
        defaults.functions.updateVolume = new standard_volume_update_function();
        //        result.functions.update_migration_bias = null;
        defaults.functions.updatePhenotype = new O2based();
        //        result.functions.custom_cell_rule = null;
        defaults.functions.updateVelocity = new standard_update_cell_velocity();
        //        result.functions.add_cell_basement_membrane_interactions = null;
        //        result.functions.calculate_distance_to_membrane = null;
        //        result.functions.set_orientation = null;
        //        cell_defaults.functions.plot_agent_SVG = standard_agent_SVG;
        //        cell_defaults.functions.plot_agent_legend = standard_agent_legend;

        // add the standard death models to the default phenotype. 
        defaults.phenotype.death.addDeathModel( 0.00319 / 60.0, apoptosis, apoptosis_parameters );
        // MCF10A, to get a 2% apoptotic index 
        defaults.phenotype.death.addDeathModel( 0.0, necrosis, necrosis_parameters );

        // set up the default phenotype (to be consistent with the default functions)
        defaults.phenotype.cycle = Ki67_advanced;//result.functions.cycleModel;

        // set molecular defaults 

        // new March 2022 : make sure Cell_Interactions and Cell_Transformations are appropriately sized. Same on motility. 
        // The Cell_Definitions constructor doesn't catch these for the cell_defaults 
        //        result.phenotype.cell_interactions.sync_to_cell_definitions();
        //        result.phenotype.cell_transformations.sync_to_cell_definitions();
        //        result.phenotype.mechanics.sync_to_cell_definitions();
        return defaults;
    }

    static void create_standard_cycle_and_death_models() throws Exception
    {
        create_standard_cell_cycle_models();
        create_standard_cell_death_models();
    }

    static boolean create_standard_cell_cycle_models() throws Exception
    {
        if( PhysiCell_standard_cycle_models_initialized )
            return false;

        create_ki67_models();
        create_live_model();
        create_cytometry_cycle_models();
        create_cycling_quiescent_model();

        PhysiCell_standard_cycle_models_initialized = true;
        if( PhysiCell_standard_death_models_initialized )
            PhysiCell_standard_models_initialized = true;
        return true;
    }

    static boolean create_standard_cell_death_models() throws Exception
    {
        if( PhysiCell_standard_death_models_initialized )
            return false;

        create_standard_apoptosis_model();
        create_standard_necrosis_model();

        PhysiCell_standard_death_models_initialized = true;
        if( PhysiCell_standard_cycle_models_initialized )
        {
            PhysiCell_standard_models_initialized = true;
        }
        return true;
    }

    public static void create_standard_apoptosis_model() throws Exception
    {
        // set default parameters for apoptosis
        apoptosis_parameters.time_units = "min";
        apoptosis_parameters.cytoplasmic_biomass_change_rate = 1.0 / 60.0;
        apoptosis_parameters.nuclear_biomass_change_rate = 0.35 / 60.0;
        apoptosis_parameters.unlysed_fluid_change_rate = 3.0 / 60.0;
        apoptosis_parameters.lysed_fluid_change_rate = 0.0;
        apoptosis_parameters.calcification_rate = 0;
        apoptosis_parameters.relative_rupture_volume = 2.0;

        // set up the apoptosis model
        apoptosis = new CycleModel();
        apoptosis.name = "Apoptosis";
        apoptosis.code = PhysiCellConstants.apoptosis_death_model;

        // add the main phase for this model, make sure it 
        // triggers the appropriate entry function, and note that 
        // it should trigger cell removal at its end 
        apoptosis.add_phase( PhysiCellConstants.apoptotic, "Apoptotic" );
        apoptosis.phases.get( 0 ).entryFunction = new PhaseEntry.standard_apoptosis_entry_function();
        apoptosis.phases.get( 0 ).removalAtExit = true;

        // add an empty junk debris phase for this model 
        apoptosis.add_phase( PhysiCellConstants.debris, "Debris" );

        // Add a link between these phases. Set the cell to be removed upon this transition. (So the "debris" phase should never be entered). 
        apoptosis.add_phase_link( 0, 1, null );
        apoptosis.setBasicTransitionRate( 0, 1, 1.0 / ( 8.6 * 60.0 ) );
        // Use the deterministic model, where this phase has fixed duration
        apoptosis.phase_link( 0, 1 ).fixedDuration = true;
    }

    public static void create_standard_necrosis_model() throws Exception
    {
        // set default parameters for necrosis
        necrosis_parameters.time_units = "min";
        necrosis_parameters.cytoplasmic_biomass_change_rate = 0.0032 / 60.0;
        necrosis_parameters.nuclear_biomass_change_rate = 0.013 / 60.0;
        necrosis_parameters.unlysed_fluid_change_rate = 0.67 / 60.0;
        necrosis_parameters.lysed_fluid_change_rate = 0.050 / 60.0;
        necrosis_parameters.calcification_rate = 0.0042 / 60.0;
        necrosis_parameters.relative_rupture_volume = 2.0;

        // set up the necrosis model 
        necrosis = new CycleModel();
        necrosis.name = "Necrosis";
        necrosis.code = PhysiCellConstants.necrosis_death_model;

        necrosis.add_phase( PhysiCellConstants.necrotic_swelling, "Necrotic (swelling)" );
        necrosis.phases.get( 0 ).entryFunction = new PhaseEntry.standard_necrosis_entry_function();

        necrosis.add_phase( PhysiCellConstants.necrotic_lysed, "Necrotic (lysed)" );
        necrosis.phases.get( 1 ).entryFunction = new PhaseEntry.standard_lysis_entry_function();
        necrosis.phases.get( 1 ).removalAtExit = true;

        // add an empty junk debris phase for this model 
        necrosis.add_phase( PhysiCellConstants.debris, "Debris" );

        necrosis.add_phase_link( 0, 1, new PhaseArrest.standard_necrosis_arrest_function() );
        necrosis.add_phase_link( 1, 2, null );

        necrosis.setBasicTransitionRate( 0, 1, 9E9 ); // set high so it's always evaluating against the "arrest"
        necrosis.setBasicTransitionRate( 1, 2, 1.0 / ( 60.0 * 24.0 * 60.0 ) ); // 60 days max  

        // Deterministically remove the necrotic cell if it has been 60 days
        necrosis.phase_link( 1, 2 ).fixedDuration = true;
    }

    static void create_ki67_models() throws Exception
    {
        createBasicKi67();
        createAdvancedKi67();
    }

    public static CycleModel createAdvancedKi67() throws Exception
    {
        Ki67_advanced = new CycleModel();
        Ki67_advanced.code = PhysiCellConstants.advanced_Ki67_cycle_model;
        Ki67_advanced.name = "Ki67 (advanced)";
        Ki67_advanced.data.timeUnits = "min";

        Ki67_advanced.add_phase( PhysiCellConstants.Ki67_negative, "Ki67-" );
        Ki67_advanced.add_phase( PhysiCellConstants.Ki67_positive_premitotic, "Ki67+ (premitotic)" );
        Ki67_advanced.add_phase( PhysiCellConstants.Ki67_positive_postmitotic, "Ki67+ (postmitotic)" );

        Ki67_advanced.phases.get( 1 ).divisionAtExit = true;

        Ki67_advanced.add_phase_link( 0, 1, null ); // - to +
        Ki67_advanced.add_phase_link( 1, 2, null ); // + (pre-mitotic) to + (post-mitotic) 
        Ki67_advanced.add_phase_link( 2, 0, null ); // + to - 

        Ki67_advanced.phase_link( 1, 2 ).fixedDuration = true;
        Ki67_advanced.phase_link( 2, 0 ).fixedDuration = true;

        Ki67_advanced.setBasicTransitionRate( 0, 1, 1.0 / ( 3.62 * 60.0 ) );// MCF10A cells ~3.62 hours in Ki67- in this fitted model
        Ki67_advanced.setBasicTransitionRate( 1, 2, 1.0 / ( 13.0 * 60.0 ) );
        Ki67_advanced.setBasicTransitionRate( 2, 0, 1.0 / ( 2.5 * 60.0 ) );

        Ki67_advanced.phases.get( 0 ).entryFunction = null; // standard_Ki67_negative_phase_entry_function;
        Ki67_advanced.phases.get( 1 ).entryFunction = new PhaseEntry.standard_Ki67_positive_phase_entry_function();
        return Ki67_advanced;
    }

    public static CycleModel createBasicKi67() throws Exception
    {
        Ki67_basic = new CycleModel();
        Ki67_basic.code = PhysiCellConstants.basic_Ki67_cycle_model;
        Ki67_basic.name = "Ki67 (basic)";
        Ki67_basic.data.timeUnits = "min";

        Ki67_basic.add_phase( PhysiCellConstants.Ki67_negative, "Ki67-" );
        Ki67_basic.add_phase( PhysiCellConstants.Ki67_positive, "Ki67+" );
        Ki67_basic.phases.get( 1 ).divisionAtExit = true;

        Ki67_basic.add_phase_link( 0, 1, null ); // - to +
        Ki67_basic.add_phase_link( 1, 0, null ); // + to - 

        Ki67_basic.setBasicTransitionRate( 0, 1, 1.0 / ( 4.59 * 60.0 ) ); // MCF10A cells are ~4.59 hours in Ki67- state
        Ki67_basic.setBasicTransitionRate( 1, 0, 1.0 / ( 15.5 * 60.0 ) );// length of Ki67+ states in advanced model 
        Ki67_basic.phase_link( 1, 0 ).fixedDuration = true;

        Ki67_basic.phases.get( 0 ).entryFunction = null; // standard_Ki67_negative_phase_entry_function;
        Ki67_basic.phases.get( 1 ).entryFunction = new PhaseEntry.standard_Ki67_positive_phase_entry_function();
        return Ki67_basic;
    }

    static void create_live_model() throws Exception
    {
        live = new CycleModel();
        live.code = PhysiCellConstants.live_cells_cycle_model;
        live.name = "Live";
        live.data.timeUnits = "min";
        live.add_phase( PhysiCellConstants.live, "Live" );
        live.phases.get( 0 ).divisionAtExit = true;
        live.add_phase_link( 0, 0, null );
        live.setBasicTransitionRate( 0, 0, 0.0432 / 60.0 ); // MCF10A have ~0.04 1/hr net birth rate
        live.phases.get( 0 ).entryFunction = new PhaseEntry.StandardLivePhaseEntry();
    }

    static void create_cytometry_cycle_models() throws Exception
    {
        // basic one first 
        flow_cytometry_cycle_model.code = PhysiCellConstants.flow_cytometry_cycle_model;
        flow_cytometry_cycle_model.name = "Flow cytometry model (basic)";

        flow_cytometry_cycle_model.data.timeUnits = "min";

        flow_cytometry_cycle_model.add_phase( PhysiCellConstants.G0G1_phase, "G0/G1" );
        flow_cytometry_cycle_model.add_phase( PhysiCellConstants.S_phase, "S" );
        flow_cytometry_cycle_model.add_phase( PhysiCellConstants.G2M_phase, "G2/M" );

        flow_cytometry_cycle_model.phases.get( 2 ).divisionAtExit = true;

        flow_cytometry_cycle_model.add_phase_link( 0, 1, null ); // G0/G1 to S
        flow_cytometry_cycle_model.add_phase_link( 1, 2, null ); // S to G2/M
        flow_cytometry_cycle_model.add_phase_link( 2, 0, null ); // G2/M to G0/G1 

        // need reference values! https://www.ncbi.nlm.nih.gov/books/NBK9876/
        flow_cytometry_cycle_model.setBasicTransitionRate( 0, 1, 0.00324 );// 5.15 hours in G0/G1 by fitting 
        flow_cytometry_cycle_model.setBasicTransitionRate( 1, 2, 0.00208 );// 8 hours in S
        flow_cytometry_cycle_model.setBasicTransitionRate( 2, 0, 0.00333 );// 5 hours in G2/M 

        flow_cytometry_cycle_model.phases.get( 0 ).entryFunction = null; //  ;
        flow_cytometry_cycle_model.phases.get( 1 ).entryFunction = new PhaseEntry.SPhaseEntry();// S_phase_entry_function; // Double nuclear volume ;
        flow_cytometry_cycle_model.phases.get( 2 ).entryFunction = null;

        // expanded flow cytometry model 
        flow_cytometry_separated_cycle_model.code = PhysiCellConstants.flow_cytometry_separated_cycle_model;
        flow_cytometry_separated_cycle_model.name = "Flow cytometry model (separated)";

        flow_cytometry_separated_cycle_model.data.timeUnits = "min";

        flow_cytometry_separated_cycle_model.add_phase( PhysiCellConstants.G0G1_phase, "G0/G1" );
        flow_cytometry_separated_cycle_model.add_phase( PhysiCellConstants.S_phase, "S" );
        flow_cytometry_separated_cycle_model.add_phase( PhysiCellConstants.G2_phase, "G2" );
        flow_cytometry_separated_cycle_model.add_phase( PhysiCellConstants.M_phase, "M" );

        flow_cytometry_separated_cycle_model.phases.get( 3 ).divisionAtExit = true;

        flow_cytometry_separated_cycle_model.add_phase_link( 0, 1, null ); // G0/G1 to S
        flow_cytometry_separated_cycle_model.add_phase_link( 1, 2, null ); // S to G2
        flow_cytometry_separated_cycle_model.add_phase_link( 2, 3, null ); // G2 to M 
        flow_cytometry_separated_cycle_model.add_phase_link( 3, 0, null ); // M to G0/G1 

        // need reference values! 
        flow_cytometry_separated_cycle_model.setBasicTransitionRate( 0, 1, 0.00335 );// 4.98 hours in G0/G1
        flow_cytometry_separated_cycle_model.setBasicTransitionRate( 1, 2, 0.00208 );// 8 hours in S 
        flow_cytometry_separated_cycle_model.setBasicTransitionRate( 2, 3, 0.00417 );// 4 hours in G2
        flow_cytometry_separated_cycle_model.setBasicTransitionRate( 3, 0, 0.0167 );// 1 hour in M 

        flow_cytometry_separated_cycle_model.phases.get( 0 ).entryFunction = null; //  ;
        flow_cytometry_separated_cycle_model.phases.get( 1 ).entryFunction = new PhaseEntry.SPhaseEntry(); // Double nuclear volume ;
        flow_cytometry_separated_cycle_model.phases.get( 2 ).entryFunction = null;
        flow_cytometry_separated_cycle_model.phases.get( 3 ).entryFunction = null;
    }

    static void create_cycling_quiescent_model() throws Exception
    {
        // Ki67_basic: 
        cycling_quiescent.code = PhysiCellConstants.cycling_quiescent_model;
        cycling_quiescent.name = "Cycling-Quiescent model";
        cycling_quiescent.data.timeUnits = "min";
        cycling_quiescent.add_phase( PhysiCellConstants.quiescent, "Quiescent" );
        cycling_quiescent.add_phase( PhysiCellConstants.cycling, "Cycling" );
        cycling_quiescent.phases.get( 1 ).divisionAtExit = true;
        cycling_quiescent.add_phase_link( 0, 1, null ); // Q to C
        cycling_quiescent.add_phase_link( 1, 0, null ); // C to Q 
        cycling_quiescent.setBasicTransitionRate( 0, 1, 1.0 / ( 4.59 * 60.0 ) ); // MCF10A cells are ~4.59 hours in Ki67- state
        cycling_quiescent.setBasicTransitionRate( 1, 0, 1.0 / ( 15.5 * 60.0 ) );// length of Ki67+ states in advanced model 
        cycling_quiescent.phase_link( 1, 0 ).fixedDuration = true;
        cycling_quiescent.phases.get( 0 ).entryFunction = null;
        cycling_quiescent.phases.get( 1 ).entryFunction = new PhaseEntry.StandardCyclingEntry();//standard_cycling_entry_function;
    }


    public static class standard_volume_update_function implements CellFunctions.volume_update_function
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            Volume v = phenotype.volume;
            v.fluid += dt * v.fluid_change_rate * ( v.target_fluid_fraction * v.total - v.fluid );
            if( v.fluid < 0.0 )
                v.fluid = 0.0;

            v.nuclear_fluid = ( v.nuclear / ( v.total + 1e-16 ) ) * ( v.fluid );
            v.cytoplasmic_fluid = v.fluid - v.nuclear_fluid;

            v.nuclear_solid += dt * v.nuclear_biomass_change_rate * ( v.target_solid_nuclear - v.nuclear_solid );
            if( v.nuclear_solid < 0.0 )
                v.nuclear_solid = 0.0;

            v.target_solid_cytoplasmic = v.target_cytoplasmic_to_nuclear_ratio * v.target_solid_nuclear;// phenotype.volume.cytoplasmic_to_nuclear_fraction * 

            v.cytoplasmic_solid += dt * v.cytoplasmic_biomass_change_rate * ( v.target_solid_cytoplasmic - v.cytoplasmic_solid );
            if( v.cytoplasmic_solid < 0.0 )
                v.cytoplasmic_solid = 0.0;

            v.solid = v.nuclear_solid + v.cytoplasmic_solid;

            v.nuclear = v.nuclear_solid + v.nuclear_fluid;
            v.cytoplasmic = v.cytoplasmic_solid + v.cytoplasmic_fluid;

            v.calcified_fraction += dt * v.calcification_rate * ( 1 - v.calcified_fraction );

            v.total = v.cytoplasmic + v.nuclear;

            v.fluid_fraction = v.fluid / ( 1e-16 + v.total );

            phenotype.geometry.update( pCell, phenotype, dt );
        }
    }

    public static class standard_update_cell_velocity implements CellFunctions.update_velocity
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            if( pCell.functions.add_cell_basement_membrane_interactions != null )
            {
                pCell.functions.add_cell_basement_membrane_interactions.execute( pCell, phenotype, dt );
            }


            pCell.state.simplePressure = 0.0;
            pCell.state.neighbors.clear(); // new 1.8.0

            //First check the neighbors in my current voxel
            //        std::vector<Cell*>::iterator neighbor;
            //        std::vector<Cell*>::iterator end = pCell.get_container().agent_grid[pCell.get_current_mechanics_voxel_index()].end();
            //        for(neighbor = pCell.get_container().agent_grid[pCell.get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
            //        {
            //            pCell.add_potentials(*neighbor);
            //        }
            CellContainer container = pCell.get_container();
            Set<Cell> neighbors = container.agent_grid.get( pCell.get_current_mechanics_voxel_index() );
            for( Cell neighbor : neighbors )
            {
                pCell.addPotentials( neighbor );
            }
            //        std::vector<int>::iterator neighbor_voxel_index;
            //        std::vector<int>::iterator neighbor_voxel_index_end = 
            //            pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()].end();               
            //        for( neighbor_voxel_index = 
            //            pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()].begin();
            //            neighbor_voxel_index != neighbor_voxel_index_end; 
            //            ++neighbor_voxel_index )
            //        {
            //            if(!is_neighbor_voxel(pCell, pCell.get_container().underlying_mesh.voxels[pCell.get_current_mechanics_voxel_index()].center, pCell.get_container().underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
            //                continue;
            //            end = pCell.get_container().agent_grid[*neighbor_voxel_index].end();
            //            for(neighbor = pCell.get_container().agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
            //            {
            //                pCell.add_potentials(neighbor);
            //            }
            //        }
            int voxelIndex = pCell.get_current_mechanics_voxel_index();
            double[] center = container.underlying_mesh.voxels[voxelIndex].center;
            int[] neighborVoxels = container.underlying_mesh.moore_connected_voxel_indices[voxelIndex];
            for( int neighbor_voxel_index : neighborVoxels )
            {
                if( !Cell.isNeighborVoxel( pCell, center, container.underlying_mesh.voxels[neighbor_voxel_index].center,
                        neighbor_voxel_index ) )
                    continue;
                for( Cell neighbor : container.agent_grid.get( neighbor_voxel_index ) )
                {
                    pCell.addPotentials( neighbor );
                }
            }
            pCell.updateMotilityVector( dt );
            VectorUtil.sum( pCell.velocity, phenotype.motility.motilityVector );
        }
    }

    public static void dynamic_spring_attachments(Cell pCell, Phenotype phenotype, double dt)
    {
        // check for detachments 
        double detachment_probability = phenotype.mechanics.detachmentRate * dt;

        for( Cell pTest : pCell.state.springAttachments )
        {
            if( PhysiCellUtilities.UniformRandom() <= detachment_probability )
            {
                Cell.detachCellsAsSpring( pCell, pTest );
            }
        }
        //        for( int j = 0; j < pCell.state.spring_attachments.size(); j++ )
        //        {
        //            Cell pTest = pCell.state.spring_attachments[j];
        //            if( PhysiCellUtilities.UniformRandom() <= detachment_probability )
        //            {
        //                detach_cells_as_spring( pCell, pTest );
        //            }
        //        }

        // check if I have max number of attachments 
        if( pCell.state.springAttachments.size() >= phenotype.mechanics.maxAttachments )
        {
            return;
        }

        // check for new attachments; 
        double attachment_probability = phenotype.mechanics.attachmentRate * dt;
        boolean done = false;
        int j = 0;

        for( Cell pTest : pCell.state.neighbors )
        {
            if( pTest.state.springAttachments.size() < pTest.phenotype.mechanics.maxAttachments )
            {
                // std::string search_string = "adhesive affinity to " + pTest.type_name; 
                // double affinity = get_single_behavior( pCell , search_string );
                double affinity = phenotype.mechanics.cell_adhesion_affinity( pTest.typeName );

                double prob = attachment_probability * affinity;
                if( PhysiCellUtilities.UniformRandom() <= prob )
                {
                    // attempt the attachment. testing for prior connection is already automated 
                    Cell.attachCellsAsSpring( pCell, pTest );
                    if( pCell.state.springAttachments.size() >= phenotype.mechanics.maxAttachments )
                    {
                        done = true;
                        break;
                    }
                }
            }
        }
        //        while( done == false && j < pCell.state.neighbors.size() )
        //        {
        //            Cell pTest = pCell.state.neighbors[j];
        //            if( pTest.state.spring_attachments.size() < pTest.phenotype.mechanics.maximum_number_of_attachments )
        //            {
        //                // std::string search_string = "adhesive affinity to " + pTest.type_name; 
        //                // double affinity = get_single_behavior( pCell , search_string );
        //                double affinity = phenotype.mechanics.cell_adhesion_affinity( pTest.type_name );
        //
        //                double prob = attachment_probability * affinity;
        //                if( PhysiCellUtilities.UniformRandom() <= prob )
        //                {
        //                    // attempt the attachment. testing for prior connection is already automated 
        //                    attach_cells_as_spring( pCell, pTest );
        //                    if( pCell.state.spring_attachments.size() >= phenotype.mechanics.maximum_number_of_attachments )
        //                    {
        //                        done = true;
        //                    }
        //                }
        //            }
        //            j++;
        //        }
    }

    public static class O2based extends UpdatePhenotype
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
        {
            // supported cycle models:
            // advanced_Ki67_cycle_model= 0;
            // basic_Ki67_cycle_model=1
            // live_cells_cycle_model = 5; 

            if( phenotype.death.dead == true )
                return;

            // set up shortcuts to find the Q and K(1) phases (assuming Ki67 basic or advanced model)
            boolean indices_initiated = false;
            int start_phase_index = 0; // Q_phase_index; 
            int end_phase_index = 0; // K_phase_index;
            int necrosis_index = 0;

            int oxygen_substrate_index = pCell.getMicroenvironment().findDensityIndex( "oxygen" );

            if( !indices_initiated )
            {
                // Ki67 models
                if( phenotype.cycle.code == PhysiCellConstants.advanced_Ki67_cycle_model
                        || phenotype.cycle.code == PhysiCellConstants.basic_Ki67_cycle_model )
                {
                    start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_negative );
                    necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

                    if( phenotype.cycle.code == PhysiCellConstants.basic_Ki67_cycle_model )
                    {
                        end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_positive );
                        indices_initiated = true;
                    }
                    if( phenotype.cycle.code == PhysiCellConstants.advanced_Ki67_cycle_model )
                    {
                        end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_positive_premitotic );
                        indices_initiated = true;
                    }
                }

                // live model 
                if( phenotype.cycle.code == PhysiCellConstants.live_cells_cycle_model )
                {
                    start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.live );
                    necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
                    end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.live );
                    indices_initiated = true;
                }

                // cytometry models 
                if( phenotype.cycle.code == PhysiCellConstants.flow_cytometry_cycle_model
                        || phenotype.cycle.code == PhysiCellConstants.flow_cytometry_separated_cycle_model )
                {
                    start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.G0G1_phase );
                    necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
                    end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.S_phase );
                    indices_initiated = true;
                }

                if( phenotype.cycle.code == PhysiCellConstants.cycling_quiescent_model )
                {
                    start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.quiescent );
                    necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
                    end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.cycling );
                    indices_initiated = true;
                }

            }

            // don't continue if we never "figured out" the current cycle model. 
            if( !indices_initiated )
                return;

            // sample the microenvironment to get the pO2 value 
            double pO2 = ( pCell.nearest_density_vector() )[oxygen_substrate_index]; // PhysiCellConstants.oxygen_index]; 
            int n = pCell.phenotype.cycle.data.currentPhaseIndex;

            // this multiplier is for linear interpolation of the oxygen value 
            double multiplier = 1.0;
            if( pO2 < pCell.parameters.o2_proliferation_saturation )
            {
                multiplier = ( pO2 - pCell.parameters.o2_proliferation_threshold )
                        / ( pCell.parameters.o2_proliferation_saturation - pCell.parameters.o2_proliferation_threshold );
            }
            if( pO2 < pCell.parameters.o2_proliferation_threshold )
            {
                multiplier = 0.0;
            }

            // now, update the appropriate cycle transition rate 
            //            pCell.parameters.pReference_live_phenotype.cycle.data.transition_rate( start_phase_index, end_phase_index )\
            //            double base = pCell.parameters.pReference_live_phenotype.cycle.data.transition_rate( start_phase_index, end_phase_index );
            //            System.out.println( pCell.toString() + " " + multiplier );
            phenotype.cycle.data.modifyTransitionRate( start_phase_index, end_phase_index, multiplier );

            // Update necrosis rate
            multiplier = 0.0;
            if( pO2 < pCell.parameters.o2_necrosis_threshold )
            {
                multiplier = ( pCell.parameters.o2_necrosis_threshold - pO2 )
                        / ( pCell.parameters.o2_necrosis_threshold - pCell.parameters.o2_necrosis_max );
            }
            if( pO2 < pCell.parameters.o2_necrosis_max )
            {
                multiplier = 1.0;
            }

            // now, update the necrosis rate 
            pCell.phenotype.death.rates.set( necrosis_index, multiplier * pCell.parameters.max_necrosis_rate );
            // check for deterministic necrosis 
            if( pCell.parameters.necrosis_type == PhysiCellConstants.deterministic_necrosis && multiplier > 1e-16 )
            {
                pCell.phenotype.death.rates.set( necrosis_index, 9e99 );
            }
        }
    }

    public static void standard_elastic_contact_function(Cell pC1, Phenotype p1, Cell pC2, Phenotype p2, double dt)
    {
        if( pC1.position.length != 3 || pC2.position.length != 3 )
            return;

        double[] displacement = VectorUtil.newDiff( pC2.position, pC1.position );
        // update May 2022 - effective adhesion 
        int ii = CellDefinition.getCellDefinitionIndex( pC1.type );
        int jj = CellDefinition.getCellDefinitionIndex( pC2.type );
        double adhesion_ii = pC1.phenotype.mechanics.attachmentElasticConstant * pC1.phenotype.mechanics.cellAdhesionAffinities[jj];
        double adhesion_jj = pC2.phenotype.mechanics.attachmentElasticConstant * pC2.phenotype.mechanics.cellAdhesionAffinities[ii];
        double effective_attachment_elastic_constant = Math.sqrt( adhesion_ii * adhesion_jj );
        VectorUtil.axpy( pC1.velocity, effective_attachment_elastic_constant, displacement );
    }

    void standard_elastic_contact_function_confluent_rest_length(Cell pC1, Phenotype p1, Cell pC2, Phenotype p2, double dt)
    {
        if( pC1.position.length != 3 || pC2.position.length != 3 )
        {
            return;
        }

        double[] displacement = VectorUtil.newDiff( pC2.position, pC1.position );

        // update May 2022 - effective adhesion 
        int ii = CellDefinition.getCellDefinitionIndex( pC1.type );
        int jj = CellDefinition.getCellDefinitionIndex( pC2.type );

        double adhesion_ii = pC1.phenotype.mechanics.attachmentElasticConstant * pC1.phenotype.mechanics.cellAdhesionAffinities[jj];
        double adhesion_jj = pC2.phenotype.mechanics.attachmentElasticConstant * pC2.phenotype.mechanics.cellAdhesionAffinities[ii];

        double effective_attachment_elastic_constant = Math.sqrt( adhesion_ii * adhesion_jj );
        // axpy( &(pC1.velocity) , effective_attachment_elastic_constant , displacement ); 

        // have the adhesion strength taper away at this rest lenght
        // set the rest length = confluent cell-cell spacing 
        // 
        double rest_length = ( p1.geometry.radius + p2.geometry.radius ) * 0.9523809523809523;

        double strength = ( VectorUtil.norm( displacement ) - rest_length ) * effective_attachment_elastic_constant;
        VectorUtil.normalize( displacement );
        VectorUtil.axpy( pC1.velocity, strength, displacement );
    }


    public static void evaluate_interactions(Cell pCell, Phenotype phenotype, double dt)
    {
        if( pCell.functions.contact_function == null )
        {
            return;
        }

        List<Cell> attached = new ArrayList<>( pCell.state.attachedCells );
        for( Cell cell : attached )
        {
            pCell.functions.contact_function.execute( pCell, phenotype, cell, cell.phenotype, dt );
        }

        //        for( int n = 0; n < pCell.state.attached_cells.size(); n++ )
        //        {
        //            pCell.functions.contact_function( pCell, phenotype, pCell.state.attached_cells[n], pCell.state.attached_cells[n].phenotype,
        //                    dt );
        //        }
    }

    public static void standard_cell_cell_interactions(Cell pCell, Phenotype phenotype, double dt)
    {
        if( phenotype.death.dead == true )
        {
            return;
        }

        //        Cell pTarget = null; 
        int type = -1;
        String type_name = "none";
        double probability = 0.0;

        boolean attacked = false;
        boolean phagocytosed = false;
        boolean fused = false;

        //        for( int n=0; n < pCell.state.neighbors.size(); n++ )
        //        {
        //            pTarget = pCell.state.neighbors[n]; 
        for( Cell pTarget : pCell.state.neighbors )
        {
            type = pTarget.type;
            type_name = pTarget.typeName;

            if( pTarget.phenotype.volume.total < 1e-15 )
            {
                break;
            }

            if( pTarget.phenotype.death.dead == true )
            {
                // dead phagocytosis 
                probability = phenotype.cellInteractions.deadPhagocytosisRate * dt;
                if( PhysiCellUtilities.UniformRandom() < probability )
                {
                    pCell.ingest_cell( pTarget );
                }
            }
            else
            {
                // live phagocytosis
                // assume you can only phagocytose one at a time for now 
                probability = phenotype.cellInteractions.getLivePhagocytosisRate( type_name ) * dt; // s[type] * dt;  
                if( PhysiCellUtilities.UniformRandom() < probability && !phagocytosed )
                {
                    pCell.ingest_cell( pTarget );
                    phagocytosed = true;
                }

                // attack 
                // assume you can only attack one cell at a time 
                // probability = phenotype.cell_interactions.attack_rate(type_name)*dt; // s[type] * dt;  

                double attack_ij = phenotype.cellInteractions.getAttackRate( type_name );
                double immunogenicity_ji = pTarget.phenotype.cellInteractions.getImmunogenicity( pCell.typeName );

                probability = attack_ij * immunogenicity_ji * dt;

                if( PhysiCellUtilities.UniformRandom() < probability && !attacked )
                {
                    pCell.attackCell( pTarget, dt );
                    attacked = true;
                }

                // fusion 
                // assume you can only fuse once cell at a time 
                probability = phenotype.cellInteractions.getFusionRate( type_name ) * dt; // s[type] * dt;  
                if( PhysiCellUtilities.UniformRandom() < probability && !fused )
                {
                    pCell.fuseCell( pTarget );
                    fused = true;
                }
            }
        }
    }

    public static void standard_cell_transformations(Cell pCell, Phenotype phenotype, double dt)
    {
        if( phenotype.death.dead == true )
        {
            return;
        }

        double probability = 0.0;
        for( int i = 0; i < phenotype.cellTransformations.transformation_rates.length; i++ )
        {
            probability = phenotype.cellTransformations.transformation_rates[i] * dt;
            if( PhysiCellUtilities.UniformRandom() <= probability )
            {
                // std::cout << "Transforming from " << pCell.type_name << " to " << cell_definitions_by_index[i].name << std::endl; 
                pCell.convert_to_cell_definition( CellDefinition.getCellDefinition( i ) );
                return;
            }
        }
    }

    public static class chemotaxis_function implements update_migration_bias
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // bias direction is gradient for the indicated substrate 
            phenotype.motility.migrationBiasDirection = pCell.nearest_gradient( phenotype.motility.chemotaxisIndex ).clone();
            // move up or down gradient based on this direction 
            VectorUtil.prod( phenotype.motility.migrationBiasDirection, phenotype.motility.chemotaxisDirection );
            VectorUtil.normalize( phenotype.motility.migrationBiasDirection );
        }
    }

    public static class advanced_chemotaxis_function_normalized implements update_migration_bias
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // We'll work directly on the migration bias direction 
            double[] pVec = phenotype.motility.migrationBiasDirection;
            // reset to zero. use memset to be faster??
            for( int i = 0; i < 3; i++ )
                pVec[i] = 0;

            // a place to put each gradient prior to normalizing it 
            double[] temp = new double[3];
            // weighted combination of the gradients 
            for( int i = 0; i < phenotype.motility.chemotacticSensitivities.length; i++ )
            {
                // get and normalize ith gradient 
                temp = pCell.nearest_gradient( i );
                VectorUtil.normalize( temp );
                VectorUtil.axpy( pVec, phenotype.motility.chemotacticSensitivities[i], temp );
            }
            // normalize that 
            VectorUtil.normalize( pVec );
        }
    }

    public static class advanced_chemotaxis_function implements update_migration_bias
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            // We'll work directly on the migration bias direction 
            double[] pVec = phenotype.motility.migrationBiasDirection;
            // reset to zero. use memset to be faster??

            for( int i = 0; i < 3; i++ )
                pVec[i] = 0;
            //            pVec = new double[3];

            // weighted combination of the gradients 
            for( int i = 0; i < phenotype.motility.chemotacticSensitivities.length; i++ )
            {
                // get and normalize ith gradient 
                VectorUtil.axpy( pVec, phenotype.motility.chemotacticSensitivities[i], pCell.nearest_gradient( i ) );
            }
            // normalize that 
            VectorUtil.normalize( pVec );
        }
    }

    public static class standard_elastic_contact_function implements contact_function
    {
        @Override
        public void execute(Cell pC1, Phenotype p1, Cell pC2, Phenotype p2, double dt)
        {
            if( pC1.position.length != 3 || pC2.position.length != 3 )
            {
                return;
            }

            double[] displacement = VectorUtil.newDiff( pC2.position, pC1.position );
            //        std::vector<double> displacement = pC2.position;
            //        displacement -= pC1.position; 

            // update May 2022 - effective adhesion 
            int ii = CellDefinition.getCellDefinitionIndex( pC1.type );
            int jj = CellDefinition.getCellDefinitionIndex( pC2.type );

            double adhesion_ii = pC1.phenotype.mechanics.attachmentElasticConstant * pC1.phenotype.mechanics.cellAdhesionAffinities[jj];
            double adhesion_jj = pC2.phenotype.mechanics.attachmentElasticConstant * pC2.phenotype.mechanics.cellAdhesionAffinities[ii];

            double effective_attachment_elastic_constant = Math.sqrt( adhesion_ii * adhesion_jj );

            // axpy( &(pC1.velocity) , p1.mechanics.attachment_elastic_constant , displacement ); 
            VectorUtil.axpy( pC1.velocity, effective_attachment_elastic_constant, displacement );
        }
    }

    public static class standard_domain_edge_avoidance_interactions implements add_cell_basement_membrane_interactions
    {
        public void execute(Cell pCell, Phenotype phenotype, double dt)
        {
            if( pCell.functions.calculate_distance_to_membrane == null )
            {
                pCell.functions.calculate_distance_to_membrane = new distance_to_domain_edge();
            }
            phenotype.mechanics.cellBMRepulsionStrength = 100;

            double max_interactive_distance = phenotype.mechanics.relMaxAdhesionDistance * phenotype.geometry.radius;
            double distance = pCell.functions.calculate_distance_to_membrane.execute( pCell, phenotype, dt );
            //Note that the distance_to_membrane function must set displacement values (as a normal vector)

            // Repulsion from basement membrane
            double temp_r = 0;
            if( distance < phenotype.geometry.radius )
            {
                temp_r = ( 1 - distance / phenotype.geometry.radius );
                temp_r *= temp_r;
                temp_r *= phenotype.mechanics.cellBMRepulsionStrength;
            }
            if( Math.abs( temp_r ) < 1e-16 )
                return;

            VectorUtil.axpy( ( pCell.velocity ), temp_r, pCell.displacement );
        }
    }

    public static class distance_to_domain_edge implements calculate_distance_to_membrane
    {
        double tolerance = 1e-7;
        double one_over_sqrt_2 = 0.70710678118;
        double one_over_sqrt_3 = 0.57735026919;
        public double execute(Cell pCell, Phenotype phenotype, double dummy)
        {
            Microenvironment m = pCell.getMicroenvironment();
            double min_distance = 9e99;
            int nearest_boundary = -1;

            // check against xL and xU
            double temp_distance = pCell.position[0] - m.mesh.boundingBox[0];
            if( temp_distance < min_distance )
            {
                min_distance = temp_distance;
                nearest_boundary = 0;
            }
            temp_distance = m.mesh.boundingBox[3] - pCell.position[0];
            if( temp_distance < min_distance )
            {
                min_distance = temp_distance;
                nearest_boundary = 1;
            }

            // check against yL and yU
            temp_distance = pCell.position[1] - m.mesh.boundingBox[1];
            if( temp_distance < min_distance )
            {
                min_distance = temp_distance;
                nearest_boundary = 2;
            }
            temp_distance = m.mesh.boundingBox[4] - pCell.position[1];
            if( temp_distance < min_distance )
            {
                min_distance = temp_distance;
                nearest_boundary = 3;
            }

            if( m.options.simulate2D == false )
            {
                // if in 3D, check against zL and zU
                temp_distance = pCell.position[2] - m.mesh.boundingBox[2];
                if( temp_distance < min_distance )
                {
                    min_distance = temp_distance;
                    nearest_boundary = 4;
                }
                temp_distance = m.mesh.boundingBox[5] - pCell.position[2];
                if( temp_distance < min_distance )
                {
                    min_distance = temp_distance;
                    nearest_boundary = 5;
                }

                // check for 3D exceptions 

                // lines 
                if( Math.abs( ( pCell.position[0] ) - ( pCell.position[1] ) ) < tolerance
                        && Math.abs( ( pCell.position[1] ) - ( pCell.position[2] ) ) < tolerance
                        && Math.abs( ( pCell.position[0] ) - ( pCell.position[2] ) ) < tolerance )
                {
                    if( pCell.position[0] > 0 )
                    {
                        if( pCell.position[0] > 0 && pCell.position[1] > 0 )
                        {
                            pCell.displacement = new double[] { -one_over_sqrt_3, -one_over_sqrt_3, -one_over_sqrt_3};
                        }
                        if( pCell.position[0] < 0 && pCell.position[1] > 0 )
                        {
                            pCell.displacement = new double[] {one_over_sqrt_3, -one_over_sqrt_3, -one_over_sqrt_3};
                        }

                        if( pCell.position[0] > 0 && pCell.position[1] < 0 )
                        {
                            pCell.displacement = new double[] { -one_over_sqrt_3, one_over_sqrt_3, -one_over_sqrt_3};
                        }
                        if( pCell.position[0] < 0 && pCell.position[1] < 0 )
                        {
                            pCell.displacement = new double[] {one_over_sqrt_3, one_over_sqrt_3, -one_over_sqrt_3};
                        }
                    }
                    else
                    {
                        if( pCell.position[0] > 0 && pCell.position[1] > 0 )
                        {
                            pCell.displacement = new double[] { -one_over_sqrt_3, -one_over_sqrt_3, one_over_sqrt_3};
                        }
                        if( pCell.position[0] < 0 && pCell.position[1] > 0 )
                        {
                            pCell.displacement = new double[] {one_over_sqrt_3, -one_over_sqrt_3, one_over_sqrt_3};
                        }

                        if( pCell.position[0] > 0 && pCell.position[1] < 0 )
                        {
                            pCell.displacement = new double[] { -one_over_sqrt_3, one_over_sqrt_3, one_over_sqrt_3};
                        }
                        if( pCell.position[0] < 0 && pCell.position[1] < 0 )
                        {
                            pCell.displacement = new double[] {one_over_sqrt_3, one_over_sqrt_3, one_over_sqrt_3};
                        }
                    }
                    return min_distance;
                }

                // planes - let's not worry for today 

            }
            else
            {
                // check for 2D  exceptions 

                if( Math.abs( ( pCell.position[0] ) - ( pCell.position[1] ) ) < tolerance )
                {
                    if( pCell.position[0] > 0 && pCell.position[1] > 0 )
                    {
                        pCell.displacement = new double[] { -one_over_sqrt_2, -one_over_sqrt_2, 0};
                    }
                    if( pCell.position[0] < 0 && pCell.position[1] > 0 )
                    {
                        pCell.displacement = new double[] {one_over_sqrt_2, -one_over_sqrt_2, 0};
                    }

                    if( pCell.position[0] > 0 && pCell.position[1] < 0 )
                    {
                        pCell.displacement = new double[] { -one_over_sqrt_2, one_over_sqrt_2, 0};
                    }
                    if( pCell.position[0] < 0 && pCell.position[1] < 0 )
                    {
                        pCell.displacement = new double[] {one_over_sqrt_2, one_over_sqrt_2, 0};
                    }
                    return min_distance;
                }
            }

            // no exceptions 
            switch( nearest_boundary )
            {
                case 0:
                    pCell.displacement = new double[] {1, 0, 0};
                    return min_distance;
                case 1:
                    pCell.displacement = new double[] { -1, 0, 0};
                    return min_distance;
                case 2:
                    pCell.displacement = new double[] {0, 1, 0};
                    return min_distance;
                case 3:
                    pCell.displacement = new double[] {0, -1, 0};
                    return min_distance;
                case 4:
                    pCell.displacement = new double[] {0, 0, 1};
                    return min_distance;
                case 5:
                    pCell.displacement = new double[] {0, 0, -1};
                    return min_distance;
                default:
                    pCell.displacement = new double[] {0, 0, 0};
                    return 9e99;
            }
            //            pCell.displacement = new double[] {0, 0, 0};
            //            return 9e99;
        }
    }
}