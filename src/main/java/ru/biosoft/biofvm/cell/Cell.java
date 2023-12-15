package ru.biosoft.biofvm.cell;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ru.biosoft.biofvm.BasicAgent;
import ru.biosoft.biofvm.Microenvironment;
import ru.biosoft.biofvm.VectorUtil;
import ru.biosoft.biofvm.cell.CellFunctions.instantiate_cell;

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
public class Cell extends BasicAgent
{
    public static Map<String, Integer> cell_definition_indices_by_name = new HashMap<>();
    public static Map<Integer, Integer> cell_definition_indices_by_type = new HashMap<>();
    public static List<CellDefinition> cell_definitions_by_index = new ArrayList<>();
    public static Map<Integer, CellDefinition> cell_definitions_by_type = new HashMap<>();

    CellContainer container;
    int current_mechanics_voxel_index;
    int updated_current_mechanics_voxel_index; // keeps the updated voxel index for later adjusting of current voxel index

    String type_name;

    CustomCellData custom_data;// = new CustomCellData();
    CellParameters parameters;// = new CellParameters();
    CellFunctions functions;// = new CellFunctions();

    CellState state = new CellState();
    public Phenotype phenotype;// = new Phenotype();
    public boolean is_out_of_domain;
    boolean is_movable;
    double[] displacement; // this should be moved to state, or made private

    public Cell()
    {
        // use the cell defaults; 
        type = StandardModels.cellDefaults.type;
        type_name = StandardModels.cellDefaults.name;
        custom_data = StandardModels.cellDefaults.custom_data;
        parameters = StandardModels.cellDefaults.parameters.clone();
        functions = StandardModels.cellDefaults.functions;
        phenotype = StandardModels.cellDefaults.phenotype.clone();
        phenotype.molecular.sync_to_cell( this );

        // cell state should be fine by the default constructor 
        current_mechanics_voxel_index = -1;
        updated_current_mechanics_voxel_index = 0;

        is_movable = true;
        is_out_of_domain = false;
        displacement = new double[3];//.resize( 3, 0.0 ); // state? 

        assign_orientation();
        container = null;
        set_total_volume( phenotype.volume.total );
    }

    @Override
    public CellContainer get_container()
    {
        if( container == null )
            container = (CellContainer)getMicroenvironment().agent_container;
        return container;
    }

    void flag_for_division()
    {
        get_container().flag_cell_for_division( this );
    }

    void flag_for_removal()
    {
        get_container().flag_cell_for_removal( this );
    }

    @Override
    public void set_total_volume(double volume)
    {
        super.set_total_volume( volume );

        // If the new volume is significantly different than the 
        // current total volume, adjust all the sub-volumes 
        // proportionally. 

        if( Math.abs( phenotype.volume.total - volume ) > 1e-16 )
        {
            double ratio = volume / ( phenotype.volume.total + 1e-16 );
            phenotype.volume.multiply_by_ratio( ratio );
        }

        phenotype.geometry.update( this, phenotype, 0.0 );

        // Here the current mechanics voxel index may not be initialized, when position is still unknown. 
        if( get_current_mechanics_voxel_index() >= 0 )
        {
            if( get_container().max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < phenotype.geometry.radius
                    * phenotype.mechanics.relative_maximum_adhesion_distance )
            {
                // get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()]= phenotype.geometry.radius*parameters.max_interaction_distance_factor;
                get_container().max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
                        * phenotype.mechanics.relative_maximum_adhesion_distance;
            }
        }
    }

    //    void start_death(int death_model_index)
    //    {
    //    // set the death data struture to the indicated death model 
    //    phenotype.death.trigger_death( death_model_index ); 
    //    // change the cycle model to the current death model 
    //    phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() ); 
    //        
    //    // turn off secretion, and reduce uptake by a factor of 10 
    //    phenotype.secretion.set_all_secretion_to_zero();
    //    phenotype.secretion.scale_all_uptake_by_factor( 0.10 );
    //        
    //    // turn off motility.
    //    phenotype.motility.is_motile = false; 
    //    phenotype.motility.motility_vector.assign( 3, 0.0 ); 
    //    functions.update_migration_bias = NULL;
    //        
    //    // make sure to run the death entry function 
    //    if( phenotype.cycle.current_phase().entry_function )
    //    {
    //        phenotype.cycle.current_phase().entry_function( this, phenotype, 0.0 ); 
    //    }
    //    
    //    return; 
    //    }
    //    
    void assign_orientation()
    {
        state.orientation = new double[3];//.orientation.resize(3,0.0);
        //TODO: check if function should be used at all
        //        if( functions.set_orientation != null )
        //        {
        //            functions.set_orientation( this, phenotype, 0.0 );
        //        }
        //        else
        //        {
        //assign a random unit vector
        double theta = Math.random() * 6.28318530717959; //rand*2*pi
        double z = 2 * Math.random() - 1;
        double temp = Math.sqrt( 1 - z * z );
        state.orientation[0] = temp * Math.cos( theta );
        state.orientation[1] = temp * Math.sin( theta );
        state.orientation[2] = z;
        //        }
    }

    @Override
    public int get_current_mechanics_voxel_index()
    {
        return current_mechanics_voxel_index;
    }

    public static Cell createCell()
    {
        return createCell( null );
    }

    public static Cell createCell(instantiate_cell custom_instantiate)
    {
        Cell pNew = custom_instantiate == null ? new Cell() : custom_instantiate.execute();
        BasicAgent.allBasicAgents.add( pNew );
        pNew.index = BasicAgent.allBasicAgents.size() - 1;
        // new usability enhancements in May 2017 
        if( Microenvironment.get_default_microenvironment() != null )
            pNew.registerMicroenvironment( Microenvironment.get_default_microenvironment() );

        // All the phenotype and other data structures are already set 
        // by virtue of the default Cell constructor. 
        pNew.set_total_volume( pNew.phenotype.volume.total );
        return pNew;
    }

    public void die()
    {
        delete_cell( this );
    }

    void update_position(double dt)
    {
        // BioFVM Basic_Agent::update_position(dt) returns without doing anything. 
        // So we remove this to avoid any future surprises. 
        // 
        // Basic_Agent::update_position(dt);

        // use Adams-Bashforth 
        double d1 = 1.5 * dt;
        double d2 = -0.5 * dt;
        //        boolean constants_defined = false; 
        //        if( constants_defined == false )
        //        {
        //            d1 = dt; 
        //            d1 *= 1.5; 
        //            d2 = dt; 
        //            d2 *= -0.5; 
        //            constants_defined = true; 
        //        }

        // new AUgust 2017
        //        if( MicroenvironmentOptions.default_microenvironment_options.simulate_2D == true )
        if( Microenvironment.get_default_microenvironment().options.simulate_2D )
        {
            velocity[2] = 0.0;
        }

        //        std::vector<double> old_position(position);

        VectorUtil.axpy( position, d1, velocity );
        VectorUtil.axpy( position, d2, previous_velocity );
        // overwrite previous_velocity for future use 
        // if(sqrt(dist(old_position, position))>3* phenotype.geometry.radius)
        // std::cout<<sqrt(dist(old_position, position))<<"old_position: "<<old_position<<", new position: "<< position<<", velocity: "<<velocity<<", previous_velocity: "<< previous_velocity<<std::endl;

        previous_velocity = velocity;

        velocity[0] = 0;
        velocity[1] = 0;
        velocity[2] = 0;
        if( get_container().underlying_mesh.isPositionValid( position[0], position[1], position[2] ) )
        {
            updated_current_mechanics_voxel_index = get_container().underlying_mesh.nearest_voxel_index( position );
        }
        else
        {
            updated_current_mechanics_voxel_index = -1;

            is_out_of_domain = true;
            isActive = false;
            is_movable = false;
        }
        return;
    }

    void delete_cell(int index)
    {
        //  std::cout << __FUNCTION__ << " " << (*all_cells)[index] 
        //  << " " << (*all_cells)[index]->type_name << std::endl; 

        Cell pDeleteMe = (Cell)BasicAgent.allBasicAgents.get( index );
        //        System.out.println( "Died " );
        // release any attached cells (as of 1.7.2 release)
        pDeleteMe.remove_all_attached_cells();
        // 1.11.0 
        pDeleteMe.remove_all_spring_attachments();

        // released internalized substrates (as of 1.5.x releases)
        pDeleteMe.release_internalized_substrates();

        // performance goal: don't delete in the middle -- very expensive reallocation
        // alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

        // move last item to index location 
        //        BasicAgent.allBasicAgents.remove( index );
        BasicAgent last = BasicAgent.allBasicAgents.get( BasicAgent.allBasicAgents.size() - 1 );
        last.index = index;
        BasicAgent.allBasicAgents.set( index, last );
        BasicAgent.allBasicAgents.remove( BasicAgent.allBasicAgents.size() - 1 );
        //        (*all_cells)[ (*all_cells).size()-1 ]->index=index;
        //        (*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];
        //        // shrink the vector
        //        (*all_cells).pop_back();    

        // deregister agent in from the agent container
        pDeleteMe.get_container().remove_agent( pDeleteMe );
        // de-allocate (delete) the cell; 
        //        delete pDeleteMe; 
    }

    void delete_cell(Cell pDelete)
    {
        delete_cell( pDelete.index );
        return;
    }

    public Cell divide()
    {
        //commented in original code
        // phenotype.flagged_for_division = false; 
        // phenotype.flagged_for_removal = false; 

        // make sure ot remove adhesions 
        remove_all_attached_cells();
        remove_all_spring_attachments();

        // version 1.10.3: 
        // conserved quantitites in custom data aer divided in half
        // so that each daughter cell gets half of the original ;
        for( int nn = 0; nn < custom_data.variables.size(); nn++ )
        {
            if( custom_data.variables.get( nn ).conserved_quantity )
            {
                custom_data.variables.get( nn ).value *= 0.5;
            }
        }
        for( int nn = 0; nn < custom_data.vector_variables.size(); nn++ )
        {
            if( custom_data.vector_variables.get( nn ).conserved_quantity )
            {
                VectorUtil.prod( custom_data.vector_variables.get( nn ).value, 0.5 );
            }
        }

        Cell child = createCell( functions.instantiate_cell );
        child.copy_data( this );
        child.copy_function_pointers( this );
        child.parameters = parameters;

        // evenly divide internalized substrates 
        // if these are not actively tracked, they are zero anyway 
        VectorUtil.prod( internalizedSubstrates, 0.5 );
        child.internalizedSubstrates = internalizedSubstrates.clone();

        // The following is already performed by create_cell(). JULY 2017 ***
        // child->register_microenvironment( get_microenvironment() );

        // randomly place the new agent close to me, accounting for orientation and 
        // polarity (if assigned)

        // May 30, 2020: 
        // Set cell_division_orientation = LegacyRandomOnUnitSphere to 
        // reproduce this code 
        /*
        double temp_angle = 6.28318530717959*UniformRandom();
        double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
        
        double radius= phenotype.geometry.radius;
        std::vector<double> rand_vec (3, 0.0);
        
        rand_vec[0]= cos( temp_angle ) * sin( temp_phi );
        rand_vec[1]= sin( temp_angle ) * sin( temp_phi );
        rand_vec[2]= cos( temp_phi );
        
        rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+ 
            rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;
        
        if( norm(rand_vec) < 1e-16 )
        {
            std::cout<<"************ERROR********************"<<std::endl;
        }
        normalize( &rand_vec ); 
        rand_vec *= radius; // multiply direction times the displacement 
        */

        double[] rand_vec = PhysiCellUtilities.UniformOnUnitSphere();//cell_division_orientation();

        double multiplier = phenotype.geometry.polarity
                * ( rand_vec[0] * state.orientation[0] + rand_vec[1] * state.orientation[1] + rand_vec[2] * state.orientation[2] );

        VectorUtil.naxpy( rand_vec, multiplier, state.orientation );
        //        rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+ 
        //            rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;   
        //        rand_vec *= phenotype.geometry.radius;
        VectorUtil.prod( rand_vec, phenotype.geometry.radius );

        child.assignPosition( position[0] + rand_vec[0], position[1] + rand_vec[1], position[2] + rand_vec[2] );

        //change my position to keep the center of mass intact 
        // and then see if I need to update my voxel index
        double negative_one_half = -0.5;
        VectorUtil.axpy( position, negative_one_half, rand_vec ); // position = position - 0.5*rand_vec; 

        //If this cell has been moved outside of the boundaries, mark it as such.
        //(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)
        if( !get_container().underlying_mesh.isPositionValid( position[0], position[1], position[2] ) )
        {
            is_out_of_domain = true;
            isActive = false;
            is_movable = false;
        }

        update_voxel_in_container();
        phenotype.volume.divide();
        child.phenotype.volume.divide();
        child.set_total_volume( child.phenotype.volume.total );//TODO: check
        set_total_volume( phenotype.volume.total );
        child.phenotype = phenotype.clone();

        if( child.phenotype.intracellular != null )
        {
            child.phenotype.intracellular.start();
            child.phenotype.intracellular.inherit( this );
        }
        // #ifdef ADDON_PHYSIDFBA
        //  child->fba_model = this->fba_model;
        // #endif
        // changes for new phenotyp March 2022
        state.damage = 0.0;
        state.total_attack_time = 0;
        child.state.damage = 0.0;
        child.state.total_attack_time = 0.0;
        return child;
    }

    void remove_all_attached_cells()
    {
        // remove self from any attached cell's list. 
        //            for( int i = 0; i < state.attached_cells.length ; i++ )
        //            {
        //                state.attached_cells[i]->detach_cell( this ); 
        //            }

        for( Cell cell : state.attached_cells )
            cell.detach_cell( this );
        state.attached_cells.clear(); // clear my list 
    }

    void remove_all_spring_attachments()
    {
        // remove self from any attached cell's list. 
        //            for( int i = 0; i < state.spring_attachments.size(); i++ )
        //            {
        //                state.spring_attachments[i].detach_cell_as_spring( this );
        //            }
        for( Cell cell : state.spring_attachments )
            cell.detach_cell_as_spring( this );
        state.spring_attachments.clear(); // clear my list 
    }

    void attach_cells(Cell pCell_1, Cell pCell_2)
    {
        pCell_1.attach_cell( pCell_2 );
        pCell_2.attach_cell( pCell_1 );
    }

    public static void attach_cells_as_spring(Cell pCell_1, Cell pCell_2)
    {
        pCell_1.attach_cell_as_spring( pCell_2 );
        pCell_2.attach_cell_as_spring( pCell_1 );
    }

    void detach_cells(Cell pCell_1, Cell pCell_2)
    {
        pCell_1.detach_cell( pCell_2 );
        pCell_2.detach_cell( pCell_1 );
    }

    public static void detach_cells_as_spring(Cell pCell_1, Cell pCell_2)
    {
        pCell_1.detach_cell_as_spring( pCell_2 );
        pCell_2.detach_cell_as_spring( pCell_1 );
    }

    void attach_cell(Cell pAddMe)
    {
        //        #pragma omp critical
        //            bool already_attached = false; 
        //            for( int i=0 ; i < state.attached_cells.size() ; i++ )
        //            {
        //                if( state.attached_cells[i] == pAddMe )
        //                { already_attached = true; }
        //            }
        //            if( already_attached == false )
        //            { state.attached_cells.push_back( pAddMe ); }
        state.attached_cells.add( pAddMe );
        // pAddMe->attach_cell( this ); 
    }

    void attach_cell_as_spring(Cell pAddMe)
    {
        //        #pragma omp critical
            //            bool already_attached = false;
            //            for( int i = 0; i < state.spring_attachments.size(); i++ )
            //            {
            //                if( state.spring_attachments[i] == pAddMe )
            //                {
            //                    already_attached = true;
            //                }
            //            }
            //            if( already_attached == false )
            //            {
            //                state.spring_attachments.push_back( pAddMe );
            //            }
            state.spring_attachments.add( pAddMe );
        // pAddMe->attach_cell( this );  
    }

    void detach_cell(Cell pRemoveMe)
    {
        //        #pragma omp critical
            //            bool found = false;
            //            int i = 0;
            //            while( !found && i < state.attached_cells.size() )
            //            {
            //                // if pRemoveMe is in the cell's list, remove it
            //                if( state.attached_cells[i] == pRemoveMe )
            //                {
            //                    int n = state.attached_cells.size();
            //                    // copy last entry to current position 
            //                    state.attached_cells[i] = state.attached_cells[n - 1];
            //                    // shrink by one 
            //                    state.attached_cells.pop_back();
            //                    found = true;
            //                }
            //                i++;
            //            }
            state.attached_cells.remove( pRemoveMe );
    }

    void detach_cell_as_spring(Cell pRemoveMe)
    {
        //        #pragma omp critical
        //        {
        //            boolean found = false;
        //            int i = 0;
        //            while( !found && i < state.spring_attachments.size() )
        //            {
        //                // if pRemoveMe is in the cell's list, remove it
        //                if( state.spring_attachments[i] == pRemoveMe )
        //                {
        //                    int n = state.spring_attachments.size();
        //                    // copy last entry to current position 
        //                    state.spring_attachments[i] = state.spring_attachments[n - 1];
        //                    // shrink by one 
        //                    state.spring_attachments.pop_back();
        //                    found = true;
        //                }
        //                i++;
        //            }
        state.spring_attachments.remove( pRemoveMe );
    }

    public void copy_data(Cell copy_me)
    {
        // phenotype=copyMe->phenotype; //it is taken care in set_phenotype
        type = copy_me.type;
        type_name = copy_me.type_name;

        custom_data = copy_me.custom_data;
        parameters = copy_me.parameters.clone();

        velocity = copy_me.velocity.clone();
        // expected_phenotype = copy_me. expected_phenotype; //it is taken care in set_phenotype
        sourceSinkTemp1 = copy_me.sourceSinkTemp1.clone();
        sourceSinkTemp2 = copy_me.sourceSinkTemp2.clone();
        //        cell_source_sink_solver_temp1 = std::vector<double>(copy_me.cell_source_sink_solver_temp1);
        //        cell_source_sink_solver_temp2 = std::vector<double>(copy_me.cell_source_sink_solver_temp2);
    }

    void copy_function_pointers(Cell copy_me)
    {
        functions = copy_me.functions;
    }

    void update_voxel_in_container()
    {
        // call the method from BioFVM_basic_agent to update microenvironment's voxel index
        updateVoxelIndex();
        // int temp_current_voxel_index;
        // Check to see if we need to remove agents that are pushed out of boundary
        // if(!get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2])) 

        if( updated_current_mechanics_voxel_index == -1 )// updated_current_mechanics_voxel_index is updated in update_position
        {
            // check if this agent has a valid voxel index, if so, remove it from previous voxel
            if( get_current_mechanics_voxel_index() >= 0 )
            {
                {
                    get_container().remove_agent_from_voxel( this, get_current_mechanics_voxel_index() );
                }
            }
            {
                get_container().add_agent_to_outer_voxel( this );
            }
            // std::cout<<"cell out of boundary..."<< __LINE__<<" "<<ID<<std::endl;
            current_mechanics_voxel_index = -1;
            is_out_of_domain = true;
            isActive = false;
            return;
        }

        // temp_current_voxel_index= get_current_mechanics_voxel_index();
        // updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );

        // update mesh indices (if needed)
        if( updated_current_mechanics_voxel_index != get_current_mechanics_voxel_index() )
        {
            {
                container.remove_agent_from_voxel( this, get_current_mechanics_voxel_index() );
                container.add_agent_to_voxel( this, updated_current_mechanics_voxel_index );
            }
            current_mechanics_voxel_index = updated_current_mechanics_voxel_index;
        }
    }

    @Override
    public boolean assignPosition(double[] new_position)
    {
        return assignPosition( new_position[0], new_position[1], new_position[2] );
    }

    void setPrevious_velocity(double xV, double yV, double zV)
    {
        previous_velocity[0] = xV;
        previous_velocity[1] = yV;
        previous_velocity[2] = zV;
    }

    @Override
    public boolean assignPosition(double x, double y, double z)
    {
        position[0] = x;
        position[1] = y;
        position[2] = z;

        // update microenvironment current voxel index
        updateVoxelIndex();
        // update current_mechanics_voxel_index
        current_mechanics_voxel_index = get_container().underlying_mesh.nearest_voxel_index( position );

        // Since it is most likely our first position, we update the max_cell_interactive_distance_in_voxel
        // which was not initialized at cell creation
        if( get_container().max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < phenotype.geometry.radius
                * phenotype.mechanics.relative_maximum_adhesion_distance )
        {
            // get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()]= phenotype.geometry.radius*parameters.max_interaction_distance_factor;
            get_container().max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
                    * phenotype.mechanics.relative_maximum_adhesion_distance;
        }

        get_container().register_agent( this );

        if( !get_container().underlying_mesh.isPositionValid( x, y, z ) )
        {
            is_out_of_domain = true;
            isActive = false;
            is_movable = false;
            return false;
        }
        return true;
    }

    void start_death(int death_model_index)
    {
        // set the death data struture to the indicated death model 
        phenotype.death.trigger_death( death_model_index );
        // change the cycle model to the current death model 
        phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() );

        // turn off secretion, and reduce uptake by a factor of 10 
        phenotype.secretion.set_all_secretion_to_zero();
        phenotype.secretion.scale_all_uptake_by_factor( 0.10 );

        // turn off motility.
        phenotype.motility.is_motile = false;
        phenotype.motility.motility_vector = new double[3];//.assign( 3, 0.0 ); 
        functions.update_migration_bias = null;

        // make sure to run the death entry function 
        if( phenotype.cycle.current_phase().entryFunction != null )
        {
            phenotype.cycle.current_phase().entryFunction.execute( this, phenotype, 0.0 );
        }
    }

    void add_potentials(Cell other_agent)
    {
        // if( this->ID == other_agent->ID )
        if( this == other_agent )
            return;

        /*
        // new April 2022: don't interact with cells with 0 volume 
        // does not seem to really help 
        if( other_agent->phenotype.volume.total < 1e-15 )
        { std::cout << "zero size cell in mechanics!" << std::endl; return; }
        */
        // 12 uniform neighbors at a close packing distance, after dividing out all constants
        double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
        // 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

        double distance = 0;
        for( int i = 0; i < 3; i++ )
        {
            displacement[i] = position[i] - ( other_agent ).position[i];
            distance += displacement[i] * displacement[i];
        }
        // Make sure that the distance is not zero

        distance = Math.max( Math.sqrt( distance ), 0.00001 );

        //Repulsive
        double R = phenotype.geometry.radius + ( other_agent ).phenotype.geometry.radius;

        double RN = phenotype.geometry.nuclear_radius + ( other_agent ).phenotype.geometry.nuclear_radius;
        double temp_r, c;
        if( distance > R )
        {
            temp_r = 0;
        }
        // else if( distance < RN ) 
        // {
        // double M = 1.0; 
        // c = 1.0 - RN/R; 
        // c *= c; 
        // c -= M; 
        // temp_r = ( c*distance/RN  + M  ); 
        // }
        else
        {
            temp_r = 1 - distance / R;
            temp_r *= temp_r; // (1-d/R)^2 
            // add the relative pressure contribution 
            state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
        }

        // August 2017 - back to the original if both have same coefficient 

        double effective_repulsion = Math
                .sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent.phenotype.mechanics.cell_cell_repulsion_strength );
        temp_r *= effective_repulsion;

        // temp_r *= phenotype.mechanics.cell_cell_repulsion_strength; // original 
        //////////////////////////////////////////////////////////////////

        // Adhesive
        //double max_interactive_distance = parameters.max_interaction_distance_factor * phenotype.geometry.radius + 
        //  (*other_agent).parameters.max_interaction_distance_factor * (*other_agent).phenotype.geometry.radius;

        double max_interactive_distance = phenotype.mechanics.relative_maximum_adhesion_distance * phenotype.geometry.radius
                + ( other_agent ).phenotype.mechanics.relative_maximum_adhesion_distance * ( other_agent ).phenotype.geometry.radius;

        if( distance < max_interactive_distance )
        {
            // double temp_a = 1 - distance/max_interactive_distance; 
            double temp_a = -distance; // -d
            temp_a /= max_interactive_distance; // -d/S
            temp_a += 1.0; // 1 - d/S 
            temp_a *= temp_a; // (1-d/S)^2 
            // temp_a *= phenotype.mechanics.cell_cell_adhesion_strength; // original 

            // August 2017 - back to the original if both have same coefficient 
            // May 2022 - back to oriinal if both affinities are 1
            int ii = find_cell_definition_index( this.type );
            int jj = find_cell_definition_index( other_agent.type );

            double adhesion_ii = phenotype.mechanics.cell_cell_adhesion_strength * phenotype.mechanics.cell_adhesion_affinities[jj];
            double adhesion_jj = other_agent.phenotype.mechanics.cell_cell_adhesion_strength
                    * other_agent.phenotype.mechanics.cell_adhesion_affinities[ii];

            // double effective_adhesion = sqrt( phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
            double effective_adhesion = Math.sqrt( adhesion_ii * adhesion_jj );
            temp_a *= effective_adhesion;

            temp_r -= temp_a;

            state.neighbors.add( other_agent );
            //            state.neighbors.push_back(other_agent); // move here in 1.10.2 so non-adhesive cells also added. 
        }
        /////////////////////////////////////////////////////////////////
        if( Math.abs( temp_r ) < 1e-16 )
            return;

        temp_r /= distance;
        // for( int i = 0 ; i < 3 ; i++ ) 
        // {
        //  velocity[i] += displacement[i] * temp_r; 
        // }
        VectorUtil.axpy( velocity, temp_r, displacement );
        // state.neighbors.push_back(other_agent); // new 1.8.0
    }

    int find_cell_definition_index(String search_string)
    {
        Integer result = cell_definition_indices_by_name.get( search_string );
        return result != null ? result : -1;
        //        auto search = cell_definition_indices_by_name.find( search_string );
        //        // safety first! 
        //        if( search != cell_definition_indices_by_name.end() )
        //        {
        //            // if the target is found, set the appropriate rate 
        //            return search -> second;
        //        }
        //        return -1;
    }

    public static int find_cell_definition_index(int search_type)
    {
        Integer result = cell_definition_indices_by_type.get( search_type );
        return result != null ? result : -1;
        //        auto search = cell_definition_indices_by_type.find( search_type );
        //        // safety first! 
        //        if( search != cell_definition_indices_by_type.end() )
        //        {
        //            // if the target is found, set the appropriate rate 
        //            return search -> second;
        //        }
        //        return -1;
    }

    public void update_motility_vector(double dt_)
    {
        if( phenotype.motility.is_motile == false )
        {
            phenotype.motility.motility_vector = new double[3];//.assign( 3, 0.0 ); 
            return;
        }

        if( PhysiCellUtilities.UniformRandom() < dt_ / phenotype.motility.persistence_time || phenotype.motility.persistence_time < dt_ )
        {
            /*
            // choose a uniformly random unit vector 
            double temp_angle = 6.28318530717959*UniformRandom();
            double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
            
            double sin_phi = sin(temp_phi);
            double cos_phi = cos(temp_phi);
            
            if( phenotype.motility.restrict_to_2D == true )
            { 
                sin_phi = 1.0; 
                cos_phi = 0.0;
            }
            
            std::vector<double> randvec; 
            randvec.resize(3,sin_phi); 
            
            randvec[0] *= cos( temp_angle ); // cos(theta)*sin(phi)
            randvec[1] *= sin( temp_angle ); // sin(theta)*sin(phi)
            randvec[2] = cos_phi; //  cos(phi)
            */
            //            std::vector<double> randvec(3,0.0);
            double[] randvec;
            if( phenotype.motility.restrict_to_2D == true )
            {
                randvec = PhysiCellUtilities.UniformOnUnitCircle();
            }
            else
            {
                randvec = PhysiCellUtilities.UniformOnUnitSphere();
            }

            // if the update_bias_vector function is set, use it  
            if( functions.update_migration_bias != null )
            {
                functions.update_migration_bias.execute( this, phenotype, dt_ );
            }

            //            phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction; // motiltiy = bias_vector
            //            phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector 
            phenotype.motility.motility_vector = VectorUtil.newProd( phenotype.motility.migration_bias_direction,
                    phenotype.motility.migration_bias );

            double one_minus_bias = 1.0 - phenotype.motility.migration_bias;

            VectorUtil.axpy( ( phenotype.motility.motility_vector ), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector

            VectorUtil.normalize( ( phenotype.motility.motility_vector ) );

            VectorUtil.prod( phenotype.motility.motility_vector, phenotype.motility.migration_speed );
            //            phenotype.motility.motility_vector *= phenotype.motility.migration_speed; 
        }
    }

    public static boolean is_neighbor_voxel(Cell pCell, double[] my_voxel_center, double[] other_voxel_center, int other_voxel_index)
    {
        double max_interactive_distance = pCell.phenotype.mechanics.relative_maximum_adhesion_distance * pCell.phenotype.geometry.radius
                + pCell.get_container().max_cell_interactive_distance_in_voxel[other_voxel_index];

        int comparing_dimension = -1, comparing_dimension2 = -1;
        if( my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[1] == other_voxel_center[1] )
        {
            comparing_dimension = 2;
        }
        else if( my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[2] == other_voxel_center[2] )
        {
            comparing_dimension = 1;
        }
        else if( my_voxel_center[1] == other_voxel_center[1] && my_voxel_center[2] == other_voxel_center[2] )
        {
            comparing_dimension = 0;
        }

        if( comparing_dimension != -1 )
        { //then it is an immediate neighbor (through side faces)
            double surface_coord = 0.5 * ( my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension] );
            if( Math.abs( pCell.position[comparing_dimension] - surface_coord ) > max_interactive_distance )
            {
                return false;
            }
            return true;
        }
        comparing_dimension = -1;

        if( my_voxel_center[0] == other_voxel_center[0] )
        {
            comparing_dimension = 1;
            comparing_dimension2 = 2;
        }
        else if( my_voxel_center[1] == other_voxel_center[1] )
        {
            comparing_dimension = 0;
            comparing_dimension2 = 2;
        }
        else if( my_voxel_center[2] == other_voxel_center[2] )
        {
            comparing_dimension = 0;
            comparing_dimension2 = 1;
        }
        if( comparing_dimension != -1 )
        {
            double line_coord1 = 0.5 * ( my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension] );
            double line_coord2 = 0.5 * ( my_voxel_center[comparing_dimension2] + other_voxel_center[comparing_dimension2] );
            double distance_squared = Math.pow( pCell.position[comparing_dimension] - line_coord1, 2 )
                    + Math.pow( pCell.position[comparing_dimension2] - line_coord2, 2 );
            if( distance_squared > max_interactive_distance * max_interactive_distance )
                return false;
            return true;
        }
        //        double[] corner_point= 0.5*(my_voxel_center+other_voxel_center);
        double[] corner_point = VectorUtil.newSum( my_voxel_center, other_voxel_center );
        VectorUtil.prod( corner_point, 0.5 );
        double distance_squared = ( corner_point[0] - pCell.position[0] ) * ( corner_point[0] - pCell.position[0] )
                + ( corner_point[1] - pCell.position[1] ) * ( corner_point[1] - pCell.position[1] )
                + ( corner_point[2] - pCell.position[2] ) * ( corner_point[2] - pCell.position[2] );
        if( distance_squared > max_interactive_distance * max_interactive_distance )
            return false;
        return true;
    }

    double[] nearest_density_vector()
    {
        return getMicroenvironment().nearest_density_vector( this.currentVoxelIndex );//current_voxel_index );
    }

    void ingest_cell(Cell pCell_to_eat)
    {
        // don't ingest self 
        if( pCell_to_eat == this )
            return;

        // don't ingest a cell that's already ingested 
        if( pCell_to_eat.phenotype.volume.total < 1e-15 )
            return;

        // make this thread safe 
        //        #pragma omp critical
        {
            /*
            if( pCell_to_eat.phenotype.death.dead == true )
            { std::cout << this.type_name << " (" << this << ")" << " eats dead " << pCell_to_eat.type_name << " (" << pCell_to_eat 
                << ") of size " << pCell_to_eat.phenotype.volume.total << std::endl; }
            else
            { std::cout << this.type_name << " (" << this << ")" << " eats live " << pCell_to_eat.type_name << " (" << pCell_to_eat 
                << ") of size " << pCell_to_eat.phenotype.volume.total << std::endl; }
            */

            // mark it as dead 
            pCell_to_eat.phenotype.death.dead = true;
            // set secretion and uptake to zero 
            pCell_to_eat.phenotype.secretion.set_all_secretion_to_zero();
            pCell_to_eat.phenotype.secretion.set_all_uptake_to_zero();

            // deactivate all custom function 
            pCell_to_eat.functions.custom_cell_rule = null;
            pCell_to_eat.functions.updatePhenotype = null;
            pCell_to_eat.functions.contact_function = null;

            // should set volume fuction to NULL too! 
            pCell_to_eat.functions.volume_update_function = null;

            // set cell as unmovable and non-secreting 
            pCell_to_eat.is_movable = false;
            pCell_to_eat.isActive = false;

            // absorb all the volume(s)

            // absorb fluid volume (all into the cytoplasm) 
            phenotype.volume.cytoplasmic_fluid += pCell_to_eat.phenotype.volume.fluid;
            pCell_to_eat.phenotype.volume.cytoplasmic_fluid = 0.0;

            // absorb nuclear and cyto solid volume (into the cytoplasm) 
            phenotype.volume.cytoplasmic_solid += pCell_to_eat.phenotype.volume.cytoplasmic_solid;
            pCell_to_eat.phenotype.volume.cytoplasmic_solid = 0.0;

            phenotype.volume.cytoplasmic_solid += pCell_to_eat.phenotype.volume.nuclear_solid;
            pCell_to_eat.phenotype.volume.nuclear_solid = 0.0;

            // consistency calculations 

            phenotype.volume.fluid = phenotype.volume.nuclear_fluid + phenotype.volume.cytoplasmic_fluid;
            pCell_to_eat.phenotype.volume.fluid = 0.0;

            phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + phenotype.volume.nuclear_solid;
            pCell_to_eat.phenotype.volume.solid = 0.0;

            // no change to nuclear volume (initially) 
            pCell_to_eat.phenotype.volume.nuclear = 0.0;
            pCell_to_eat.phenotype.volume.nuclear_fluid = 0.0;

            phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_solid + phenotype.volume.cytoplasmic_fluid;
            pCell_to_eat.phenotype.volume.cytoplasmic = 0.0;

            phenotype.volume.total = phenotype.volume.nuclear + phenotype.volume.cytoplasmic;
            pCell_to_eat.phenotype.volume.total = 0.0;

            phenotype.volume.fluid_fraction = phenotype.volume.fluid / ( phenotype.volume.total + 1e-16 );
            pCell_to_eat.phenotype.volume.fluid_fraction = 0.0;

            phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / ( phenotype.volume.nuclear_solid + 1e-16 );

            // update corresponding BioFVM parameters (self-consistency) 
            set_total_volume( phenotype.volume.total );
            pCell_to_eat.set_total_volume( 0.0 );

            // absorb the internalized substrates 

            // multiply by the fraction that is supposed to be ingested (for each substrate) 
            //            *(pCell_to_eat.internalized_substrates) *= 
            //                *(pCell_to_eat.fraction_transferred_when_ingested); // 
            VectorUtil.prod( pCell_to_eat.internalizedSubstrates, pCell_to_eat.fraction_transferred_when_ingested );

            //            *internalized_substrates += *(pCell_to_eat.internalized_substrates);
            VectorUtil.sum( internalizedSubstrates, pCell_to_eat.internalizedSubstrates );
            //            int n_substrates = internalizedSubstrates.length;
            pCell_to_eat.internalizedSubstrates = new double[internalizedSubstrates.length];//assign( n_substrates , 0.0 );    

            // trigger removal from the simulation 
            // pCell_to_eat.die(); // I don't think this is safe if it's in an OpenMP loop 

            // flag it for removal 
            // pCell_to_eat.flag_for_removal(); 

            // remove all adhesions 
            // pCell_to_eat.remove_all_attached_cells();

        }

        // things that have their own thread safety 
        pCell_to_eat.flag_for_removal();
        pCell_to_eat.remove_all_attached_cells();
        pCell_to_eat.remove_all_spring_attachments();

        return;
    }

    void attack_cell(Cell pCell_to_attack, double dt)
    {
        // don't attack self 
        if( pCell_to_attack == this )
        {
            return;
        }

        // don't attack a dead or tiny cell 
        if( pCell_to_attack.phenotype.death.dead == true || pCell_to_attack.phenotype.volume.total < 1e-15 )
        {
            return;
        }

        // make this thread safe 
        //        #pragma omp critical
        {
            // std::cout << this.type_name << " attacks " << pCell_to_attack.type_name << std::endl;
            // 
            pCell_to_attack.state.damage += phenotype.cell_interactions.damage_rate * dt;
            pCell_to_attack.state.total_attack_time += dt;
        }
        return;
    }



    void fuse_cell(Cell pCell_to_fuse)
    {
        // don't ingest a cell that's already fused or fuse self 
        if( pCell_to_fuse.phenotype.volume.total < 1e-15 || this == pCell_to_fuse )
        {
            return;
        }

        // make this thread safe 
        //        #pragma omp critical
        {

            // set new position at center of volume 
            // x_new = (vol_B * x_B + vol_S * x_S ) / (vol_B + vol_S )

            //            std::vector<double> new_position = position; // x_B
            //            new_position *= phenotype.volume.total; // vol_B * x_B 
            double[] new_position = VectorUtil.newProd( position, phenotype.volume.total );
            double total_volume = phenotype.volume.total;
            total_volume += pCell_to_fuse.phenotype.volume.total;

            VectorUtil.axpy( new_position, pCell_to_fuse.phenotype.volume.total, pCell_to_fuse.position ); // vol_B*x_B + vol_S*x_S
            //            new_position /= total_volume; // (vol_B*x_B+vol_S*x_S)/(vol_B+vol_S);
            VectorUtil.div( new_position, total_volume );

            double xL = Microenvironment.get_default_microenvironment().mesh.bounding_box[0];
            double xU = Microenvironment.get_default_microenvironment().mesh.bounding_box[3];

            double yL = Microenvironment.get_default_microenvironment().mesh.bounding_box[1];
            double yU = Microenvironment.get_default_microenvironment().mesh.bounding_box[4];

            double zL = Microenvironment.get_default_microenvironment().mesh.bounding_box[2];
            double zU = Microenvironment.get_default_microenvironment().mesh.bounding_box[5];

            if( new_position[0] < xL || new_position[0] > xU || new_position[1] < yL || new_position[1] > yU || new_position[2] < zL
                    || new_position[2] > zU )
            {
                System.out.println( "cell fusion at " + new_position + " violates domain bounds" );
                System.out.println( Microenvironment.get_default_microenvironment().mesh.bounding_box );
                //                std::cout << "cell fusion at " << new_position << " violates domain bounds" << std::endl; 
                //                std::cout << get_default_microenvironment().mesh.bounding_box << std::endl << std::endl; 
            }
            position = new_position;
            update_voxel_in_container();

            // set number of nuclei 

            state.number_of_nuclei += pCell_to_fuse.state.number_of_nuclei;

            // absorb all the volume(s)

            // absorb fluid volume (all into the cytoplasm) 
            phenotype.volume.cytoplasmic_fluid += pCell_to_fuse.phenotype.volume.cytoplasmic_fluid;
            pCell_to_fuse.phenotype.volume.cytoplasmic_fluid = 0.0;

            phenotype.volume.nuclear_fluid += pCell_to_fuse.phenotype.volume.nuclear_fluid;
            pCell_to_fuse.phenotype.volume.nuclear_fluid = 0.0;

            // absorb nuclear and cyto solid volume (into the cytoplasm) 
            phenotype.volume.cytoplasmic_solid += pCell_to_fuse.phenotype.volume.cytoplasmic_solid;
            pCell_to_fuse.phenotype.volume.cytoplasmic_solid = 0.0;

            phenotype.volume.nuclear_solid += pCell_to_fuse.phenotype.volume.nuclear_solid;
            pCell_to_fuse.phenotype.volume.nuclear_solid = 0.0;

            // consistency calculations 

            phenotype.volume.fluid = phenotype.volume.nuclear_fluid + phenotype.volume.cytoplasmic_fluid;
            pCell_to_fuse.phenotype.volume.fluid = 0.0;

            phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + phenotype.volume.nuclear_solid;
            pCell_to_fuse.phenotype.volume.solid = 0.0;

            phenotype.volume.nuclear = phenotype.volume.nuclear_fluid + phenotype.volume.nuclear_solid;
            pCell_to_fuse.phenotype.volume.nuclear = 0.0;

            phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_fluid + phenotype.volume.cytoplasmic_solid;
            pCell_to_fuse.phenotype.volume.cytoplasmic = 0.0;

            phenotype.volume.total = phenotype.volume.nuclear + phenotype.volume.cytoplasmic;
            pCell_to_fuse.phenotype.volume.total = 0.0;

            phenotype.volume.fluid_fraction = phenotype.volume.fluid / ( phenotype.volume.total + 1e-16 );
            pCell_to_fuse.phenotype.volume.fluid_fraction = 0.0;

            phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / ( phenotype.volume.nuclear_solid + 1e-16 );

            // update corresponding BioFVM parameters (self-consistency) 
            set_total_volume( phenotype.volume.total );
            pCell_to_fuse.set_total_volume( 0.0 );

            // absorb the internalized substrates 

            //            *internalized_substrates += *(pCell_to_fuse.internalized_substrates); 
            //            static int n_substrates = internalized_substrates.size(); 
            //            pCell_to_fuse.internalized_substrates.assign( n_substrates , 0.0 );   
            VectorUtil.sum( internalizedSubstrates, pCell_to_fuse.internalizedSubstrates );
            pCell_to_fuse.internalizedSubstrates = new double[internalizedSubstrates.length];
            // set target volume(s)

            phenotype.volume.target_solid_cytoplasmic += pCell_to_fuse.phenotype.volume.target_solid_cytoplasmic;
            phenotype.volume.target_solid_nuclear += pCell_to_fuse.phenotype.volume.target_solid_nuclear;

            // trigger removal from the simulation 
            // pCell_to_eat.die(); // I don't think this is safe if it's in an OpenMP loop 

            // flag it for removal 
            // pCell_to_eat.flag_for_removal(); 
            // mark it as dead 
            pCell_to_fuse.phenotype.death.dead = true;
            // set secretion and uptake to zero 
            pCell_to_fuse.phenotype.secretion.set_all_secretion_to_zero();
            pCell_to_fuse.phenotype.secretion.set_all_uptake_to_zero();

            // deactivate all custom function 
            pCell_to_fuse.functions.custom_cell_rule = null;
            pCell_to_fuse.functions.updatePhenotype = null;
            pCell_to_fuse.functions.contact_function = null;
            pCell_to_fuse.functions.volume_update_function = null;

            // remove all adhesions 
            // pCell_to_eat.remove_all_attached_cells();

            // set cell as unmovable and non-secreting 
            pCell_to_fuse.is_movable = false;
            pCell_to_fuse.isActive = false;

        }

        // things that have their own thread safety 
        pCell_to_fuse.flag_for_removal();
        pCell_to_fuse.remove_all_attached_cells();
        pCell_to_fuse.remove_all_spring_attachments();
    }

    public void convert_to_cell_definition(CellDefinition cd)
    {
        // use the cell defaults; 
        type = cd.type;
        type_name = cd.name;

        custom_data = cd.custom_data;
        parameters = cd.parameters;
        functions = cd.functions;

        phenotype = cd.phenotype;
        // is_movable = true;
        // is_out_of_domain = false;

        // displacement.resize(3,0.0); // state? 

        assign_orientation();

        set_total_volume( phenotype.volume.total );
    }

    public void advance_bundled_phenotype_functions(double dt_)
    {
        // New March 2022
        // perform transformations 
        //        System.out.println( "Step: " + dt_ );
        StandardModels.standard_cell_transformations( this, this.phenotype, dt_ );

        // New March 2023 in Version 1.12.0 
        // call the rules-based code to update the phenotype 
        //        if( PhysiCellSettings.rules_enabled )
        //        { 
        //            apply_ruleset( this ); TODO: later
        //        }
        //        if( SignalBehavior.get_single_signal( this, "necrotic" ) > 0.5 )
        //        {
        //            double rupture = this.phenotype.volume.rupture_volume;
        //            double volume = this.phenotype.volume.total;
        //            if( volume > rupture )
        //            {
        //                System.out.println( volume + " vs " + this.phenotype.volume.rupture_volume + " dead: "
        //                        + SignalBehavior.get_single_signal( this, "dead" ) );
        //                System.out.println( this.phenotype.cycle.current_phase_index() + " " + this.phenotype.cycle.pCycle_Model.name );
        //                //                std::cout << this->phenotype.volume.total << " vs " << this->phenotype.volume.rupture_volume << 
        //                //                " dead: " << get_single_signal( this, "dead") <<    std::endl; 
        //                //                std::cout << this->phenotype.cycle.current_phase_index() << " " 
        //                //                << this->phenotype.cycle.pCycle_Model->name << std::endl; 
        //            }
        //
        //        }

        //  if( functions.update_phenotype )
        //  { functions.update_phenotype( this , phenotype , dt_ ); }

        // call the custom code to update the phenotype 
        if( functions.updatePhenotype != null )
        {
            functions.updatePhenotype.execute( this, phenotype, dt_ );
        }

        // update volume 
        if( functions.volume_update_function != null )
        {
            functions.volume_update_function.execute( this, phenotype, dt_ );

            // The following line is needed in every volume 
            // regulation method (it sets BioFVM total_volume)

            set_total_volume( phenotype.volume.total );
        }

        // update geometry
        phenotype.geometry.update( this, phenotype, dt_ );

        // check for new death events 
        if( phenotype.death.check_for_death( dt_ ) == true )
        {
            // if so, change the cycle model to the current death model 
            phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() );

            // also, turn off motility.

            phenotype.motility.is_motile = false;
            phenotype.motility.motility_vector = new double[3];//.assign( 3, 0.0 ); 
            functions.update_migration_bias = null;

            // turn off secretion, and reduce uptake by a factor of 10 
            phenotype.secretion.set_all_secretion_to_zero();
            phenotype.secretion.scale_all_uptake_by_factor( 0.10 );

            // make sure to run the death entry function 
            if( phenotype.cycle.current_phase().entryFunction != null )
            {
                phenotype.cycle.current_phase().entryFunction.execute( this, phenotype, dt_ );
            }
        }

        // advance cycle model (for both cell cycle and death cycle models)
        phenotype.cycle.advance_cycle( this, phenotype, dt_ );
        if( phenotype.flagged_for_removal )
        {
            flag_for_removal();
            phenotype.flagged_for_removal = false;
        }
        if( phenotype.flagged_for_division )
        {
            flag_for_division();
            phenotype.flagged_for_division = false;
        }
    }

    @Override
    public double getRadius()
    {
        return this.phenotype.geometry.radius;
    }
}