package ru.biosoft.physicell.core;

import java.util.Set;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.CellFunctions.instantiate_cell;

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
    public static double[] cell_division_orientation;
    CellContainer container;
    int currentMechanicsVoxelIndex;
    int updated_current_mechanics_voxel_index; // keeps the updated voxel index for later adjusting of current voxel index

    public String typeName;

    CellDefinition definition;
    public CustomCellData custom_data;// = new CustomCellData();
    CellParameters parameters;// = new CellParameters();
    public CellFunctions functions;// = new CellFunctions();

    public CellState state = new CellState();
    public Phenotype phenotype;// = new Phenotype();
    public boolean isOutOfDomain;
    boolean isMovable;
    public double[] displacement; // this should be moved to state, or made private

    public Cell(CellDefinition cd, Microenvironment m)
    {
        super( m );
        type = cd.type;
        typeName = cd.name;
        custom_data = cd.custom_data.clone();
        parameters = cd.parameters.clone();
        functions = cd.functions.clone();
        phenotype = cd.phenotype.clone();
        phenotype.molecular.sync( this );
        this.definition = cd;
        // cell state should be fine by the default constructor 
        currentMechanicsVoxelIndex = -1;
        updated_current_mechanics_voxel_index = 0;

        isMovable = true;
        isOutOfDomain = false;
        displacement = new double[3];// state? 

        assignOrientation();
        container = null;
        setTotalVolume( phenotype.volume.total );
    }

    @Override
    public CellContainer get_container()
    {
        if( container == null )
            container = (CellContainer)getMicroenvironment().agentContainer;
        return container;
    }

    void flagForDivision()
    {
        get_container().flag_cell_for_division( this );
    }

    void flagForRemoval()
    {
        get_container().flag_cell_for_removal( this );
    }

    @Override
    public void setTotalVolume(double volume)
    {
        super.setTotalVolume( volume );

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
                    * phenotype.mechanics.relMaxAdhesionDistance )
            {
                // get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()]= phenotype.geometry.radius*parameters.max_interaction_distance_factor;
                get_container().max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
                        * phenotype.mechanics.relMaxAdhesionDistance;
            }
        }
    }

    void assignOrientation()
    {
        state.orientation = new double[3];
        //TODO: check if function should be used at all
        //        if( functions.set_orientation != null )
        //        {
        //            functions.set_orientation( this, phenotype, 0.0 );
        //        }
        //        else
        //        {
        //assign a random unit vector
        double theta = PhysiCellUtilities.UniformRandom() * 6.28318530717959; //rand*2*pi
        double z = 2 * PhysiCellUtilities.UniformRandom() - 1;
        double temp = Math.sqrt( 1 - z * z );
        state.orientation[0] = temp * Math.cos( theta );
        state.orientation[1] = temp * Math.sin( theta );
        state.orientation[2] = z;
        //        }
    }

    @Override
    public int get_current_mechanics_voxel_index()
    {
        return currentMechanicsVoxelIndex;
    }

    public static Cell createCell(CellDefinition cd, Microenvironment m, double[] position)
    {
        return createCell( null, cd, m, position );
    }

    public static Cell createCell(instantiate_cell custom_instantiate, CellDefinition cd, Microenvironment m, double[] position)
    {
        Cell pNew = custom_instantiate == null ? new Cell( cd, m ) : custom_instantiate.execute();
        pNew.index = m.getAgentsCount();
        pNew.registerMicroenvironment( m );
        // All the phenotype and other data structures are already set by virtue of the default Cell constructor. 
        pNew.setTotalVolume( pNew.phenotype.volume.total );
        pNew.assignPosition( position );
        return pNew;
    }

    public void die()
    {
        deleteCell( this );
    }

    @Override
    public void updatePosition(double dt)
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
        if( microenvironment.options.simulate2D )
        {
            velocity[2] = 0.0;
        }

        //        std::vector<double> old_position(position);

        VectorUtil.axpy( position, d1, velocity );
        VectorUtil.axpy( position, d2, previous_velocity );
        // overwrite previous_velocity for future use 
        // if(sqrt(dist(old_position, position))>3* phenotype.geometry.radius)
        // std::cout<<sqrt(dist(old_position, position))<<"old_position: "<<old_position<<", new position: "<< position<<", velocity: "<<velocity<<", previous_velocity: "<< previous_velocity<<std::endl;

        previous_velocity = velocity.clone();

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

            isOutOfDomain = true;
            isActive = false;
            isMovable = false;
        }
    }

    public static void deleteCell(Cell cell)
    {
        //  std::cout << __FUNCTION__ << " " << (*all_cells)[index] 
        //  << " " << (*all_cells)[index]->type_name << std::endl; 

        //        Cell pDeleteMe = (Cell)BasicAgent.allBasicAgents.get( index );
        //        System.out.println( "Died " );
        // release any attached cells (as of 1.7.2 release)
        cell.removeAllAttachedCells();
        // 1.11.0 
        cell.removeAllSpringAttachments();

        // released internalized substrates (as of 1.5.x releases)
        cell.release_internalized_substrates();

        // performance goal: don't delete in the middle -- very expensive reallocation
        // alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

        // move last item to index location 
        //        BasicAgent.allBasicAgents.remove( index );

        //        BasicAgent last = BasicAgent.allBasicAgents.get( BasicAgent.allBasicAgents.size() - 1 );
        //        last.index = index;
        //        BasicAgent.allBasicAgents.set( index, last );
        //        BasicAgent.allBasicAgents.remove( BasicAgent.allBasicAgents.size() - 1 );
        //        (*all_cells)[ (*all_cells).size()-1 ]->index=index;
        //        (*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];
        //        // shrink the vector
        //        (*all_cells).pop_back();    

        // deregister agent in from the agent container
        cell.getMicroenvironment().removeAgent( cell );
        cell.get_container().remove_agent( cell );
        // de-allocate (delete) the cell; 
        //        delete pDeleteMe; 
    }

    public Cell divide()
    {
        //        System.out.println( toString() + " divided" );
        //commented in original code
        // phenotype.flagged_for_division = false; 
        // phenotype.flagged_for_removal = false; 

        // make sure ot remove adhesions 
        removeAllAttachedCells();
        removeAllSpringAttachments();

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
        for( int nn = 0; nn < custom_data.vectorVariables.size(); nn++ )
        {
            if( custom_data.vectorVariables.get( nn ).conserved_quantity )
            {
                VectorUtil.prod( custom_data.vectorVariables.get( nn ).value, 0.5 );
            }
        }

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
        double[] pos = VectorUtil.newSum( position, rand_vec );
        //        System.out.println( type_name + "\tdivided\t" + position[0] + "\t" + position[1] + "\t" + position[2] + "\t->\t" + pos[0] + "\t"
        //                + pos[1] + "\t" + pos[2] );
        Cell child = createCell( functions.instantiate_cell, definition, microenvironment, pos );
        child.copyData( this );
        child.copyFunctionPointers( this );
        child.parameters = parameters;

        // evenly divide internalized substrates 
        // if these are not actively tracked, they are zero anyway 
        VectorUtil.prod( internalizedSubstrates, 0.5 );
        child.internalizedSubstrates = internalizedSubstrates.clone();

        //change my position to keep the center of mass intact 
        // and then see if I need to update my voxel index
        double negative_one_half = -0.5;
        VectorUtil.axpy( position, negative_one_half, rand_vec ); // position = position - 0.5*rand_vec; 

        //If this cell has been moved outside of the boundaries, mark it as such.
        //(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)
        if( !get_container().underlying_mesh.isPositionValid( position[0], position[1], position[2] ) )
        {
            isOutOfDomain = true;
            isActive = false;
            isMovable = false;
        }

        updateVoxelInContainer();
        phenotype.volume.divide();
        child.phenotype.volume.divide();
        child.setTotalVolume( child.phenotype.volume.total );//TODO: check
        setTotalVolume( phenotype.volume.total );
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
        state.totalAttackTime = 0;
        child.state.damage = 0.0;
        child.state.totalAttackTime = 0.0;
        return child;
    }

    public void removeAllAttachedCells()
    {
        for( Cell cell : state.attachedCells )
            cell.detachCell( this );
        state.attachedCells.clear();
    }

    public void removeAllSpringAttachments()
    {
        for( Cell cell : state.springAttachments )
            cell.detachCellAsSpring( this );
        state.springAttachments.clear();
    }

    public static void attachcCells(Cell pCell_1, Cell pCell_2)
    {
        pCell_1.attachCell( pCell_2 );
        pCell_2.attachCell( pCell_1 );
    }

    public static void attachCellsAsSpring(Cell pCell_1, Cell pCell_2)
    {
        pCell_1.attachCellAsSpring( pCell_2 );
        pCell_2.attachCellAsSpring( pCell_1 );
    }

	public static void detachCells(Cell pCell_1, Cell pCell_2)
    {
        pCell_1.detachCell( pCell_2 );
        pCell_2.detachCell( pCell_1 );
    }

    public static void detachCellsAsSpring(Cell pCell_1, Cell pCell_2)
    {
        pCell_1.detachCellAsSpring( pCell_2 );
        pCell_2.detachCellAsSpring( pCell_1 );
    }

    void attachCell(Cell pAddMe)
    {
        state.attachedCells.add( pAddMe );
    }

    void attachCellAsSpring(Cell pAddMe)
    {
        state.springAttachments.add( pAddMe );
    }

    void detachCell(Cell pRemoveMe)
    {
        state.attachedCells.remove( pRemoveMe );
    }

    void detachCellAsSpring(Cell pRemoveMe)
    {

        state.springAttachments.remove( pRemoveMe );
    }

    public void copyData(Cell copy_me)
    {
        // phenotype=copyMe->phenotype; //it is taken care in set_phenotype
        type = copy_me.type;
        typeName = copy_me.typeName;

        custom_data = copy_me.custom_data.clone();
        parameters = copy_me.parameters.clone();

        velocity = copy_me.velocity.clone();
        // expected_phenotype = copy_me. expected_phenotype; //it is taken care in set_phenotype
        sourceSinkTemp1 = copy_me.sourceSinkTemp1.clone();
        sourceSinkTemp2 = copy_me.sourceSinkTemp2.clone();
    }

    void copyFunctionPointers(Cell copy_me)
    {
        functions = copy_me.functions.clone();
    }

    void updateVoxelInContainer()
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
            currentMechanicsVoxelIndex = -1;
            isOutOfDomain = true;
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
            currentMechanicsVoxelIndex = updated_current_mechanics_voxel_index;
        }
    }

    @Override
    public boolean assignPosition(double[] new_position)
    {
        return assignPosition( new_position[0], new_position[1], new_position[2] );
    }

    public void setPreviousVelocity(double xV, double yV, double zV)
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
        currentMechanicsVoxelIndex = get_container().underlying_mesh.nearest_voxel_index( position );

        // Since it is most likely our first position, we update the max_cell_interactive_distance_in_voxel
        // which was not initialized at cell creation
        if( get_container().max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < phenotype.geometry.radius
                * phenotype.mechanics.relMaxAdhesionDistance )
        {
            // get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()]= phenotype.geometry.radius*parameters.max_interaction_distance_factor;
            get_container().max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
                    * phenotype.mechanics.relMaxAdhesionDistance;
        }

        get_container().register_agent( this );

        if( !get_container().underlying_mesh.isPositionValid( x, y, z ) )
        {
            isOutOfDomain = true;
            isActive = false;
            isMovable = false;
            return false;
        }
        return true;
    }

    public void startDeath(int deathModelIndex)
    {
        // set the death data struture to the indicated death model 
        phenotype.death.triggerDeath( deathModelIndex );
        // change the cycle model to the current death model 
        phenotype.cycle = phenotype.death.currentModel();

        // turn off secretion, and reduce uptake by a factor of 10 
        phenotype.secretion.setSecretionToZero();
        phenotype.secretion.scaleUptake( 0.10 );

        // turn off motility.
        phenotype.motility.isMotile = false;
        phenotype.motility.motilityVector = new double[3];//.assign( 3, 0.0 ); 
        functions.update_migration_bias = null;

        // make sure to run the death entry function 
        if( phenotype.cycle.data.currentPhase().entryFunction != null )
        {
            phenotype.cycle.data.currentPhase().entryFunction.execute( this, phenotype, 0.0 );
        }
    }

    void addPotentials(Cell other)
    {
        if( this.ID == other.ID )
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
            displacement[i] = position[i] - other.position[i];
            distance += displacement[i] * displacement[i];
        }
        // Make sure that the distance is not zero

        distance = Math.max( Math.sqrt( distance ), 0.00001 );
        //        double distance = Math.max( VectorUtil.dist( position, other.position ), 0.00001 );
        //Repulsive
        double R = phenotype.geometry.radius + other.phenotype.geometry.radius;
        //        double RN = phenotype.geometry.nuclear_radius + other.phenotype.geometry.nuclear_radius;
        double temp_r;
        if( distance > R )
            temp_r = 0;
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
            state.simplePressure += ( temp_r / simple_pressure_scale ); // New July 2017 
        }

        Mechanics m1 = phenotype.mechanics;
        Mechanics m2 = other.phenotype.mechanics;
        // August 2017 - back to the original if both have same coefficient 
        double effective_repulsion = Math.sqrt( m1.cellCellRepulsionStrength * m2.cellCellRepulsionStrength );
        temp_r *= effective_repulsion;

        // Adhesive
        double max_interactive_distance = m1.relMaxAdhesionDistance * phenotype.geometry.radius
                + m2.relMaxAdhesionDistance * other.phenotype.geometry.radius;

        if( distance < max_interactive_distance )
        {
            // double temp_a = 1 - distance/max_interactive_distance; 
            //            double temp_a = -distance; // -d
            //            temp_a /= max_interactive_distance; // -d/S
            //            temp_a += 1.0; // 1 - d/S
            double temp_a = 1 - distance / max_interactive_distance;
            temp_a *= temp_a; // (1-d/S)^2 
            // temp_a *= phenotype.mechanics.cell_cell_adhesion_strength; // original 

            // August 2017 - back to the original if both have same coefficient 
            // May 2022 - back to original if both affinities are 1
            int ii = CellDefinition.getCellDefinitionIndex( type );//CellDefinition.getCellDefinition( type ).;
            int jj = CellDefinition.getCellDefinitionIndex( other.type );//find_cell_definition_index( other.type );

            double adhesion_ii = m1.cellCellAdhesionStrength * m1.cellAdhesionAffinities[jj];
            double adhesion_jj = m2.cellCellAdhesionStrength * m2.cellAdhesionAffinities[ii];

            // double effective_adhesion = sqrt( phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
            double effective_adhesion = Math.sqrt( adhesion_ii * adhesion_jj );
            temp_a *= effective_adhesion;
            temp_r -= temp_a;
            state.neighbors.add( other );// move here in 1.10.2 so non-adhesive cells also added. 
            //            state.neighbors.push_back(other_agent); 
        }
        if( Math.abs( temp_r ) < 1e-16 )
            return;
        temp_r /= distance;
        VectorUtil.axpy( velocity, temp_r, displacement );
    }

    public void updateMotilityVector(double dt_)
    {
        if( phenotype.motility.isMotile == false )
        {
            phenotype.motility.motilityVector = new double[3];//.assign( 3, 0.0 ); 
            return;
        }

        if( PhysiCellUtilities.UniformRandom() < dt_ / phenotype.motility.persistenceTime || phenotype.motility.persistenceTime < dt_ )
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
            if( phenotype.motility.restrictTo2D == true )
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
            phenotype.motility.motilityVector = VectorUtil.newProd( phenotype.motility.migrationBiasDirection,
                    phenotype.motility.migrationBias );

            double one_minus_bias = 1.0 - phenotype.motility.migrationBias;

            VectorUtil.axpy( ( phenotype.motility.motilityVector ), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector

            VectorUtil.normalize( ( phenotype.motility.motilityVector ) );

            VectorUtil.prod( phenotype.motility.motilityVector, phenotype.motility.migrationSpeed );
            //            phenotype.motility.motility_vector *= phenotype.motility.migration_speed; 
        }
    }

    public static boolean isNeighborVoxel(Cell pCell, double[] my_voxel_center, double[] other_voxel_center, int other_voxel_index)
    {
        double max_interactive_distance = pCell.phenotype.mechanics.relMaxAdhesionDistance * pCell.phenotype.geometry.radius
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

    public double[] nearest_density_vector()
    {
        if( currentVoxelIndex == -1 )
            return new double[] { -1};
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
            pCell_to_eat.phenotype.secretion.setSecretionToZero();
            pCell_to_eat.phenotype.secretion.setUptakeToZero();

            // deactivate all custom function 
            pCell_to_eat.functions.custom_cell_rule = null;
            pCell_to_eat.functions.updatePhenotype = null;
            pCell_to_eat.functions.contact_function = null;

            // should set volume fuction to NULL too! 
            pCell_to_eat.functions.updateVolume = null;

            // set cell as unmovable and non-secreting 
            pCell_to_eat.isMovable = false;
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
            setTotalVolume( phenotype.volume.total );
            pCell_to_eat.setTotalVolume( 0.0 );

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
        pCell_to_eat.flagForRemoval();
        pCell_to_eat.removeAllAttachedCells();
        pCell_to_eat.removeAllSpringAttachments();

        return;
    }

    void attackCell(Cell pCellToAttack, double dt)
    {
        if( pCellToAttack == this )
            return; // don't attack self 

        if( pCellToAttack.phenotype.death.dead == true || pCellToAttack.phenotype.volume.total < 1e-15 )
            return; // don't attack a dead or tiny cell 

        // make this thread safe 
        //        #pragma omp critical
        {
            // std::cout << this.type_name << " attacks " << pCell_to_attack.type_name << std::endl;
            pCellToAttack.state.damage += phenotype.cellInteractions.damageRate * dt;
            pCellToAttack.state.totalAttackTime += dt;
        }
    }

    void fuseCell(Cell pCellToFuse)
    {

        if( pCellToFuse.phenotype.volume.total < 1e-15 || this == pCellToFuse )
            return; // don't ingest a cell that's already fused or fuse self 

        // make this thread safe 
        //        #pragma omp critical
        {

            // set new position at center of volume 
            // x_new = (vol_B * x_B + vol_S * x_S ) / (vol_B + vol_S )

            //            std::vector<double> new_position = position; // x_B
            //            new_position *= phenotype.volume.total; // vol_B * x_B 
            double[] new_position = VectorUtil.newProd( position, phenotype.volume.total );
            double total_volume = phenotype.volume.total;
            total_volume += pCellToFuse.phenotype.volume.total;

            VectorUtil.axpy( new_position, pCellToFuse.phenotype.volume.total, pCellToFuse.position ); // vol_B*x_B + vol_S*x_S
            //            new_position /= total_volume; // (vol_B*x_B+vol_S*x_S)/(vol_B+vol_S);
            VectorUtil.div( new_position, total_volume );

            double xL = microenvironment.mesh.boundingBox[0];
            double xU = microenvironment.mesh.boundingBox[3];

            double yL = microenvironment.mesh.boundingBox[1];
            double yU = microenvironment.mesh.boundingBox[4];

            double zL = microenvironment.mesh.boundingBox[2];
            double zU = microenvironment.mesh.boundingBox[5];

            if( new_position[0] < xL || new_position[0] > xU || new_position[1] < yL || new_position[1] > yU || new_position[2] < zL
                    || new_position[2] > zU )
            {
                System.out.println( "cell fusion at " + new_position + " violates domain bounds" );
                System.out.println( microenvironment.mesh.boundingBox );
                //                std::cout << "cell fusion at " << new_position << " violates domain bounds" << std::endl; 
                //                std::cout << get_default_microenvironment().mesh.bounding_box << std::endl << std::endl; 
            }
            position = new_position;
            updateVoxelInContainer();

            // set number of nuclei 

            state.numberNuclei += pCellToFuse.state.numberNuclei;

            // absorb all the volume(s)

            // absorb fluid volume (all into the cytoplasm) 
            phenotype.volume.cytoplasmic_fluid += pCellToFuse.phenotype.volume.cytoplasmic_fluid;
            pCellToFuse.phenotype.volume.cytoplasmic_fluid = 0.0;

            phenotype.volume.nuclear_fluid += pCellToFuse.phenotype.volume.nuclear_fluid;
            pCellToFuse.phenotype.volume.nuclear_fluid = 0.0;

            // absorb nuclear and cyto solid volume (into the cytoplasm) 
            phenotype.volume.cytoplasmic_solid += pCellToFuse.phenotype.volume.cytoplasmic_solid;
            pCellToFuse.phenotype.volume.cytoplasmic_solid = 0.0;

            phenotype.volume.nuclear_solid += pCellToFuse.phenotype.volume.nuclear_solid;
            pCellToFuse.phenotype.volume.nuclear_solid = 0.0;

            // consistency calculations 

            phenotype.volume.fluid = phenotype.volume.nuclear_fluid + phenotype.volume.cytoplasmic_fluid;
            pCellToFuse.phenotype.volume.fluid = 0.0;

            phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + phenotype.volume.nuclear_solid;
            pCellToFuse.phenotype.volume.solid = 0.0;

            phenotype.volume.nuclear = phenotype.volume.nuclear_fluid + phenotype.volume.nuclear_solid;
            pCellToFuse.phenotype.volume.nuclear = 0.0;

            phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_fluid + phenotype.volume.cytoplasmic_solid;
            pCellToFuse.phenotype.volume.cytoplasmic = 0.0;

            phenotype.volume.total = phenotype.volume.nuclear + phenotype.volume.cytoplasmic;
            pCellToFuse.phenotype.volume.total = 0.0;

            phenotype.volume.fluid_fraction = phenotype.volume.fluid / ( phenotype.volume.total + 1e-16 );
            pCellToFuse.phenotype.volume.fluid_fraction = 0.0;

            phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / ( phenotype.volume.nuclear_solid + 1e-16 );

            // update corresponding BioFVM parameters (self-consistency) 
            setTotalVolume( phenotype.volume.total );
            pCellToFuse.setTotalVolume( 0.0 );

            // absorb the internalized substrates 

            //            *internalized_substrates += *(pCell_to_fuse.internalized_substrates); 
            //            static int n_substrates = internalized_substrates.size(); 
            //            pCell_to_fuse.internalized_substrates.assign( n_substrates , 0.0 );   
            VectorUtil.sum( internalizedSubstrates, pCellToFuse.internalizedSubstrates );
            pCellToFuse.internalizedSubstrates = new double[internalizedSubstrates.length];
            // set target volume(s)

            phenotype.volume.target_solid_cytoplasmic += pCellToFuse.phenotype.volume.target_solid_cytoplasmic;
            phenotype.volume.target_solid_nuclear += pCellToFuse.phenotype.volume.target_solid_nuclear;

            // trigger removal from the simulation 
            // pCell_to_eat.die(); // I don't think this is safe if it's in an OpenMP loop 

            // flag it for removal 
            // pCell_to_eat.flag_for_removal(); 
            // mark it as dead 
            pCellToFuse.phenotype.death.dead = true;
            // set secretion and uptake to zero 
            pCellToFuse.phenotype.secretion.setSecretionToZero();
            pCellToFuse.phenotype.secretion.setUptakeToZero();

            // deactivate all custom function 
            pCellToFuse.functions.custom_cell_rule = null;
            pCellToFuse.functions.updatePhenotype = null;
            pCellToFuse.functions.contact_function = null;
            pCellToFuse.functions.updateVolume = null;

            // remove all adhesions 
            // pCell_to_eat.remove_all_attached_cells();

            // set cell as unmovable and non-secreting 
            pCellToFuse.isMovable = false;
            pCellToFuse.isActive = false;

        }

        // things that have their own thread safety 
        pCellToFuse.flagForRemoval();
        pCellToFuse.removeAllAttachedCells();
        pCellToFuse.removeAllSpringAttachments();
    }

    public void convert_to_cell_definition(CellDefinition cd)
    {
        // use the cell defaults; 
        type = cd.type;
        typeName = cd.name;
        custom_data = cd.custom_data.clone();
        parameters = cd.parameters.clone();
        functions = cd.functions.clone();
        phenotype = cd.phenotype.clone();
        // is_movable = true;
        // is_out_of_domain = false;
        // displacement.resize(3,0.0); // state? 
        assignOrientation();
        setTotalVolume( phenotype.volume.total );
    }

    public void advanceBundledPhenotype(double dt) throws Exception
    {
        // New March 2022
        // perform transformations 
        //        System.out.println( "Step: " + dt_ );
        StandardModels.standard_cell_transformations( this, this.phenotype, dt );

        // New March 2023 in Version 1.12.0 
        // call the rules-based code to update the phenotype 
        if( PhysiCellSettings.rules_enabled )
        {
            Rules.apply_ruleset( this );
        }
        //        if( SignalBehavior.getSingleSignal( this, "necrotic" ) > 0.5 )
        //        {
        //            double rupture = this.phenotype.volume.rupture_volume;
        //            double volume = this.phenotype.volume.total;
        //            if( volume > rupture )
        //            {
        //                System.out.println( volume + " vs " + this.phenotype.volume.rupture_volume + " dead: "
        //                        + SignalBehavior.getSingleSignal( this, "dead" ) );
        //                System.out.println( this.phenotype.cycle.currentPhase().name + " " + this.phenotype.cycle.name );
        //            }
        //
        //        }

        // call the custom code to update the phenotype 
        if( functions.updatePhenotype != null )
        {
            functions.updatePhenotype.execute( this, phenotype, dt );
        }

        // update volume 
        if( functions.updateVolume != null )
        {
            functions.updateVolume.execute( this, phenotype, dt );
            // The following line is needed in every volume regulation method (it sets BioFVM total_volume)
            setTotalVolume( phenotype.volume.total );
        }

        // update geometry
        phenotype.geometry.update( this, phenotype, dt );

        // check for new death events 
        if( phenotype.death.checkForDeath( dt ) )
        {
            // if so, change the cycle model to the current death model 
            phenotype.cycle = phenotype.death.currentModel();

            // also, turn off motility.
            phenotype.motility.isMotile = false;
            phenotype.motility.motilityVector = new double[3];
            functions.update_migration_bias = null;

            // turn off secretion, and reduce uptake by a factor of 10 
            phenotype.secretion.setSecretionToZero();
            phenotype.secretion.scaleUptake( 0.10 );

            // make sure to run the death entry function 
            if( phenotype.cycle.currentPhase().entryFunction != null )
            {
                phenotype.cycle.currentPhase().entryFunction.execute( this, phenotype, dt );
            }
        }

        // advance cycle model (for both cell cycle and death cycle models)
        phenotype.cycle.advance( this, phenotype, dt );
        if( phenotype.flaggedForRemoval )
        {
            flagForRemoval();
            phenotype.flaggedForRemoval = false;
        }
        if( phenotype.flaggedForDivision )
        {
            flagForDivision();
            phenotype.flaggedForDivision = false;
        }
    }

    @Override
    public double getRadius()
    {
        return this.phenotype.geometry.radius;
    }

    //TODO: move somewhere else
    static boolean cell_definitions_by_name_constructed = false;

    @Deprecated
    public double get_total_volume()
    {
        return phenotype.volume.total;
    }

    public CellContainer getContainer()
    {
        if( container == null )
        {
            container = (CellContainer)getMicroenvironment().agentContainer;
        }

        return container;
    }

    public Set<Cell> cells_in_my_container()
    {
        return getContainer().agent_grid.get( get_current_mechanics_voxel_index() );
    }

    public double[] nearest_gradient(int substrate_index)
    {
        return getMicroenvironment().gradient_vector( currentVoxelIndex )[substrate_index];
    }

    @Override
    public String toString()
    {
        return "Cell " + ID + ", type: " + typeName + " (" + type + ") in phase " + phenotype.cycle.name + " at "
                + VectorUtil.print( position, 2 );
    }
}