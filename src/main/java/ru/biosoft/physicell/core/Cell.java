package ru.biosoft.physicell.core;

import java.util.HashSet;
import java.util.Set;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.CellFunctions.Instantiator;
import ru.biosoft.physicell.core.standard.StandardModels;

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
    public static int calls = 0;
    public static int badCalls = 0;

    public static double[] cell_division_orientation;
    CellContainer container;
    int currentMechanicsVoxelIndex;
    int updated_current_mechanics_voxel_index; // keeps the updated voxel index for later adjusting of current voxel index

    public String typeName;

    CellDefinition definition;
    public CustomCellData customData;// = new CustomCellData();
    public CellParameters parameters;// = new CellParameters();
    public CellFunctions functions;// = new CellFunctions();

    public CellState state = new CellState();
    public Phenotype phenotype;// = new Phenotype();
    public boolean isOutOfDomain;
    boolean isMovable;
    public double[] displacement; // this should be moved to state, or made private

    public Cell(CellDefinition cd, Model model)
    {
        super( model );
        type = cd.type;
        typeName = cd.name;
        customData = cd.custom_data.clone();
        parameters = cd.parameters.clone();
        functions = cd.functions.clone();
        phenotype = cd.phenotype.clone();
        phenotype.molecular.sync( this );
        this.definition = cd;
        // cell state should be fine by the default constructor 
        currentMechanicsVoxelIndex = -1;
        updated_current_mechanics_voxel_index = 0;

        isMovable = cd.isMovable;
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
        get_container().flagDivision( this );
    }

    void flagForRemoval()
    {
        get_container().flagRemoval( this );
    }

    @Override
    public void setTotalVolume(double volume)
    {
        super.setTotalVolume( volume );

        // If the new volume is significantly different than the current total volume, adjust all the sub-volumes proportionally. 
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
        double theta = model.rng.UniformRandom() * 6.28318530717959; //rand*2*pi
        double z = 2 * model.rng.UniformRandom() - 1;
        double temp = Math.sqrt( 1 - z * z );
        state.orientation[0] = temp * Math.cos( theta );
        state.orientation[1] = temp * Math.sin( theta );
        state.orientation[2] = z;
    }

    @Override
    public int get_current_mechanics_voxel_index()
    {
        return currentMechanicsVoxelIndex;
    }

    public static Cell createCell(CellDefinition cd, Model model, double[] position) throws Exception
    {
        Instantiator instantiator = cd.functions.instantiator;
        return createCell( instantiator, cd, model, position );
    }

    public static Cell createCell(Instantiator instantiator, CellDefinition cd, Model model, double[] position) throws Exception
    {
        Cell pNew = instantiator == null ? new Cell( cd, model ) : instantiator.execute( cd, model );
        pNew.index = model.getMicroenvironment().getAgentsCount();
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
        if( m.options.simulate2D )
            velocity[2] = 0.0;

        VectorUtil.axpy( position, 1.5 * dt, velocity );
        VectorUtil.axpy( position, -0.5 * dt, prevVelocity );
        prevVelocity = velocity.clone();
        VectorUtil.zero( velocity );

        if( get_container().mesh.isPositionValid( position[0], position[1], position[2] ) )
        {
            updated_current_mechanics_voxel_index = get_container().mesh.nearestVoxelIndex( position );
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
        cell.removeAllAttachedCells();
        cell.removeAllSpringAttachments();
        cell.releaseSubstrates();
        cell.getMicroenvironment().removeAgent( cell );
        cell.get_container().remove_agent( cell );
    }

    public Cell divide() throws Exception
    {
        removeAllAttachedCells();
        removeAllSpringAttachments();

        // conserved quantitites in custom data aer divided in half so that each daughter cell gets half of the original ;
        for( int nn = 0; nn < customData.variables.size(); nn++ )
        {
            if( customData.variables.get( nn ).conserved_quantity )
            {
                customData.variables.get( nn ).value *= 0.5;
            }
        }
        for( int nn = 0; nn < customData.vectorVariables.size(); nn++ )
        {
            if( customData.vectorVariables.get( nn ).conserved_quantity )
            {
                VectorUtil.prod( customData.vectorVariables.get( nn ).value, 0.5 );
            }
        }

        double[] rand_vec = PhysiCellUtilities.UniformOnUnitSphere( model );//cell_division_orientation();
        if( getMicroenvironment().options.simulate2D )
            rand_vec[2] = 0;
        double multiplier = phenotype.geometry.polarity
                * ( rand_vec[0] * state.orientation[0] + rand_vec[1] * state.orientation[1] + rand_vec[2] * state.orientation[2] );

        VectorUtil.naxpy( rand_vec, multiplier, state.orientation );
        VectorUtil.prod( rand_vec, phenotype.geometry.radius );
        double[] pos = VectorUtil.newSum( position, rand_vec );
        Cell child = createCell( functions.instantiator, definition, getModel(), pos );
        child.copyData( this );
        child.copyFunctionPointers( this );
        child.parameters = parameters;

        // evenly divide internalized substrates if these are not actively tracked, they are zero anyway 
        VectorUtil.prod( internalizedSubstrates, 0.5 );
        child.internalizedSubstrates = internalizedSubstrates.clone();

        //change my position to keep the center of mass intact and then see if I need to update my voxel index
        double negative_one_half = -0.5;
        VectorUtil.axpy( position, negative_one_half, rand_vec ); // position = position - 0.5*rand_vec; 

        //If this cell has been moved outside of the boundaries, mark it as such.
        //(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)
        if( !get_container().mesh.isPositionValid( position[0], position[1], position[2] ) )
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

//        state.damage = 0.0;
        state.totalAttackTime = 0;
//        child.state.damage = 0.0;
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

    public static void attachcCells(Cell cell1, Cell cell2)
    {
        cell1.attachCell( cell2 );
        cell2.attachCell( cell1 );
    }

    public static void attachCellsAsSpring(Cell cell1, Cell cell2)
    {
        cell1.attachCellAsSpring( cell2 );
        cell2.attachCellAsSpring( cell1 );
    }

    public static void detachCells(Cell cell1, Cell cell2)
    {
        cell1.detachCell( cell2 );
        cell2.detachCell( cell1 );
    }

    public static void detachCellsAsSpring(Cell cell1, Cell cell2)
    {
        cell1.detachCellAsSpring( cell2 );
        cell2.detachCellAsSpring( cell1 );
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

    public void copyData(Cell copyMe)
    {
        // phenotype=copyMe->phenotype; //it is taken care in set_phenotype
        type = copyMe.type;
        typeName = copyMe.typeName;

        customData = copyMe.customData.clone();
        parameters = copyMe.parameters.clone();

        velocity = copyMe.velocity.clone();
        // expected_phenotype = copy_me. expected_phenotype; //it is taken care in set_phenotype
        sourceSinkTemp1 = copyMe.sourceSinkTemp1.clone();
        sourceSinkTemp2 = copyMe.sourceSinkTemp2.clone();
    }

    void copyFunctionPointers(Cell copyMe)
    {
        functions = copyMe.functions.clone();
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
    public boolean assignPosition(double[] position)
    {
        return assignPosition( position[0], position[1], position[2] );
    }

    public void setPreviousVelocity(double xV, double yV, double zV)
    {
        prevVelocity[0] = xV;
        prevVelocity[1] = yV;
        prevVelocity[2] = zV;
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
        currentMechanicsVoxelIndex = get_container().mesh.nearestVoxelIndex( position );

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

        if( !get_container().mesh.isPositionValid( x, y, z ) )
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
        functions.updateMigration = null;

        // make sure to run the death entry function 
        if( phenotype.cycle.data.currentPhase().entryFunction != null )
        {
            phenotype.cycle.data.currentPhase().entryFunction.execute( this, phenotype, 0.0 );
        }
    }


    
 
    
    public double addPotentials(Cell other)
    {
        if( this.ID == other.ID )
            return 0;

        calls++;
        if( other.isOutOfDomain || !other.isActive )
            badCalls++;

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
            double temp_a = 1 - distance / max_interactive_distance;
            temp_a *= temp_a; // (1-d/S)^2 

            // August 2017 - back to the original if both have same coefficient 
            // May 2022 - back to original if both affinities are 1
            int ii = getModel().getCellDefinitionIndex( type );//CellDefinition.getCellDefinition( type ).;
            int jj = getModel().getCellDefinitionIndex( other.type );//find_cell_definition_index( other.type );

            double adhesion_ii = m1.cellCellAdhesionStrength * m1.cellAdhesionAffinities[jj];
            double adhesion_jj = m2.cellCellAdhesionStrength * m2.cellAdhesionAffinities[ii];

            // double effective_adhesion = sqrt( phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
            double effective_adhesion = Math.sqrt( adhesion_ii * adhesion_jj );
            temp_a *= effective_adhesion;
            temp_r -= temp_a;
            state.neighbors.add( other );// move here in 1.10.2 so non-adhesive cells also added. 
        }
        if( Math.abs( temp_r ) < 1e-16 )
            return 0;
        temp_r /= distance;
        VectorUtil.axpy( velocity, temp_r, displacement );
        return temp_r;
    }

    public void updateMotilityVector(double dt) throws Exception
    {
        if( phenotype.motility.isMotile == false )
        {
            phenotype.motility.motilityVector = new double[3];//.assign( 3, 0.0 ); 
            return;
        }

        if( phenotype.motility.persistenceTime < dt || model.rng.checkRandom( dt / phenotype.motility.persistenceTime ) )
        {
            double[] randvec;
            if( phenotype.motility.restrictTo2D )
            {
                randvec = PhysiCellUtilities.UniformOnUnitCircle( model );
            }
            else
            {
                randvec = PhysiCellUtilities.UniformOnUnitCircle( model );
            }

            // if the update_bias_vector function is set, use it  
            if( functions.updateMigration != null )
            {
                functions.updateMigration.execute( this, phenotype, dt );
            }

            phenotype.motility.motilityVector = VectorUtil.newProd( phenotype.motility.migrationBiasDirection,
                    phenotype.motility.migrationBias );

            double one_minus_bias = 1.0 - phenotype.motility.migrationBias;
            VectorUtil.axpy( ( phenotype.motility.motilityVector ), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector
            VectorUtil.normalize( ( phenotype.motility.motilityVector ) );
            VectorUtil.prod( phenotype.motility.motilityVector, phenotype.motility.migrationSpeed );
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
        double[] cornerPoint = VectorUtil.newSum( my_voxel_center, other_voxel_center );
        VectorUtil.prod( cornerPoint, 0.5 );
        double distance_squared = ( cornerPoint[0] - pCell.position[0] ) * ( cornerPoint[0] - pCell.position[0] )
                + ( cornerPoint[1] - pCell.position[1] ) * ( cornerPoint[1] - pCell.position[1] )
                + ( cornerPoint[2] - pCell.position[2] ) * ( cornerPoint[2] - pCell.position[2] );
        if( distance_squared > max_interactive_distance * max_interactive_distance )
            return false;
        return true;
    }

    public double[] nearest_density_vector()
    {
        if( currentVoxelIndex == -1 )
            return new double[] { -1};
        return getMicroenvironment().nearestDensity( this.currentVoxelIndex );//current_voxel_index );
    }

    public void ingestCell(Cell cellToEat)
    {
        // don't ingest self 
        if( cellToEat == this )
            return;

        // don't ingest a cell that's already ingested 
        if( cellToEat.phenotype.volume.total < 1e-15 )
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
            cellToEat.phenotype.death.dead = true;
            // set secretion and uptake to zero 
            cellToEat.phenotype.secretion.setSecretionToZero();
            cellToEat.phenotype.secretion.setUptakeToZero();

            // deactivate all custom function 
            cellToEat.functions.customCellRule = null;
            cellToEat.functions.updatePhenotype = null;
            cellToEat.functions.contact = null;

            // should set volume fuction to NULL too! 
            cellToEat.functions.updateVolume = null;

            // set cell as unmovable and non-secreting 
            cellToEat.isMovable = false;
            cellToEat.isActive = false;

            // absorb all the volume(s)

            // absorb fluid volume (all into the cytoplasm) 
            phenotype.volume.cytoplasmic_fluid += cellToEat.phenotype.volume.fluid;
            cellToEat.phenotype.volume.cytoplasmic_fluid = 0.0;

            // absorb nuclear and cyto solid volume (into the cytoplasm) 
            phenotype.volume.cytoplasmic_solid += cellToEat.phenotype.volume.cytoplasmic_solid;
            cellToEat.phenotype.volume.cytoplasmic_solid = 0.0;

            phenotype.volume.cytoplasmic_solid += cellToEat.phenotype.volume.nuclear_solid;
            cellToEat.phenotype.volume.nuclear_solid = 0.0;

            // consistency calculations 

            phenotype.volume.fluid = phenotype.volume.nuclear_fluid + phenotype.volume.cytoplasmic_fluid;
            cellToEat.phenotype.volume.fluid = 0.0;

            phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + phenotype.volume.nuclear_solid;
            cellToEat.phenotype.volume.solid = 0.0;

            // no change to nuclear volume (initially) 
            cellToEat.phenotype.volume.nuclear = 0.0;
            cellToEat.phenotype.volume.nuclear_fluid = 0.0;

            phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_solid + phenotype.volume.cytoplasmic_fluid;
            cellToEat.phenotype.volume.cytoplasmic = 0.0;

            phenotype.volume.total = phenotype.volume.nuclear + phenotype.volume.cytoplasmic;
            cellToEat.phenotype.volume.total = 0.0;

            phenotype.volume.fluid_fraction = phenotype.volume.fluid / ( phenotype.volume.total + 1e-16 );
            cellToEat.phenotype.volume.fluid_fraction = 0.0;

            phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / ( phenotype.volume.nuclear_solid + 1e-16 );

            // update corresponding BioFVM parameters (self-consistency) 
            setTotalVolume( phenotype.volume.total );
            cellToEat.setTotalVolume( 0.0 );

            // absorb the internalized substrates 

            // multiply by the fraction that is supposed to be ingested (for each substrate) 
            //            *(pCell_to_eat.internalized_substrates) *= 
            //                *(pCell_to_eat.fraction_transferred_when_ingested); // 
            VectorUtil.prod( cellToEat.internalizedSubstrates, cellToEat.fractionTransferredIngested );

            //            *internalized_substrates += *(pCell_to_eat.internalized_substrates);
            VectorUtil.sum( internalizedSubstrates, cellToEat.internalizedSubstrates );
            //            int n_substrates = internalizedSubstrates.length;
            cellToEat.internalizedSubstrates = new double[internalizedSubstrates.length];//assign( n_substrates , 0.0 );    

            // trigger removal from the simulation 
            // pCell_to_eat.die(); // I don't think this is safe if it's in an OpenMP loop 

            // flag it for removal 
            // pCell_to_eat.flag_for_removal(); 

            // remove all adhesions 
            // pCell_to_eat.remove_all_attached_cells();

        }

        // things that have their own thread safety 
        cellToEat.flagForRemoval();
        cellToEat.removeAllAttachedCells();
        cellToEat.removeAllSpringAttachments();
    }

    public void attackCell(Cell pCellToAttack, double dt)
    {
        if( pCellToAttack == this )
            return; // don't attack self 

        if( pCellToAttack.phenotype.death.dead == true || pCellToAttack.phenotype.volume.total < 1e-15 )
            return; // don't attack a dead or tiny cell 

        // make this thread safe 
        //        #pragma omp critical
        {
            // std::cout << this.type_name << " attacks " << pCell_to_attack.type_name << std::endl;
            pCellToAttack.phenotype.cellIntegrity.damage += phenotype.cellInteractions.damageRate * dt;
            pCellToAttack.state.totalAttackTime += dt;

            phenotype.cellInteractions.total_damage_delivered += phenotype.cellInteractions.damageRate * dt; 
        }
    }

    public void fuseCell(Cell pCellToFuse)
    {

        if( pCellToFuse.phenotype.volume.total < 1e-15 || this == pCellToFuse )
            return; // don't ingest a cell that's already fused or fuse self 

        // make this thread safe 
        {

            // set new position at center of volume 
            // x_new = (vol_B * x_B + vol_S * x_S ) / (vol_B + vol_S )
            double[] newPosition = VectorUtil.newProd( position, phenotype.volume.total );
            double total_volume = phenotype.volume.total;
            total_volume += pCellToFuse.phenotype.volume.total;

            VectorUtil.axpy( newPosition, pCellToFuse.phenotype.volume.total, pCellToFuse.position ); // vol_B*x_B + vol_S*x_S
            VectorUtil.div( newPosition, total_volume );

            double xL = m.mesh.boundingBox[0];
            double xU = m.mesh.boundingBox[3];

            double yL = m.mesh.boundingBox[1];
            double yU = m.mesh.boundingBox[4];

            double zL = m.mesh.boundingBox[2];
            double zU = m.mesh.boundingBox[5];

            if( newPosition[0] < xL || newPosition[0] > xU || newPosition[1] < yL || newPosition[1] > yU || newPosition[2] < zL
                    || newPosition[2] > zU )
            {
                System.out.println( "cell fusion at " + newPosition + " violates domain bounds" );
                System.out.println( m.mesh.boundingBox );
            }
            position = newPosition;
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
            pCellToFuse.functions.customCellRule = null;
            pCellToFuse.functions.updatePhenotype = null;
            pCellToFuse.functions.contact = null;
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

    public void lyseCell()
    {
        // don't lyse a cell that's already lysed 
        if( phenotype.volume.total < 1e-15 )
            return;

        // flag for removal 
        flagForRemoval(); // should be safe now 

        // mark it as dead 
        phenotype.death.dead = true;

        // set secretion and uptake to zero 
        phenotype.secretion.setSecretionToZero();
        phenotype.secretion.setUptakeToZero();

        // deactivate all custom function 
        functions.customCellRule = null;
        functions.updatePhenotype = null;
        functions.contact = null;

        // remove all adhesions 

        removeAllAttachedCells();

        // set volume to zero 
        setTotalVolume( 0.0 );

        // set cell as unmovable and non-secreting 
        isMovable = false;
        isActive = false;
    }

    public void convert(CellDefinition cd)
    {
        type = cd.type;
        typeName = cd.name;
        customData = cd.custom_data.clone();
        parameters = cd.parameters.clone();
        functions = cd.functions.clone();
        phenotype = cd.phenotype.clone();
        // is_movable = true;
        // is_out_of_domain = false;
        // displacement.resize(3,0.0); // state? 
        assignOrientation();
        setTotalVolume( phenotype.volume.total );
    }

    public void advanceBundledPhenotype(double dt, boolean rulesEnabled) throws Exception
    {
        // New March 2022
        // perform transformations 
        StandardModels.standard_cell_transformations( this, this.phenotype, dt );

        // New March 2023 in Version 1.12.0 
        // call the rules-based code to update the phenotype 
        if( rulesEnabled )
            model.rules.applyRuleset( this );
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

        // update integrity 
        phenotype.cellIntegrity.advanceDamage( dt ); 
        
        // check for new death events 
        if( phenotype.death.checkForDeath( model.getRNG(), dt ) )
        {
            // if so, change the cycle model to the current death model 
            phenotype.cycle = phenotype.death.currentModel();
            phenotype.motility.disable();
            functions.updateMigration = null;

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
        return getContainer().agentGrid.get( get_current_mechanics_voxel_index() );
    }

    public double[] nearest_gradient(int substrate_index)
    {
        return getMicroenvironment().getGradient( currentVoxelIndex )[substrate_index];
    }

    public double[] nearestGradient(String substrate)
    {
        int index = getMicroenvironment().findDensityIndex( substrate );
        return getMicroenvironment().getGradient( currentVoxelIndex )[index];
    }

    public Set<Cell> nearby_interacting_cells()
    {
        return find_nearby_interacting_cells( this );
    }


    Set<Cell> find_nearby_interacting_cells(Cell pCell)
    {
        Set<Cell> neighbors = new HashSet<>();

        for( Cell neighbor : pCell.get_container().agentGrid.get( pCell.get_current_mechanics_voxel_index() ) )
        {
            double distance = VectorUtil.dist( neighbor.position, pCell.position );
            if( distance <= pCell.phenotype.mechanics.relMaxAdhesionDistance * pCell.phenotype.geometry.radius
                    + ( neighbor ).phenotype.mechanics.relMaxAdhesionDistance * ( neighbor ).phenotype.geometry.radius
                    && ( neighbor ) != pCell )
            {
                neighbors.add( neighbor );
            }
        }

        for( int neighbor_voxel_index : pCell.get_container().mesh.moore_connected_voxel_indices[pCell
                .get_current_mechanics_voxel_index()] )
        {
            if( !Cell.isNeighborVoxel( pCell, pCell.get_container().mesh.voxels[pCell.get_current_mechanics_voxel_index()].center,
                    pCell.get_container().mesh.voxels[neighbor_voxel_index].center, neighbor_voxel_index ) )
                continue;

            for( Cell neighbor : pCell.get_container().agentGrid.get( neighbor_voxel_index ) )
            {
                double distance = VectorUtil.dist( neighbor.position, pCell.position );
                if( distance <= pCell.phenotype.mechanics.relMaxAdhesionDistance * pCell.phenotype.geometry.radius
                        + ( neighbor ).phenotype.mechanics.relMaxAdhesionDistance * ( neighbor ).phenotype.geometry.radius
                        && ( neighbor ) != pCell )
                {
                    neighbors.add( neighbor );
                }
            }
        }
        return neighbors;
    }

    @Override
    public String toString()
    {
        return "Cell " + ID + ", type: " + typeName + " (" + type + ") in phase " + phenotype.cycle.name + " at "
                + VectorUtil.print( position, 2 );
    }
}