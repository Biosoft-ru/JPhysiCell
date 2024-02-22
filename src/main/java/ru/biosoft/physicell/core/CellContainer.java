package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ru.biosoft.physicell.biofvm.AgentContainer;
import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.CartesianMesh;
import ru.biosoft.physicell.biofvm.Microenvironment;
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
public class CellContainer extends AgentContainer
{
    Set<Cell> cellsReadyToDivide; // the index of agents ready to divide
    Set<Cell> cellsReadyToDie;
    int boundary_condition_for_pushed_out_agents; // what to do with pushed out cells
    boolean initialzed = false;

    public CartesianMesh underlying_mesh = new CartesianMesh();
    double[] max_cell_interactive_distance_in_voxel;
    public int numDivisionsCurStep = 0;
    public int numDeathsCurStep = 0;

    double last_diffusion_time = 0.0;
    double lastCellCycleTime = 0.0;
    double last_mechanics_time = 0.0;

    public List<Set<Cell>> agent_grid;
    List<Set<Cell>> agents_in_outer_voxels;

    public CellContainer()
    {
        boundary_condition_for_pushed_out_agents = PhysiCellConstants.default_boundary_condition_for_pushed_out_agents;
        cellsReadyToDivide = new HashSet<>();
        cellsReadyToDie = new HashSet<>();
    }

    void initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, double voxel_size)
    {
        initialize( x_start, x_end, y_start, y_end, z_start, z_end, voxel_size, voxel_size, voxel_size );
    }

    void initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, double dx, double dy,
            double dz)
    {
        //        allCells = BasicAgent.allBasicAgents;
        boundary_condition_for_pushed_out_agents = PhysiCellConstants.default_boundary_condition_for_pushed_out_agents;
        cellsReadyToDivide = new HashSet<>();
        cellsReadyToDie = new HashSet<>();
        underlying_mesh.resize( x_start, x_end, y_start, y_end, z_start, z_end, dx, dy, dz );
        max_cell_interactive_distance_in_voxel = new double[underlying_mesh.voxels.length];
        agent_grid = new ArrayList<>();//[underlying_mesh.voxels.length];
        for( int i = 0; i < underlying_mesh.voxels.length; i++ )
            agent_grid.add( new HashSet<Cell>() );
        max_cell_interactive_distance_in_voxel = new double[underlying_mesh.voxels.length];
        agents_in_outer_voxels = new ArrayList<>( 6 );
        for( int i = 0; i < 6; i++ )
            agents_in_outer_voxels.add( new HashSet<Cell>() );
        //            underlying_mesh.resize(x_start, x_end, y_start, y_end, z_start, z_end , dx, dy, dz);
        //            agent_grid.resize(underlying_mesh.voxels.size());
        //            max_cell_interactive_distance_in_voxel.resize(underlying_mesh.voxels.size(), 0.0);
        //            agents_in_outer_voxels.resize(6);

    }

    public void updateAllCells(Microenvironment m, double t, double dt) throws Exception
    {
        updateAllCells( m, t, dt, dt, dt );
    }

    public void updateAllCells(Microenvironment m, double t, double phenotypeDT, double mechanicsDT, double diffusionDT)
            throws Exception
    {
        // secretions and uptakes. Syncing with BioFVM is automated. 
        //            #pragma omp parallel for 
        Set<Cell> cells = m.getAgents( Cell.class );
        for( Cell cell : cells )
        {

            if( !cell.isOutOfDomain )
            {
                cell.phenotype.secretion.advance( cell, cell.phenotype, diffusionDT );
            }
        }
        //if it is the time for running cell cycle, do it!
        double timeSinceLastCycle = t - lastCellCycleTime;
        double phenotypeDTtolerance = 0.001 * phenotypeDT;
        double mechanicsDTtolerance = 0.001 * mechanicsDT;

        // intracellular update. called for every diffusion_dt, but actually depends on the intracellular_dt of each cell (as it can be noisy)
        //            #pragma omp parallel for 
        for( Cell cell : cells )
        {
            if( !cell.isOutOfDomain && initialzed )
            {
                if( cell.phenotype.intracellular != null && cell.phenotype.intracellular.need_update() )
                {
                    if( ( cell.functions.pre_update_intracellular != null ) )
                        cell.functions.pre_update_intracellular.execute( cell, cell.phenotype, diffusionDT );

                    cell.phenotype.intracellular.update( cell, cell.phenotype, diffusionDT );

                    if( cell.functions.post_update_intracellular != null )
                        cell.functions.post_update_intracellular.execute( cell, cell.phenotype, diffusionDT );
                }
            }
        }

        if( Math.abs( timeSinceLastCycle - phenotypeDT ) < phenotypeDTtolerance || !initialzed )
        {
            // Reset the max_radius in each voxel. It will be filled in set_total_volume
            // It might be better if we calculate it before mechanics each time 
            // std::fill(max_cell_interactive_distance_in_voxel.begin(), max_cell_interactive_distance_in_voxel.end(), 0.0);

            if( !initialzed )
            {
                timeSinceLastCycle = phenotypeDT;
            }

            // new as of 1.2.1 -- bundles cell phenotype parameter update, volume update, geometry update, 
            // checking for death, and advancing the cell cycle. Not motility, though. (that's in mechanics)
            //                #pragma omp parallel for 
            for( Cell cell : cells )
            {
                if( !cell.isOutOfDomain )
                {
                    cell.advanceBundledPhenotype( timeSinceLastCycle );
                }
            }

            // process divides / removes 
            for( Cell cell : cellsReadyToDivide )
                cell.divide();

            for( Cell cell : cellsReadyToDie )
                cell.die();

            numDivisionsCurStep += cellsReadyToDivide.size();
            numDeathsCurStep += cellsReadyToDie.size();

            cellsReadyToDie.clear();
            cellsReadyToDivide.clear();
            lastCellCycleTime = t;
        }

        double time_since_last_mechanics = t - last_mechanics_time;

        if( Math.abs( time_since_last_mechanics - mechanicsDT ) < mechanicsDTtolerance || !initialzed )
        {
            if( !initialzed )
            {
                time_since_last_mechanics = mechanicsDT;
            }
            // new February 2018  if we need gradients, compute them
            if( m.options.calculate_gradients )
            {
                m.computeAllGradientVectors();
            }

            //TOD: commented by now
            // end of new in Feb 2018 
            // perform interactions -- new in June 2020 
            //                #pragma omp parallel for 
            for( Cell cell : cells )
            {
                if( cell.functions.contact != null && !cell.isOutOfDomain )
                {
                    StandardModels.evaluate_interactions( cell, cell.phenotype, time_since_last_mechanics );
                }
            }
            // perform custom computations 
            //                #pragma omp parallel for 
            for( Cell cell : cells )
            {
                if( cell.functions.customCellRule != null && !cell.isOutOfDomain )
                {
                    cell.functions.customCellRule.execute( cell, cell.phenotype, time_since_last_mechanics );
                }
            }
            // update velocities 
            //                #pragma omp parallel for 
            for( Cell cell : cells )
            {
                if( cell.functions.updateVelocity != null && !cell.isOutOfDomain && cell.isMovable )
                {
                    cell.functions.updateVelocity.execute( cell, cell.phenotype, time_since_last_mechanics );
                }
            }
            // new March 2023: 
            // dynamic spring attachments, followed by built-in springs
            if( !PhysiCellSettings.disable_automated_spring_adhesions )
            {
                //                    #pragma omp parallel for 
                for( Cell cell : cells )
                {
                    StandardModels.dynamic_spring_attachments( cell, cell.phenotype, time_since_last_mechanics );
                }
                //                    #pragma omp parallel for 
                for( Cell cell : cells )
                {
                    if( cell.isMovable )
                    {
                        for( Cell pC1 : cell.state.springAttachments )
                        {
                            StandardModels.standard_elastic_contact_function( cell, cell.phenotype, pC1, pC1.phenotype,
                                    time_since_last_mechanics );
                        }
                    }
                }
            }

            // new March 2022: 
            // run standard interactions (phagocytosis, attack, fusion) here 
            //                #pragma omp parallel for 
            for( Cell cell : cells )
            {
                StandardModels.standard_cell_cell_interactions( cell, cell.phenotype, time_since_last_mechanics );
            }
            // super-critical to performance! clear the "dummy" cells from phagocytosis / fusion
            // otherwise, comptuational cost increases at polynomial rate VERY fast, as O(10,000) 
            // dummy cells of size zero are left ot interact mechanically, etc. 
            if( cellsReadyToDie.size() > 0 )
            {
                /*
                std::cout << "\tClearing dummy cells from phagocytosis and fusion events ... " << std::endl; 
                std::cout << "\t\tClearing " << cells_ready_to_die.size() << " cells ... " << std::endl; 
                // there might be a lot of "dummy" cells ready for removal. Let's do it.        
                */
                //                    for( int i=0; i < cells_ready_to_die.size(); i++ )
                //                    { cells_ready_to_diei].die(); }

                for( Cell cell : cellsReadyToDie )
                    cell.die();
                cellsReadyToDie.clear();
            }
            // update positions         
            //                #pragma omp parallel for 
            for( Cell cell : cells )
            {
                if( !cell.isOutOfDomain && cell.isMovable )
                {
                    cell.updatePosition( time_since_last_mechanics );
                }
            }

            // When somebody reviews this code, let's add proper braces for clarity!!! 
            // Update cell indices in the container
            for( Cell cell : cells )
            {
                if( !cell.isOutOfDomain && cell.isMovable )
                    cell.updateVoxelInContainer();
            }
            last_mechanics_time = t;
        }
        initialzed = true;
    }

    @Override
    public void register_agent(BasicAgent agent)
    {
        int index = ( (Cell)agent ).get_current_mechanics_voxel_index();
        agent_grid.get( index ).add( (Cell)agent );
        //        agent_grid[agent.get_current_mechanics_voxel_index()].push_back( agent );
    }


    @Override
    public void remove_agent(BasicAgent agent)
    {
        remove_agent_from_voxel( agent, agent.get_current_mechanics_voxel_index() );
        return;
    }

    @Override
    public void add_agent_to_outer_voxel(BasicAgent agent)
    {
        int escaping_face = find_escaping_face_index( agent );
        agents_in_outer_voxels.get( escaping_face ).add( (Cell)agent );
        ( (Cell)agent ).isOutOfDomain = true;
        return;
    }
    //
    @Override
    public void remove_agent_from_voxel(BasicAgent agent, int voxel_index)
    {
        if( voxel_index < 0 )
        {
            return;
        }
        agent_grid.get( voxel_index ).remove( agent );
        //        int delete_index = 0;
        //        while( agent_grid[voxel_index][delete_index] != agent )
        //        {
        //            delete_index++;
        //        }
        //
        //        // move last item to index location
        //        agent_grid[voxel_index][delete_index] = agent_grid[voxel_index][agent_grid[voxel_index].size() - 1];
        //        // shrink the vector
        //        agent_grid[voxel_index].pop_back();
        //
        //        return;
    }
    //
    @Override
    public void add_agent_to_voxel(BasicAgent agent, int voxel_index)
    {
        agent_grid.get( voxel_index ).add( (Cell)agent );//.push_back( agent );
    }
    //
    boolean contain_any_cell(int voxel_index)
    {
        // Let's replace this with clearer statements. 
        return agent_grid.get( voxel_index ).size() == 0 ? false : true;
    }

    int find_escaping_face_index(BasicAgent agent)
    {
        if( agent.position[0] <= agent.get_container().underlying_mesh.boundingBox[PhysiCellConstants.mesh_min_x_index] )
        {
            return PhysiCellConstants.mesh_lx_face_index;
        }
        if( agent.position[0] >= agent.get_container().underlying_mesh.boundingBox[PhysiCellConstants.mesh_max_x_index] )
        {
            return PhysiCellConstants.mesh_ux_face_index;
        }
        if( agent.position[1] <= agent.get_container().underlying_mesh.boundingBox[PhysiCellConstants.mesh_min_y_index] )
        {
            return PhysiCellConstants.mesh_ly_face_index;
        }
        if( agent.position[1] >= agent.get_container().underlying_mesh.boundingBox[PhysiCellConstants.mesh_max_y_index] )
        {
            return PhysiCellConstants.mesh_uy_face_index;
        }
        if( agent.position[2] <= agent.get_container().underlying_mesh.boundingBox[PhysiCellConstants.mesh_min_z_index] )
        {
            return PhysiCellConstants.mesh_lz_face_index;
        }
        if( agent.position[2] >= agent.get_container().underlying_mesh.boundingBox[PhysiCellConstants.mesh_max_z_index] )
        {
            return PhysiCellConstants.mesh_uz_face_index;
        }
        return -1;
    }

    void flag_cell_for_division(Cell pCell)
    {
        cellsReadyToDivide.add( pCell );
    }

    void flag_cell_for_removal(Cell pCell)
    {
        cellsReadyToDie.add( pCell );
    }

    public static CellContainer createCellContainer(Microenvironment m, double mechanicsVoxelSize)
    {
        CellContainer cellContainer = new CellContainer();
        cellContainer.initialize( m.mesh.boundingBox[0], m.mesh.boundingBox[3], m.mesh.boundingBox[1], m.mesh.boundingBox[4],
                m.mesh.boundingBox[2], m.mesh.boundingBox[5], mechanicsVoxelSize );
        m.agentContainer = (AgentContainer)cellContainer;
        return cellContainer;
    }
}
