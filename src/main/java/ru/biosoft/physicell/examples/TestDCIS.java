package ru.biosoft.physicell.examples;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.calculate_distance_to_membrane;
import ru.biosoft.physicell.core.Output;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.ui.Visualizer.Section;

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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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
public class TestDCIS
{
    static double o2_conc = 7.1930;
    static double duct_radius = 158.75;
    private static String resultPath = "C:/Users/Damag/BIOFVM/DCIS/POV";
    //    private static String resultName = "result";
    private static int zSlice = 100;

    static double dt = 0.01; // reaction-diffusion time-step
    static double mechanics_dt = 0.1;
    static double cell_cycle_dt = 6;
    static double t_output_interval = 60; // 1.0; 
    static double t_max = 365 * 24 * 60;
    private static List<Visualizer> resultListeners = new ArrayList<>();
    static
    {
        //        for( int i = 0; i < 20; i++ )
        //        {
        //            addListener( "Result_Y_" + i, Section.Y, -200 + i * 20 );
        //            addListener( "Result_Z_" + i, Section.Z, -200 + i * 20 );
        //        }
        //        for (int i=0; i<60; i++)
        //            addListener( "Result_X" + i, Section.X, -200 + i * 20 );
        addListener( "Result_X", Section.X, 0 );
        addListener( "Result_Y", Section.Y, 0 );
        addListener( "Result_Z", Section.Z, 0 );
    }

    private static void addListener(String name, Section sec, int slice)
    {
        Visualizer visualizer3 = new Visualizer( resultPath, name, sec, slice );
        visualizer3.setSaveImage( false );
        //        visualizer3.setDrawDensity( false );
        visualizer3.setMaxDensity( o2_conc * 2 );
        resultListeners.add( visualizer3 );
    }

    public static void main(String[] args) throws Exception
    {
        double dx = 20;
        double dy = 20;
        double dz = 20;
        // figure out the bounding box 
        double[] bounding_box = new double[6];
        bounding_box[PhysiCellConstants.mesh_min_x_index] = -200;
        bounding_box[PhysiCellConstants.mesh_max_x_index] = 1000;
        bounding_box[PhysiCellConstants.mesh_min_y_index] = -200;
        bounding_box[PhysiCellConstants.mesh_max_y_index] = 200;
        bounding_box[PhysiCellConstants.mesh_min_z_index] = -200;
        bounding_box[PhysiCellConstants.mesh_max_z_index] = 200;
        // create a microenvironment
        Microenvironment m = new Microenvironment( "substrate scale", "minutes", "microns" );
        // add a microenvironment for simulating substrates 	
        m.setDensity( 0, "oxygen", "mmHg" );
        // microenvironment.add_density( "glucose" , "dimensionless" );

        m.resizeSpace( bounding_box[0], bounding_box[3], bounding_box[1], bounding_box[4], bounding_box[2], bounding_box[5], dx, dy, dz );

        // Cell_Container 
        double mechanics_voxel_size = 30;
        CellContainer.createCellContainer( m, mechanics_voxel_size );

        for( int n = 0; n < m.number_of_voxels(); n++ )
            m.getDensity( n )[0] = o2_conc;

        // register substrates properties 
        m.diffusion_coefficients[0] = 1.0e5; // microns^2 / min 
        m.decay_rates[0] = 0.1;

        CellDefinition cd = StandardModels.createDefaultCellDefinition( "tumor cell", m );
        // set default cell cycle model 
        cd.phenotype.cycle = StandardModels.Ki67_advanced;
        // set default_cell_functions; 
        cd.functions.updatePhenotype = new StandardModels.update_cell_and_death_parameters_O2_based();

        int Q_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_negative );
        int oxygen_substrate_index = m.findDensityIndex( "oxygen" );
        int K1_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_positive_premitotic );
        int K2_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_positive_postmitotic );
        int apoptosis_model_index = cd.phenotype.death.find_death_model_index( "apoptosis" );
        int necrosis_model_index = cd.phenotype.death.find_death_model_index( "necrosis" );
        // cells apoptose after about 7 days 
        cd.phenotype.death.rates.set( apoptosis_model_index, 1.0 / ( 7.0 * 24.0 * 60.0 ) );
        // initially no necrosis 
        cd.phenotype.death.rates.set( necrosis_model_index, 0.0 );

        // make sure the cells uptake oxygen at the right rate 
        cd.phenotype.secretion.uptakeRates[oxygen_substrate_index] = 10;

        // update transition times 
        cd.phenotype.cycle.data.setTransitionRate( Q_index, K1_index, 1.0 / ( 8.5 * 60.0 ) );//  transition_rate(Q_index,K1_index) = 1.0 / ( 8.5 * 60.0 ); 
        cd.phenotype.cycle.data.setTransitionRate( K1_index, K2_index, 1.0 / ( 13.0 * 60.0 ) );
        cd.phenotype.cycle.data.setTransitionRate( K2_index, Q_index, 1.0 / ( 2.5 * 60.0 ) );

        // let's make necrotic cells survive 6 hours in minimal oxygen conditions  
        cd.parameters.max_necrosis_rate = 1.0 / ( 6.0 * 60.0 );


        cd.functions.calculate_distance_to_membrane = new distance_to_membrane_duct();
        CellDefinition.registerCellDefinition( cd );
        double cell_radius = 10;
        double sphere_radius = duct_radius - 10;

        List<double[]> positions = createSphere( cell_radius, sphere_radius );
        //add Dirichlet node for all the voxels located outside of the duct
        //	std::vector<double> dirichlet_o2( 1 , o2_conc );
        double[] dirichlet_o2 = new double[] {o2_conc};
        for( int i = 0; i < m.number_of_voxels(); i++ )
        {
            if( m.voxels( i ).center[0] >= 0 )
            {
                if( Math.sqrt( m.voxels( i ).center[1] * m.voxels( i ).center[1]
                        + m.voxels( i ).center[2] * m.voxels( i ).center[2] ) > duct_radius )
                    m.add_dirichlet_node( i, dirichlet_o2 );
            }
            else
            {
                if( VectorUtil.dist( m.voxels( i ).center, new double[] {0.0, 0.0, 0.0} ) > duct_radius )
                    m.add_dirichlet_node( i, dirichlet_o2 );
            }
        }

        for( double[] pos : positions )
        {
            if( pos[0] > 0 )
                continue;
            Cell cell = Cell.createCell( cd, m, pos );
            cell.phenotype.cycle.data.currentPhaseIndex = Q_index;
            if( cell.phenotype.cycle.currentPhase().entryFunction != null )
                cell.phenotype.cycle.currentPhase().entryFunction.execute( cell, cell.phenotype, dt );
            // pCell.parameters.necrosis_type= PhysiCell_constants::deterministic_necrosis;
        }
        System.out.println( "Cells created: " + m.getAgentsCount() );

        for( BasicAgent agent : m.getAgents() )
            agent.setUptakeConstants( dt );

        for( Visualizer listener : resultListeners )
            listener.init();

        double t = 0.0;
        double t_next_output_time = 0;
        int counter = 0;
        while( t < t_max )
        {
            if( Math.abs( t - t_next_output_time ) < 0.0001 )
            {
                t_next_output_time += t_output_interval;
                for( Visualizer listener : resultListeners )
                    listener.saveResult( m, t );
                Output.writePov( m.getAgents( Cell.class ), t, 1000.0, resultPath + "/result_" + counter + ".pov" );
                counter++;
            }
            m.simulate_cell_sources_and_sinks( dt );
            m.simulate_diffusion_decay( dt );
            ( (CellContainer)m.agentContainer ).updateAllCells( m, t, cell_cycle_dt, mechanics_dt, dt );
            t += dt;
        }
        for( Visualizer listener : resultListeners )
            listener.finish();
    }

    void log_output(double t, int output_index, Microenvironment m)
    {
        int num_new_cells = 0;
        int num_deaths = 0;

        StringBuilder sb = new StringBuilder();
        sb.append( "time: \n" );
        Set<BasicAgent> agents = m.getAgents();
        num_new_cells = t == 0 ? m.getAgentsCount() : ( (CellContainer)m.agentContainer ).num_divisions_in_current_step;
        num_deaths = ( (CellContainer)m.agentContainer ).num_deaths_in_current_step;
        sb.append( "total number of agents (newly born, deaths): " + agents.size() + "(" + num_new_cells + ", " + num_deaths + ")\n" );
        sb.append( t + "\t" + agents.size() + "\t" + num_new_cells + "\t" + num_deaths + "\n" );//+BioFVM::stopwatch_value()+ std::endl; 

        ( (CellContainer)m.agentContainer ).num_divisions_in_current_step = 0;
        ( (CellContainer)m.agentContainer ).num_deaths_in_current_step = 0;
    }

    static List<double[]> createSphere(double cell_radius, double sphere_radius)
    {
        List<double[]> result = new ArrayList<>();
        int xc = 0, yc = 0, zc = 0;
        double x_spacing = cell_radius * Math.sqrt( 3 );
        double y_spacing = cell_radius * 2;
        double z_spacing = cell_radius * Math.sqrt( 3 );

        for( double z = -sphere_radius; z < sphere_radius; z += z_spacing, zc++ )
            for( double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++ )
                for( double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++ )
                {
                    double[] tempPoint = new double[3];
                    tempPoint[0] = x + ( zc % 2 ) * 0.5 * cell_radius;
                    tempPoint[1] = y + ( xc % 2 ) * cell_radius;
                    tempPoint[2] = z;

                    if( Math.sqrt( VectorUtil.norm_squared( tempPoint ) ) < sphere_radius )
                        result.add( tempPoint );
                }
        return result;
    }


    public static class distance_to_membrane_duct implements calculate_distance_to_membrane
    {
        @Override
        public double execute(Cell pCell, Phenotype phenotype, double dummy)
        {
            double epsillon = 1e-7;
            //Note that this function assumes that duct cap center is located at <0, 0, 0>
            if( pCell.position[0] >= 0 ) // Cell is within the cylinder part of the duct
            {
                double distance_to_x_axis = Math.sqrt( pCell.position[1] * pCell.position[1] + pCell.position[2] * pCell.position[2] );
                distance_to_x_axis = Math.max( distance_to_x_axis, epsillon ); // prevents division by zero
                pCell.displacement[0] = 0;
                pCell.displacement[1] = -pCell.position[1] / distance_to_x_axis;
                pCell.displacement[2] = -pCell.position[2] / distance_to_x_axis;
                return Math.abs( duct_radius - distance_to_x_axis );
            }

            // Cell is inside the cap of the duct
            double distance_to_origin = VectorUtil.dist( pCell.position, new double[] {0.0, 0.0, 0.0} ); // distance to the origin 
            distance_to_origin = Math.max( distance_to_origin, epsillon ); // prevents division by zero
            pCell.displacement[0] = -pCell.position[0] / distance_to_origin;
            pCell.displacement[1] = -pCell.position[1] / distance_to_origin;
            pCell.displacement[2] = -pCell.position[2] / distance_to_origin;
            return Math.abs( duct_radius - distance_to_origin );
        }
    }
}