package ru.biosoft.physicell.examples;

import java.util.ArrayList;
import java.util.List;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.standard.O2based;
import ru.biosoft.physicell.core.standard.StandardModels;
import ru.biosoft.physicell.ui.Visualizer2D;
import ru.biosoft.physicell.ui.Visualizer2D.Section;

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
public class TestHDS
{
    static double o2_conc = 38.0; //check this value to make sure it has support from literature
    private static String resultPath = "C:/Users/Damag/BIOFVM/HDS";
    //    private static String resultName = "result";

    private static List<Visualizer2D> resultListeners = new ArrayList<>();
    static
    {
        addListener( "Result_X", Section.X, 0 );
        addListener( "Result_Y", Section.Y, 0 );
        addListener( "Result_Z", Section.Z, 0 );
    }

    private static void addListener(String name, Section sec, int slice)
    {
        Visualizer2D visualizer = Visualizer2D.createWithGIF( resultPath, name, sec, slice );
        visualizer.setSaveImage( false );
        //        visualizer3.setDrawDensity( false );
        visualizer.setMaxDensity( o2_conc * 2 );
        resultListeners.add( visualizer );
    }


    static List<double[]> create_sphere(double cell_radius, double sphere_radius)
    {
        List<double[]> cells = new ArrayList<>();
        int xc = 0, zc = 0;
        double x_spacing = cell_radius * Math.sqrt( 3 );
        double y_spacing = cell_radius * 2;
        double z_spacing = cell_radius * Math.sqrt( 3 );

        //        std::vector<double> tempPoint(3,0.0);
        // std::vector<double> cylinder_center(3,0.0);

        for( double z = -sphere_radius; z < sphere_radius; z += z_spacing, zc++ )
            for( double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++ )
                for( double y = -sphere_radius; y < sphere_radius; y += y_spacing )
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
        return cells;

    }

    public static void main(String[] args) throws Exception
    {
        double t = 0.0;
        double dt = 0.01; // reaction-diffusion time-step
        double mechanics_dt = 0.1;
        double cell_cycle_dt = 6;

        double t_output_interval = 60.0; // 1.0; 
        double t_max = 365 * 24 * 60;
        double t_next_output_time = 0;
//        int next_output_index = 0;

        double dx;
        double dy;
        double dz;

        // figure out the bounding box 
        double[] bounding_box = new double[6];//( 6, 0.0 );
        bounding_box[PhysiCellConstants.mesh_min_x_index] = -1000;
        bounding_box[PhysiCellConstants.mesh_max_x_index] = 1000;
        bounding_box[PhysiCellConstants.mesh_min_y_index] = -1000;
        bounding_box[PhysiCellConstants.mesh_max_y_index] = 1000;
        bounding_box[PhysiCellConstants.mesh_min_z_index] = -1000;
        bounding_box[PhysiCellConstants.mesh_max_z_index] = 1000;
        dx = 20;
        dy = 20;
        dz = 20;


        // create a microenvironment
        Microenvironment microenvironment = new Microenvironment( "substrate scale", "minutes", "microns" );
        Model model = new Model( microenvironment );
        microenvironment.setDensity( 0, "oxygen", "mmHg", 1.0e5, 0.1 );
        microenvironment.resizeSpace( bounding_box[0], bounding_box[3], bounding_box[1], bounding_box[4], bounding_box[2], bounding_box[5],
                dx, dy, dz );

        // Cell_Container 
        double mechanics_voxel_size = 30;
        CellContainer.createCellContainer( microenvironment, mechanics_voxel_size );

        for( int n = 0; n < microenvironment.numberVoxels(); n++ )
            microenvironment.getDensity( n )[0] = o2_conc;

        // register substrates properties 
        //        microenvironment.diffusion_coefficients[0] = 1.0e5; // microns^2 / min 
        //        microenvironment.decay_rates[0] = 0.1;

        CellDefinition cd = StandardModels.createFromDefault( "tumor cell", 0, microenvironment );

        cd.phenotype.cycle = StandardModels.Ki67_advanced;
        // set default_cell_functions; 
        cd.functions.updatePhenotype = new O2based();
        //        cd.phenotype.secretion.sync_to_microenvironment( microenvironment );
        //        cd.phenotype.sync_to_functions( cd.functions );

        int Q_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_negative );
        int oxygen_substrate_index = microenvironment.findDensityIndex( "oxygen" );
        int K1_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_positive_premitotic );
        int K2_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_positive_postmitotic );
        int apoptosis_model_index = cd.phenotype.death.findDeathModelIndex( "apoptosis" );
        int necrosis_model_index = cd.phenotype.death.findDeathModelIndex( "necrosis" );
        // cells apoptose after about 7 days 
        cd.phenotype.death.rates.set( apoptosis_model_index, 1.0 / ( 7.0 * 24.0 * 60.0 ) );
        // initially no necrosis 
        cd.phenotype.death.rates.set( necrosis_model_index, 0.0 );

        // make sure the cells uptake oxygen at the right rate 
        cd.phenotype.secretion.uptakeRates[oxygen_substrate_index] = 10;

        // update transition times 
        cd.phenotype.cycle.data.setTransitionRate( Q_index, K1_index, 1.0 / ( 8.5 * 60.0 ) );
        cd.phenotype.cycle.data.setTransitionRate( K1_index, K2_index, 1.0 / ( 13.0 * 60.0 ) );
        cd.phenotype.cycle.data.setTransitionRate( K2_index, Q_index, 1.0 / ( 2.5 * 60.0 ) );

        // let's make necrotic cells survive 6 hours in minimal oxygen conditions  
        cd.parameters.max_necrosis_rate = 1.0 / ( 6.0 * 60.0 );
        model.registerCellDefinition( cd );

        double cell_radius = 10;
        double sphere_radius = 150;
        // std::cout << __FILE__ << " custom " << __LINE__ << std::endl; 
        //        std::vector<std::vector<double>> cell_positions;
        List<double[]> cell_positions = create_sphere( cell_radius, sphere_radius );

        //add Dirichlet node for all the voxels located outside of the duct
        double[] dirichlet_o2 = new double[] {o2_conc};

        double min_x = microenvironment.mesh.boundingBox[0];
        double max_x = microenvironment.mesh.boundingBox[3];
        double min_y = microenvironment.mesh.boundingBox[1];
        double max_y = microenvironment.mesh.boundingBox[4];
        double min_z = microenvironment.mesh.boundingBox[2];
        double max_z = microenvironment.mesh.boundingBox[5];
        double strip_width = 40;

        for( int i = 0; i < microenvironment.numberVoxels(); i++ )
        {
            if( Math.abs( max_x - microenvironment.voxels( i ).center[0] ) < strip_width
                    || Math.abs( microenvironment.voxels( i ).center[0] - min_x ) < strip_width
                    || Math.abs( max_y - microenvironment.voxels( i ).center[1] ) < strip_width
                    || Math.abs( microenvironment.voxels( i ).center[1] - min_y ) < strip_width
                    || Math.abs( max_z - microenvironment.voxels( i ).center[2] ) < strip_width
                    || Math.abs( microenvironment.voxels( i ).center[2] - min_z ) < strip_width )
            {
                microenvironment.addDirichletNode( i, dirichlet_o2 );
            }
        }

        for( double[] pos : cell_positions )
        {
            if( pos[0] > 0 )
                continue;
            Cell pCell = Cell.createCell( cd, model, pos );
            pCell.phenotype.cycle.data.currentPhaseIndex = Q_index;
            if( pCell.phenotype.cycle.currentPhase().entryFunction != null )
                pCell.phenotype.cycle.currentPhase().entryFunction.execute( pCell, pCell.phenotype, dt );
            // pCell->parameters.necrosis_type= PhysiCellConstants.deterministic_necrosis;
        }

        for( BasicAgent agent : microenvironment.getAgents() )
        {
            agent.setUptakeConstants( dt );
        }

        //        std::cout << (*all_cells).size() <<" agents created successfully." <<std::endl;

        //        std::vector<double> position (3, 0.0);
        //        position[0]=0;
        //        position[1]=0;
        //        position[2]=0;

        //        int output_index =0; 
        //        BioFVM::RUNTIME_TIC();
        //        BioFVM::TIC();

        //        std::ofstream report_file ("report_spheroid.txt");
        //        report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
        //        try 
        //        {       

        for( Visualizer2D listener : resultListeners )
            listener.init();
        while( t < t_max )
        {
            if( Math.abs( t - t_next_output_time ) < 0.0001 )
            {
                //                    log_output(t, output_index, microenvironment, report_file);

                for( Visualizer2D listener : resultListeners )
                    listener.saveResult( microenvironment, t );
                t_next_output_time += t_output_interval;
            }
            microenvironment.simulateSourcesSinks( dt );
            microenvironment.simulateDiffusionDecay( dt );
            ( (CellContainer)microenvironment.agentContainer ).updateAllCells( model, t, cell_cycle_dt, mechanics_dt, dt );
            t += dt;
            //                output_index++;
        }
        for( Visualizer2D listener : resultListeners )
            listener.finish();
        //            log_output(t, output_index, microenvironment, report_file);
        //            report_file.close();
        //        }
        //        catch( const std::exception& e ) { // reference to the base of a polymorphic object
        //            std::cout << e.what(); // information from length_error printed
        //        }
        //        return 0; 
    }

}
