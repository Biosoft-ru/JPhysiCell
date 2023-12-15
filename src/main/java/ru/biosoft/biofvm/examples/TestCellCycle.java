package ru.biosoft.biofvm.examples;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.text.DecimalFormat;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;

import ru.biosoft.biofvm.BasicAgent;
import ru.biosoft.biofvm.ConstantCoefficientsLOD3D;
import ru.biosoft.biofvm.Microenvironment;
import ru.biosoft.biofvm.VectorUtil;
import ru.biosoft.biofvm.Visualizer;
import ru.biosoft.biofvm.cell.Cell;
import ru.biosoft.biofvm.cell.CellContainer;
import ru.biosoft.biofvm.cell.CellDefinition;
import ru.biosoft.biofvm.cell.Phenotype;
import ru.biosoft.biofvm.cell.PhysiCellConstants;
import ru.biosoft.biofvm.cell.StandardModels;

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
public class TestCellCycle
{
    private static double o2Сonc = 6.06;
    private static DecimalFormat format = new DecimalFormat( "##.##" );
    private static String resultPath = "C:/Users/Damag/BIOFVM/CellCycle";

    //visualizer settings
    private static Visualizer visualizer = new Visualizer();
    static
    {
        visualizer.setDrawDensity( false );
        visualizer.setColorPhase( "Ki67-", Color.gray );
        visualizer.setColorPhase( "Ki67+ (premitotic)", Color.blue );
        visualizer.setColorPhase( "Ki67+ (postmitotic)", Color.cyan );
        visualizer.setColorPhase( "Apoptotic", Color.red );
    }

    //Cell types
    private static int num_ki67_q = 0;//100;//20;
    private static int num_ki67_positive_pre = 1;//100;//20;
    private static int num_ki67_positive_post = 0;//100;//20;
    private static int num_apoptotic = 0;//100;//20;

    private static int size = 1000;
    private static int zSlice = 100;

    //Time options
    private static double tMax = 60 * 24 * 6;
    private static double dt = 1;
    private static double outputInterval = 60;

    private static boolean outputReport = true;
    private static boolean outputImage = true;
    private static boolean outputGIF = true;

    public static String format(double val)
    {
        return format.format( val );
    }

    public static void main(String ... argv) throws Exception
    {
        double sideLength = size;// 2000;
        double dx = 20;
        double dy = 20;
        double dz = 20;

        Microenvironment microenvironment = new Microenvironment();
        microenvironment.name = "substrate scale";
        microenvironment.setDensity( 0, "oxygen", "mmHg" );
        microenvironment.resizeSpace( 0, sideLength, 0, sideLength, 0, sideLength, dx, dy, dz );
        microenvironment.spatialUnits = "microns";
        microenvironment.timeUnits = "minutes";
        microenvironment.mesh.units = "microns";

        // Cell_Container 
        double mechanicsVoxelSize = 30;
        CellContainer.create_cell_container_for_microenvironment( microenvironment, mechanicsVoxelSize );
        for( int n = 0; n < microenvironment.number_of_voxels(); n++ )
            microenvironment.density_vector( n )[0] = o2Сonc;
        microenvironment.setSolver( new ConstantCoefficientsLOD3D() );

        StandardModels.initialize_default_cell_definition();
        CellDefinition cellDefaults = StandardModels.cellDefaults;
        cellDefaults.type = 0;
        cellDefaults.name = "tumor cell";
        cellDefaults.functions.cycle_model = StandardModels.Ki67_advanced; // set default cell cycle model 
        cellDefaults.functions.updatePhenotype = new StandardModels.update_cell_and_death_parameters_O2_based(); // set default_cell_functions; 
        cellDefaults.functions.updateVelocity = (pCell, phenotype, dt) -> { // disable cell's movement
            return;
        };
        Phenotype defaultPhenotype = cellDefaults.phenotype;
        defaultPhenotype.secretion.sync_to_microenvironment( microenvironment );
        cellDefaults.phenotype.sync_to_functions( cellDefaults.functions );
        // first find index for a few key variables. 
        int apoptosisModelIndex = defaultPhenotype.death.find_death_model_index( "Apoptosis" );
        int necrosisModelIndex = defaultPhenotype.death.find_death_model_index( "Necrosis" );
        int oxygenSubstrateIndex = microenvironment.find_density_index( "oxygen" );

        int K1_index = cellDefaults.functions.cycle_model.find_phase_index( PhysiCellConstants.Ki67_positive_premitotic );
        int K2_index = cellDefaults.functions.cycle_model.find_phase_index( PhysiCellConstants.Ki67_positive_postmitotic );
        int Q_index = cellDefaults.functions.cycle_model.find_phase_index( PhysiCellConstants.Ki67_negative );
        int A_index = cellDefaults.functions.cycle_model.find_phase_index( PhysiCellConstants.apoptotic );
        //        int N_index = cell_defaults.functions.cycle_model.find_phase_index( PhysiCellConstants.necrotic_swelling );

        // cells apoptose after about 7 days 
        defaultPhenotype.death.rates.set( apoptosisModelIndex, 1.0 / ( 7.0 * 24.0 * 60.0 ) );

        // initially no necrosis 
        defaultPhenotype.death.rates.set( necrosisModelIndex, 0.0 );

        // make sure the cells uptake oxygen at the right rate 
        defaultPhenotype.secretion.uptake_rates[oxygenSubstrateIndex] = 0;

        // cells leave the Q phase and enter the K1 phase after 5 hours 
        defaultPhenotype.cycle.data.setTransitionRate( Q_index, K1_index, 1.0 / ( 5.0 * 60.0 ) );

        // let's make necrotic cells survive 6 hours in minimal oxygen conditions  
        cellDefaults.parameters.max_necrosis_rate = 1.0 / ( 6.0 * 60.0 );

        int total = num_ki67_positive_pre + num_ki67_positive_post + num_ki67_q + num_apoptotic;
        //        double T1 = 13 * 60;
        //        double T2 = 2.5 * 60;
        //        double TQ = 74.35 * 60;
        //        double TA = 8.6 * 60;
        int phaseIndex;
        for( int i = 0; i < total; i++ )
        {
            double[] tempPosition = VectorUtil.random( 3, 0, size );
            tempPosition[2] = zSlice;//keep z at slice
            Cell cell = Cell.createCell();
            cell.registerMicroenvironment( microenvironment );
            cell.tag = "Source";
            cell.assignPosition( tempPosition );
            if( i < num_ki67_positive_pre )
            {
                phaseIndex = K1_index;
                //                T = T1;
            }
            else if( i < num_ki67_positive_pre + num_ki67_positive_post )
            {
                phaseIndex = K2_index;
                //                T = T2;
            }
            else if( i < num_ki67_positive_pre + num_ki67_positive_post + num_apoptotic )
            {
                phaseIndex = A_index;
                //                T = TA;
                cell.phenotype.death.trigger_death( apoptosisModelIndex );
                cell.phenotype.cycle.sync_to_cycle_model( cell.phenotype.death.current_model() );
            }
            else
            {
                phaseIndex = Q_index;
                //                T = T/Q;
            }
            cell.phenotype.cycle.data.current_phase_index = phaseIndex;
            if( cell.phenotype.cycle.current_phase().entryFunction != null )
                cell.phenotype.cycle.current_phase().entryFunction.execute( cell, cell.phenotype, dt );
        }

        for( BasicAgent agent : BasicAgent.allBasicAgents )
            agent.setUptakeConstants( dt );

        ImageWriter writer = null;
        ImageOutputStream ios = null;
        if( outputImage && outputGIF )
        {
            writer = ImageIO.getImageWritersByFormatName( "GIF" ).next();
            ios = ImageIO.createImageOutputStream( new File( resultPath + "/result.gif" ) );
            writer.setOutput( ios );
            writer.prepareWriteSequence( null );
        }

        //Simulation starts
        double t = 0.0;
        double tOutputInterval = Math.max( outputInterval, dt );
        double tNextOutputTime = 0;
        while( t < tMax )
        {
            if( Math.abs( t - tNextOutputTime ) < 0.0001 )
            {
                if( outputReport )
                    writTestReport( t );

                if( outputImage )
                {
                    BufferedImage image = visualizer.draw( microenvironment, zSlice, (int)t, resultPath + "/figure_" + (int)t + ".png" );
                    if( outputGIF )
                        writer.writeToSequence( new IIOImage( image, null, null ), writer.getDefaultWriteParam() );
                }
                tNextOutputTime += tOutputInterval;
                System.out.println( t + ", cell count " + BasicAgent.allBasicAgents.size() );
            }
            ( (CellContainer)microenvironment.agent_container ).update_all_cells( t, dt, dt, dt );
            t += dt;
        }

        if( outputImage && outputGIF )
        {
            writer.endWriteSequence();
            writer.reset();
            writer.dispose();
            ios.flush();
            ios.close();
        }
    }

    public static void writTestReport(double timepoint)
    {
        String filename = resultPath + "/cells_" + (int)timepoint + ".txt";
        File f = new File( filename );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( f ) ))
        {
            bw.write( "\tID\tx\ty\tz\tradius\tphenotype\telapsed_time\n" );
            for( int i = 0; i < Cell.allBasicAgents.size(); i++ )
            {
                Cell cell = (Cell)Cell.allBasicAgents.get( i );
                String phenotypeCode = cell.phenotype.cycle.current_phase().name;
                bw.write( i + "\t" + cell.ID + "\t" + format( cell.position[0] ) + "\t" + format( cell.position[1] ) + "\t"
                        + format( cell.position[1] ) + "\t" + format( cell.position[2] ) + "\t" + format( cell.phenotype.geometry.radius )
                        + "\t" + phenotypeCode + "\t" + cell.phenotype.cycle.data.elapsed_time_in_phase + "\n" );
            }
        }
        catch( Exception ex )
        {

        }
    }
}