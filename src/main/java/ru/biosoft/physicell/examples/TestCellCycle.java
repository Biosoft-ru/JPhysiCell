package ru.biosoft.physicell.examples;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.text.DecimalFormat;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.standard.O2based;
import ru.biosoft.physicell.core.standard.StandardModels;
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
public class TestCellCycle
{
    private static double o2Сonc = 6.06;
    private static DecimalFormat format = new DecimalFormat( "##.##" );
    private static String resultPath = "C:/Users/Damag/BIOFVM/CellCycle";
    private static String resultName = "CellCycle4";
    private static int zSlice = 100;
    //visualizer settings
    private static Visualizer visualizer = Visualizer.createWithGIF( resultPath, resultName, Section.Z, zSlice );
    static
    {
        visualizer.setDrawDensity( false );
        visualizer.setColorPhase( "Ki67-", Color.gray );
        visualizer.setColorPhase( "Ki67+ (premitotic)", Color.blue );
        visualizer.setColorPhase( "Ki67+ (postmitotic)", new Color( 0, 180, 0 ) );
        visualizer.setColorPhase( "Apoptotic", Color.red );
    }

    //Cell types
    private static int num_ki67_q = 100;//100;//20;
    private static int num_ki67_positive_pre = 100;//100;//100;//20;
    private static int num_ki67_positive_post = 100;//100;//100;//20;
    private static int num_apoptotic = 100;//100;//20;

    private static int size = 1000;
    private static int cellSize = 20;

    //Time options
    private static double tMax = 60 * 24 * 6;
    private static double dt = 1;
    private static double outputInterval = 60;

    private static boolean outputReport = true;

    public static String format(double val)
    {
        return format.format( val );
    }

    public static void main(String ... argv) throws Exception
    {
        Microenvironment m = new Microenvironment( "substrate scale", size, cellSize, "minutes", "microns" );
        Model model = new Model( m );
        m.setDensity( 0, "oxygen", "mmHg", 0, 0 );
        for( int n = 0; n < m.numberVoxels(); n++ )
            m.getDensity( n )[0] = o2Сonc;

        CellContainer.createCellContainer( m, 30 );

        CellDefinition cd = StandardModels.createFromDefault( "tumor cell", 0, m );
        model.registerCellDefinition( cd );
        cd.phenotype.cycle = StandardModels.Ki67_advanced; // set default cell cycle model 
        cd.functions.updatePhenotype = new O2based(); // set default_cell_functions; 
        cd.functions.updateVelocity = null;
        Phenotype defaultPhenotype = cd.phenotype;

        // first find index for a few key variables. 
        int apoptosisModelIndex = defaultPhenotype.death.findDeathModelIndex( "Apoptosis" );
        int necrosisModelIndex = defaultPhenotype.death.findDeathModelIndex( "Necrosis" );
        int oxygenSubstrateIndex = m.findDensityIndex( "oxygen" );

        int K1_index = cd.phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_positive_premitotic );
        int K2_index = cd.phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_positive_postmitotic );
        int Q_index = cd.phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_negative );
        //        int A_index = cd.phenotype.cycle.findPhaseIndex( PhysiCellConstants.apoptotic );
        //        int N_index = cell_defaults.functions.cycle_model.find_phase_index( PhysiCellConstants.necrotic_swelling );

        // cells apoptose after about 7 days 
        defaultPhenotype.death.rates.set( apoptosisModelIndex, 1.0 / ( 7.0 * 24.0 * 60.0 ) );

        // initially no necrosis 
        defaultPhenotype.death.rates.set( necrosisModelIndex, 0.0 );

        // make sure the cells uptake oxygen at the right rate 
        defaultPhenotype.secretion.uptakeRates[oxygenSubstrateIndex] = 0;

        // cells leave the Q phase and enter the K1 phase after 5 hours 
        defaultPhenotype.cycle.data.setTransitionRate( Q_index, K1_index, 1.0 / ( 5.0 * 60.0 ) );

        // let's make necrotic cells survive 6 hours in minimal oxygen conditions  
        cd.parameters.max_necrosis_rate = 1.0 / ( 6.0 * 60.0 );

        int total = num_ki67_positive_pre + num_ki67_positive_post + num_ki67_q + num_apoptotic;
        //        double T1 = 13 * 60;
        //        double T2 = 2.5 * 60;
        //        double TQ = 74.35 * 60;
        //        double TA = 8.6 * 60;

        for( int i = 0; i < total; i++ )
        {
            int phaseIndex = 0;
            double[] tempPosition = VectorUtil.random( model.getRNG(), 3, 0, size );
            tempPosition[2] = zSlice;//keep z at slice
            Cell cell = Cell.createCell( cd, model, tempPosition );
            cell.tag = "Source";
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
                //                phaseIndex = A_index;
                //                T = TA;
                cell.phenotype.death.triggerDeath( apoptosisModelIndex );
                cell.phenotype.cycle = cell.phenotype.death.currentModel();
            }
            else
            {
                phaseIndex = Q_index;
                //                T = T/Q;
            }
            cell.phenotype.cycle.data.currentPhaseIndex = phaseIndex;
            if( cell.phenotype.cycle.currentPhase().entryFunction != null )
                cell.phenotype.cycle.currentPhase().entryFunction.execute( cell, cell.phenotype, dt );
        }

        for( BasicAgent agent : m.getAgents() )
            agent.setUptakeConstants( dt );

        visualizer.init();

        //Simulation starts
        double t = 0.0;
        double tOutputInterval = Math.max( outputInterval, dt );
        double tNextOutputTime = 0;
        while( t < tMax )
        {
            if( Math.abs( t - tNextOutputTime ) < 0.0001 )
            {
                if( outputReport )
                    writTestReport( m, t );

                visualizer.saveResult( m, t );
                tNextOutputTime += tOutputInterval;
                System.out.println( t + ", cell count " + m.getAgentsCount() );
            }
            ( (CellContainer)m.agentContainer ).updateAllCells( m, t, dt, dt, dt );
            t += dt;
        }
        visualizer.finish();
    }

    public static void writTestReport(Microenvironment m, double timepoint)
    {
        String filename = resultPath + "/cells_" + (int)timepoint + ".txt";
        File f = new File( filename );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( f ) ))
        {
            bw.write( "\tID\tx\ty\tz\tradius\tphenotype\telapsed_time\n" );
            for( BasicAgent agent : m.getAgents() )
            {
                Cell cell = (Cell)agent;
                String phenotypeCode = cell.phenotype.cycle.currentPhase().name;
                bw.write( cell.ID + "\t" + format( cell.position[0] ) + "\t" + format( cell.position[1] ) + "\t"
                        + format( cell.position[1] ) + "\t" + format( cell.position[2] ) + "\t" + format( cell.phenotype.geometry.radius )
                        + "\t" + phenotypeCode + "\t" + cell.phenotype.cycle.data.elapsedTimePhase + "\n" );
            }
        }
        catch( Exception ex )
        {

        }
    }
}