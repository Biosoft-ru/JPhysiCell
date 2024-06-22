package ru.biosoft.physicell.examples;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.text.DecimalFormat;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.standard.O2based;
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
public class TestVolume
{
    private static DecimalFormat format = new DecimalFormat( "##.##" );
    static double outputInterval = 60.0;
    static double dt = 5;
    static double o2Conc = 5.01;

    static String APOPTOSIS = "Apoptosis";
    static String NECROSIS = "Necrosis";

    private static String resultPath = "C:/Users/Damag/BIOFVM/Volume/";

    public static void main(String[] argv) throws Exception
    {
        run( 2402, APOPTOSIS, resultPath + "/Apoptosis.txt" );
        run( 2402, NECROSIS, resultPath + "/Necrosis.txt" );
        run( 4000, "", resultPath + "/Default.txt" );
    }

    public static void run(double tMax, String type, String name) throws Exception
    {
        Microenvironment m = new Microenvironment( "substrate scale", 2000, 20, "minutes", "microns" );
        Model model = new Model( m );
        CellContainer.createCellContainer( m, 30 );

        m.setDensity( 0, "oxygen", "mmHg", 0, 0 );
        for( int n = 0; n < m.numberVoxels(); n++ )
            m.getDensity( n )[0] = o2Conc;

        model.clearCellDefinitions();
        CellDefinition cd = StandardModels.createFromDefault( "tumor cell", 0, m );
        model.registerCellDefinition( cd );
        cd.phenotype.cycle = StandardModels.createAdvancedKi67();
        cd.functions.updatePhenotype = new O2based();
        //cell_defaults.functions.volume_update_function = standard_volume_update_function;

        // first find index for a few key variables. 
        int apoptosisModelIndex = cd.phenotype.death.findDeathModelIndex( "Apoptosis" );
        int necrosisModelIndex = cd.phenotype.death.findDeathModelIndex( "Necrosis" );

        int K1_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_positive_premitotic );
        int K2_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_positive_postmitotic );
        int Q_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.Ki67_negative );
        //        int A_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.apoptotic );
        //        int N_index = StandardModels.Ki67_advanced.findPhaseIndex( PhysiCellConstants.necrotic_swelling );

        Cell cell = Cell.createCell( cd, model, new double[] {500, 500, 500} );
        if( type.equals( APOPTOSIS ) )
        {
            //            cell.phenotype.cycle.data.currentPhaseIndex = A_index;
            cell.phenotype.death.triggerDeath( apoptosisModelIndex );
            cell.phenotype.cycle = cell.phenotype.death.currentModel();
        }
        else if( type.equals( NECROSIS ) )
        {
            //            cell.phenotype.cycle.data.currentPhaseIndex = N_index;
            cell.phenotype.death.triggerDeath( necrosisModelIndex );
            cell.phenotype.cycle = cell.phenotype.death.currentModel();
        }
        else
        {
            cell.phenotype.cycle.data.currentPhaseIndex = K1_index;
            cell.phenotype.death.rates.set( apoptosisModelIndex, 0.0 ); // disable apoptosis
            cell.phenotype.cycle.data.setTransitionRate( Q_index, K1_index, 1e9 ); // set Q duration to a large value
        }
        cell.phenotype.cycle.currentPhase().entryFunction.execute( cell, cell.phenotype, dt );
        for( BasicAgent agent : m.getAgents() )
            agent.setUptakeConstants( dt );

        //Simulation starts
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( new File( name ) ) ))
        {
            bw.append( "Total\tFluid\tNuclear\tCytoplasmatic\n" );
            double t = 0.0;
            double nextOutputTime = 0;
            while( t < tMax )
            {
                if( Math.abs( t - nextOutputTime ) < 0.001 )
                {
                    String report = cell.ID + "\t" + cell.phenotype.cycle.currentPhase().name + "\t"
                            + format.format( cell.get_total_volume() ) + "\t" + format.format( cell.phenotype.volume.fluid ) + "\t"
                            + format.format( cell.phenotype.volume.nuclear_solid ) + "\t"
                            + format.format( cell.phenotype.volume.cytoplasmic_solid ) + "\n";
                    //                    String report = cell.ID + "\t" + m.get( cell.currentVoxelIndex )[0];
                    bw.append( report );
                    nextOutputTime += outputInterval;
                }
                if( m.getAgentsCount() > 1 )
                {
                    Cell toDelete = null;
                    for( BasicAgent agent : m.getAgents() )
                    {
                        if( agent.ID == cell.ID )
                            continue;
                        else
                            toDelete = (Cell)agent;
                    }
                    Cell.deleteCell( toDelete );
                    bw.append( "New cell deleted \n" );
                }
                ( (CellContainer)m.agentContainer ).updateAllCells( m, t, dt, dt, dt );
                t += dt;
            }
        }
        m.getAgents().clear();
    }
}