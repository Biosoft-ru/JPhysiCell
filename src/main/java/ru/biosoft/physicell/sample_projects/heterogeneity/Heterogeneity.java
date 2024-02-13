package ru.biosoft.physicell.sample_projects.heterogeneity;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.ui.Visualizer;

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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

public class Heterogeneity
{
    private static final String CUSTOM_ONCOPROTEIN = "custom:oncoprotein";

    public static void init(Model model) throws Exception
    {
        PhysiCellUtilities.setSeed( model.getParameterInt( "random_seed" ) );
        SignalBehavior.setup_signal_behavior_dictionaries( model.getMicroenvironment() );
        createCellTypes( model );
        setupTissue( model );
        printSummary( model.getMicroenvironment(), CUSTOM_ONCOPROTEIN );
        for( Visualizer visualizer : model.getVisualizers() )
        {
            visualizer.setAgentVisualizer( new HeterogeneityVisualizer( model ) );
        }
    }

    static void createCellTypes(Model model)
    {
        CellDefinition pCD = CellDefinition.getCellDefinition( "cancer cell" );
        pCD.functions.updatePhenotype = new TumorPhenotype();
        pCD.parameters.o2_proliferation_saturation = 38;
        pCD.parameters.o2_reference = 38;
    }

    static void setupTissue(Model model) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        double xMin = m.mesh.boundingBox[0];
        double yMin = m.mesh.boundingBox[1];
        double zMin = m.mesh.boundingBox[2];
        double xMax = m.mesh.boundingBox[3];
        double yMax = m.mesh.boundingBox[4];
        double zMax = m.mesh.boundingBox[5];

        if( m.options.simulate2D == true )
        {
            zMin = 0.0;
            zMax = 0.0;
        }

        // create some of each type of cell 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            System.out.println( "Placing cells of type " + cd.name + " ... " );
            for( int n = 0; n < model.getParameterInt( "number_of_cells" ); n++ )
            {
                double[] position = new double[3];
                position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
                position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
                position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
                Cell.createCell( cd, m, position );
            }
        }

        CellDefinition pCD = CellDefinition.getCellDefinition( "cancer cell" );
        double cellRadius = pCD.phenotype.geometry.radius;
        double cellSpacing = 0.95 * 2.0 * cellRadius;
        double tumorRadius = model.getParameterDouble( "tumor_radius" ); // 250.0; 
        double x = 0.0;
        double xOuter = tumorRadius;
        double y = 0.0;
        double pMean = model.getParameterDouble( "oncoprotein_mean" );
        double pSD = model.getParameterDouble( "oncoprotein_sd" );
        double pMin = model.getParameterDouble( "oncoprotein_min" );
        double pMax = model.getParameterDouble( "oncoprotein_max" );

        int n = 0;
        while( y < tumorRadius )
        {
            x = 0.0;
            if( n % 2 == 1 )
            {
                x = 0.5 * cellSpacing;
            }
            xOuter = Math.sqrt( tumorRadius * tumorRadius - y * y );

            while( x < xOuter )
            {
                Cell pCell = Cell.createCell( pCD, m, new double[] {x, y, 0.0} ); // tumor cell 
                double p = PhysiCellUtilities.NormalRandom( pMean, pSD );
                p = PhysiCellUtilities.restrict( p, pMin, pMax );
                SignalBehavior.setSingleBehavior( pCell, CUSTOM_ONCOPROTEIN, p );

                if( Math.abs( y ) > 0.01 )
                {
                    Cell pCell2 = Cell.createCell( pCD, m, new double[] {x, -y, 0.0} ); // tumor cell 
                    p = PhysiCellUtilities.NormalRandom( pMean, pSD );
                    p = PhysiCellUtilities.restrict( p, pMin, pMax );
                    SignalBehavior.setSingleBehavior( pCell2, CUSTOM_ONCOPROTEIN, p );
                }
                if( Math.abs( x ) > 0.01 )
                {
                    Cell pCell3 = Cell.createCell( pCD, m, new double[] { -x, y, 0} ); // tumor cell 
                    p = PhysiCellUtilities.NormalRandom( pMean, pSD );
                    p = PhysiCellUtilities.restrict( p, pMin, pMax );
                    SignalBehavior.setSingleBehavior( pCell3, CUSTOM_ONCOPROTEIN, p );

                    if( Math.abs( y ) > 0.01 )
                    {
                        Cell pCell4 = Cell.createCell( pCD, m, new double[] { -x, -y, 0} ); // tumor cell
                        p = PhysiCellUtilities.NormalRandom( pMean, pSD );
                        p = PhysiCellUtilities.restrict( p, pMin, pMax );
                        SignalBehavior.setSingleBehavior( pCell4, CUSTOM_ONCOPROTEIN, p );
                    }
                }
                x += cellSpacing;
            }
            y += cellSpacing * Math.sqrt( 3.0 ) / 2.0;
            n++;
        }
    }

    static void printSummary(Microenvironment m, String parameter)
    {
        double sum = 0.0;
        double min = 9e9;
        double max = -9e9;
        for( Cell cell : m.getAgents( Cell.class ) )
        {
            double r = SignalBehavior.getSingleSignal( cell, parameter );
            sum += r;
            min = Math.min( r, min );
            max = Math.max( r, max );
        }
        int size = m.getAgentsCount();
        double mean = sum / ( size + 1e-15 );
        // compute standard deviation 
        sum = 0.0;
        for( Cell cell : m.getAgents( Cell.class ) )
        {
            double r = SignalBehavior.getSingleSignal( cell, parameter );
            sum += ( r - mean ) * ( r - mean );
        }
        double sd = Math.sqrt( sum / ( size - 1.0 + 1e-15 ) );

        System.out.println( "Oncoprotein summary: " );
        System.out.println( "===================" );
        System.out.println( "mean: " + mean );
        System.out.println( "standard deviation: " + sd );
        System.out.println( "[min max]: [" + min + " " + max + "]" );
    }
}