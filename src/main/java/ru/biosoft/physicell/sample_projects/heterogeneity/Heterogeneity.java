package ru.biosoft.physicell.sample_projects.heterogeneity;

import java.awt.Color;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.update_phenotype;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.ui.AgentVisualizer;
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
    public static void init(Model model) throws Exception
    {
        SignalBehavior.setup_signal_behavior_dictionaries( model.getMicroenvironment() );
        create_cell_types( model );
        setup_tissue( model );
        for( Visualizer visualizer : model.getVisualizers() )
        {
            visualizer.setAgentVisualizer( new HeterogeneityVisualizer( model ) );
        }
    }

    static void create_cell_types(Model model)
    {
        PhysiCellUtilities.setSeed( model.getParameterInt( "random_seed" ) );
        CellDefinition pCD = CellDefinition.getCellDefinition( "cancer cell" );
        pCD.functions.updatePhenotype = new tumor_cell_phenotype_with_oncoprotein();
        pCD.parameters.o2_proliferation_saturation = 38;
        pCD.parameters.o2_reference = 38;
    }

    static void setup_tissue(Model model) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Zmin = m.mesh.boundingBox[2];
        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];
        double Zmax = m.mesh.boundingBox[5];

        if( m.options.simulate_2D == true )
        {
            Zmin = 0.0;
            Zmax = 0.0;
        }

        // create some of each type of cell 
        for( CellDefinition cd : CellDefinition.getCellDefinitions() )
        {
            System.out.println( "Placing cells of type " + cd.name + " ... " );
            for( int n = 0; n < model.getParameterInt( "number_of_cells" ); n++ )
            {
                double[] position = new double[3];
                position[0] = PhysiCellUtilities.UniformRandom( Xmin, Xmax );
                position[1] = PhysiCellUtilities.UniformRandom( Ymin, Ymax );
                position[2] = PhysiCellUtilities.UniformRandom( Zmin, Zmax );
                Cell.createCell( cd, m, position );
            }
        }

        CellDefinition pCD = CellDefinition.getCellDefinition( "cancer cell" );
        double cell_radius = pCD.phenotype.geometry.radius;
        double cell_spacing = 0.95 * 2.0 * cell_radius;
        double tumor_radius = model.getParameterDouble( "tumor_radius" ); // 250.0; 
        double x = 0.0;
        double x_outer = tumor_radius;
        double y = 0.0;
        double p_mean = model.getParameterDouble( "oncoprotein_mean" );
        double p_sd = model.getParameterDouble( "oncoprotein_sd" );
        double p_min = model.getParameterDouble( "oncoprotein_min" );
        double p_max = model.getParameterDouble( "oncoprotein_max" );

        int n = 0;
        //        Cell pCell = Cell.createCell( pCD, m, new double[] {x, y, 0.0} ); // tumor cell 
        //        double p = PhysiCellUtilities.NormalRandom( p_mean, p_sd );
        //        p = PhysiCellUtilities.restrict( p, p_min, p_max );
        //        SignalBehavior.setSingleBehavior( pCell, "custom:oncoprotein", p );
        while( y < tumor_radius )
        {
            x = 0.0;
            if( n % 2 == 1 )
            {
                x = 0.5 * cell_spacing;
            }
            x_outer = Math.sqrt( tumor_radius * tumor_radius - y * y );

            while( x < x_outer )
            {
                Cell pCell = Cell.createCell( pCD, m, new double[] {x, y, 0.0} ); // tumor cell 
                double p = PhysiCellUtilities.NormalRandom( p_mean, p_sd );
                p = PhysiCellUtilities.restrict( p, p_min, p_max );
                SignalBehavior.setSingleBehavior( pCell, "custom:oncoprotein", p );

                if( Math.abs( y ) > 0.01 )
                {
                    pCell = Cell.createCell( pCD, m, new double[] {x, -y, 0.0} ); // tumor cell 
                    p = PhysiCellUtilities.NormalRandom( p_mean, p_sd );
                    p = PhysiCellUtilities.restrict( p, p_min, p_max );
                    SignalBehavior.setSingleBehavior( pCell, "custom:oncoprotein", p );
                }
                if( Math.abs( x ) > 0.01 )
                {
                    pCell = Cell.createCell( pCD, m, new double[] { -x, y, 0} ); // tumor cell 
                    p = PhysiCellUtilities.NormalRandom( p_mean, p_sd );
                    p = PhysiCellUtilities.restrict( p, p_min, p_max );
                    SignalBehavior.setSingleBehavior( pCell, "custom:oncoprotein", p );

                    if( Math.abs( y ) > 0.01 )
                    {
                        pCell = Cell.createCell( pCD, m, new double[] { -x, -y, 0} ); // tumor cell
                        p = PhysiCellUtilities.NormalRandom( p_mean, p_sd );
                        p = PhysiCellUtilities.restrict( p, p_min, p_max );
                        SignalBehavior.setSingleBehavior( pCell, "custom:oncoprotein", p );
                    }
                }
                x += cell_spacing;
            }
            y += cell_spacing * Math.sqrt( 3.0 ) / 2.0;
            n++;
        }

        double sum = 0.0;
        double min = 9e9;
        double max = -9e9;
        for( Cell cell : m.getAgents( Cell.class ) )
        {
            double r = SignalBehavior.get_single_signal( cell, "custom:oncoprotein" );
            sum += r;
            r = PhysiCellUtilities.restrict( r, min, max ); //TODO: check
        }
        int size = m.getAgentsCount();
        double mean = sum / ( size + 1e-15 );
        // compute standard deviation 
        sum = 0.0;
        for( Cell cell : m.getAgents( Cell.class ) )
        {
            double r = SignalBehavior.get_single_signal( cell, "custom:oncoprotein" );
            sum += ( r - mean ) * ( r - mean );
        }
        double standard_deviation = Math.sqrt( sum / ( size - 1.0 + 1e-15 ) );

        System.out.println( "Oncoprotein summary: " );
        System.out.println( "===================" );
        System.out.println( "mean: " + mean );
        System.out.println( "standard deviation: " + standard_deviation );
        System.out.println( "[min max]: [" + min + " " + max + "]" );
    }

    public static class HeterogeneityVisualizer extends AgentVisualizer
    {
        double p_min;
        double p_max;
        public HeterogeneityVisualizer(Model model)
        {
            p_min = model.getParameterDouble( "oncoprotein_min" );
            p_max = model.getParameterDouble( "oncoprotein_max" );
        }

        @Override
        public Color findColor(Cell pCell)
        {
            double p = SignalBehavior.get_single_signal( pCell, "custom:oncoprotein" );

            // immune are black
            Color output = Color.black;

            if( pCell.type == 1 )
            {
                return output;
            }

            // live cells are green, but shaded by oncoprotein value 
            if( pCell.phenotype.death.dead == false )
            {
                int oncoprotein = (int)Math.round( ( 1.0 / ( p_max - p_min ) ) * ( p - p_min ) * 255.0 );
                output = new Color( (int) ( oncoprotein / p_max ), (int) ( oncoprotein / p_max ), (int) ( ( 255 - oncoprotein ) / p_max ) );
                //		char szTempString [128];
                //		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
                //		output[0].assign( szTempString );
                //		output[1].assign( szTempString );
                //
                //		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/p_max) , (int)round(output[0][1]/p_max) , (int)round(output[0][2]/p_max) );
                //		output[2].assign( szTempString );
                return output;
            }

            // if not, dead colors 
            if( SignalBehavior.get_single_signal( pCell, "apoptotic" ) > 0.5 )
            {
                output = new Color( 125, 0, 0 );
                //		output[0] = "rgb(255,0,0)";
                //		output[2] = "rgb(125,0,0)";
            }

            // Necrotic - Brown
            if( SignalBehavior.get_single_signal( pCell, "necrotic" ) > 0.5 )
            {
                output = new Color( 139, 69, 19 );
                //		output[0] = "rgb(250,138,38)";
                //		output[2] = "rgb(139,69,19)";
            }
            return output;
        }
    }

    public static class tumor_cell_phenotype_with_oncoprotein extends update_phenotype
    {
        @Override
        public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
        {
            new StandardModels.update_cell_and_death_parameters_O2_based().execute( pCell, phenotype, dt );

            // if cell is dead, don't bother with future phenotype changes. 
            if( SignalBehavior.get_single_signal( pCell, "dead" ) > 0.5 )
            {
                pCell.functions.updatePhenotype = null;
                return;
            }
            // multiply proliferation rate by the oncoprotein 
            double cycle_rate = SignalBehavior.get_single_behavior( pCell, "cycle entry" );
            //            System.out.println( "BEFORE " + SignalBehavior.get_single_behavior( pCell, "cycle entry" ) );
            cycle_rate *= SignalBehavior.get_single_signal( pCell, "custom:oncoprotein" );
            SignalBehavior.setSingleBehavior( pCell, "cycle entry", cycle_rate );
            //            System.out.println( "AFTER " + SignalBehavior.get_single_behavior( pCell, "cycle entry" ) );
        }
    }
}