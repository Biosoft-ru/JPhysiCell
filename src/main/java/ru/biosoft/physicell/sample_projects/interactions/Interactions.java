package ru.biosoft.physicell.sample_projects.interactions;
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

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;

public class Interactions
{

    public static void init(Model model)
    {
        createCellTypes( model);
        setupTissue( model );
    }

    static void createCellTypes(Model model)
    {
        // set the random seed 
        PhysiCellUtilities.setSeed( model.getParameterInt( "random_seed" ) );
        // set up bacteria 
        CellDefinition pCD = CellDefinition.getCellDefinition( "bacteria" );
        pCD.functions.updatePhenotype = new BacteriaPhenotype();
        // pCD.functions.update_migration_bias = advanced_chemotaxis_function; 
        // pCD.phenotype.motility.chemotactic_sensitivity( "resource" ) = 1; 
        // pCD.phenotype.motility.chemotactic_sensitivity( "quorum" ) = 0.1; 

        // set up blood vessels 
        pCD = CellDefinition.getCellDefinition( "blood vessel" );
        pCD.functions.updatePhenotype = null;
        pCD.isMovable = false;

        // set up stem cells 
        pCD = CellDefinition.getCellDefinition( "stem" );
        pCD.functions.updatePhenotype = new StemPhenotype( model );
        // pCD.phenotype.cell_transformations.transformation_rate("differentiated") = 0.0001; 

        // set up differentiated cells 
        pCD = CellDefinition.getCellDefinition( "differentiated" );
        pCD.functions.updatePhenotype = new DifferentiatedPhenotype();

        // set up macrophages 
        pCD = CellDefinition.getCellDefinition( "macrophage" );
        // pCD.phenotype.cell_interactions.dead_phagocytosis_rate = 0.05; 
        pCD.functions.updatePhenotype = new MacrophagePhenotype();
        // pCD.functions.update_migration_bias = advanced_chemotaxis_function; 
        // pCD.phenotype.motility.chemotactic_sensitivity( "debris" ) = 0.1; 
        // pCD.phenotype.motility.chemotactic_sensitivity( "quorum" ) = 1; 


        // set up CD8+ T cells 
        pCD = CellDefinition.getCellDefinition( "CD8+ T cell" );
        pCD.functions.updatePhenotype = new CD8TcellPhenotype();
        // pCD.phenotype.cell_interactions.attack_rate("bacteria") = 0.05; 

        // set up neutrophil  
        pCD = CellDefinition.getCellDefinition( "neutrophil" );
        pCD.functions.updatePhenotype = new NeutrophilPhenotype();
        // pCD.phenotype.cell_interactions.live_phagocytosis_rate("bacteria") = 0.05; 

    }

    static void setupTissue(Model model)
    {
        Microenvironment microenvironment = model.getMicroenvironment();

        double xMin = microenvironment.mesh.boundingBox[0];
        double yMin = microenvironment.mesh.boundingBox[1];
        double zMin = microenvironment.mesh.boundingBox[2];

        double xMax = microenvironment.mesh.boundingBox[3];
        double yMax = microenvironment.mesh.boundingBox[4];
        double zMax = microenvironment.mesh.boundingBox[5];

        if( microenvironment.options.simulate2D == true )
        {
            zMin = 0.0;
            zMax = 0.0;
        }

        // create some of each type of cell 
        for( CellDefinition pCD : CellDefinition.getCellDefinitions() )
        {
            int num_cells = model.getParameterInt( "number_of_cells" );
            if( num_cells > 0 )
            {
                System.out.println( "Placing cells of type " + pCD.name + " ... " );
                for( int n = 0; n < num_cells; n++ )
                {
                    double[] position = {0, 0, 0};
                    position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
                    position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
                    position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
                    Cell.createCell( pCD, microenvironment, position );
                }
            }
        }

        // parameter-based placement 
        // bacteria 
        CellDefinition pCD = CellDefinition.getCellDefinition( "bacteria" );
        System.out.println( "Placing cells of type " + pCD.name + " ... " );
        ;
        for( int n = 0; n < model.getParameterInt( "number_of_bacteria" ); n++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
            position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
            position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
            Cell.createCell( pCD, microenvironment, position );
        }

        // blood vessels 
        pCD = CellDefinition.getCellDefinition( "blood vessel" );
        System.out.println( "Placing cells of type " + pCD.name + " ... " );
        ;
        for( int n = 0; n < model.getParameterInt( "number_of_blood_vessels" ); n++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
            position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
            position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
            Cell.createCell( pCD, microenvironment, position );
        }

        // stem cells 
        pCD = CellDefinition.getCellDefinition( "stem" );
        System.out.println( "Placing cells of type " + pCD.name + " ... " );
        ;
        for( int n = 0; n < model.getParameterInt( "number_of_stem_cells" ); n++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
            position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
            position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
            Cell.createCell( pCD, microenvironment, position );
        }

        // differentiated cells 
        pCD = CellDefinition.getCellDefinition( "differentiated" );
        System.out.println( "Placing cells of type " + pCD.name + " ... " );
        for( int n = 0; n < model.getParameterInt( "number_of_differentiated_cells" ); n++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
            position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
            position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
            Cell.createCell( pCD, microenvironment, position );
        }

        // macrophages 
        pCD = CellDefinition.getCellDefinition( "macrophage" );
        System.out.println( "Placing cells of type " + pCD.name + " ... " );
        ;
        for( int n = 0; n < model.getParameterInt( "number_of_macrophages" ); n++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
            position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
            position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
            Cell.createCell( pCD, microenvironment, position );
        }

        // neutrophils  
        pCD = CellDefinition.getCellDefinition( "neutrophil" );
        System.out.println( "Placing cells of type " + pCD.name + " ... " );
        for( int n = 0; n < model.getParameterInt( "number_of_neutrophils" ); n++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
            position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
            position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
            Cell.createCell( pCD, microenvironment, position );
        }

        // CD8+ T cells   
        pCD = CellDefinition.getCellDefinition( "CD8+ T cell" );
        System.out.println( "Placing cells of type " + pCD.name + " ... " );
        ;
        for( int n = 0; n < model.getParameterInt( "number_of_CD8T_cells" ); n++ )
        {
            double[] position = {0, 0, 0};
            position[0] = PhysiCellUtilities.UniformRandom( xMin, xMax );
            position[1] = PhysiCellUtilities.UniformRandom( yMin, yMax );
            position[2] = PhysiCellUtilities.UniformRandom( zMin, zMax );
            Cell.createCell( pCD, microenvironment, position );
        }
    }
}