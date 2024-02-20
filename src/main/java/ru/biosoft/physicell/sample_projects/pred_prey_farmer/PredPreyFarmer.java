package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import java.util.ArrayList;
import java.util.List;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;

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
public class PredPreyFarmer
{
    public static void init(Model model)
    {
        createCellTypes( model );
        setupTissue( model );
        model.getVisualizers().forEach( v -> v.setAgentVisualizer( new PPFVisualizer() ) );
    }

    static void createCellTypes(Model model)
    {
        PhysiCellUtilities.setSeed( model.getParameterInt( "random_seed" ) );

        //	cell_defaults.functions.update_phenotype = phenotype_function; 
        //	cell_defaults.functions.custom_cell_rule = custom_function; 

        CellDefinition pFarmerDef = CellDefinition.getCellDefinition( "farmer" );
        CellDefinition pPreyDef = CellDefinition.getCellDefinition( "prey" );
        CellDefinition pPredDef = CellDefinition.getCellDefinition( "predator" );

        pFarmerDef.functions.custom_cell_rule = new FarmerCustomRule();
        pFarmerDef.functions.updatePhenotype = null;
        pFarmerDef.functions.update_migration_bias = null;
        // pFarmerDef.functions.update_phenotype = prey_phenotype_function; 
        // pFarmerDef.functions.update_migration_bias = prey_motility_function; 

        pPreyDef.functions.custom_cell_rule = new PreyCustomRule();
        pPreyDef.functions.updatePhenotype = new PreyPhenotype();
        pPreyDef.functions.update_migration_bias = new PreyMotility();

        pPredDef.functions.custom_cell_rule = new PredatorCustomRule();
        pPredDef.functions.updatePhenotype = new PredatorPhenotype();
        pPredDef.functions.update_migration_bias = new PredatorMotility();
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

        if( microenvironment.options.simulate2D )
        {
            zMin = 0.0;
            zMax = 0.0;
        }

        //        double Xrange = Xmax - Xmin;
        //        double Yrange = Ymax - Ymin;
        //        double Zrange = Zmax - Zmin;

        double[] box = new double[] {xMin, yMin, zMin, xMax, yMax, zMax};
        CellDefinition pCD = CellDefinition.getCellDefinition( "farmer" );
        PhysiCellUtilities.placeInBox( box, pCD, model.getParameterInt( "number_of_farmers" ), microenvironment );

        //        for( int n = 0; n < model.getParameterInt( "number_of_farmers" ); n++ )
        //        {
        //            double[] position = {0, 0, 0};
        //            position[0] = Xmin + PhysiCellUtilities.UniformRandom() * Xrange;
        //            position[1] = Ymin + PhysiCellUtilities.UniformRandom() * Yrange;
        //            position[2] = Zmin + PhysiCellUtilities.UniformRandom() * Zrange;
        //            Cell.createCell( pCD, microenvironment, position );
        //        }

        // place prey 
        pCD = CellDefinition.getCellDefinition( "prey" );
        PhysiCellUtilities.placeInBox( box, pCD, model.getParameterInt( "number_of_prey" ), microenvironment );
        //        for( int n = 0; n < model.getParameterInt( "number_of_prey" ); n++ )
        //        {
        //            double[] position = {0, 0, 0};
        //            position[0] = Xmin + PhysiCellUtilities.UniformRandom() * Xrange;
        //            position[1] = Ymin + PhysiCellUtilities.UniformRandom() * Yrange;
        //            position[2] = Zmin + PhysiCellUtilities.UniformRandom() * Zrange;
        //            Cell.createCell( pCD, microenvironment, position );
        //        }

        // place predators 
        pCD = CellDefinition.getCellDefinition( "predator" );
        PhysiCellUtilities.placeInBox( box, pCD, model.getParameterInt( "number_of_predators" ), microenvironment );
        //        for( int n = 0; n < model.getParameterInt( "number_of_predators" ); n++ )
        //        {
        //            double[] position = {0, 0, 0};
        //            position[0] = Xmin + PhysiCellUtilities.UniformRandom() * Xrange;
        //            position[1] = Ymin + PhysiCellUtilities.UniformRandom() * Yrange;
        //            position[2] = Zmin + PhysiCellUtilities.UniformRandom() * Zrange;
        //            Cell.createCell( pCD, microenvironment, position );
        //        }
    }

    public static List<Cell> get_possible_neighbors(Cell pCell)
    {
        List<Cell> neighbors = new ArrayList<>();
        for( Cell neighbor : pCell.get_container().agent_grid.get( pCell.get_current_mechanics_voxel_index() ) )
        {
            neighbors.add( neighbor );
        }
        for( int ind : pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()] )
        {
            for( Cell neighbor : pCell.get_container().agent_grid.get( ind ) )
            {
                neighbors.add( neighbor );
            }
        }
        return neighbors;
    }
}