package ru.biosoft.physicell.sample_projects.celltypes3;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
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
public class Celltype3
{
    public static void init(Model model) throws Exception
    {
        PhysiCellUtilities.setSeed( model.getParameterInt( "random_seed" ) );
        SignalBehavior.setupDictionaries( model.getMicroenvironment() );
        createCellTypes( model );
        setupTissue( model );
        for( Visualizer visualizer : model.getVisualizers() )
        {
            visualizer.setAgentVisualizer( new RegularAgentVisualizer( model ) );
        }
    }

    /**
     * Modifies CellDefinitions
     */
    static void createCellTypes(Model model) throws Exception
    {
        CellDefinition.getCellDefinition( "A" ).functions.updatePhenotype = new A_phenotype( model );
        CellDefinition.getCellDefinition( "B" ).functions.updatePhenotype = new B_phenotype( model );
        CellDefinition.getCellDefinition( "C" ).functions.updatePhenotype = new C_phenotype( model );
    }

    /**
     * Creates all cells for the model
     */
    static void setupTissue(Model model)
    {
        Microenvironment m = model.getMicroenvironment();
        double xMin = m.mesh.boundingBox[0];
        double yMin = m.mesh.boundingBox[1];
        double zMin = m.mesh.boundingBox[2];
        double xMax = m.mesh.boundingBox[3];
        double yMax = m.mesh.boundingBox[4];
        double zMax = m.mesh.boundingBox[5];

        double maxRadius = model.getParameterDouble( "max_distance_from_origin" );
        xMax = Math.min( xMax, maxRadius );
        xMin = Math.max( xMin, -maxRadius );
        yMax = Math.min( yMax, maxRadius );
        yMin = Math.max( yMin, -maxRadius );
        zMax = Math.min( zMax, maxRadius );
        zMin = Math.max( zMin, -maxRadius );

        if( m.options.simulate2D )
        {
            zMin = 0.0;
            zMax = 0.0;
        }

        double[] range = new double[] {xMin, yMin, zMin, xMax, yMax, zMin, zMax};

        CellDefinition aCD = CellDefinition.getCellDefinition( "A" );
        CellDefinition bCD = CellDefinition.getCellDefinition( "B" );
        CellDefinition cCD = CellDefinition.getCellDefinition( "C" );

        int number = model.getParameterInt( "number_of_A" );
        placeInRadius( aCD, m, number, range, maxRadius );

        number = model.getParameterInt( "number_of_B" );
        placeInRadius( bCD, m, number, range, maxRadius );

        number = model.getParameterInt( "number_of_C" );
        placeInRadius( cCD, m, number, range, maxRadius );

        for( Cell cell : m.getAgents( Cell.class ) )
        {
            for( int k = 0; k < cell.phenotype.death.rates.size(); k++ )
            {
                cell.phenotype.death.rates.set( k, 0.0 );
            }
        }
    }

    /**
     * Places cells in given range but not exceeding given maxRadius 
     * @param cd - CellDefinition for cells
     * @param m - microenvironment
     * @param number - number of cells to place
     * @param range - bounding box: xMin, yMin, zMin, xMax, yMax, zMax
     * @param maxRadius - maximum radius
     */
    private static void placeInRadius(CellDefinition cd, Microenvironment m, int number, double[] range, double maxRadius)
    {
        for( int n = 0; n < number; n++ )
        {
            double[] position = {0, 0, 0};
            double r = maxRadius + 1;
            while( r > maxRadius )
            {
                position[0] = PhysiCellUtilities.UniformRandom( range[0], range[3] );
                position[1] = PhysiCellUtilities.UniformRandom( range[1], range[4] );
                position[2] = PhysiCellUtilities.UniformRandom( range[2], range[5] );
                r = VectorUtil.norm( position );
            }
            Cell.createCell( cd, m, position );
        }
    }
}