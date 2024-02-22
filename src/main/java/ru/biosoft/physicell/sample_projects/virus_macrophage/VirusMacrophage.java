package ru.biosoft.physicell.sample_projects.virus_macrophage;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
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

public class VirusMacrophage
{


    public static void init(Model model) throws Exception
    {
        PhysiCellUtilities.setSeed( model.getParameterInt( "random_seed" ) );
        createCellTypes( model );
        setupTissue( model );
        model.getVisualizers().forEach( v -> v.setAgentVisualizer( new VirusVisualizer( model ) ) );
    }

    static void createCellTypes(Model model) throws Exception
    {
        Microenvironment microenvironment = model.getMicroenvironment();

        // first find index for a few key variables. 
        int virus_index = microenvironment.findDensityIndex( "virus" );
        int nInterferon = microenvironment.findDensityIndex( "interferon" );

        CellDefinition pEpithelial = CellDefinition.getCellDefinition( "epithelial cell" );
        CellDefinition pMacrophage = CellDefinition.getCellDefinition( "macrophage" );

        pEpithelial.functions.updatePhenotype = new Epithelial();

        pEpithelial.phenotype.molecular.fraction_released_at_death[virus_index] = model.getParameterDouble( "fraction_released_at_death" );
        pEpithelial.phenotype.molecular.fraction_transferred_when_ingested[virus_index] = model
                .getParameterDouble( "fraction_transferred_when_ingested" );
        /*		
        	pEpithelial.phenotype.molecular.fraction_released_at_death[ nInterferon ] = 0;
        	pEpithelial.phenotype.molecular.fraction_transferred_when_ingested[ nInterferon ] = 0; 		
        */
        pMacrophage.phenotype.mechanics.cellCellAdhesionStrength *= model.getParameterDouble( "macrophage_relative_adhesion" );
        pMacrophage.phenotype.molecular.fraction_released_at_death[virus_index] = 0.0;
        pMacrophage.phenotype.molecular.fraction_transferred_when_ingested[virus_index] = 0.0;

        pMacrophage.functions.updatePhenotype = new Macrophage();
        pMacrophage.functions.customCellRule = new AvoidBoundaries();
    }

    static void setupTissue(Model model)
    {
        Microenvironment microenvironment = model.getMicroenvironment();
        int nVirus = microenvironment.findDensityIndex( "virus" );
        // create some cells near the origin

        double length_x = microenvironment.mesh.boundingBox[3] - microenvironment.mesh.boundingBox[0];

        double length_y = microenvironment.mesh.boundingBox[4] - microenvironment.mesh.boundingBox[1];

        int number_of_infected_cells = model.getParameterInt( "number_of_infected_cells" );

        CellDefinition pCD = CellDefinition.getCellDefinition( "epithelial cell" );

        for( int n = 0; n < number_of_infected_cells; n++ )
        {
            double x = microenvironment.mesh.boundingBox[0] + PhysiCellUtilities.UniformRandom() * length_x;
            double y = microenvironment.mesh.boundingBox[1] + PhysiCellUtilities.UniformRandom() * length_y;
            Cell pC = Cell.createCell( pCD, microenvironment, new double[] {x, y, 0.0} );
            pC.phenotype.molecular.internalized_total_substrates[nVirus] = 1;
        }

        int number_of_uninfected_cells = model.getParameterInt( "number_of_uninfected_cells" );

        for( int n = 0; n < number_of_uninfected_cells; n++ )
        {
            double x = microenvironment.mesh.boundingBox[0] + PhysiCellUtilities.UniformRandom() * length_x;
            double y = microenvironment.mesh.boundingBox[1] + PhysiCellUtilities.UniformRandom() * length_y;
            Cell.createCell( pCD, microenvironment, new double[] {x, y, 0.0} );
        }

        pCD = CellDefinition.getCellDefinition( "macrophage" );

        for( int n = 0; n < model.getParameterInt( "number_of_macrophages" ); n++ )
        {
            double x = microenvironment.mesh.boundingBox[0] + PhysiCellUtilities.UniformRandom() * length_x;
            double y = microenvironment.mesh.boundingBox[1] + PhysiCellUtilities.UniformRandom() * length_y;
            Cell.createCell( pCD, microenvironment, new double[] {x, y, 0.0} );
        }
    }

    public static double[] integrate_total_substrates(Microenvironment microenvironment)
    {
        // start with 0 vector 
        //	std::vector<double> out( microenvironment.number_of_densities() , 0.0 ); 
        //        List<Double> out = new ArrayList<>();
        double[] out = new double[microenvironment.numberDensities()];
        // integrate extracellular substrates 
        for( int n = 0; n < microenvironment.numberVoxels(); n++ )
        {
            // out = out + microenvironment(n) * dV(n) 
            VectorUtil.axpy( out, microenvironment.mesh.voxels[n].volume, microenvironment.get( n ) );
        }

        // inte
        for( Cell cell : microenvironment.getAgents( Cell.class ) )//int n = 0; n < ( all_cells ).size(); n++ )
        {
            VectorUtil.sum( out, cell.phenotype.molecular.internalized_total_substrates );
        }

        return out;
    }

}