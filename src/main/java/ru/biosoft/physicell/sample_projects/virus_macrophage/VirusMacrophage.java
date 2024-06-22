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

public class VirusMacrophage extends Model
{

    @Override
    public void init() throws Exception
    {
        super.init();
        PhysiCellUtilities.setSeed( getParameterInt( "random_seed" ) );
        createCellTypes();
        setupTissue();
        getVisualizers().forEach( v -> v.setAgentVisualizer( new VirusVisualizer( this ) ) );
    }

    void createCellTypes() throws Exception
    {

        // first find index for a few key variables. 
        int virus_index = m.findDensityIndex( "virus" );
        int nInterferon = m.findDensityIndex( "interferon" );

        CellDefinition pEpithelial = getCellDefinition( "epithelial cell" );
        CellDefinition pMacrophage = getCellDefinition( "macrophage" );

        pEpithelial.functions.updatePhenotype = new Epithelial();

        pEpithelial.phenotype.molecular.fractionReleasedDeath[virus_index] = getParameterDouble( "fraction_released_at_death" );
        pEpithelial.phenotype.molecular.fractionTransferredIngested[virus_index] = getParameterDouble(
                "fraction_transferred_when_ingested" );
        /*		
        	pEpithelial.phenotype.molecular.fraction_released_at_death[ nInterferon ] = 0;
        	pEpithelial.phenotype.molecular.fraction_transferred_when_ingested[ nInterferon ] = 0; 		
        */
        pMacrophage.phenotype.mechanics.cellCellAdhesionStrength *= getParameterDouble( "macrophage_relative_adhesion" );
        pMacrophage.phenotype.molecular.fractionReleasedDeath[virus_index] = 0.0;
        pMacrophage.phenotype.molecular.fractionTransferredIngested[virus_index] = 0.0;

        pMacrophage.functions.updatePhenotype = new Macrophage();
        pMacrophage.functions.customCellRule = new AvoidBoundaries();
    }

    void setupTissue()
    {
        int nVirus = m.findDensityIndex( "virus" );
        PhysiCellUtilities.place2D( this, "epithelial cell", getParameterInt( "number_of_infected_cells" ) );
        for( Cell cell : m.getAgents( Cell.class ) )
            cell.phenotype.molecular.internSubstrates[nVirus] = 1;

        PhysiCellUtilities.place2D( this, "epithelial cell", getParameterInt( "number_of_uninfected_cells" ) );
        PhysiCellUtilities.place2D( this, "macrophage", getParameterInt( "number_of_macrophages" ) );

        //        for (Cell cell)
        // create some cells near the origin
        //        CellDefinition pCD = CellDefinition.getCellDefinition( "epithelial cell" );
        //        for( int n = 0; n < getParameterInt( "number_of_infected_cells" ); n++ )
        //        {
        //            PhysiCellUtilities.place( m, getInitialPath(), n );
        //            double x = PhysiCellUtilities.UniformRandom( m.mesh.boundingBox[0], m.mesh.boundingBox[3] );
        //            double y = PhysiCellUtilities.UniformRandom( m.mesh.boundingBox[1], m.mesh.boundingBox[4] );
        //            Cell pC = Cell.createCell( pCD, m, new double[] {x, y, 0.0} );
        //            pC.phenotype.molecular.internSubstrates[nVirus] = 1;
        //        }
        //
        //        for( int n = 0; n < getParameterInt( "number_of_uninfected_cells" ); n++ )
        //        {
        //            double x = PhysiCellUtilities.UniformRandom( m.mesh.boundingBox[0], m.mesh.boundingBox[3] );
        //            double y = PhysiCellUtilities.UniformRandom( m.mesh.boundingBox[1], m.mesh.boundingBox[4] );
        //            Cell.createCell( pCD, m, new double[] {x, y, 0.0} );
        //        }
        //
        //        pCD = CellDefinition.getCellDefinition( "macrophage" );
        //        for( int n = 0; n < getParameterInt( "number_of_macrophages" ); n++ )
        //        {
        //            double x = PhysiCellUtilities.UniformRandom( m.mesh.boundingBox[0], m.mesh.boundingBox[3] );
        //            double y = PhysiCellUtilities.UniformRandom( m.mesh.boundingBox[1], m.mesh.boundingBox[4] );
        //            Cell.createCell( pCD, m, new double[] {x, y, 0.0} );
        //        }
    }

    public static double[] integrateTotalSubstrates(Microenvironment microenvironment)
    {
        double[] out = new double[microenvironment.numberDensities()];
        for( int n = 0; n < microenvironment.numberVoxels(); n++ )
        {
            VectorUtil.axpy( out, microenvironment.mesh.voxels[n].volume, microenvironment.get( n ) );
        }

        for( Cell cell : microenvironment.getAgents( Cell.class ) )
        {
            VectorUtil.sum( out, cell.phenotype.molecular.internSubstrates );
        }
        return out;
    }

}