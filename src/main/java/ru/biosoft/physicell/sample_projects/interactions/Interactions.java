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

import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;

public class Interactions extends Model
{

    @Override
    public void init() throws Exception
    {
        super.init();
        setSeed( getParameterInt( "random_seed" ) );
        createCellTypes();
        setupTissue();
    }

    void createCellTypes()
    {
        CellDefinition pCD = getCellDefinition( "bacteria" );
        pCD.functions.updatePhenotype = new BacteriaPhenotype();
        // pCD.functions.update_migration_bias = advanced_chemotaxis_function; 
        // pCD.phenotype.motility.chemotactic_sensitivity( "resource" ) = 1; 
        // pCD.phenotype.motility.chemotactic_sensitivity( "quorum" ) = 0.1; 

        // set up blood vessels 
        pCD = getCellDefinition( "blood vessel" );
        pCD.functions.updatePhenotype = null;
        pCD.isMovable = false;

        // set up stem cells 
        pCD = getCellDefinition( "stem" );
        pCD.functions.updatePhenotype = new StemPhenotype( );
        // pCD.phenotype.cell_transformations.transformation_rate("differentiated") = 0.0001; 

        // set up differentiated cells 
        pCD = getCellDefinition( "differentiated" );
        pCD.functions.updatePhenotype = new DifferentiatedPhenotype();

        // set up macrophages 
        pCD = getCellDefinition( "macrophage" );
        // pCD.phenotype.cell_interactions.dead_phagocytosis_rate = 0.05; 
        pCD.functions.updatePhenotype = new MacrophagePhenotype();
        // pCD.functions.update_migration_bias = advanced_chemotaxis_function; 
        // pCD.phenotype.motility.chemotactic_sensitivity( "debris" ) = 0.1; 
        // pCD.phenotype.motility.chemotactic_sensitivity( "quorum" ) = 1; 


        // set up CD8+ T cells 
        pCD = getCellDefinition( "CD8+ T cell" );
        pCD.functions.updatePhenotype = new CD8TcellPhenotype();
        // pCD.phenotype.cell_interactions.attack_rate("bacteria") = 0.05; 

        // set up neutrophil  
        pCD = getCellDefinition( "neutrophil" );
        pCD.functions.updatePhenotype = new NeutrophilPhenotype();
        // pCD.phenotype.cell_interactions.live_phagocytosis_rate("bacteria") = 0.05; 

    }

    void setupTissue() throws Exception
    {
        PhysiCellUtilities.place( this, "bacteria", getParameterInt( "number_of_bacteria" ) );
        PhysiCellUtilities.place( this, "blood vessel", getParameterInt( "number_of_blood_vessels" ) );
        PhysiCellUtilities.place( this, "stem", getParameterInt( "number_of_stem_cells" ) );
        PhysiCellUtilities.place( this, "differentiated", getParameterInt( "number_of_differentiated_cells" ) );
        PhysiCellUtilities.place( this, "macrophage", getParameterInt( "number_of_macrophages" ) );
        PhysiCellUtilities.place( this, "neutrophil", getParameterInt( "number_of_neutrophils" ) );
        PhysiCellUtilities.place( this, "CD8+ T cell", getParameterInt( "number_of_CD8T_cells" ) );
    }
}