package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.standard.StandardModels;
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

public class CancerBiorobots extends Model
{

    @Override
    public void init() throws Exception
    {
        super.init();
        signals.setupDictionaries( this );
        createCellTypes();
        setupTissue();
        addEvent( new TherapyEvent( this ) );
        for( Visualizer visualizer : getVisualizers() )
        {
            visualizer.setAgentVisualizer( new CancerBiorobotsVisualizer( this ) );
        }
    }

    private void createCellTypes()
    {
        setSeed( getParameterInt( "random_seed" ) );

        //cancer cell
        CellDefinition cd = getCellDefinition( "cancer cell" );
        cd.functions.updatePhenotype = new TumorPhenotype( this );
        cd.parameters.o2_proliferation_saturation = 38.0;
        cd.parameters.o2_reference = 38.0;

        // cargo cells
        cd = getCellDefinition( "cargo cell" );
        // figure out mechanics parameters
        cd.phenotype.mechanics.relMaxAttachmentDistance = cd.custom_data.get( "max_attachment_distance" ) / cd.phenotype.geometry.radius;
        cd.phenotype.mechanics.relDetachmentDistance = cd.custom_data.get( "max_elastic_displacement" ) / cd.phenotype.geometry.radius;
        cd.phenotype.mechanics.attachmentElasticConstant = cd.custom_data.get( "elastic_coefficient" );

        // set functions
        cd.functions.updatePhenotype = new CargoPhenotype( this );
        cd.functions.customCellRule = new CargoCellRule( this );
        cd.functions.contact = new BiorobotsContact();
        cd.functions.updateMigration = null;

        // worker cells
        cd = getCellDefinition( "worker cell" );
        cd.phenotype.mechanics.relMaxAttachmentDistance = cd.custom_data.get( "max_attachment_distance" ) / cd.phenotype.geometry.radius;
        cd.phenotype.mechanics.relDetachmentDistance = cd.custom_data.get( "max_elastic_displacement" ) / cd.phenotype.geometry.radius;
        cd.phenotype.mechanics.attachmentElasticConstant = cd.custom_data.get( "elastic_coefficient" );
        cd.functions.updatePhenotype = null; // worker_cell_rule;
        cd.functions.customCellRule = new WorkerCellRule( this );
        cd.functions.contact = new BiorobotsContact();
    }

    private void setupTissue() throws Exception
    {
        CellDefinition defaults = StandardModels.getDefaultCellDefinition();
        double cellRadius = defaults.phenotype.geometry.radius;
        double cellSpacing = 0.95 * 2.0 * cellRadius;
        double tumorRadius = getParameterDouble( "tumor_radius" ); // 200.0;

        CellDefinition cd = getCellDefinition( "cancer cell" );
        double x = 0.0;
        double xOuter = tumorRadius;
        double y = 0.0;

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
                Cell.createCell( cd, this, new double[] {x, y, 0.0} ); // tumor cell

                if( Math.abs( y ) > 0.01 )
                    Cell.createCell( cd, this, new double[] {x, -y, 0.0} ); // tumor cell			

                if( Math.abs( x ) > 0.01 )
                {
                    Cell.createCell( cd, this, new double[] { -x, y, 0.0} );

                    if( Math.abs( y ) > 0.01 )
                        Cell.createCell( cd, this, new double[] { -x, -y, 0.0} );
                }
                x += cellSpacing;
            }

            y += cellSpacing * Math.sqrt( 3.0 ) / 2.0;
            n++;
        }
    }
}