package ru.biosoft.physicell.sample_projects.cancer_metabolism;

import java.util.ArrayList;
import java.util.List;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.Secretion;
import ru.biosoft.physicell.core.standard.StandardModels;
import ru.biosoft.physicell.core.standard.StandardVolumeUpdate;
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
public class CancerMetabolism extends Model
{
    @Override
    public void init() throws Exception
    {
        super.init();
        createCellTypes();
        setupTissue();
        for( Visualizer visualizer : getVisualizers() )
            visualizer.setAgentVisualizer( new CancerMetabolismVisualizer() );
    }

    void createCellTypes() throws Exception
    {
        CellDefinition cellDefaults = StandardModels.getDefaultCellDefinition();
        cellDefaults.parameters.o2_proliferation_saturation = 38.0;
        cellDefaults.parameters.o2_reference = 38.0;

        cellDefaults.functions.updatePhenotype = null;//update_cell;
        cellDefaults.functions.updateVolume = new StandardVolumeUpdate();
        cellDefaults.functions.updateVelocity = null;
        cellDefaults.functions.updateMigration = null;
        cellDefaults.functions.customCellRule = null;
    }

    List<double[]> createSpherePositions(double cell_radius, double sphere_radius)
    {
        List<double[]> cells = new ArrayList<>();
        int xc = 0, yc = 0, zc = 0;
        double x_spacing = cell_radius * Math.sqrt( 3 );
        double y_spacing = cell_radius * 2;
        double z_spacing = cell_radius * Math.sqrt( 3 );

        //        for( double z = -sphere_radius; z < sphere_radius; z += z_spacing, zc++ )
        //        {
            for( double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++ )
            {
                for( double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++ )
                {
                    double[] position = new double[3];
                    position[0] = x + ( zc % 2 ) * 0.5 * cell_radius;
                    position[1] = y + ( xc % 2 ) * cell_radius;
                    position[2] = 0;

                    if( Math.sqrt( VectorUtil.norm_squared( position ) ) < sphere_radius )
                        cells.add( position );
                }
            }
            //        }
        return cells;

    }

    void setupTissue() throws Exception
    {
        CellDefinition defaults = StandardModels.getDefaultCellDefinition();
        double cellRadius = defaults.phenotype.geometry.radius;
        double tumorRadius = getParameterDouble( "tumor_radius" ); // 250.0; 
        List<double[]> positions = createSpherePositions( cellRadius, tumorRadius );

        CellDefinition cd = getCellDefinition( "metabolic cell" );
        //        Cell.createCell( cd, m, new double[] {0, 0, 0} );
        for( int i = 0; i < positions.size(); i++ )
            Cell.createCell( cd, this, positions.get( i ) );
    }

    public void updateIntracellular() throws Exception
    {
        for( Cell cell : m.getAgents( Cell.class ) )
            updateCell( cell, cell.phenotype, diffusion_dt );
    }

    void updateCell(Cell cell, Phenotype phenotype, double dt)
    {
        phenotype.intracellular.update( cell, phenotype, dt );
    }

    void setup_default_metabolic_model()
    {
        return;
    }

    void anuclear_volume_model(Cell pCell, Phenotype phenotype, double dt)
    {
        return;
    }

    void metabolic_cell_phenotype(Cell pCell, Phenotype phenotype, double dt)
    {
        // if cell is dead, don't bother with future phenotype changes.
        if( phenotype.death.dead )
        {
            pCell.functions.updatePhenotype = null;
            return;
        }

        // update the transition rate according to growth rate?
        int cycle_start_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
        int cycle_end_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );

        //static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" );
        //phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
    }

    @Override
    public String getReportHeader()
    {
        return "ID\tX\tY\tZ\tvoxel\toxygen\tglucose\tacetate\toxygen_conc\tglucose_conc\tacetate_conc\tcycle";
    }

    @Override
    public String getReport(Cell cell) throws Exception
    {
        Microenvironment m = cell.getMicroenvironment();
        int oxygen_idx = m.findDensityIndex( "oxygen" );
        int glucose_idx = m.findDensityIndex( "glucose" );
        int lactate_idx = m.findDensityIndex( "lactate" );

        Secretion secretion = cell.phenotype.secretion;
        double uptakeOxygen = secretion.uptakeRates[oxygen_idx];
        double uptakeGlucose = secretion.uptakeRates[glucose_idx];
        double uptakeAcetate = secretion.uptakeRates[lactate_idx];

        int voxelIndex = cell.currentVoxelIndex;
        double[] density = m.nearestDensity( voxelIndex );

        return "\n" + cell.ID + "\t" + cell.position[0] + "\t" + cell.position[1] + "\t" + cell.position[2] + "\t" + voxelIndex + "\t"
                + uptakeOxygen + "\t" + uptakeGlucose + "\t" + uptakeAcetate + "\t" + density[oxygen_idx] + "\t" + density[glucose_idx]
                + "\t" + density[lactate_idx] + "\t" + cell.phenotype.cycle.name;
    }
}