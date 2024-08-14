package ru.biosoft.physicell.sample_projects.ode_energy;

import java.util.ArrayList;
import java.util.List;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Intracellular;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.standard.StandardUpdateVelocity;
import ru.biosoft.physicell.core.standard.StandardVolumeUpdate;

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

public class OdeEnergy extends Model
{

    @Override
    public void init() throws Exception
    {
        super.init();
        setSeed( getParameterInt( "random_seed" ) );
        createCellTypes();
        setupTissue();
        getVisualizers().forEach( v -> v.setAgentVisualizer( new EnergyVisualizer() ) );
    }

    public void createCellTypes() throws Exception
    {
        CellDefinition defaults = getCellDefinition( "default" );
        defaults.functions.updateVolume = new StandardVolumeUpdate();
        defaults.functions.updateVelocity = new StandardUpdateVelocity();
        defaults.functions.updateMigration = null;
        defaults.functions.updatePhenotype = null; // update_cell_and_death_parameters_O2_based; 
        defaults.functions.customCellRule = null;
        defaults.functions.membraneInteraction = null;
        defaults.functions.membraneDistanceCalculator = null;
        signals.setupDictionaries( this );
    }

    void setupTissue() throws Exception
    {
        int oxygenIndex = m.findDensityIndex( "oxygen" );
        int glucoseIndex = m.findDensityIndex( "glucose" );
        int lactateIndex = m.findDensityIndex( "lactate" );

        CellDefinition cd = getCellDefinition( "default" );
        double cellRadius = cd.phenotype.geometry.radius;
        double initialTumorRadius = 100;

        List<double[]> positions = createCirclePositions( cellRadius, initialTumorRadius );
        for( int i = 0; i < positions.size(); i++ )
        {
            Cell pCell = Cell.createCell( cd, this, positions.get( i ) );
            signals.setSingleBehavior( pCell, "custom:intra_oxy", getParameterDouble( "initial_internal_oxygen" ) );
            signals.setSingleBehavior( pCell, "custom:intra_glu", getParameterDouble( "initial_internal_glucose" ) );
            signals.setSingleBehavior( pCell, "custom:intra_lac", getParameterDouble( "initial_internal_lactate" ) );
            signals.setSingleBehavior( pCell, "custom:intra_energy", getParameterDouble( "initial_energy" ) );
            double cellVolume = pCell.phenotype.volume.total;
            double[] substrates = pCell.phenotype.molecular.internSubstrates;
            substrates[oxygenIndex] = signals.getSingleSignal( pCell, "custom:intra_oxy" ) * cellVolume;
            substrates[glucoseIndex] = signals.getSingleSignal( pCell, "custom:intra_glu" ) * cellVolume;
            substrates[lactateIndex] = signals.getSingleSignal( pCell, "custom:intra_lac" ) * cellVolume;
            pCell.phenotype.intracellular.start();
            pCell.phenotype.intracellular.setParameterValue( "Energy", signals.getSingleSignal( pCell, "custom:intra_energy" ) );
        }
    }

    @Override
    public void updateIntracellular() throws Exception
    {
        int oxygenIndex = m.findDensityIndex( "oxygen" );
        int glucoseIndex = m.findDensityIndex( "glucose" );
        int lactateIndex = m.findDensityIndex( "lactate" );

        m.getAgents( Cell.class ).parallelStream().filter( cell -> !cell.isOutOfDomain ).forEach( cell -> {
            try
            {
                double cellVolume = cell.phenotype.volume.total;

                // Get Intracellular Concentrations
                double cellOxygen = signals.getSingleSignal( cell, "intracellular oxygen" );
                double cellGlucose = signals.getSingleSignal( cell, "intracellular glucose" );
                double cellLactate = signals.getSingleSignal( cell, "intracellular lactate" );

                Intracellular intra = cell.phenotype.intracellular;
                // Update SBML 
                intra.setParameterValue( "Oxygen", cellOxygen );
                intra.setParameterValue( "Glucose", cellGlucose );
                intra.setParameterValue( "Lactate", cellLactate );

                // SBML Simulation
                intra.step();

                // Phenotype Simulation
                intra.updatePhenotypeParameters( m, cell.phenotype );

                double[] substrates = cell.phenotype.molecular.internSubstrates;
                // Internalized Chemical Update After SBML Simulation
                substrates[oxygenIndex] = intra.getParameterValue( "Oxygen" ) * cellVolume;
                substrates[glucoseIndex] = intra.getParameterValue( "Glucose" ) * cellVolume;
                substrates[lactateIndex] = intra.getParameterValue( "Lactate" ) * cellVolume;

                //Save custom data
                signals.setSingleBehavior( cell, "custom:intra_oxy", intra.getParameterValue( "Oxygen" ) );
                signals.setSingleBehavior( cell, "custom:intra_glu", intra.getParameterValue( "Glucose" ) );
                signals.setSingleBehavior( cell, "custom:intra_lac", intra.getParameterValue( "Lactate" ) );
                signals.setSingleBehavior( cell, "custom:intra_energy", intra.getParameterValue( "Energy" ) );
            }
            catch( Exception ex )
            {
                ex.printStackTrace();
            }
        } );
    }

    public static List<double[]> createCirclePositions(double cellRadius, double sphereRadius)
    {
        List<double[]> result = new ArrayList<>();
        int xc = 0;
        double xSpacing = cellRadius * Math.sqrt( 3 );
        double ySpacing = cellRadius * Math.sqrt( 3 );

        for( double x = -sphereRadius; x < sphereRadius; x += xSpacing, xc++ )
        {
            for( double y = -sphereRadius; y < sphereRadius; y += ySpacing )
            {
                double[] tempPoint = new double[3];
                tempPoint[1] = y + ( xc % 2 ) * cellRadius;
                tempPoint[0] = x;
                tempPoint[2] = 0;
                if( Math.sqrt( VectorUtil.norm_squared( tempPoint ) ) < sphereRadius )
                {
                    result.add( tempPoint );
                }
            }
        }
        return result;
    }

    @Override
    public String getReport(Cell cell) throws Exception
    {
        return "\n" + cell.ID + "\t" + cell.isOutOfDomain + "\t" + VectorUtil.print( cell.position ) + "\t" + cell.phenotype.cycle.name
                + "\t" + cell.phenotype.cycle.currentPhase().name + "\t" + signals.getSingleSignal( cell, "intracellular lactate" ) + "\t"
                + signals.getSingleSignal( cell, "intracellular glucose" ) + "\t" + signals.getSingleSignal( cell, "intracellular oxygen" )
                + "\t"
                + cell.phenotype.intracellular.getParameterValue( "Energy" );
    }

    @Override
    public String getReportHeader()
    {
        return "ID\tLactate\tGlucose\tOxygen\tEnergy";
    }

    @Override
    public String getLogInfo() throws Exception
    {
        Cell cell = m.getAgents( Cell.class ).iterator().next();
        String addon = signals.getSingleSignal( cell, "intracellular lactate" ) + "\t"
                + signals.getSingleSignal( cell, "intracellular glucose" ) + "\t" + signals.getSingleSignal( cell, "intracellular oxygen" )
                + "\t"
                + cell.phenotype.intracellular.getParameterValue( "Energy" );

        return PhysiCellUtilities.getCurrentTime() + "\tElapsed\t" + ( System.currentTimeMillis() - startTime ) / 1000 + "\tTime:\t"
                + (int)Math.round( curTime ) + "\tCells\t" + m.getAgentsCount() + "\t" + addon;
    }

}