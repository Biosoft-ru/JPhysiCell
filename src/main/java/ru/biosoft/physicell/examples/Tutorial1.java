package ru.biosoft.physicell.examples;

import ru.biosoft.physicell.biofvm.BasicAgent;
import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.ui.Visualizer;
import ru.biosoft.physicell.ui.Visualizer.Section;

/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.0.4) [1]        #
#                                                                           #
# [1] A. Ghaffarizaeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient  #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics, 2015. DOI: 10.1093/bioinformatics/btv730               #
#############################################################################
#                                                                           #
# Copyright 2015 Paul Macklin and the BioFVM Project                        #
#                                                                           #
# Licensed under the Apache License, Version 2.0 (the "License");           #
# you may not use this file except in compliance with the License.          #
# You may obtain a copy of the License at                                   #
#                                                                           #
#    http://www.apache.org/licenses/LICENSE-2.0                             #
#                                                                           #
# Unless required by applicable law or agreed to in writing, software       #
# distributed under the License is distributed on an "AS IS" BASIS,         #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #
# See the License for the specific language governing permissions and       #
# limitations under the License.                                            #
#############################################################################
*/
public class Tutorial1
{
    //Grid options
    private static int size = 1000; // x size = y size = z size
    private static int cellSize = 10; // x size = y size = z size

    //Substrate options
    private static int diffusionCoefficients = 1000;// microns^2 / min 
    private static double decayRate = 0.01;
    private static double saturationDensity = 1;

    //Agents options
    private static int sourceAgentsNumber = 500;
    private static int sinkAgentsNumber = 500;
    private static double agentRadius = 5;
    private static boolean generateAtSlice = true; //create agents only at slice we watch
    private static int agentSecretionRate = 10;
    private static double agentUptakeRate = 0.8;

    //Simulation options
    private static double tMax = 10;
    private static double dt = 0.01;
    private static String resultPath = "C:/Users/Damag/BIOFVM/Tutorial1";
    private static int zSlice = 500; //watch slice in the middle (for images)
    private static double outputInterval = 0.1; // time delta between images
    private static boolean outputTables = false;

    static Visualizer visualizer = new Visualizer( resultPath, "Z5003", Section.Z, zSlice );
    static
    {
        visualizer.setSaveImage( false );
        visualizer.setSaveGIF( true );
    }

    public static void main(String ... args) throws Exception
    {
        double timeStart = System.currentTimeMillis();
        Microenvironment m = new Microenvironment( "substrate scale", size, cellSize, "minutes", "microns" );
        //        m.setDensity( 0, "substrate1", "dimensionless" );
        m.setDensity( 0, "substrate1", "dimensionless", diffusionCoefficients, decayRate );

        //set initial density - more in the middle, gradually decrease
        double[] center = PhysiCellUtilities.getCenter( m.mesh );
        double stddevSquared = size / cellSize;
        stddevSquared *= -stddevSquared;
        for( int i = 0; i < m.number_of_voxels(); i++ )
        {
            double[] displacement = VectorUtil.newDiff( m.voxels( i ).center, center );
            double coeff = VectorUtil.norm_squared( displacement ) / stddevSquared;
            m.getDensity( i )[0] = Math.exp( coeff );
        }

        // register substrates properties 
        //        m.diffusion_coefficients[0] = diffusionCoefficients; // microns^2 / min 
        //        m.decay_rates[0] = decayRate;
        m.options.track_internalized_substrates_in_each_agent = false;

        //populate model
        for( int i = 0; i < sourceAgentsNumber; i++ )
        {
            //substrate source agents
            double[] pos = VectorUtil.random( 3, 0, size );
            if( generateAtSlice )
                pos[2] = zSlice; //keep agents at z slice

            BasicAgent agentSource = BasicAgent.createBasicAgent( m );
            agentSource.registerMicroenvironment( m );
            agentSource.assignPosition( pos );
            agentSource.setRadius( agentRadius );
            agentSource.secretionRates[0] = agentSecretionRate;
            agentSource.saturationDensities[0] = saturationDensity;
            agentSource.setUptakeConstants( dt );
            agentSource.setTag( "Source" );
        }

        for( int i = 0; i < sinkAgentsNumber; i++ )
        {
            double[] pos = VectorUtil.random( 3, 0, size );
            if( generateAtSlice )
                pos[2] = zSlice; //keep agents at z slice 

            BasicAgent agentSink = BasicAgent.createBasicAgent( m );
            agentSink.registerMicroenvironment( m );
            agentSink.assignPosition( pos );
            agentSink.setRadius( agentRadius );
            agentSink.uptakeRates[0] = agentUptakeRate;
            agentSink.setUptakeConstants( dt );
            agentSink.setTag( "Sink" );
        }
        //visualization of all layers
        //        for( int z = 0; z < 100; z++ )
        //            SimpleVisualizer.draw( m, z, 1000, 1000, resultPath + "/init_slice" + z + ".png" );

        //Simulation starts
        visualizer.init();

        if( outputTables )
            m.writeDensity( resultPath + "/initial.txt" );

        visualizer.saveResult( m, 0 );
        double t = 0.0;
        double nextOutputTime = outputInterval;
        double timeSimulationStart = System.currentTimeMillis();
        int counter = 0;
        while( t < tMax )
        {
            m.simulate_cell_sources_and_sinks( dt );
            m.simulate_diffusion_decay( dt );
            t += dt;
            counter++;
            if( Math.abs( t - nextOutputTime ) < 0.0001 )
            {
                if( outputTables )
                    m.writeDensity( resultPath + "/step_" + counter + ".txt" );
                visualizer.saveResult( m, t );
                nextOutputTime += outputInterval;
            }
        }
        visualizer.finish();
        double totalTime = ( System.currentTimeMillis() - timeStart ) / 1000;
        double simulationTime = ( System.currentTimeMillis() - timeSimulationStart ) / 1000;
        System.out.println( "Done! Elapsed time: " + totalTime + " s. Simulation time: " + simulationTime + " s." );
    }
}